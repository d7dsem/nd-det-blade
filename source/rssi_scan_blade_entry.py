#!/usr/bin/env python3
"""RSSI sweep for bladeRF 2.0 using AD9361 hardware (symbol) RSSI.

CLI:
  --freq-spans "[{flo~fhi,step,wdt},...]"  Hz, suffixes k/M/G; wdt = sample rate (+ HW bandwidth)
  --dwell <sec>                            minimum time on each frequency, counted from the end
                                           of the retune (retune time not included); RSSI results
                                           are read as often as they update until dwell is used up
                                           (checked after each read -> at least one read)
  [--rssi-n <int> | --rssi-dur <sec>]      optional: on-board (AD9361) RSSI averaging length in
                                           RX samples, or its duration -> nearest representable N.
                                           Not given or not valid -> board default is used.
  -g, --gain <dB>                          manual RX gain (default 42)
  --warmup {short,pass,none}               before measuring: short = a few retunes + reads per
                                           span, pass = one full silent pass, none (default short)
  --metrics <file>                         performance metrics (TSV) for analysis: one row per
                                           frequency per pass, plus '# event' lines with timestamps
                                           of startup steps, individual startup retunes and span
                                           reconfiguration (read with comment='#' to skip them)
  --gui                                    Tk window: continuous sweeps, plot after each pass

Data -> stdout, info/warnings -> stderr.
Ctrl+C: stops after the current RSSI read / libbladeRF call and closes the device cleanly;
a second Ctrl+C forces exit.
"""
import argparse
import bisect
import gc
from itertools import combinations_with_replacement
import math
from dataclasses import dataclass
import os
import queue
import re
import signal
import sys
import threading
import time

import numpy as np

import bladerf
from bladerf import _bladerf
# Low-level handles for the hot path: bladerf.get_rfic_rssi() allocates two cffi buffers
# and a namedtuple on every call; we call libbladeRF directly with preallocated buffers.
from bladerf._bladerf import ffi, libbladeRF, _check_error

CALIBRATION_READS = 50
CALIBRATION_TUNES = 5
WARMUP_TUNES = 10            # short warmup: retunes per span
WARMUP_READS = 20            # short warmup: reads after the retunes, per span
LATE_GAP_FACTOR = 1.5        # read gap > LATE_GAP_FACTOR * period counts as a late read (assumption)

RETUNE_ANOMALY_FACTOR = 3.0  # retune > factor * median retune of the pass = anomaly (assumption)
SLOW_CALL_S = 5e-3           # RSSI read call longer than this = anomaly (assumption)

METRIC_COLS = ("pass", "t_ms", "freq_MHz", "wdt_MHz", "reconf_ms", "jump_MHz", "retune_ms", "set_ms",
               "readback_ms", "first_read_ms", "reads", "mean_gap_ms", "max_gap_ms", "late_reads",
               "max_call_ms", "sum_call_ms", "calls_gt1ms", "calls_gt5ms", "max_sched_late_ms",
               "gc_count", "gc_max_ms", "used_ms")
# Pause after the retune has finished, before the first read, in hardware RSSI periods:
# one full RSSI window + margin (assumption).
SETTLE_PERIODS = 2

# AD9361 RSSI registers (no-OS ad9361 driver, ad9361_rssi_setup()).
# Averaging length = up to 4 segments of 2^k RX samples (k = 0..14), weights sum to 255.
REG_RSSI_DUR_01 = 0x150      # [7:4] duration exp of segment 1, [3:0] segment 0
REG_RSSI_DUR_23 = 0x151      # [7:4] segment 3, [3:0] segment 2
REG_RSSI_WEIGHT_0 = 0x152    # 0x152..0x155: weights of segments 0..3
REG_RSSI_WAIT = 0x157        # wait between measurements, units of 4 samples
REG_RSSI_CONFIG = 0x158      # bit0: default RSSI meas mode (single power-of-two segment)
RSSI_CFG_DEFAULT_MODE = 0x01
RSSI_SEG_MAX_EXP = 14
RSSI_MAX_SEGMENTS = 4

_SUFFIX = {"": 1.0, "k": 1e3, "K": 1e3, "M": 1e6, "G": 1e9}


def info(msg):
    print(msg, file=sys.stderr, flush=True)


HZ_TO_MHZ = 1e-6  # single scale for frequency display (text and plot axis)


def fmt_freq(hz):
    """Single formatter for every frequency-like value (Hz) in text output -> MHz."""
    return f"{hz * HZ_TO_MHZ:6.2f}"


# Output table: one spec for header and rows -> columns always aligned
TABLE_COLS = (("freq_MHz", 9), ("wdt_MHz", 8), ("rssi_dB", 8), ("reads", 6), ("std_dB", 7))


def fmt_row(values):
    return "".join(f"{v:>{w}}" for v, (_, w) in zip(values, TABLE_COLS))


class ConfigError(Exception):
    pass


# ---------------------------------------------------------------- CLI parsing

def parse_hz(text):
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)\s*([kKMG]?)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError(f"bad frequency value: {text!r}")
    return float(m.group(1)) * _SUFFIX[m.group(2)]


def parse_spans(text):
    body = text.strip()
    if body.startswith("[") and body.endswith("]"):
        body = body[1:-1]
    items = re.findall(r"\{([^{}]*)\}", body)
    if not items or re.sub(r"\{[^{}]*\}", "", body).strip(" ,"):
        raise argparse.ArgumentTypeError(f"bad --freq-spans: {text!r}")

    spans = []
    for item in items:
        parts = [p.strip() for p in item.split(",")]
        if len(parts) != 3 or "~" not in parts[0]:
            raise argparse.ArgumentTypeError(f"span must be {{flo~fhi,step,wdt}}: {{{item}}}")
        lo, hi = (parse_hz(x) for x in parts[0].split("~", 1))
        step, wdt = parse_hz(parts[1]), parse_hz(parts[2])
        if hi < lo or step <= 0 or wdt <= 0:
            raise argparse.ArgumentTypeError(f"bad span values: {{{item}}}")
        count = int(math.floor((hi - lo) / step + 1e-9)) + 1
        spans.append((lo + step * np.arange(count, dtype=np.float64), wdt))
    return spans


def positive_float(text):
    v = float(text)
    if v <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return v


def parse_args():
    p = argparse.ArgumentParser(description="bladeRF 2.0 hardware RSSI sweep")
    p.add_argument("--freq-spans", required=True, type=parse_spans,
                   help='e.g. "[{2.40G~2.48G,5M,10M},{900M~930M,1M,2M}]"')
    p.add_argument("--dwell", required=True, type=positive_float,
                   help="minimum time per frequency after the retune, s (RSSI is read until it is used up)")
    g = p.add_mutually_exclusive_group(required=False)
    g.add_argument("--rssi-n", type=int,
                   help="on-board RSSI averaging length, RX samples (optional)")
    g.add_argument("--rssi-dur", type=float,
                   help="on-board RSSI averaging duration, s (optional)")
    p.add_argument("-g", "--gain", type=int, default=42, help="manual RX gain, dB (default 42)")
    p.add_argument("--gui", action="store_true",
                   help="Tk window, continuous sweeps until the window is closed")
    p.add_argument("--warmup", choices=("short", "pass", "none"), default="short",
                   help="warmup before measuring (default short)")
    p.add_argument("--metrics", metavar="FILE",
                   help="write per-frequency performance metrics (TSV) to FILE")
    return p.parse_args()


# ------------------------------------------------------------------ timing

def wait_until(t_target):
    """Coarse sleep, then spin on perf_counter for sub-ms accuracy."""
    while True:
        left = t_target - time.perf_counter()
        if left <= 0:
            return
        if left > 2e-3:
            time.sleep(left - 1.5e-3)


def db_to_lin(db):
    return 10.0 ** (db / 10.0)


def lin_to_db(lin):
    return 10.0 * math.log10(lin)


# ------------------------------------------------- hardware RSSI averaging

def _build_rssi_table():
    """All averaging lengths the AD9361 can represent: sums of <= 4 powers of two (2^0..2^14).
    Returns (sorted totals, {total: exponents with the fewest segments})."""
    best = {}
    for k in range(1, RSSI_MAX_SEGMENTS + 1):
        for exps in combinations_with_replacement(range(RSSI_SEG_MAX_EXP, -1, -1), k):
            total = sum(1 << e for e in exps)
            if total not in best or len(exps) < len(best[total]):
                best[total] = exps
    return sorted(best), best


RSSI_TOTALS, RSSI_SEGMENTS = _build_rssi_table()
RSSI_MAX_SAMPLES = RSSI_TOTALS[-1]


def nearest_rssi_samples(n):
    """Nearest representable averaging length (samples) to n, and its segment exponents."""
    i = bisect.bisect_left(RSSI_TOTALS, n)
    total = min(RSSI_TOTALS[max(i - 1, 0):i + 1], key=lambda t: (abs(t - n), t))
    return total, RSSI_SEGMENTS[total]


def rssi_register_values(exps):
    """Register values for the given segment exponents, as the AD9361 driver computes them."""
    exps = sorted(exps, reverse=True)
    total = sum(1 << e for e in exps)
    w = [max(1, int(255 * (1 << e) / total + 0.5)) for e in exps]  # every used segment counts
    # weights must sum to exactly 255; the driver corrects the LAST (smallest) weight,
    # which can zero a 1-sample segment -> we correct the largest one instead
    w[0] -= sum(w) - 255
    d = exps + [0] * (RSSI_MAX_SEGMENTS - len(exps))
    w += [0] * (RSSI_MAX_SEGMENTS - len(w))
    regs = {REG_RSSI_DUR_01: (d[1] << 4) | d[0], REG_RSSI_DUR_23: (d[3] << 4) | d[2]}
    regs.update({REG_RSSI_WEIGHT_0 + i: w[i] for i in range(RSSI_MAX_SEGMENTS)})
    return regs


def rssi_samples_from_registers(dur01, dur23, weights):
    """Effective averaging length (samples) currently programmed: sum of segments with
    non-zero weight (a zero-weight segment does not contribute to the average)."""
    d = (dur01 & 0xF, dur01 >> 4, dur23 & 0xF, dur23 >> 4)
    return sum(1 << d[i] for i in range(RSSI_MAX_SEGMENTS) if weights[i])


# ------------------------------------------------------------ device state

RETUNE_ORDER = ("sample_rate", "bandwidth", "frequency")  # order in which changes are applied


@dataclass
class DeviceState:
    """Last REQUESTED RX configuration applied to the device; None = unknown.

    Requested (not actual) values are stored on purpose: the device returns quantized
    actual values (e.g. 1491999998 Hz for 1492000000 Hz), so comparing a request with
    an actual value would report a change every time."""
    sample_rate: float | None = None
    bandwidth: float | None = None
    frequency: float | None = None

    def reset(self):
        self.sample_rate = self.bandwidth = self.frequency = None


DEVICE_STATE = DeviceState()


def get_retune_params(current, required):
    """Return [(name, value), ...] for parameters in `required` (dict) whose value
    differs from `current` (DeviceState), in RETUNE_ORDER. Empty list if nothing changes."""
    return [(name, required[name]) for name in RETUNE_ORDER
            if name in required and getattr(current, name) != required[name]]


# ------------------------------------------------------------------ scanner

class Scanner:
    """All bladeRF calls live here and must run in one thread."""

    def __init__(self, args, stop=None):
        self.args = args
        self.stop = stop if stop is not None else threading.Event()  # cooperative stop (Ctrl+C, GUI close)
        self.spans = args.freq_spans
        self.dev = None
        self.ch_id = bladerf.CHANNEL_RX(0)
        self.enabled = False
        self.gain_log = []       # (label, gain reported by the device), in the order they were read
        self.t_read = None       # time of one RSSI read over USB, s
        self.t_tune = None       # time of one retune, s
        self.period = None       # read period for the current span, s
        self.settle = None       # pause after retune for the current span, s
        self.hw_avg = None       # current on-board RSSI averaging length, samples
        self.rssi_cfg_sr = None  # sample rate the RSSI averaging was configured for
        self._reported = set()   # RSSI configs already written to the dump (no repeats per pass)
        self._quiet = False      # True during warmup: no dump, no info
        self.pass_no = 0
        # anomaly totals over the whole run, printed by close()
        self.anom = {"slow_retune": 0, "single_read_freq": 0, "slow_calls": 0, "late_reads": 0,
                     "wakeup_delay_freq": 0, "gc_runs": 0, "reconfigs": 0}
        self.anom_max = {"retune_ms": 0.0, "call_ms": 0.0, "reconf_ms": 0.0}
        self.t_origin = time.perf_counter()  # time base for all timestamps (t_ms)
        self._t_set = self._t_readback = 0.0  # split of the last frequency change
        # GC activity (Python), updated by a gc.callbacks hook; per-frequency deltas in metrics
        self._gc_n = 0
        self._gc_max = 0.0
        self._gc_t = 0.0
        self._metrics = open(args.metrics, "w", encoding="utf-8") if args.metrics else None
        if self._metrics:
            self._metrics.write("# event\tt_ms\tname\tdur_ms\n")
            self._metrics.write("\t".join(METRIC_COLS) + "\n")
        self.actual = {}         # actual values reported by the device for the current config
        # per-frequency stats of raw reads in dB (Welford), reset by _reset_read_stats()
        self._rd_n = 0
        self._rd_mean = 0.0
        self._rd_m2 = 0.0
        self._pre = ffi.new("int32_t *")   # preallocated output buffers for RSSI reads
        self._sym = ffi.new("int32_t *")
        self._u8 = ffi.new("uint8_t *")    # preallocated buffer for RFIC register reads

    def open(self):
        a = self.args
        gc.callbacks.append(self._gc_hook)
        t = time.perf_counter()
        self.dev = bladerf.BladeRF()
        ch = self.dev.Channel(self.ch_id)
        self.ch = ch
        self._event("open_device", t)

        if self.dev.get_rfic_rssi(self.ch_id) is None:
            raise ConfigError("RFIC RSSI is not supported by this board")

        t = time.perf_counter()
        ch.gain_mode = _bladerf.GainMode.Manual
        ch.gain = a.gain
        self._event("set_gain", t)
        self._read_gain("right after the request")

        DEVICE_STATE.reset()  # device state is unknown right after open
        self.rssi_cfg_sr = None
        first_freqs, first_wdt = self.spans[0]
        t = time.perf_counter()
        self.apply(sample_rate=first_wdt, bandwidth=first_wdt, frequency=first_freqs[0])
        self._event("initial_config", t)

        t = time.perf_counter()
        self.dev.enable_module(self.ch_id, True)
        self.enabled = True
        self._event("enable_rx", t)
        self._read_gain("start, tuned, RX enabled")

        t0 = time.perf_counter()
        for _ in range(CALIBRATION_READS):
            self._read_rssi()
        self.t_read = (time.perf_counter() - t0) / CALIBRATION_READS
        self._event("calib_reads", t0)

        self.t_tune = 0.0  # provisional, measured after the warmup
        t = time.perf_counter()
        self._warmup()
        self._event(f"warmup_{a.warmup}", t)

        # Retune time, measured over a few points of the first span
        # (the first probe point is the current frequency and is skipped by apply())
        probe = first_freqs[:: max(1, len(first_freqs) // CALIBRATION_TUNES)][:CALIBRATION_TUNES + 1]
        t0 = time.perf_counter()
        n_tunes = 0
        for f in probe:
            t = time.perf_counter()
            if self.apply(frequency=f):
                n_tunes += 1
                self._event(f"calib_tune {fmt_freq(f).strip()} MHz", t)
        self.t_tune = (time.perf_counter() - t0) / n_tunes if n_tunes else 0.0
        info(f"RSSI read: {self.t_read * 1e3:.3f} ms/call, retune: {self.t_tune * 1e3:.3f} ms")

    def _event(self, name, t_start, t_end=None):
        """Timestamped event with duration -> '# event' line in the metrics file."""
        t_end = time.perf_counter() if t_end is None else t_end
        if self._metrics:
            self._metrics.write(f"# event\t{(t_start - self.t_origin) * 1e3:.3f}\t{name}\t"
                                f"{(t_end - t_start) * 1e3:.3f}\n")

    def _gc_hook(self, phase, _info):
        if phase == "start":
            self._gc_t = time.perf_counter()
        else:
            dt = time.perf_counter() - self._gc_t
            self._gc_n += 1
            if dt > self._gc_max:
                self._gc_max = dt

    def _warmup(self):
        mode = self.args.warmup
        if mode == "none":
            return
        t0 = time.perf_counter()
        self._quiet = True
        try:
            if mode == "pass":
                self.sweep(self.stop)
            else:
                for freqs, wdt in self.spans:
                    self.apply(sample_rate=wdt, bandwidth=wdt)
                    probe = freqs[:: max(1, len(freqs) // WARMUP_TUNES)][:WARMUP_TUNES]
                    for f in probe:
                        t = time.perf_counter()
                        if self.apply(frequency=f):
                            self._event(f"warmup_tune {fmt_freq(f).strip()} MHz", t)
                    for _ in range(WARMUP_READS):
                        self._read_rssi()
                    if self.stop.is_set():
                        return
        finally:
            self._quiet = False
        # back to the first span; RSSI averaging gets configured (and reported) by the sweep
        first_freqs, first_wdt = self.spans[0]
        self.apply(sample_rate=first_wdt, bandwidth=first_wdt, frequency=first_freqs[0])
        self.rssi_cfg_sr = None
        self._reported.clear()
        info(f"warmup ({mode}): {(time.perf_counter() - t0) * 1e3:.1f} ms")

    def apply(self, **required):
        """Change only the parameters that differ from DEVICE_STATE.
        Returns the list of applied (name, value) changes."""
        changes = get_retune_params(DEVICE_STATE, required)
        for name, value in changes:
            if name == "sample_rate":
                self.actual[name] = self.dev.set_sample_rate(self.ch_id, value)
                # libbladeRF re-programs the RSSI registers to its defaults on a sample
                # rate change -> our averaging setup must be applied again
                self.rssi_cfg_sr = None
            elif name == "bandwidth":
                self.actual[name] = self.dev.set_bandwidth(self.ch_id, value)
            elif name == "frequency":
                t = time.perf_counter()
                self.ch.frequency = int(round(value))
                t_mid = time.perf_counter()
                self.actual[name] = self.ch.frequency
                self._t_set, self._t_readback = t_mid - t, time.perf_counter() - t_mid
            setattr(DEVICE_STATE, name, value)
        return changes

    # ---------------------------------------------------- RFIC registers

    def _get_reg(self, addr):
        ret = libbladeRF.bladerf_get_rfic_register(self.dev.dev[0], addr, self._u8)
        if ret:
            _check_error(ret)
        return self._u8[0]

    def _set_reg(self, addr, val):
        ret = libbladeRF.bladerf_set_rfic_register(self.dev.dev[0], addr, val)
        if ret:
            _check_error(ret)

    def _requested_rssi_samples(self, sr):
        """(samples, description, warning) requested by the user for sample rate sr;
        samples is None if not requested or not valid (-> board default)."""
        a = self.args
        if a.rssi_n is not None:
            n, what = a.rssi_n, f"--rssi-n {a.rssi_n}"
        elif a.rssi_dur is not None:
            n, what = int(round(a.rssi_dur * sr)), f"--rssi-dur {a.rssi_dur} s at {fmt_freq(sr).strip()} MHz"
        else:
            return None, "board default", None
        if not 1 <= n <= RSSI_MAX_SAMPLES:
            return (None, f"board default: {what} not valid",
                    f"WARNING: {what} -> {n} samples, outside 1..{RSSI_MAX_SAMPLES}; using board default")
        return n, what, None

    def _configure_rssi(self):
        """Program on-board RSSI averaging for the current sample rate (or keep the board
        default), read the result back from the registers, set read period and settle time,
        and write the parameters to the dump (stdout)."""
        sr = self.actual["sample_rate"]
        n_req, what, warning = self._requested_rssi_samples(sr)
        if n_req is not None:
            total, exps = nearest_rssi_samples(n_req)
            for addr, val in rssi_register_values(exps).items():
                self._set_reg(addr, val)
            cfg = self._get_reg(REG_RSSI_CONFIG)
            cfg = (cfg | RSSI_CFG_DEFAULT_MODE) if len(exps) == 1 else (cfg & ~RSSI_CFG_DEFAULT_MODE)
            self._set_reg(REG_RSSI_CONFIG, cfg)
            source = f"{what} -> {n_req} samples -> nearest representable"
        else:
            source = what

        # read back what is actually programmed
        weights = [self._get_reg(REG_RSSI_WEIGHT_0 + i) for i in range(RSSI_MAX_SEGMENTS)]
        self.hw_avg = rssi_samples_from_registers(self._get_reg(REG_RSSI_DUR_01),
                                                  self._get_reg(REG_RSSI_DUR_23), weights)
        wait = self._get_reg(REG_RSSI_WAIT) * 4
        hw_period = (self.hw_avg + wait) / sr
        self.period = max(hw_period, self.t_read)
        self.settle = SETTLE_PERIODS * self.period
        self.rssi_cfg_sr = sr

        # libbladeRF resets the registers on every sample rate change, so with spans of
        # different width this runs every pass; report each distinct config only once
        key = (sr, self.hw_avg, wait, source)
        if self._quiet or key in self._reported:
            return
        self._reported.add(key)
        if warning:
            info(warning)
        print(f"# wdt {fmt_freq(sr).strip()} MHz: RSSI hw averaging {self.hw_avg} samples "
              f"= {self.hw_avg / sr * 1e3:.3f} ms ({source}); wait {wait} samples; "
              f"read period {self.period * 1e3:.3f} ms", flush=True)

        dwell = self.args.dwell
        # dwell is counted from the end of the retune; reads at settle end + i*period,
        # check after each read -> first m with settle + (m-1)*period >= dwell
        rest = dwell - self.settle
        k = 1 if rest <= 0 else 1 + math.ceil(rest / self.period)
        info(f"wdt {fmt_freq(sr).strip()} MHz: dwell {dwell * 1e3:.3f} ms -> ~{k} reads per frequency"
             + (" (dwell shorter than settle)" if k == 1 else ""))

    def anomaly_summary(self):
        a, m = self.anom, self.anom_max
        return (f"ANOMALIES over {self.pass_no} pass(es): "
                f"slow retunes (>{RETUNE_ANOMALY_FACTOR:g}x median) {a['slow_retune']} (max {m['retune_ms']:.1f} ms); "
                f"frequencies with <2 reads {a['single_read_freq']}; "
                f"read calls >{SLOW_CALL_S * 1e3:g} ms {a['slow_calls']} (max {m['call_ms']:.1f} ms); "
                f"late reads {a['late_reads']}; "
                f"frequencies with wake-up delay > period {a['wakeup_delay_freq']}; "
                f"GC runs {a['gc_runs']}; "
                f"span reconfigs {a['reconfigs']} (max {m['reconf_ms']:.1f} ms)")

    def _read_gain(self, label):
        """READ-ONLY: ask the device for its gain right now and report it (stderr). Nothing is
        applied or corrected here. Returns the value or None if it cannot be read; every
        reading is kept in self.gain_log for the final summary."""
        want = self.args.gain
        try:
            got = self.ch.gain
        except Exception as e:
            info(f"WARNING: gain [{label}]: cannot be read: {e}")
            self.gain_log.append((label, None))
            return None
        self.gain_log.append((label, got))
        if got == want:
            info(f"gain [{label}]: device reports {got} dB (requested {want} dB)")
        else:
            info(f"WARNING: gain [{label}]: device reports {got} dB, requested {want} dB")
        return got

    def gain_summary(self):
        parts = "; ".join(f"[{label}] {'?' if v is None else v}" for label, v in self.gain_log)
        return f"# gain_dB: requested {self.args.gain}; as reported by the device: {parts}"

    def close(self):
        if self.pass_no:
            info(self.anomaly_summary())
        DEVICE_STATE.reset()
        if self._gc_hook in gc.callbacks:
            gc.callbacks.remove(self._gc_hook)
        if self._metrics:
            self._metrics.close()
            self._metrics = None
        if self.dev is None:
            return
        if self.enabled:
            try:
                self.dev.enable_module(self.ch_id, False)
            except Exception:
                pass
        self.dev.close()
        self.dev = None

    def header(self):
        # column header is printed by sweep() right above the rows of each span,
        # after the span's configuration lines
        return f"# gain_dB requested {self.args.gain}  dwell_ms={self.args.dwell * 1e3:.3f}"

    # ----------------------------------------------------------- reading

    def _read_rssi(self):
        """Symbol RSSI, dB. No allocations: preallocated cffi buffers, direct libbladeRF call."""
        ret = libbladeRF.bladerf_get_rfic_rssi(self.dev.dev[0], self.ch_id, self._pre, self._sym)
        if ret:
            _check_error(ret)
        return self._sym[0]

    def _reset_read_stats(self):
        self._rd_n = 0
        self._rd_mean = 0.0
        self._rd_m2 = 0.0

    def _read_std_db(self):
        """Sample standard deviation (ddof=1) of all reads on the frequency, dB; NaN if < 2 reads."""
        return math.sqrt(self._rd_m2 / (self._rd_n - 1)) if self._rd_n > 1 else math.nan

    def sweep(self, stop=None):
        """One pass over all spans. Returns a list (one item per span) of dicts of
        numpy arrays: freq (Hz), mean (dB, power mean), std (dB), n (reads),
        plus performance metrics: retune, first_read, mean_gap, max_gap, max_call, used (s),
        late (reads); or None if stopped.

        On each frequency: retune, settle, then read the hardware RSSI once per hardware
        averaging period until dwell has elapsed since the END of the retune; the check is
        done after each read, so there is always at least one. If a read starts late, the
        schedule restarts from that read (no back-to-back catch-up reads of the same RSSI
        value). Mean is in linear power, std over the reads in dB. No per-read allocations."""
        result = []
        dwell = self.args.dwell
        t_pass = time.perf_counter()
        total_reads = 0
        if not self._quiet:
            self.pass_no += 1
        prev_f = DEVICE_STATE.frequency
        for freqs, wdt in self.spans:
            t_rc = time.perf_counter()
            changed = self.apply(sample_rate=wdt, bandwidth=wdt)
            if changed and not self._quiet:
                self._event(f"span_config wdt {fmt_freq(wdt).strip()} MHz", t_rc)
            sr, bw = self.actual["sample_rate"], self.actual["bandwidth"]
            if changed and bw != sr and not self._quiet:
                info(f"WARNING: wdt {fmt_freq(sr).strip()} MHz: actual bandwidth "
                     f"{fmt_freq(bw).strip()} MHz differs")
            if self.rssi_cfg_sr != sr:
                t = time.perf_counter()
                self._configure_rssi()
                if not self._quiet:
                    self._event(f"rssi_config wdt {fmt_freq(sr).strip()} MHz", t)
                changed = changed or [("rssi_config", None)]
            # reconfiguration time of this span (0 if nothing had to change)
            reconf = time.perf_counter() - t_rc if changed else 0.0
            if not self._quiet:
                print(fmt_row(name for name, _ in TABLE_COLS), flush=True)
            period = self.period
            late_gap = LATE_GAP_FACTOR * period
            k = len(freqs)
            out = {"freq": np.empty(k), "mean": np.empty(k), "std": np.empty(k),
                   "n": np.empty(k, dtype=np.int64)}
            out["reconf"] = np.zeros(k)
            out["reconf"][0] = reconf          # attributed to the first frequency of the span
            for key in ("t", "jump", "retune", "set", "readback", "first_read", "mean_gap",
                        "max_gap", "max_call", "sum_call", "sched_late", "gc_max", "used"):
                out[key] = np.empty(k)
            for key in ("late", "gt1", "gt5", "gc_n"):
                out[key] = np.empty(k, dtype=np.int64)
            for i in range(k):
                if stop is not None and stop.is_set():
                    return None
                gc_n0 = self._gc_n
                self._gc_max = 0.0
                self._t_set = self._t_readback = 0.0
                t0 = time.perf_counter()
                self.apply(frequency=freqs[i])
                t_tuned = time.perf_counter()
                fr = self.actual["frequency"]
                t_next = t_tuned + self.settle
                sum_lin = 0.0
                self._reset_read_stats()
                t_first = t_prev = t_end = 0.0
                max_gap = max_call = sum_call = sched_late = 0.0
                late = gt1 = gt5 = 0
                while True:
                    t_ready = t_next if t_next > t_end else t_end  # earliest possible start
                    wait_until(t_next)
                    t_start = time.perf_counter()
                    v = self._read_rssi()
                    t_end = time.perf_counter()
                    call = t_end - t_start
                    sum_call += call
                    if call > 1e-3:
                        gt1 += 1
                        if call > SLOW_CALL_S:
                            gt5 += 1
                    # wake-up delay: start later than both the schedule and the end of the
                    # previous call -> time lost outside libbladeRF (Python/OS), not a slow call
                    if t_start - t_ready > sched_late:
                        sched_late = t_start - t_ready
                    if self._rd_n == 0:
                        t_first = t_start
                    else:
                        gap = t_start - t_prev
                        if gap > max_gap:
                            max_gap = gap
                        if gap > late_gap:
                            late += 1
                    if call > max_call:
                        max_call = call
                    t_prev = t_start
                    # next read one period after this one's start: on time -> same as
                    # t_next + period; late -> schedule restarts, no catch-up burst
                    t_next = max(t_next, t_start) + period
                    sum_lin += db_to_lin(v)
                    self._rd_n += 1
                    d = v - self._rd_mean
                    self._rd_mean += d / self._rd_n
                    self._rd_m2 += d * (v - self._rd_mean)
                    if t_end - t_tuned >= dwell:   # dwell counted after the retune
                        break
                    if stop is not None and stop.is_set():
                        return None
                n = self._rd_n
                mean_db, std_db = lin_to_db(sum_lin / n), self._read_std_db()
                out["freq"][i], out["mean"][i], out["std"][i], out["n"][i] = fr, mean_db, std_db, n
                out["t"][i] = t0 - self.t_origin
                out["jump"][i] = abs(freqs[i] - prev_f) if prev_f is not None else math.nan
                prev_f = freqs[i]
                out["retune"][i] = t_tuned - t0
                out["set"][i], out["readback"][i] = self._t_set, self._t_readback
                out["sum_call"][i], out["sched_late"][i] = sum_call, sched_late
                out["gt1"][i], out["gt5"][i] = gt1, gt5
                out["gc_n"][i], out["gc_max"][i] = self._gc_n - gc_n0, self._gc_max
                out["first_read"][i] = t_first - t0
                out["mean_gap"][i] = (t_prev - t_first) / (n - 1) if n > 1 else math.nan
                out["max_gap"][i] = max_gap if n > 1 else math.nan
                out["max_call"][i] = max_call
                out["used"][i] = t_end - t0
                out["late"][i] = late
                total_reads += n
                if not self._quiet:
                    std_txt = "n/a" if math.isnan(std_db) else f"{std_db:.2f}"
                    print(fmt_row((fmt_freq(fr), fmt_freq(sr), f"{mean_db:.2f}", n, std_txt)), flush=True)
            out["sr"] = sr
            out["period"] = period
            result.append(out)
        if not self._quiet:
            n_freq = sum(len(r["freq"]) for r in result)
            info(f"pass {self.pass_no}: {n_freq} frequencies, "
                 f"{(time.perf_counter() - t_pass) * 1e3 / n_freq:.2f} ms/frequency, "
                 f"{total_reads / n_freq:.1f} reads/frequency")
            self._report_perf(result)
        return result

    def _report_perf(self, result):
        """Pass summary to stderr, per-frequency metrics to the --metrics file (outside hot path)."""
        self._read_gain(f"end of pass {self.pass_no}")
        cat = lambda key: np.concatenate([r[key] for r in result])  # noqa: E731
        retune, setf, rb = cat("retune") * 1e3, cat("set") * 1e3, cat("readback") * 1e3
        call, gap, sched = cat("max_call") * 1e3, cat("max_gap") * 1e3, cat("sched_late") * 1e3
        late, gt5, gc_n, gc_max = cat("late"), cat("gt5"), cat("gc_n"), cat("gc_max") * 1e3
        reads, reconf = cat("n"), cat("reconf") * 1e3
        period_ms = np.concatenate([np.full(len(r["freq"]), r["period"]) for r in result]) * 1e3
        gap_ok = gap[~np.isnan(gap)]
        p = lambda a: f"{np.median(a):.2f}/{np.percentile(a, 90):.2f}/{a.max():.2f}"  # noqa: E731
        rc = reconf[reconf > 0]
        info(f"perf pass {self.pass_no} (median/p90/max, ms): retune {p(retune)} "
             f"[set {p(setf)}, readback {p(rb)}]; read call max {p(call)}; "
             f"read gap max {gap_ok.max() if gap_ok.size else math.nan:.2f}; "
             f"wake-up delay max {sched.max():.3f}; "
             + (f"reconfig {rc.size}x, max {rc.max():.2f}" if rc.size else "reconfig none"))
        info(f"perf pass {self.pass_no}: late reads {int(late.sum())} on {int(np.count_nonzero(late))} "
             f"frequencies; calls >5 ms {int(gt5.sum())}; GC runs {int(gc_n.sum())}, "
             f"longest GC pause {gc_max.max():.3f} ms")

        # anomaly totals for the end-of-run summary
        a, m = self.anom, self.anom_max
        a["slow_retune"] += int(np.count_nonzero(retune > RETUNE_ANOMALY_FACTOR * np.median(retune)))
        a["single_read_freq"] += int(np.count_nonzero(reads < 2))
        a["slow_calls"] += int(gt5.sum())
        a["late_reads"] += int(late.sum())
        a["wakeup_delay_freq"] += int(np.count_nonzero(sched > period_ms))
        a["gc_runs"] += int(gc_n.sum())
        a["reconfigs"] += int(rc.size)
        m["retune_ms"] = max(m["retune_ms"], float(retune.max()))
        m["call_ms"] = max(m["call_ms"], float(call.max()))
        m["reconf_ms"] = max(m["reconf_ms"], float(rc.max()) if rc.size else 0.0)
        if not self._metrics:
            return
        f = self._metrics
        ms = lambda v: f"{v * 1e3:.3f}"  # noqa: E731
        for r in result:
            wdt_txt = fmt_freq(r["sr"]).strip()
            for i in range(len(r["freq"])):
                jump = r["jump"][i]
                f.write("\t".join((
                    str(self.pass_no), ms(r["t"][i]), fmt_freq(r["freq"][i]).strip(), wdt_txt,
                    ms(r["reconf"][i]), "nan" if math.isnan(jump) else f"{jump * HZ_TO_MHZ:.2f}",
                    ms(r["retune"][i]), ms(r["set"][i]), ms(r["readback"][i]), ms(r["first_read"][i]),
                    str(r["n"][i]), ms(r["mean_gap"][i]), ms(r["max_gap"][i]), str(r["late"][i]),
                    ms(r["max_call"][i]), ms(r["sum_call"][i]), str(r["gt1"][i]), str(r["gt5"][i]),
                    ms(r["sched_late"][i]), str(r["gc_n"][i]), ms(r["gc_max"][i]), ms(r["used"][i]),
                )) + "\n")
        f.flush()


# ---------------------------------------------------------------------- CLI

def run_cli(args, stop):
    sc = Scanner(args, stop)
    try:
        sc.open()
        if stop.is_set():
            info("interrupted")
            return 130
        print(sc.header())
        if sc.sweep(stop) is None:
            info("interrupted")
            return 130
        return 0
    except ConfigError as e:
        info(f"ERROR: {e}")
        return 2
    finally:
        sc.close()


# ---------------------------------------------------------------------- GUI

def run_gui(args, stop):
    import tkinter as tk
    from tkinter import messagebox

    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure

    q = queue.Queue()
    sc = Scanner(args, stop)

    def worker():
        try:
            sc.open()
            q.put(("ready", None))
            print(sc.header(), flush=True)
            k = 0
            while not stop.is_set():
                k += 1
                print(f"# pass {k}", flush=True)
                res = sc.sweep(stop)
                if res is None:
                    break
                q.put(("pass", (k, res)))
        except Exception as e:  # report to GUI thread
            q.put(("error", e))
        finally:
            sc.close()
            q.put(("closed", None))

    root = tk.Tk()
    root.title("bladeRF RSSI sweep")
    n_spans = len(args.freq_spans)
    fig = Figure(figsize=(10, max(3.0, 2.6 * n_spans)), dpi=100, layout="constrained")
    axes = [fig.add_subplot(n_spans, 1, i + 1) for i in range(n_spans)]
    for ax in axes:
        ax.set_ylabel("RSSI, dB")
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("Frequency, MHz")
    canvas = FigureCanvasTkAgg(fig, master=root)
    NavigationToolbar2Tk(canvas, root)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    status = tk.StringVar(value="opening device...")
    tk.Label(root, textvariable=status, anchor="w").pack(fill=tk.X)

    # per axis: drawn artists and the limits we set automatically last time
    state = {"artists": [[] for _ in axes], "auto_lims": [None] * n_spans, "closing": False}

    def draw(k, res):
        for i, (ax, r) in enumerate(zip(axes, res)):
            auto = state["auto_lims"][i]
            user_zoomed = auto is not None and (ax.get_xlim(), ax.get_ylim()) != auto
            x = r["freq"] * HZ_TO_MHZ  # vectorized, same scale as fmt_freq
            sd = np.nan_to_num(r["std"])  # no std (single read) -> zero-width band
            lo, hi = r["mean"] - sd, r["mean"] + sd
            arts = state["artists"][i]
            if arts:  # update existing artists in place, no re-creation
                band, ln = arts
                band.set_data(x, lo, hi)
                ln.set_data(x, r["mean"])
            else:
                band = ax.fill_between(x, lo, hi, color="C0", alpha=0.2, linewidth=0)
                (ln,) = ax.plot(x, r["mean"], color="C0", marker="." if len(x) < 200 else None)
                state["artists"][i] = [band, ln]
            if not user_zoomed:
                dx = max((x[-1] - x[0]) * 0.02, 0.5)
                ax.set_xlim(x[0] - dx, x[-1] + dx)
                ax.set_ylim(lo.min() - 2, hi.max() + 2)
                state["auto_lims"][i] = (ax.get_xlim(), ax.get_ylim())
        fig.suptitle(f"pass {k}   gain (device, last reading)={sc.gain_log[-1][1] if sc.gain_log else '?'} dB   dwell={args.dwell * 1e3:.1f} ms   "
                     "(line: power mean, band: ±std of reads)")
        canvas.draw_idle()
        status.set(f"pass {k} done, {time.strftime('%H:%M:%S')}")

    def poll():
        # Ctrl+C in the console sets `stop` (signal handler runs in this, the main, thread
        # between Tk callbacks) -> close the window the same way as the close button
        if stop.is_set() and not state["closing"]:
            on_close()
        try:
            while True:
                kind, payload = q.get_nowait()
                if kind == "ready":
                    status.set("sweeping...")
                elif kind == "pass":
                    draw(*payload)
                elif kind == "error":
                    if not state["closing"]:
                        msg = str(payload)
                        info(f"ERROR: {msg}")
                        status.set(f"ERROR: {msg}")
                        messagebox.showerror("RSSI sweep", msg)
                elif kind == "closed":
                    if state["closing"]:
                        root.destroy()
                        return
        except queue.Empty:
            pass
        root.after(100, poll)

    def on_close():
        state["closing"] = True
        status.set("stopping...")
        stop.set()
        if not th.is_alive():
            root.destroy()

    th = threading.Thread(target=worker, daemon=True)
    th.start()
    root.protocol("WM_DELETE_WINDOW", on_close)
    root.after(100, poll)
    root.mainloop()
    return 0


def install_sigint(stop):
    """First Ctrl+C: cooperative stop (current read / libbladeRF call finishes, device is closed
    cleanly). Second Ctrl+C: immediate exit without cleanup."""
    def handler(_sig, _frame):
        if stop.is_set():
            os._exit(130)
        stop.set()
        info("Ctrl+C: stopping after the current operation (press again to force exit)")
    signal.signal(signal.SIGINT, handler)


def main():
    args = parse_args()
    stop = threading.Event()
    install_sigint(stop)
    return run_gui(args, stop) if args.gui else run_cli(args, stop)


if __name__ == "__main__":
    try:
      main()
    except Exception as e:
        info(f"FATAL: {e}")
        traceback = getattr(e, "__traceback__", None)
        if traceback:
            import traceback as tb
            tb.print_tb(traceback)
        sys.exit(1)