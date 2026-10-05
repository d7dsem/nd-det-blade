#!/usr/bin/env python3
"""bladeRF 2.0: AD9361 hardware (register) RSSI vs RSSI estimated from the spectrum of an IQ capture.

Flow
  1. Register RSSI sweep over --freq-spans (same mechanics as rssi_scan.py): on every frequency the
     hardware RSSI is read until --dwell has elapsed since the END of the retune; mean (power) and
     std in dB are kept.
  2. One IQ capture of the band flo..fhi that contains all spans (lowest flo and highest fhi as typed in
     --freq-spans). Its width must not exceed the maximum capture (--max-bw, 56 MHz for bladeRF 2.0),
     otherwise the script refuses to run before anything is done. The capture sample rate is the
     width of this range (not lower than the minimum sample rate). There is NO margin for the
     filter roll-off and no marking of edge points: a point whose RSSI band [f - wdt/2, f + wdt/2]
     reaches beyond the captured band is integrated over the part that lies inside it. If edge
     distortions matter, set a narrower range.
  3. Welch spectrum of the capture. FFT size: --df XOR --fft-n. Number of averaged windows: all windows
     that cover --dwell (an integer count, may exceed dwell by up to one window).
     --avr-n / --avr-dur take effect only when given: not shorter than dwell -> the capture is
     lengthened to them; shorter than dwell -> conflict: the script stops, or, if started from a
     terminal, offers to resolve it in a dialog.
  4. RSSI estimate for every grid frequency = spectrum power inside a band around it, compared with the
     register value (raw difference, and difference after removing a constant offset).

CLI
  --freq-spans "[{flo~fhi,step,wdt},...]"  Hz, suffixes k/M/G; wdt = sample rate (+ HW bandwidth);
                                           the RSSI band of a point f is f +- wdt/2
  --dwell <sec>                            time on each frequency, counted after the retune
  -g, --gain <dB>                          manual RX gain (default 42)
  [--rssi-n <int> | --rssi-dur <sec>]      optional on-board RSSI averaging (samples / duration)
  --df <Hz> | --fft-n <int>                FFT resolution (guide) / size; nearest 2^a*3^b*5^c*7^d is used
  [--avr-n <int> | --avr-dur <sec>]        optional Welch averaging (windows / duration)
  --window {rect,hann,hamming,blackman}    default rect (no multiplication at all)
  --overlap <0..1)                         window overlap, default 0
  --rssi-bw <Hz>                           band integrated per frequency (default: the span's wdt)
  --max-bw <Hz>                            maximum capture, default 56e6; compared with the width of
                                           the range flo..fhi (see step 2)
  --max-capture-mb <MB>                    memory limit for the capture, default 1024
  --warmup {short,pass,none}               before measuring, default short

Data -> stdout, info/warnings -> stderr. Ctrl+C stops cleanly; a second Ctrl+C forces exit.
"""
import argparse
import bisect
from dataclasses import dataclass
from itertools import combinations_with_replacement
import math
import os
import re
import signal
import sys
import threading
import time

import numpy as np

import bladerf
from bladerf import _bladerf
# Low-level handles for the hot path (no per-call allocations): direct libbladeRF calls
from bladerf._bladerf import ffi, libbladeRF, _check_error

# ---------------------------------------------------------------- constants

# bladeRF 2.0 limits (libbladeRF 2025.10, fpga_common/include/bladerf2_common.h)
SR_MIN = 520834.0            # base sample rate range, Hz
SR_MAX = 61440000.0
BW_MAX_DEFAULT = 56000000.0  # maximum analog bandwidth = maximum capture
FS_FULL_SCALE = 2048.0       # SC16_Q11: full scale per component, counts

CALIBRATION_READS = 50
CALIBRATION_TUNES = 5
WARMUP_TUNES = 10            # short warmup: retunes per span
WARMUP_READS = 20            # short warmup: reads after the retunes, per span
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

# Streaming (assumptions: sizes taken from typical libbladeRF usage, not tuned on hardware)
SYNC_BUFFERS = 32
SYNC_BUF_SAMPLES = 65536     # multiple of 8192
SYNC_TRANSFERS = 16          # must be < SYNC_BUFFERS
SYNC_TIMEOUT_MS = 10000
SYNC_CHUNK = 1 << 20         # samples per sync_rx call
DISCARD_MIN = 131072         # samples dropped after the stream starts (assumption)
DISCARD_S = 0.01             # ... or this much time, whichever is more (assumption)
WELCH_CHUNK_BYTES = 8 << 20  # work buffer size for one batch of windows
MIN_BINS_PER_BAND = 4        # fewer FFT bins inside the RSSI band -> coarse-estimate warning (assumption)

_SUFFIX = {"": 1.0, "k": 1e3, "K": 1e3, "M": 1e6, "G": 1e9}
HZ_TO_MHZ = 1e-6


def info(msg):
    print(msg, file=sys.stderr, flush=True)


def fmt_freq(hz, digits=2):
    """Single formatter for every frequency-like value (Hz) in text output -> MHz."""
    return f"{hz * HZ_TO_MHZ:{digits + 4}.{digits}f}"


# Output tables: one spec per table for header and rows -> columns always aligned
REG_COLS = (("freq_MHz", 9), ("wdt_MHz", 8), ("rssi_dB", 8), ("reads", 6), ("std_dB", 7))
CMP_COLS = (("freq_MHz", 9), ("wdt_MHz", 8), ("reg_dB", 8), ("fft_dBFS", 10), ("diff_dB", 8))


def fmt_row(cols, values):
    return "".join(f"{v:>{w}}" for v, (_, w) in zip(values, cols))


def fmt_header(cols):
    return fmt_row(cols, [name for name, _ in cols])


class ConfigError(Exception):
    pass


# ---------------------------------------------------------------- CLI parsing

def parse_hz(text):
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)\s*([kKMG]?)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError(f"bad frequency value: {text!r}")
    return float(m.group(1)) * _SUFFIX[m.group(2)]


class Spans(list):
    """List of (frequencies, wdt); .ranges holds (flo, fhi) of every span exactly as typed."""
    ranges: list


def parse_spans(text):
    body = text.strip()
    if body.startswith("[") and body.endswith("]"):
        body = body[1:-1]
    items = re.findall(r"\{([^{}]*)\}", body)
    if not items or re.sub(r"\{[^{}]*\}", "", body).strip(" ,"):
        raise argparse.ArgumentTypeError(f"bad --freq-spans: {text!r}")

    spans = Spans()
    spans.ranges = []
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
        spans.ranges.append((lo, hi))
    return spans


def positive_float(text):
    v = float(text)
    if v <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return v


def positive_int(text):
    v = int(text)
    if v <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return v


def fraction(text):
    v = float(text)
    if not 0.0 <= v < 1.0:
        raise argparse.ArgumentTypeError("must be in [0, 1)")
    return v


def parse_args():
    p = argparse.ArgumentParser(description="bladeRF 2.0: register RSSI vs RSSI from the spectrum of an IQ capture")
    p.add_argument("--freq-spans", required=True, type=parse_spans,
                   help='e.g. "[{2.40G~2.43G,1M,1M}]"; the RSSI band of a point f is f +- wdt/2')
    p.add_argument("--dwell", required=True, type=positive_float,
                   help="time per frequency after the retune, s")
    p.add_argument("-g", "--gain", type=int, default=42, help="manual RX gain, dB (default 42)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--rssi-n", type=int, help="on-board RSSI averaging length, RX samples (optional)")
    g.add_argument("--rssi-dur", type=float, help="on-board RSSI averaging duration, s (optional)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--df", type=parse_hz, help="FFT resolution guide, Hz (nearest fast FFT size is used)")
    g.add_argument("--fft-n", type=positive_int, help="FFT size (2^a*3^b*5^c*7^d)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--avr-n", type=positive_int, help="Welch averaging, number of windows (optional)")
    g.add_argument("--avr-dur", type=positive_float, help="Welch averaging duration guide, s (optional)")
    p.add_argument("--window", choices=("rect", "hann", "hamming", "blackman"), default="rect",
                   help="FFT window; rect = no multiplication (default)")
    p.add_argument("--overlap", type=fraction, default=0.0, help="window overlap 0..1 (default 0)")
    p.add_argument("--rssi-bw", type=parse_hz, help="band integrated per frequency, Hz (default: span wdt)")
    p.add_argument("--max-bw", type=parse_hz, default=BW_MAX_DEFAULT,
                   help="maximum capture, Hz (default 56M); the width of the range flo..fhi "
                        "must not exceed it, otherwise the script refuses to run")
    p.add_argument("--max-capture-mb", type=positive_float, default=1024.0,
                   help="memory limit for the IQ capture, MB (default 1024)")
    p.add_argument("--warmup", choices=("short", "pass", "none"), default="short",
                   help="warmup of the register sweep (default short)")
    return p.parse_args()


# ------------------------------------------------------------------ small helpers

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


# =============================================================================
# CORE (no device access): planning, Welch spectrum, band power, statistics.
# Pure functions over preallocated arrays; written so that it can be ported to C++.
# =============================================================================

# ------------------------------------------------------------------ FFT sizes

def _build_fft_sizes(limit=1 << 22):
    sizes = []
    a = 1
    while a <= limit:
        b = a
        while b <= limit:
            c = b
            while c <= limit:
                d = c
                while d <= limit:
                    sizes.append(d)
                    d *= 7
                c *= 5
            b *= 3
        a *= 2
    return np.array(sorted(set(sizes)), dtype=np.int64)


FFT_SIZES = _build_fft_sizes()


def is_fft_size(n):
    i = int(np.searchsorted(FFT_SIZES, n))
    return i < len(FFT_SIZES) and int(FFT_SIZES[i]) == n


def nearest_fft_size(target):
    """Allowed FFT size (2^a*3^b*5^c*7^d) closest to target (ties -> smaller), at least 16."""
    target = max(16, int(round(target)))
    i = int(np.searchsorted(FFT_SIZES, target))
    cand = FFT_SIZES[max(i - 1, 0):i + 1]
    return int(min(cand, key=lambda s: (abs(int(s) - target), int(s))))


# --------------------------------------------------------------------- planning

def range_of_spans(spans):
    """Lowest flo and highest fhi over all spans, Hz, exactly as typed by the user."""
    return min(r[0] for r in spans.ranges), max(r[1] for r in spans.ranges)


@dataclass
class CapturePlan:
    lo: float          # lowest frequency of the range, Hz
    hi: float          # highest, Hz
    center: float
    fs: float          # requested sample rate, Hz
    bw: float          # requested analog bandwidth, Hz


def plan_capture(lo, hi, max_bw):
    """Capture of the range lo..hi (Hz). Refuses ranges wider than the maximum capture."""
    width = hi - lo
    if width > max_bw:
        raise ConfigError(f"range {fmt_freq(lo).strip()}..{fmt_freq(hi).strip()} MHz is "
                          f"{fmt_freq(width).strip()} MHz wide, the maximum capture is "
                          f"{fmt_freq(max_bw).strip()} MHz (--max-bw); nothing was done")
    fs = min(max(width, SR_MIN), SR_MAX)
    fs = float(2 * math.ceil(fs / 2))
    return CapturePlan(lo, hi, float(round((lo + hi) / 2)), fs, min(fs, max_bw))


@dataclass
class WelchPlan:
    fs: float
    nfft: int
    hop: int
    df: float
    n_dwell: int               # windows needed to cover dwell
    n_user: int | None         # windows asked by --avr-n / --avr-dur (None if not given)
    notes: list


def welch_samples(nfft, hop, n_win):
    return (n_win - 1) * hop + nfft


def plan_welch(fs, dwell, df, fft_n, avr_n, avr_dur, overlap):
    """FFT size, hop and window counts for sample rate fs."""
    notes = []
    if fft_n is not None:
        nfft = fft_n
        if not is_fft_size(nfft):
            nfft = nearest_fft_size(fft_n)
            notes.append(f"--fft-n {fft_n} is not 2^a*3^b*5^c*7^d: using {nfft}")
    else:
        nfft = nearest_fft_size(fs / df)
        notes.append(f"--df {df:g} Hz -> FFT size {nfft}, actual df {fs / nfft:.3f} Hz")
    hop = max(1, int(round(nfft * (1.0 - overlap))))
    need = dwell * fs
    n_dwell = 1 if need <= nfft else 1 + math.ceil((need - nfft) / hop)
    n_user = None
    if avr_n is not None:
        n_user = avr_n
    elif avr_dur is not None:
        n_user = max(1, int(round((avr_dur * fs - nfft) / hop)) + 1)
        notes.append(f"--avr-dur {avr_dur:g} s -> {n_user} windows "
                     f"({welch_samples(nfft, hop, n_user) / fs * 1e3:.3f} ms)")
    return WelchPlan(fs, nfft, hop, fs / nfft, n_dwell, n_user, notes)


def avr_conflict_text(wp, dwell):
    return (f"requested averaging is shorter than dwell:\n"
            f"  --avr-*  -> {wp.n_user} windows ({welch_samples(wp.nfft, wp.hop, wp.n_user) / wp.fs * 1e3:.3f} ms)\n"
            f"  --dwell  -> {wp.n_dwell} windows ({welch_samples(wp.nfft, wp.hop, wp.n_dwell) / wp.fs * 1e3:.3f} ms)"
            f" are needed to cover {dwell * 1e3:.3f} ms")


def ask_avr_conflict(wp, dwell):
    """Dialog for the --avr-* < dwell conflict. Returns 'dwell' or 'user'; raises ConfigError otherwise."""
    text = avr_conflict_text(wp, dwell)
    if not sys.stdin.isatty():
        raise ConfigError("CONFLICT: " + text + "\n  resolve: omit --avr-*, set it >= dwell, or run from a terminal "
                          "to choose interactively; nothing was done")
    info("CONFLICT: " + text)
    info("  1) average everything captured within dwell, ignore --avr-*\n"
         "  2) keep the requested averaging, the capture will be shorter than dwell\n"
         "  3) stop")
    while True:
        sys.stderr.write("choose 1/2/3: ")
        sys.stderr.flush()
        line = sys.stdin.readline()
        if line == "":
            raise ConfigError("no answer (input closed); nothing was done")
        ans = line.strip()
        if ans == "1":
            return "dwell"
        if ans == "2":
            return "user"
        if ans == "3":
            raise ConfigError("stopped by user; nothing was done")


def choose_windows(wp, dwell, policy=None):
    """Number of windows to average and the policy used for a --avr-* conflict.
    No --avr-*: all windows covering dwell. --avr-* >= dwell: the capture is lengthened to it.
    --avr-* < dwell: conflict -> dialog / stop."""
    if wp.n_user is None:
        return wp.n_dwell, None
    if wp.n_user >= wp.n_dwell:
        return wp.n_user, "user"
    if policy is None:
        policy = ask_avr_conflict(wp, dwell)
    return (wp.n_dwell if policy == "dwell" else wp.n_user), policy


def check_capture_size(n_samples, max_mb):
    mb = n_samples * 4 / (1 << 20)
    if mb > max_mb:
        raise ConfigError(f"capture of {n_samples} samples needs {mb:.0f} MB, above --max-capture-mb "
                          f"{max_mb:g}; reduce dwell/--avr-* or raise the limit; nothing was done")


# ----------------------------------------------------------------- Welch spectrum

def make_window(name, n):
    """Periodic window (float32); None for rect: no multiplication is performed at all."""
    if name == "rect":
        return None
    fn = {"hann": np.hanning, "hamming": np.hamming, "blackman": np.blackman}[name]
    return fn(n + 1)[:-1].astype(np.float32)


class Welch:
    """Welch periodogram averaging with all work buffers allocated once."""

    def __init__(self, nfft, window, n_win):
        self.nfft = nfft
        self.win = make_window(window, nfft)
        self.sum_w2 = float(nfft) if self.win is None else float(np.sum(self.win.astype(np.float64) ** 2))
        self.chunk = max(1, min(n_win, WELCH_CHUNK_BYTES // (nfft * 8)))
        self.work = np.empty((self.chunk, nfft), np.complex64)
        self.spec = np.empty((self.chunk, nfft), np.complex64)
        self.pw = np.empty((self.chunk, nfft), np.float32)
        self.tmp = np.empty((self.chunk, nfft), np.float32)
        self.colsum = np.empty(nfft, np.float64)
        self.acc = np.zeros(nfft, np.float64)
        self.n_acc = 0

    def accumulate(self, iq, hop, n_win):
        """iq: int16 array (n_samples, 2) = I, Q. Adds the power spectra of n_win windows."""
        nfft = self.nfft
        view = np.lib.stride_tricks.sliding_window_view(iq, nfft, axis=0)[::hop][:n_win]  # (win, 2, nfft)
        if view.shape[0] < n_win:
            raise ConfigError(f"capture too short: {view.shape[0]} windows available, {n_win} needed")
        for i in range(0, n_win, self.chunk):
            m = min(self.chunk, n_win - i)
            blk = view[i:i + m]
            w, s, pw, tmp = self.work[:m], self.spec[:m], self.pw[:m], self.tmp[:m]
            np.copyto(w.real, blk[:, 0, :], casting="unsafe")   # int16 -> float32, no scaling here
            np.copyto(w.imag, blk[:, 1, :], casting="unsafe")
            if self.win is not None:
                np.multiply(w, self.win, out=w)
            try:
                np.fft.fft(w, axis=1, out=s)
            except TypeError:                                   # numpy < 2.0 has no out=
                s[...] = np.fft.fft(w, axis=1)
            np.square(s.real, out=pw)
            np.square(s.imag, out=tmp)
            pw += tmp
            np.add.reduce(pw, axis=0, dtype=np.float64, out=self.colsum)
            self.acc += self.colsum
        self.n_acc += n_win

    def power_spectrum(self):
        """Power per bin, fftshifted (bin nfft//2 = center), relative to full scale:
        a complex tone of amplitude FS_FULL_SCALE counts gives 1.0 (0 dBFS) in its bin."""
        scale = 1.0 / (self.n_acc * self.nfft * self.sum_w2 * FS_FULL_SCALE ** 2)
        return np.fft.fftshift(self.acc * scale)


def band_power_db(p_shift, fs, nfft, center, freqs, bws):
    """Power (dB re full scale) in [f - bw/2, f + bw/2] for every f, from the shifted per-bin
    power. Band edges are fractional-bin accurate (linear interpolation of the cumulative sum)."""
    df = fs / nfft
    cs = np.empty(nfft + 1)
    cs[0] = 0.0
    np.cumsum(p_shift, out=cs[1:])
    k0 = nfft // 2
    u_lo = (freqs - bws / 2 - center) / df + k0 + 0.5
    u_hi = (freqs + bws / 2 - center) / df + k0 + 0.5
    x = np.arange(nfft + 1, dtype=np.float64)
    power = np.interp(u_hi, x, cs) - np.interp(u_lo, x, cs)
    with np.errstate(divide="ignore"):
        return 10.0 * np.log10(power)


@dataclass
class CompareStats:
    n: int
    mean: float
    std: float
    offset: float     # median of (reg - fft)
    rms_res: float    # RMS of (reg - fft - offset)
    max_res: float
    corr: float


def compare_stats(reg, spec):
    d = reg - spec
    ok = np.isfinite(d)
    d, r, s = d[ok], reg[ok], spec[ok]
    if d.size == 0:
        return CompareStats(0, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan)
    off = float(np.median(d))
    res = d - off
    corr = float(np.corrcoef(r, s)[0, 1]) if d.size > 1 and np.std(r) > 0 and np.std(s) > 0 else math.nan
    return CompareStats(int(d.size), float(d.mean()), float(d.std()), off,
                        float(np.sqrt(np.mean(res * res))), float(np.max(np.abs(res))), corr)


def count_clipped(iq, limit=2047, step=1 << 20):
    """Samples at or beyond +-limit counts (ADC/format clipping), counted in chunks (no big temporaries)."""
    n = 0
    for i in range(0, iq.shape[0], step):
        c = iq[i:i + step]
        n += int(np.count_nonzero(c >= limit)) + int(np.count_nonzero(c <= -limit - 1))
    return n


# =============================================================================
# DEVICE: state tracking and the register-RSSI sweep (as in rssi_scan.py)
# =============================================================================

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

    def __init__(self, args, stop):
        self.args = args
        self.stop = stop
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
        self._reported = set()   # RSSI configs already written to the dump
        self._quiet = False      # True during warmup: no dump, no info
        self.actual = {}         # actual values reported by the device for the current config
        # per-frequency stats of raw reads in dB (Welford)
        self._rd_n = 0
        self._rd_mean = 0.0
        self._rd_m2 = 0.0
        self._pre = ffi.new("int32_t *")   # preallocated output buffers for RSSI reads
        self._sym = ffi.new("int32_t *")
        self._u8 = ffi.new("uint8_t *")    # preallocated buffer for RFIC register reads

    def open(self):
        a = self.args
        self.dev = bladerf.BladeRF()
        ch = self.dev.Channel(self.ch_id)
        self.ch = ch

        if self.dev.get_rfic_rssi(self.ch_id) is None:
            raise ConfigError("RFIC RSSI is not supported by this board")

        ch.gain_mode = _bladerf.GainMode.Manual
        ch.gain = a.gain
        self._read_gain("right after the request")

        DEVICE_STATE.reset()  # device state is unknown right after open
        self.rssi_cfg_sr = None
        first_freqs, first_wdt = self.spans[0]
        self.apply(sample_rate=first_wdt, bandwidth=first_wdt, frequency=first_freqs[0])

        self.dev.enable_module(self.ch_id, True)
        self.enabled = True
        self._read_gain("sweep start, tuned, RX enabled")

        t0 = time.perf_counter()
        for _ in range(CALIBRATION_READS):
            self._read_rssi()
        self.t_read = (time.perf_counter() - t0) / CALIBRATION_READS

        self.t_tune = 0.0
        self._warmup()

        # Retune time, measured over a few points of the first span
        # (the first probe point is the current frequency and is skipped by apply())
        probe = first_freqs[:: max(1, len(first_freqs) // CALIBRATION_TUNES)][:CALIBRATION_TUNES + 1]
        t0 = time.perf_counter()
        n_tunes = 0
        for f in probe:
            if self.apply(frequency=f):
                n_tunes += 1
        self.t_tune = (time.perf_counter() - t0) / n_tunes if n_tunes else 0.0
        info(f"RSSI read: {self.t_read * 1e3:.3f} ms/call, retune: {self.t_tune * 1e3:.3f} ms")

    def _warmup(self):
        mode = self.args.warmup
        if mode == "none":
            return
        t0 = time.perf_counter()
        self._quiet = True
        try:
            if mode == "pass":
                self.sweep(None)
            else:
                for freqs, wdt in self.spans:
                    self.apply(sample_rate=wdt, bandwidth=wdt)
                    probe = freqs[:: max(1, len(freqs) // WARMUP_TUNES)][:WARMUP_TUNES]
                    for f in probe:
                        self.apply(frequency=f)
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
                self.ch.frequency = int(round(value))
                self.actual[name] = self.ch.frequency
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
        DEVICE_STATE.reset()
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

    def sweep(self, reg):
        """One pass over all spans. If reg (dict of preallocated arrays freq, wdt, mean, std, n) is
        given, results are stored there. Returns True when complete, False if stopped.

        On each frequency: retune, settle, then read the hardware RSSI once per hardware
        averaging period until dwell has elapsed since the END of the retune; the check is
        done after each read, so there is always at least one. If a read starts late, the
        schedule restarts from that read. Mean is in linear power, std over the reads in dB."""
        dwell = self.args.dwell
        idx = 0
        for freqs, wdt in self.spans:
            changed = self.apply(sample_rate=wdt, bandwidth=wdt)
            sr, bw = self.actual["sample_rate"], self.actual["bandwidth"]
            if changed and bw != sr and not self._quiet:
                info(f"WARNING: wdt {fmt_freq(sr).strip()} MHz: actual bandwidth "
                     f"{fmt_freq(bw).strip()} MHz differs")
            if self.rssi_cfg_sr != sr:
                self._configure_rssi()
            if not self._quiet:
                print(fmt_header(REG_COLS), flush=True)
                self._read_gain(f"span wdt {fmt_freq(wdt).strip()} MHz, start")
            period = self.period
            for f in freqs:
                if self.stop.is_set():
                    return False
                self.apply(frequency=f)
                t_tuned = time.perf_counter()
                fr = self.actual["frequency"]
                t_next = t_tuned + self.settle
                sum_lin = 0.0
                self._reset_read_stats()
                while True:
                    wait_until(t_next)
                    t_start = time.perf_counter()
                    v = self._read_rssi()
                    t_end = time.perf_counter()
                    # next read one period after this one's start: on time -> t_next + period;
                    # late -> the schedule restarts, no catch-up burst of the same value
                    t_next = max(t_next, t_start) + period
                    sum_lin += db_to_lin(v)
                    self._rd_n += 1
                    d = v - self._rd_mean
                    self._rd_mean += d / self._rd_n
                    self._rd_m2 += d * (v - self._rd_mean)
                    if t_end - t_tuned >= dwell:   # dwell counted after the retune
                        break
                    if self.stop.is_set():
                        return False
                n = self._rd_n
                mean_db, std_db = lin_to_db(sum_lin / n), self._read_std_db()
                if reg is not None:
                    reg["freq"][idx], reg["wdt"][idx] = fr, wdt
                    reg["mean"][idx], reg["std"][idx], reg["n"][idx] = mean_db, std_db, n
                    idx += 1
                if not self._quiet:
                    std_txt = "n/a" if math.isnan(std_db) else f"{std_db:.2f}"
                    print(fmt_row(REG_COLS, (fmt_freq(fr), fmt_freq(wdt), f"{mean_db:.2f}", n, std_txt)),
                          flush=True)
            if not self._quiet:
                self._read_gain(f"span wdt {fmt_freq(wdt).strip()} MHz, end")
        return True

    # ----------------------------------------------------------- capture

    def tune_capture(self, cap):
        """Stop the RSSI mode and tune the receiver for the IQ capture. Returns actual (fs, bw, center)."""
        self.dev.enable_module(self.ch_id, False)
        self.enabled = False
        self.apply(sample_rate=cap.fs, bandwidth=cap.bw, frequency=cap.center)
        return self.actual["sample_rate"], self.actual["bandwidth"], self.actual["frequency"]

    def capture_into(self, iq, n_discard):
        """Fill the preallocated int16 array iq (n, 2) with IQ samples after dropping n_discard
        samples. Nothing else runs between the reads. Returns elapsed seconds, or None if stopped."""
        self.dev.sync_config(layout=_bladerf.ChannelLayout.RX_X1, fmt=_bladerf.Format.SC16_Q11,
                             num_buffers=SYNC_BUFFERS, buffer_size=SYNC_BUF_SAMPLES,
                             num_transfers=SYNC_TRANSFERS, stream_timeout=SYNC_TIMEOUT_MS)
        self.dev.enable_module(self.ch_id, True)
        self.enabled = True
        self._read_gain("capture start, tuned, RX enabled")
        n = iq.shape[0]
        left = n_discard
        while left > 0:                       # stale/transient samples go to the start of iq (overwritten later)
            k = min(left, n, SYNC_CHUNK)
            self.dev.sync_rx(iq[:k], k, SYNC_TIMEOUT_MS)
            left -= k
            if self.stop.is_set():
                return None
        t0 = time.perf_counter()
        pos = 0
        while pos < n:
            k = min(n - pos, SYNC_CHUNK)
            self.dev.sync_rx(iq[pos:pos + k], k, SYNC_TIMEOUT_MS)
            pos += k
            if self.stop.is_set():
                return None
        elapsed = time.perf_counter() - t0
        self._read_gain("capture end")
        self.dev.enable_module(self.ch_id, False)
        self.enabled = False
        return elapsed


# ---------------------------------------------------------------------- main flow

def print_table(cols, rows):
    print(fmt_header(cols), flush=True)
    for r in rows:
        print(fmt_row(cols, r))


def run(args, stop):
    spans = args.freq_spans
    dwell = args.dwell
    # ---- 1. checks that need no hardware: refuse / conflict BEFORE the sweep starts
    lo, hi = range_of_spans(spans)
    cap = plan_capture(lo, hi, args.max_bw)
    wp = plan_welch(cap.fs, dwell, args.df, args.fft_n, args.avr_n, args.avr_dur, args.overlap)
    n_win, policy = choose_windows(wp, dwell)
    check_capture_size(welch_samples(wp.nfft, wp.hop, n_win), args.max_capture_mb)
    for note in wp.notes:
        info(note)
    info(f"capture band {lo * HZ_TO_MHZ:.2f}..{hi * HZ_TO_MHZ:.2f} MHz -> center {cap.center * HZ_TO_MHZ:.2f} MHz, "
         f"fs {cap.fs * HZ_TO_MHZ:.3f} MHz, analog bandwidth {cap.bw * HZ_TO_MHZ:.3f} MHz")

    # ---- 2. register sweep
    n_points = sum(len(f) for f, _ in spans)
    reg = {k: np.empty(n_points) for k in ("freq", "wdt", "mean", "std")}
    reg["n"] = np.empty(n_points, dtype=np.int64)
    sc = Scanner(args, stop)
    try:
        sc.open()
        if stop.is_set():
            info("interrupted")
            return 130
        print(sc.header())
        t0 = time.perf_counter()
        if not sc.sweep(reg):
            info("interrupted")
            return 130
        t_sweep = time.perf_counter() - t0

        # ---- 3. one IQ capture
        fs, bw, fc = sc.tune_capture(cap)
        wp = plan_welch(fs, dwell, args.df, args.fft_n, args.avr_n, args.avr_dur, args.overlap)
        n_win, policy = choose_windows(wp, dwell, policy)
        n_samples = welch_samples(wp.nfft, wp.hop, n_win)
        check_capture_size(n_samples, args.max_capture_mb)
        iq = np.empty((n_samples, 2), dtype=np.int16)          # the only large allocation
        n_discard = max(DISCARD_MIN, int(DISCARD_S * fs))
        elapsed = sc.capture_into(iq, n_discard)
        if elapsed is None:
            info("interrupted")
            return 130
    finally:
        sc.close()

    # ---- 4. spectrum and comparison (no device needed any more)
    w = Welch(wp.nfft, args.window, n_win)
    w.accumulate(iq, wp.hop, n_win)
    p = w.power_spectrum()
    bws = np.full(n_points, args.rssi_bw) if args.rssi_bw else reg["wdt"]
    spec_db = band_power_db(p, fs, wp.nfft, fc, reg["freq"], bws)
    clipped = count_clipped(iq)

    dur_ms = n_samples / fs * 1e3
    print(sc.gain_summary())
    print(f"# register sweep took {t_sweep:.3f} s; capture: center {fc * HZ_TO_MHZ:.3f} MHz, fs {fs * HZ_TO_MHZ:.3f} MHz, "
          f"bandwidth {bw * HZ_TO_MHZ:.3f} MHz, {n_samples} samples = {dur_ms:.3f} ms "
          f"(read in {elapsed * 1e3:.1f} ms)")
    print(f"# spectrum: FFT {wp.nfft}, df {wp.df:.3f} Hz, window {args.window}, overlap {args.overlap:g}, "
          f"{n_win} windows averaged = {dur_ms:.3f} ms vs dwell {dwell * 1e3:.3f} ms")
    if wp.n_user is not None:
        info(f"averaging from --avr-*: {wp.n_user} windows, used {n_win} (policy: {policy})")
    if n_win == 1 and wp.nfft / fs > dwell:
        info("NOTE: one FFT window is longer than dwell; the capture is longer than dwell")
    if np.min(bws) / wp.df < MIN_BINS_PER_BAND:
        info(f"WARNING: RSSI band holds only {np.min(bws) / wp.df:.1f} FFT bins (< {MIN_BINS_PER_BAND}); "
             f"the spectral estimate is coarse, reduce --df / raise --fft-n")
    if clipped:
        info(f"WARNING: {clipped} clipped samples in the capture; reduce the gain, the estimate is invalid")

    diff = reg["mean"] - spec_db
    rows = ((fmt_freq(reg["freq"][i]), fmt_freq(reg["wdt"][i]), f"{reg['mean'][i]:.2f}",
             f"{spec_db[i]:.2f}", f"{diff[i]:.2f}") for i in range(n_points))
    print_table(CMP_COLS, rows)

    def stats_line(label, st):
        return (f"# {label} ({st.n} points): diff reg-fft mean {st.mean:.2f} dB, std {st.std:.2f} dB; "
                f"offset (median) {st.offset:.2f} dB; residual RMS {st.rms_res:.2f} dB, max {st.max_res:.2f} dB; "
                f"correlation {st.corr:.3f}")
    print(stats_line("all", compare_stats(reg["mean"], spec_db)))
    return 0


# ---------------------------------------------------------------------- entry

def install_sigint(stop):
    """First Ctrl+C: cooperative stop (current call finishes, device is closed cleanly).
    Second Ctrl+C: immediate exit without cleanup."""
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
    try:
        return run(args, stop)
    except ConfigError as e:
        info(f"ERROR: {e}")
        return 2


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