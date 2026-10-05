#!/usr/bin/env python3
"""bladeRF 2.0: waterfall + PSD of one wide capture, with RSSI for the channels that lie inside it.

What it does
  One capture window (center, rate) is streamed continuously. Every `row` of the waterfall is one Welch
  PSD over a capture interval (the averaging). The newest row is at the top, older rows move down.
  RSSI of every channel inside the capture band is computed from the representative PSD (the latest
  row, or the mean of the rows in the waterfall) as the power inside [f - bw/2, f + bw/2].

CLI (frequencies in Hz: 40M, 1.328G, 25k, 48e6; a unit may follow: 25 kHz)
  -c, --center <Hz>                        capture center; default: middle of the channel ranges
  -r, --rate <Hz>                          sample rate, default 40e6
  -g, --gain <dB>                          manual RX gain, default 42
  -w, --bandwidth <Hz>                     analog bandwidth, default: the sample rate (limited to 56 MHz)
  --freq-spans "[{flo~fhi,step,wdt},...]"  channels, as before: centers flo..fhi with step, width wdt
  --freqs-list "[{freq,wdt},...]"          instead of --freq-spans (one or the other): every channel given by its
                                           centre and width. Only channels that lie fully inside the capture band are used.
  --ring-dur <sec>                         seconds of raw samples kept in a ring buffer in memory, default 3
  --df <Hz> | --fft-n <int>                FFT resolution guide / size (2^a*3^b*5^c*7^d); default --fft-n 1024
  --avr-n <int> | --avr-dur <sec>          windows per row / duration of one row; default --avr-dur 0.05
  --window {rect,hann,hamming,blackman}    default rect (no multiplication)
  --overlap <0..1)                         window overlap, default 0
  --rows <int>                             waterfall depth in rows, default 100
  --nb-max-width <Hz>                      NarrowBand Det: widest signal that counts as narrowband, default 25k
  --nb-guard <frac>                        NarrowBand Det: SNR noise region = signal width x (1 + guard), default 1.25
  --nb-thr <dB>                            NarrowBand Det: threshold above the local noise floor, default 6
  --metrics <file>                         performance metrics (TSV): one row per second, '# event' and '# summary' lines
  --cap-rate <Hz>                          width of a captured band = its output sample rate, default 48k. A guide: the
                                           nearest rate that a ratio of whole numbers L/M of the capture rate gives is used
  --cap-dir <dir>                          folder of the recordings, default ./captures
  --cap-pattern <text>                     file name (without the extension), default {name}_{time}_{freq}MHz. Fields:
                                           {name} {time} {utc} {freq} {rate} {id} {n}; a format after a colon:
                                           {time:%Y-%m-%d_%H-%M-%S} {freq:.3f} {n:03d}; other text is kept; {{ }} = braces

Definitions
  Power is relative to full scale: a complex tone of 2048 counts amplitude = 0 dBFS.
  PSD is shown in dBFS/Hz. RSSI is the PSD integrated over the channel band, in dBFS.
  Gain is shown as reported by the device after the stream has started (it can differ from the request).

NarrowBand Det (the switch above the table; detection runs on the displayed PSD in both modes, only the view changes)
  Bands: bins above the local noise floor (block medians) plus the threshold; the band edges are the bins within
  10 dB of the band peak; width = bins x resolution; only bands up to --nb-max-width count as narrowband.
  Power = PSD summed over the band, dBFS. SNR = band power / noise power in the same band, dB; the noise is the
  median per-bin power over a region centred on the signal, width = signal width x (1 + guard) (the region
  includes the signal's own bins, so the SNR comes out low: by several dB when the signal fills about half of
  it, and for a signal of only a few bins the region is a few bins wide and lies in the window leakage).
  A band that the leakage of a stronger signal explains (sidelobes of the chosen window) is dropped.
  The threshold (floor + Threshold, dashed) and the floor (dotted) are drawn on the PSD: Display, Detection threshold.
  Signals are tracked by frequency. Present = found in the displayed PSD. Last seen = start time of the last PSD
  window in which the signal was found. Time is counted from samples: t = stream sample / sample rate, with t = 0 at
  the moment the device started the wideband stream at the capture centre (the first sample it delivers; the start-up
  samples that are thrown away are counted: stream sample = ring sample + discarded; rows skipped by the processing too).
  For the mean PSD the window starts at the oldest row of the waterfall.

Window: the left panel changes the launch parameters, applied by "Apply and restart";
the Display and NarrowBand Det groups apply immediately.
The groups fold and unfold by a click on their title. In Spectrum, FFT size / resolution and windows / duration are
pairs: edit either one of a pair and the other follows (the value you typed is the guide, the shown ones are reachable). Ctrl+Enter = apply and restart, Space = pause (not while
typing in a field), Esc = close the window.
Plots: mouse wheel = zoom the frequency axis around the cursor (both plots), Ctrl+wheel on the PSD = zoom its
levels, drag = pan, double-click or "Reset zoom" = back to the whole capture band.
Capture goes into a ring buffer of the last --ring-dur seconds of raw samples (memory only, nothing is written to
disk). The search reads its rows from the ring; ring sample i is stream sample i + discarded (index 0 = the first
sample after the start-up discard), time = stream sample / sample rate. If the processing falls behind by more than the ring holds, the rows that were
overwritten are skipped (counted as Dropped; hover the field for the three causes); a lag of more than a row is shown
next to the buffer in the status.
Capture of a band: Shift + left click on a plot places the band there (or moves it; during a recording this asks and
stops the capture), or pick a signal (NarrowBand Det list) or a channel (RSSI list), or drag/type a band. The panel below
the list holds centre, width (output rate), name and folder; the overlay of the band on the waterfall and the PSD is the
same object: drag its edges to resize, drag it to move (it follows the fields and the fields follow it). Start, Pause /
Resume and Stop control the recording of the band from the ring buffer; Stop closes the file. A file belongs to one band: the centre, the
width, the name, the pattern and the folder are fixed while it is being written. Another centre asked for during a
recording (typed, picked in a list, or the band dragged) is confirmed in a pop-up: the file is closed and a new one is
opened. During a pause any of them can be changed: the button then says "Start new file"("Start new file"). The recording is made from the ring buffer into <stem>.wav: a WAV file (RIFF/WAVE, PCM, the canonical
44-byte header, 2 channels of int16: I = left, Q = right, counts as the radio gives them, full scale 2048) and <stem>.json,
the sidecar with the context of the capture (radio, analysis, the signal, the band, every segment with its centre and
time, the result). The header is written first with empty sizes, refreshed every few seconds and at a pause, and finished
when the file is closed. The header holds the sample rate rounded to whole Hz (the exact one is in the sidecar); a file
is closed at 4 GiB (the WAV limit) and the recording goes on in the next one. The file name
comes from the pattern (field of the panel, --cap-pattern); an existing file is never overwritten: a number is added (or
{n} counts up). The rate is
the capture rate x L / M with small whole numbers; the actual rate is shown. Pause / lost blocks start a new segment.
Signals: raw detections of every PSD are kept per signal and the registry shows their medians; a signal is listed after
REG_MIN_HITS (3) detections.
Performance: every second the capture thread, the processing thread and the window report how they are doing (see
--metrics for the columns). The same numbers are in the tooltip of the Rows field of the status bar, and a summary
with an ANOMALIES line is printed when a run ends (stderr, and the --metrics file).
Reading a range of the ring (SampleRing.copy_range) returns a RingStatus: OK, OVERWRITTEN and LOST_WHILE_COPYING (both
mean we did not keep up; `late`), NOT_YET (the end is still being recorded: wait for it) or INVALID.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import queue
import re
import signal
import struct
import sys
import threading
import time
from dataclasses import dataclass
from enum import IntEnum
from fractions import Fraction
from tkinter import font as tkfont
import tkinter as tk
from tkinter import ttk

import numpy as np



def _find_bladerf_dll():
    """Windows: make bladeRF.dll (and the libusb-1.0.dll next to it) loadable even when the PATH of this shell lacks their
    folder. Looks in BLADERF_DLL_DIR, in the PATH entries and in Program Files\\bladeRF. Returns the folder or None."""
    if sys.platform != "win32":
        return None
    dirs = [os.environ.get("BLADERF_DLL_DIR", "")] + os.environ.get("PATH", "").split(os.pathsep)
    for var in ("ProgramFiles", "ProgramW6432", "ProgramFiles(x86)"):
        base = os.environ.get(var)
        if base:
            dirs += [os.path.join(base, "bladeRF", "x64"), os.path.join(base, "bladeRF")]
    seen = set()
    for d in dirs:
        d = d.strip().strip('"')
        if not d or d.lower() in seen or not os.path.isdir(d):
            continue
        seen.add(d.lower())
        if os.path.isfile(os.path.join(d, "bladeRF.dll")):
            os.environ["PATH"] = d + os.pathsep + os.environ.get("PATH", "")
            try:
                os.add_dll_directory(d)
            except (AttributeError, OSError):
                pass
            return d
    return None


_BLADERF_DLL_DIR = _find_bladerf_dll()
try:
    import bladerf
    from bladerf import _bladerf
except OSError as _e:                      # the DLL or one of its own dependencies cannot be loaded
    _where = _BLADERF_DLL_DIR
    _lines = [f"ERROR: the bladeRF library cannot be loaded: {_e}"]
    if _where:
        _usb = os.path.isfile(os.path.join(_where, "libusb-1.0.dll"))
        _lines.append(f"bladeRF.dll was found in {_where}; libusb-1.0.dll there: {'yes' if _usb else 'NO (put it next to bladeRF.dll)'}")
    else:
        _lines.append("bladeRF.dll was not found in BLADERF_DLL_DIR, in the PATH of this shell or in Program Files\\bladeRF.")
        _lines.append("Set BLADERF_DLL_DIR to the folder with bladeRF.dll and libusb-1.0.dll, or add it to PATH "
                      "(where.exe bladeRF.dll shows what Windows finds).")
    sys.stderr.write("\n".join(_lines) + "\n")
    sys.exit(2)

# ---------------------------------------------------------------- constants

# bladeRF 2.0 limits (libbladeRF 2025.10, fpga_common/include/bladerf2_common.h)
SR_MIN, SR_MAX = 520834.0, 61440000.0
BW_MIN, BW_MAX = 200000.0, 56000000.0
FREQ_MIN, FREQ_MAX = 70e6, 6e9
FS_FULL_SCALE = 2048.0       # SC16_Q11: full scale per component, counts

DEFAULT_RATE = 40e6
DEFAULT_GAIN = 42
DEFAULT_FFT_N = 1024
DEFAULT_AVR_DUR = 0.05
DEFAULT_ROWS = 100
ROWS_MIN, ROWS_MAX = 8, 1000
FFT_MAX = 65536
MAX_BLOCK_MB = 128           # one row of samples must not be larger than this (assumption)

# Streaming (assumptions: typical libbladeRF usage, not tuned on hardware)
SYNC_BUFFERS = 32
SYNC_BUF_SAMPLES = 65536
SYNC_TRANSFERS = 16
SYNC_TIMEOUT_MS = 1000
SYNC_CHUNK = 1 << 20
DISCARD_MIN = 131072
DISCARD_S = 0.01
REG_HISTORY = 64             # raw detections kept per signal for the aggregation
REG_MIN_HITS = 3             # a signal enters the registry list after this many detections (assumption)
RAW_LOG_CAP = 16384          # raw detections of the last PSD windows
WAV_HEADER_BYTES = 44
WAV_MAX_DATA = 4 * 1024 ** 3 - 1 - WAV_HEADER_BYTES - (1 << 20)   # sizes in a RIFF header are 32 bits: stay below, with margin
HEADER_REFRESH_S = 5.0       # the header and the sidecar of an open file are brought up to date this often
DEFAULT_CAP_PATTERN = "{name}_{time}_{freq}MHz"   # the file name before the pattern existed
CAP_PATTERN_PRESETS = ("{name}_{time}_{freq}MHz", "{name}_{time:%Y-%m-%d_%H-%M-%S}", "{time:%Y%m%d_%H%M%S}_{name}",
                       "{name}_{freq:.3f}MHz_{n:03d}", "REC_{name}_{n}")
CAP_DEFAULT_RATE = 48e3      # output sample rate = width of a captured band, default
CAP_RATIO_MAX_DEN = 4096     # largest M in the rate ratio L/M (assumption)
CAP_PASS_FRAC = 0.8          # flat part of the passband as a fraction of the output band (assumption)
CAP_OVERLAP = 32             # output samples that are cut off at each end of a block (assumption)
CAP_N_MAX = 1 << 21          # largest FFT of the band extractor
PERF_PERIOD_S = 1.0          # seconds between performance snapshots
PERF_SLOW_RATIO = 0.98       # received samples / (rate x time) below this means samples were lost (assumption)
DEFAULT_RING_DUR = 3.0       # seconds of raw samples in the ring buffer (assumption)
MAX_RING_MB = 2048           # the ring buffer must not be larger than this (assumption)
RING_MIN_ROWS = 4            # the ring holds at least this many rows
RING_CHUNK_S = 0.025         # the capture thread writes the ring in pieces of about this many seconds (assumption)
ROW_RING = 16                # PSD rows in flight between the processing thread and the window
WELCH_CHUNK_BYTES = 8 << 20
MAX_CONSECUTIVE_READ_ERRORS = 3

POLL_MS = 40
DRAW_MIN_S = 0.08            # at most ~12 redraws per second; a slow draw lowers the rate: the window keeps half its time free
TEXT_MIN_S = 0.25
STATUS_MIN_S = 0.5
LN10 = math.log(10.0)
ZOOM_STEP = 1.25             # span factor per wheel notch
ZOOM_MIN_BINS = 16           # the zoomed span never gets narrower than this many FFT bins
TEXT_INPUT_CLASSES = ("TEntry", "Entry", "Text", "TCombobox", "TButton", "TCheckbutton", "TRadiobutton")  # use Space themselves

DEFAULT_NB_MAX_WIDTH = 25e3
DEFAULT_NB_GUARD = 1.25
DEFAULT_NB_THR = 6.0
NB_TRACK_LIMIT = 256         # tracked signals; when full, the absent one seen longest ago is replaced
NB_MAX_DET = 128             # strongest detections kept per PSD
NB_EDGE_DB = 10.0            # band edges: bins within this many dB of the band peak (assumption)
NB_FLOOR_MIN_BINS = 128      # noise-floor block: at least this many bins ...
NB_FLOOR_WIDTHS = 16         # ... and at least this many maximum widths (assumption)
NB_SEG_GAP = 2               # bins of a band may be separated by gaps of up to this many bins
NB_PASSES = 3                # a band is cut out and the remaining candidates searched again: weaker neighbours
NB_LEAK_MARGIN_DB = 6.0      # a band must stand this far above the window leakage of a stronger signal (assumption)
NB_LEAK_MAX_RUNS = 400       # candidate bands compared pairwise for leakage (the strongest ones)

WINDOWS = ("rect", "hann", "hamming", "blackman")
CMAPS = ("inferno", "magma", "viridis", "cividis")
PSD_SOURCES = ("Latest row", "Mean of waterfall")

_SUFFIX = {"": 1.0, "k": 1e3, "K": 1e3, "M": 1e6, "G": 1e9}
HZ_TO_MHZ = 1e-6


def info(msg):
    print(msg, file=sys.stderr, flush=True)


def fmt_freq(hz, digits=2):
    """Single formatter for every frequency-like value (Hz) in text output -> MHz."""
    return f"{hz * HZ_TO_MHZ:{digits + 4}.{digits}f}"


def fmt_hz_short(hz, digits=9):
    """Hz with a k/M/G suffix, as typed on the command line: 40000000 -> '40M'."""
    for div, suf in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(hz) >= div:
            return f"{hz / div:.{digits}g}{suf}"
    return f"{hz:.{digits}g}"


# ---------------------------------------------------------------- parsing

_HZ_RE = re.compile(r"\s*((?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)\s*([kKMG]?)\s*(?:[hH][zZ])?\s*")


def parse_hz(text):
    """A frequency in Hz: 40M, 1.328G, 25k, 48e6, 4.8e7, 48.e6, with an optional unit and spaces (48 MHz, 25 kHz, 48e6 Hz)."""
    m = _HZ_RE.fullmatch(text)
    if not m:
        raise ValueError(f"bad frequency value: {text!r} (examples: 40M, 48e6, 1.328G, 25 kHz)")
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
        raise ValueError(f"bad channel ranges: {text!r} (example: [{{1.30G~1.35G,1M,1.25M}}])")
    spans = Spans()
    spans.ranges = []
    for item in items:
        parts = [p.strip() for p in item.split(",")]
        if len(parts) != 3 or "~" not in parts[0]:
            raise ValueError(f"a range must be {{flo~fhi,step,wdt}}: {{{item}}}")
        lo, hi = (parse_hz(x) for x in parts[0].split("~", 1))
        step, wdt = parse_hz(parts[1]), parse_hz(parts[2])
        if hi < lo or step <= 0 or wdt <= 0:
            raise ValueError(f"bad range values: {{{item}}}")
        count = int(math.floor((hi - lo) / step + 1e-9)) + 1
        spans.append((lo + step * np.arange(count, dtype=np.float64), wdt))
        spans.ranges.append((lo, hi))
    return spans


def parse_freqs_list(text):
    """[{freq,wdt},...]: every channel by its centre and width. Same container as parse_spans: one-point spans."""
    body = text.strip()
    if body.startswith("[") and body.endswith("]"):
        body = body[1:-1]
    items = re.findall(r"\{([^{}]*)\}", body)
    if not items or re.sub(r"\{[^{}]*\}", "", body).strip(" ,"):
        raise ValueError(f"bad channel list: {text!r} (example: [{{1.30G,1.25M}},{{1.31G,2M}}])")
    spans = Spans()
    spans.ranges = []
    for item in items:
        parts = [p.strip() for p in item.split(",")]
        if len(parts) != 2:
            raise ValueError(f"an entry must be {{freq,wdt}}: {{{item}}}")
        f, wdt = parse_hz(parts[0]), parse_hz(parts[1])
        if f <= 0 or wdt <= 0:
            raise ValueError(f"bad entry values: {{{item}}}")
        spans.append((np.array([f]), wdt))
        spans.ranges.append((f, f))
    return spans


# ---------------------------------------------------------------- FFT sizes

def _build_fft_sizes(limit=FFT_MAX * 4):
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


# ---------------------------------------------------------------- configuration

class ConfigError(Exception):
    """Invalid launch parameter; `field` names the control that has to be fixed."""

    def __init__(self, field, message):
        super().__init__(message)
        self.field = field


@dataclass
class Config:
    center: float | None         # None = middle of the channel ranges
    rate: float
    bandwidth: float | None      # None = the sample rate
    gain: int
    fft_mode: str                # 'df' | 'fft_n'
    fft_value: float
    avr_mode: str                # 'dur' | 'n'
    avr_value: float
    window: str
    overlap: float
    spans: Spans
    rows: int
    nb_max_width: float
    nb_guard: float
    nb_thr: float
    ring_dur: float


@dataclass
class Plan:
    """Everything derived from a Config that needs no device."""
    center: float
    bandwidth: float             # effective analog bandwidth requested from the device
    nfft: int
    hop: int
    n_win: int
    row_samples: int
    notes: list
    ch_f: np.ndarray             # all channel centers, Hz
    ch_bw: np.ndarray            # their widths, Hz
    ring_rows: int               # rows in the ring buffer


def _field_value(text, field, parser, what):
    try:
        return parser(text)
    except ValueError as e:
        raise ConfigError(field, f"{what}: {e}") from None


def ring_chunk(rate):
    """Samples per write into the ring: about RING_CHUNK_S of signal, so that the rows reach the processing with a
    delay of that order at every rate, but never fewer than one USB buffer or more than SYNC_CHUNK."""
    return int(min(SYNC_CHUNK, max(SYNC_BUF_SAMPLES, RING_CHUNK_S * rate)))


def welch_samples(nfft, hop, n_win):
    return (n_win - 1) * hop + nfft


def parse_nb(d):
    """NarrowBand Det parameters from strings: (max width Hz, guard fraction, threshold dB)."""
    mw = _field_value(d["nb_max_width"], "nb_max_width", parse_hz, "Max width")
    if mw <= 0:
        raise ConfigError("nb_max_width", "Max width must be above 0")
    try:
        guard = float(d["nb_guard"].strip())
    except ValueError:
        raise ConfigError("nb_guard", "Guard must be a number (a fraction of the signal width)") from None
    if guard < 0:
        raise ConfigError("nb_guard", "Guard must be 0 or more")
    try:
        thr = float(d["nb_thr"].strip())
    except ValueError:
        raise ConfigError("nb_thr", "Threshold must be a number of dB") from None
    if thr <= 0:
        raise ConfigError("nb_thr", "Threshold must be above 0 dB")
    return mw, guard, thr


def validate(d):
    """d: dict of strings from the command line or the panel. Returns (Config, Plan) or raises ConfigError."""
    rate = _field_value(d["rate"], "rate", parse_hz, "Rate")
    if not SR_MIN <= rate <= SR_MAX:
        raise ConfigError("rate", f"Rate {fmt_hz_short(rate)} is outside {SR_MIN * 1e-6:.2f}..{SR_MAX * 1e-6:.2f} MS/s")

    bw_txt = d["bandwidth"].strip()
    bandwidth = None if bw_txt.lower() in ("", "auto") else _field_value(bw_txt, "bandwidth", parse_hz, "Bandwidth")
    if bandwidth is not None and bandwidth <= 0:
        raise ConfigError("bandwidth", "Bandwidth must be above 0")

    try:
        gain = int(d["gain"].strip())
    except ValueError:
        raise ConfigError("gain", "Gain must be a whole number of dB") from None

    spans_txt = d["spans"].strip()
    parser = parse_freqs_list if d.get("spans_mode", "spans") == "list" else parse_spans
    try:
        spans = parser(spans_txt) if spans_txt else Spans()
    except ValueError as e:
        raise ConfigError("spans", str(e)) from None
    if not spans:
        spans.ranges = []

    c_txt = d["center"].strip()
    center = None if c_txt.lower() in ("", "auto") else _field_value(c_txt, "center", parse_hz, "Center")
    if center is None:
        if not spans.ranges:
            raise ConfigError("center", "Set the center frequency, or enter a channel range to derive it from")
        center = float(round((min(r[0] for r in spans.ranges) + max(r[1] for r in spans.ranges)) / 2))
    if not FREQ_MIN <= center <= FREQ_MAX:
        raise ConfigError("center", f"Center {fmt_hz_short(center)} is outside {FREQ_MIN * 1e-6:.0f} MHz..{FREQ_MAX * 1e-9:.0f} GHz")

    notes = []
    fft_mode, fft_txt = d["fft_mode"], d["fft_value"].strip()
    if fft_mode == "df":
        df = _field_value(fft_txt, "fft", parse_hz, "Resolution")
        if df <= 0:
            raise ConfigError("fft", "Resolution must be above 0")
        fft_value = df
        need = rate / df
        if need > FFT_MAX:
            raise ConfigError("fft", f"Resolution {fmt_hz_short(df)} Hz needs about {need:.0f} bins; the limit is {FFT_MAX}")
        nfft = nearest_fft_size(need)
        notes.append(f"resolution {fmt_hz_short(df)} Hz -> FFT size {nfft}, actual {rate / nfft:.3f} Hz")
    else:
        try:
            n = int(fft_txt)
        except ValueError:
            raise ConfigError("fft", "FFT size must be a whole number") from None
        if n < 16:
            raise ConfigError("fft", "FFT size must be at least 16")
        fft_value = float(n)
        nfft = n if is_fft_size(n) else nearest_fft_size(n)
        if nfft != n:
            notes.append(f"FFT size {n} is not 2^a*3^b*5^c*7^d: using {nfft}")
    if nfft > FFT_MAX:
        raise ConfigError("fft", f"FFT size {nfft} is above {FFT_MAX}; use a coarser resolution")

    try:
        overlap = float(d["overlap"].strip())
    except ValueError:
        raise ConfigError("overlap", "Overlap must be a number from 0 to below 1") from None
    if not 0.0 <= overlap < 1.0:
        raise ConfigError("overlap", "Overlap must be from 0 to below 1")
    hop = max(1, int(round(nfft * (1.0 - overlap))))

    avr_mode, avr_txt = d["avr_mode"], d["avr_value"].strip()
    if avr_mode == "n":
        try:
            n_win = int(avr_txt)
        except ValueError:
            raise ConfigError("avr", "Windows must be a whole number") from None
        if n_win < 1:
            raise ConfigError("avr", "Windows must be at least 1")
        avr_value = float(n_win)
    else:
        try:
            dur = float(avr_txt)
        except ValueError:
            raise ConfigError("avr", "Duration must be a number of seconds") from None
        if dur <= 0:
            raise ConfigError("avr", "Duration must be above 0")
        avr_value = dur
        n_win = max(1, int(round((dur * rate - nfft) / hop)) + 1)
        notes.append(f"averaging {dur:g} s -> {n_win} windows ({welch_samples(nfft, hop, n_win) / rate * 1e3:.2f} ms)")
    row_samples = welch_samples(nfft, hop, n_win)
    mb = row_samples * 4 / (1 << 20)
    if mb > MAX_BLOCK_MB:
        raise ConfigError("avr", f"One row needs {mb:.0f} MB of samples (limit {MAX_BLOCK_MB} MB); "
                                 f"average over fewer windows or lower the rate")

    window = d["window"]
    if window not in WINDOWS:
        raise ConfigError("window", f"Window must be one of {', '.join(WINDOWS)}")

    try:
        rows = int(d["rows"].strip())
    except ValueError:
        raise ConfigError("rows", "Rows must be a whole number") from None
    if not ROWS_MIN <= rows <= ROWS_MAX:
        raise ConfigError("rows", f"Rows must be {ROWS_MIN}..{ROWS_MAX}")

    try:
        ring_dur = float(d["ring_dur"].strip())
    except ValueError:
        raise ConfigError("ring_dur", "Buffer must be a number of seconds") from None
    if ring_dur <= 0:
        raise ConfigError("ring_dur", "Buffer must be above 0 s")
    ring_rows = max(RING_MIN_ROWS, math.ceil(ring_dur * rate / row_samples), math.ceil(4 * ring_chunk(rate) / row_samples))
    ring_mb = ring_rows * row_samples * 4 / (1 << 20)
    if ring_mb > MAX_RING_MB:
        raise ConfigError("ring_dur", f"A buffer of {ring_dur:g} s at {rate * 1e-6:.2f} MS/s needs {ring_mb:.0f} MB "
                                      f"(limit {MAX_RING_MB} MB); shorten it or lower the rate")
    notes.append(f"buffer: {ring_rows} rows = {ring_rows * row_samples / rate:.2f} s, {ring_mb:.0f} MB")

    bw_eff = min(max(bandwidth if bandwidth is not None else rate, BW_MIN), BW_MAX)
    if bandwidth is None and rate > BW_MAX:
        notes.append(f"bandwidth limited to {BW_MAX * 1e-6:.0f} MHz (the device maximum)")
    elif bandwidth is not None and bw_eff != bandwidth:
        notes.append(f"bandwidth {fmt_hz_short(bandwidth)} limited to {fmt_hz_short(bw_eff)} (device range)")

    if spans:
        ch_f = np.concatenate([f for f, _ in spans])
        ch_bw = np.concatenate([np.full(len(f), w) for f, w in spans])
    else:
        ch_f = ch_bw = np.empty(0)

    nb_max_width, nb_guard, nb_thr = parse_nb(d)
    if rate / nfft > nb_max_width:
        notes.append(f"NarrowBand Det: resolution {rate / nfft / 1e3:.3f} kHz is coarser than Max width "
                     f"{fmt_hz_short(nb_max_width)}Hz, nothing can be detected")

    cfg = Config(None if c_txt.lower() in ("", "auto") else center, rate, bandwidth, gain, fft_mode, fft_value,
                 avr_mode, avr_value, window, overlap, spans, rows, nb_max_width, nb_guard, nb_thr, ring_dur)
    plan = Plan(center, bw_eff, nfft, hop, n_win, row_samples, notes, ch_f, ch_bw, ring_rows)
    return cfg, plan


# --------------------------------------------------------------------- command line

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="bladeRF 2.0 wide-capture waterfall with channel RSSI")
    p.add_argument("-c", "--center", help="capture center, Hz (default: middle of the channel ranges)")
    p.add_argument("-r", "--rate", default=fmt_hz_short(DEFAULT_RATE), help="sample rate, Hz (default 40e6)")
    p.add_argument("-g", "--gain", default=str(DEFAULT_GAIN), help="manual RX gain, dB (default 42)")
    p.add_argument("-w", "--bandwidth", help="analog bandwidth, Hz (default: the sample rate)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--freq-spans", default="", help='channels, e.g. "[{1.30G~1.35G,1M,1.25M}]"')
    g.add_argument("--freqs-list", default="", help='channels as a list, e.g. "[{1.30G,1.25M},{1.31G,2M}]"')
    p.add_argument("--ring-dur", default=str(DEFAULT_RING_DUR),
                   help="seconds of raw samples kept in the ring buffer (default 3)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--df", help="FFT resolution guide, Hz")
    g.add_argument("--fft-n", help="FFT size")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--avr-n", help="Welch windows per row")
    g.add_argument("--avr-dur", help="duration of one row, s")
    p.add_argument("--window", choices=WINDOWS, default="rect", help="FFT window (default rect: no multiplication)")
    p.add_argument("--overlap", default="0", help="window overlap, 0..1 (default 0)")
    p.add_argument("--rows", default=str(DEFAULT_ROWS), help="waterfall rows (default 100)")
    p.add_argument("--nb-max-width", default=fmt_hz_short(DEFAULT_NB_MAX_WIDTH),
                   help="NarrowBand Det: widest signal that counts as narrowband, Hz (default 25k)")
    p.add_argument("--nb-guard", default=str(DEFAULT_NB_GUARD),
                   help="NarrowBand Det: SNR noise region = signal width x (1 + guard) (default 1.25)")
    p.add_argument("--metrics", metavar="FILE", help="write performance metrics (TSV) to FILE")
    p.add_argument("--cap-rate", default=fmt_hz_short(CAP_DEFAULT_RATE),
                   help="width of a captured band = output sample rate, Hz (default 48k)")
    p.add_argument("--cap-dir", default="captures", help="folder of the recordings (default ./captures)")
    p.add_argument("--cap-pattern", default=DEFAULT_CAP_PATTERN,
                   help="file name without the extension: {name} {time} {utc} {freq} {rate} {id} {n}, e.g. "
                        "{time:%%Y-%%m-%%d_%%H-%%M-%%S}_{name} (default {name}_{time}_{freq}MHz)")
    p.add_argument("--nb-thr", default=str(DEFAULT_NB_THR),
                   help="NarrowBand Det: threshold above the local noise floor, dB (default 6)")
    return p.parse_args(argv)


def args_to_fields(a):
    """Command-line values as the strings the panel uses, so both go through validate()."""
    if a.df is not None:
        fft_mode, fft_value = "df", a.df
    else:
        fft_mode, fft_value = "fft_n", a.fft_n if a.fft_n is not None else str(DEFAULT_FFT_N)
    if a.avr_n is not None:
        avr_mode, avr_value = "n", a.avr_n
    else:
        avr_mode, avr_value = "dur", a.avr_dur if a.avr_dur is not None else str(DEFAULT_AVR_DUR)
    return {"center": a.center if a.center is not None else "auto", "rate": a.rate,
            "bandwidth": a.bandwidth if a.bandwidth is not None else "auto", "gain": a.gain,
            "fft_mode": fft_mode, "fft_value": fft_value, "avr_mode": avr_mode, "avr_value": avr_value,
            "window": a.window, "overlap": a.overlap, "spans": a.freqs_list or a.freq_spans,
            "spans_mode": "list" if a.freqs_list else "spans", "ring_dur": a.ring_dur, "rows": a.rows,
            "nb_max_width": a.nb_max_width, "nb_guard": a.nb_guard, "nb_thr": a.nb_thr}


# =============================================================================
# CORE (no device, no GUI): Welch PSD, band power, waterfall storage.
# Preallocated arrays and vectorized numpy only; written so that it can be ported to C++.
# =============================================================================

def make_window(name, n):
    """Periodic window (float32); None for rect: no multiplication is performed at all."""
    if name == "rect":
        return None
    fn = {"hann": np.hanning, "hamming": np.hamming, "blackman": np.blackman}[name]
    return fn(n + 1)[:-1].astype(np.float32)


class Welch:
    """Welch periodogram averaging; all work buffers are allocated once."""

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

    def reset(self):
        self.acc.fill(0.0)
        self.n_acc = 0

    def accumulate(self, iq, hop, n_win):
        """iq: int16 array (n_samples, 2) = I, Q. Adds the power spectra of n_win windows."""
        nfft = self.nfft
        view = np.lib.stride_tricks.sliding_window_view(iq, nfft, axis=0)[::hop][:n_win]  # (win, 2, nfft)
        if view.shape[0] < n_win:
            raise ValueError(f"block too short: {view.shape[0]} windows available, {n_win} needed")
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

    def power_into(self, out):
        """Power per bin, fftshifted (bin nfft//2 = center), relative to full scale: a complex tone of
        amplitude FS_FULL_SCALE counts gives 1.0 (0 dBFS) in its bin. No allocation."""
        scale = 1.0 / (self.n_acc * self.nfft * self.sum_w2 * FS_FULL_SCALE ** 2)
        s = self.nfft // 2
        np.multiply(self.acc[self.nfft - s:], scale, out=out[:s])
        np.multiply(self.acc[:self.nfft - s], scale, out=out[s:])


class BandIntegrator:
    """Power in [f - bw/2, f + bw/2] for fixed channels from a per-bin power spectrum.
    Band edges are fractional-bin accurate; everything is precomputed, calls do not allocate."""

    def __init__(self, nbins, fs, fc, f_ch, bw_ch):
        n = nbins
        df = fs / n
        k0 = n // 2
        u_lo = np.clip((f_ch - bw_ch / 2 - fc) / df + k0 + 0.5, 0, n)
        u_hi = np.clip((f_ch + bw_ch / 2 - fc) / df + k0 + 0.5, 0, n)
        self.i_lo = np.minimum(np.floor(u_lo).astype(np.intp), n - 1)
        self.i_hi = np.minimum(np.floor(u_hi).astype(np.intp), n - 1)
        self.w_lo = u_lo - self.i_lo
        self.w_hi = u_hi - self.i_hi
        self.i_lo1, self.i_hi1 = self.i_lo + 1, self.i_hi + 1
        m = len(f_ch)
        self.cs = np.zeros(n + 1)
        self.t1, self.t2, self.t3, self.t4 = (np.empty(m) for _ in range(4))
        self.n = n

    def db(self, p_bin):
        """Band power in dBFS for every channel (returns a shared buffer, copy it to keep)."""
        cs, t1, t2, t3, t4 = self.cs, self.t1, self.t2, self.t3, self.t4
        np.cumsum(p_bin, out=cs[1:])
        np.take(cs, self.i_hi1, out=t1, mode="clip")
        np.take(cs, self.i_hi, out=t2, mode="clip")
        t1 -= t2
        t1 *= self.w_hi
        t1 += t2
        np.take(cs, self.i_lo1, out=t3, mode="clip")
        np.take(cs, self.i_lo, out=t4, mode="clip")
        t3 -= t4
        t3 *= self.w_lo
        t3 += t4
        t1 -= t3
        np.maximum(t1, 1e-30, out=t1)
        np.log10(t1, out=t1)
        t1 *= 10.0
        return t1


class WaterfallBuffer:
    """Rows of dB values, newest row first, with the index of the first sample of every row. A buffer of twice
    the depth makes a push O(1) amortized and the displayed array a zero-copy view."""

    def __init__(self, rows, nbins):
        self.rows, self.nbins = rows, nbins
        self.buf = np.full((2 * rows, nbins), np.nan, np.float32)
        self.sbuf = np.zeros(2 * rows, np.int64)
        self.ptr = rows
        self.count = 0
        self.tmp = np.empty((rows, nbins), np.float32)

    def push(self, row_db, start=0):
        r = self.rows
        if self.ptr == 0:
            self.buf[r + 1:2 * r] = self.buf[0:r - 1]
            self.sbuf[r + 1:2 * r] = self.sbuf[0:r - 1]
            self.ptr = r + 1
        self.ptr -= 1
        self.buf[self.ptr] = row_db
        self.sbuf[self.ptr] = start
        self.count = min(self.count + 1, r)

    def view(self):
        return self.buf[self.ptr:self.ptr + self.rows]

    def starts(self):
        """First-sample index of every row of view(), newest first."""
        return self.sbuf[self.ptr:self.ptr + self.rows]

    def mean_linear_into(self, out):
        """Mean over the valid rows of 10^(dB/10) (power density), no allocation."""
        c = self.count
        if c == 0:
            out.fill(0.0)
            return
        t = self.tmp[:c]
        np.multiply(self.view()[:c], LN10 / 10.0, out=t)
        np.exp(t, out=t)
        np.add.reduce(t, axis=0, dtype=np.float64, out=out)
        out *= 1.0 / c

    def resized(self, rows):
        nb = WaterfallBuffer(rows, self.nbins)
        keep = min(rows, self.count)
        if keep:
            nb.buf[rows - keep:rows] = self.view()[:keep]        # newest first: the newest row stays at the top
            nb.sbuf[rows - keep:rows] = self.starts()[:keep]
            nb.ptr = rows - keep
            nb.count = keep
        return nb


class NbDetector:
    """Narrowband signals in one PSD. Per-bin work buffers are allocated once; only the (small) per-signal
    arrays are created per call.
    Bands are bins above the local noise floor plus the threshold; the floor is the median of blocks of bins,
    interpolated between block centres. Edges are the bins within NB_EDGE_DB of the band peak. A band that the
    window leakage of a stronger signal explains (the sidelobes of the window used) is dropped."""

    def __init__(self, n, fs, fc, window="rect"):
        self.n, self.df, self.fc = n, fs / n, fc
        # leakage of the window: highest response at every whole-bin distance, relative to the peak
        w = make_window(window, n)
        w = np.ones(n) if w is None else w.astype(np.float64)
        spec = np.abs(np.fft.rfft(w, 4 * n))
        spec /= spec[0]
        kmax = min(n // 2, 16384) + 1
        around = 4 * np.arange(1, kmax)[:, None] + np.arange(-2, 3)[None, :]
        env = np.ones(kmax)
        env[1:] = spec[np.minimum(around, len(spec) - 1)].max(axis=1)
        self.leak_cum = np.concatenate(([0.0], np.cumsum(env ** 2)))   # sum of the leakage power over distances < k
        self.k = np.arange(n, dtype=np.float64)
        self.db_pad = np.empty(n + 1)        # PSD in dB plus a sentinel: reduceat needs indices < length
        self.floor = np.full(n, np.nan)         # noise floor of the last PSD, dB per bin
        self.thr_db = np.full(n, np.nan)        # floor + threshold: what a bin must exceed
        self.t1, self.t2 = np.empty(n), np.empty(n)
        self.mask = np.empty(n, bool)
        self.m8 = np.zeros(n + 2, np.int8)
        self.cs, self.cs2 = np.zeros(n + 1), np.zeros(n + 1)
        self.set_params(DEFAULT_NB_MAX_WIDTH, DEFAULT_NB_GUARD, DEFAULT_NB_THR)

    def set_params(self, max_width, guard, thr_db):
        self.max_width, self.guard, self.thr = max_width, guard, thr_db
        self.possible = self.df <= max_width          # a single bin must fit into the maximum width
        w_bins = max(1, math.ceil(max_width / self.df))
        blk = max(NB_FLOOR_MIN_BINS, NB_FLOOR_WIDTHS * w_bins)
        self.nblk = max(1, self.n // blk)
        self.blk = self.n // self.nblk
        self.body = self.nblk * self.blk
        pos = (self.k + 0.5) / self.blk - 0.5            # position in block-centre units
        if self.nblk >= 2:                                # linear between centres, linearly extended at both ends
            self.i0 = np.clip(np.floor(pos), 0, self.nblk - 2).astype(np.intp)
            self.i1 = self.i0 + 1
            self.wgt = pos - self.i0
        else:
            self.i0 = self.i1 = np.zeros(self.n, np.intp)
            self.wgt = np.zeros(self.n)

    def detect(self, p_bin, db):
        """p_bin: power per bin (linear), db: PSD per bin (dB), both fftshifted, length n.
        Returns (center Hz, width Hz, power dBFS, SNR dB), one entry per detected signal."""
        n, df = self.n, self.df
        none = (np.empty(0), np.empty(0), np.empty(0), np.empty(0))
        if not self.possible:
            return none
        med = np.median(db[:self.body].reshape(self.nblk, self.blk), axis=1)
        np.take(med, self.i0, out=self.t1)
        np.take(med, self.i1, out=self.t2)
        np.subtract(self.t2, self.t1, out=self.t2)
        self.t2 *= self.wgt
        np.add(self.t1, self.t2, out=self.floor)
        np.add(self.floor, self.thr, out=self.thr_db)
        np.greater(db, self.thr_db, out=self.mask)
        m = self.mask
        if n > 2:                                       # fill single-bin gaps inside a signal
            m[1:-1] |= m[:-2] & m[2:] & ~m[1:-1]
        self.db_pad[:n] = db
        self.db_pad[n] = -np.inf
        firsts, lasts, peaks = [], [], []
        for _ in range(NB_PASSES):
            self.m8[1:n + 1] = m
            d = np.diff(self.m8)
            s = np.flatnonzero(d == 1)
            e = np.flatnonzero(d == -1)
            runs = len(s)
            if runs == 0:
                break
            nrun = e - s
            idx = np.empty(2 * runs, np.intp)
            idx[0::2] = s
            idx[1::2] = e
            peak = np.maximum.reduceat(self.db_pad, idx)[0::2]
            # trim every run to the bins within NB_EDGE_DB of its peak (leakage skirts of strong signals)
            run_bins = np.flatnonzero(m)
            rep = np.repeat(np.arange(runs), nrun)
            keep = db[run_bins] >= peak[rep] - NB_EDGE_DB
            rb = run_bins[keep]
            # split the kept bins into separate bands: a gap wider than NB_SEG_GAP bins starts a new one, so two
            # signals joined only by their leakage skirts are not taken for one wide signal
            brk = np.flatnonzero(np.diff(rb) > NB_SEG_GAP + 1)
            first = rb[np.concatenate(([0], brk + 1))]
            last = rb[np.concatenate((brk, [len(rb) - 1]))]
            idx2 = np.empty(2 * len(first), np.intp)
            idx2[0::2] = first
            idx2[1::2] = last + 1
            seg_peak = np.maximum.reduceat(self.db_pad, idx2)[0::2]
            firsts.append(first)
            lasts.append(last)
            peaks.append(seg_peak)
            # take the found bands (and a bin around them) out of the candidates; the next pass finds what hid in
            # their skirts (the leakage test below drops the skirts themselves)
            cut = np.zeros(n + 1, np.int32)
            np.add.at(cut, np.maximum(first - 1, 0), 1)
            np.add.at(cut, np.minimum(last + 2, n), -1)
            m &= ~(np.cumsum(cut[:n]) > 0)
        if not firsts:
            return none
        first, last, peak = np.concatenate(firsts), np.concatenate(lasts), np.concatenate(peaks)
        nb = last - first + 1
        pk = peak
        # band power and centroid from cumulative sums
        np.cumsum(p_bin, out=self.cs[1:])
        np.multiply(p_bin, self.k, out=self.t1)
        np.cumsum(self.t1, out=self.cs2[1:])
        power = np.maximum(self.cs[last + 1] - self.cs[first], 1e-30)
        cidx = (self.cs2[last + 1] - self.cs2[first]) / power
        if len(pk) > 1:              # drop bands that the window leakage of a stronger band explains
            if len(pk) > NB_LEAK_MAX_RUNS:
                sel = np.argsort(-pk)[:NB_LEAK_MAX_RUNS]
                first, last, nb, pk, power, cidx = first[sel], last[sel], nb[sel], pk[sel], power[sel], cidx[sel]
            # per-bin power expected at every band's centre from the leakage of the stronger bands (each spread
            # evenly over its bins) plus the noise floor; a band that does not stand out of it is leakage
            K = len(self.leak_cum) - 1
            ci = np.rint(cidx).astype(np.intp)
            dn = np.maximum(np.maximum(first[None, :] - ci[:, None], ci[:, None] - last[None, :]), 0)
            lk = (power[None, :] / nb[None, :]) * (self.leak_cum[np.minimum(dn + nb[None, :], K)]
                                                   - self.leak_cum[np.minimum(dn, K)])
            lk = np.where(pk[None, :] > pk[:, None], lk, 0.0)
            floor_p = df * 10.0 ** (self.floor[np.clip(ci, 0, n - 1)] / 10.0)
            real = df * 10.0 ** (pk / 10.0) > (lk.sum(axis=1) + floor_p) * 10.0 ** (NB_LEAK_MARGIN_DB / 10.0)
            first, last, nb, pk, power, cidx = first[real], last[real], nb[real], pk[real], power[real], cidx[real]
        ok = nb * df <= self.max_width * (1.0 + 1e-9)       # only now: wide bands also explain leakage
        if not ok.any():
            return none
        first, last, nb, pk, power, cidx = first[ok], last[ok], nb[ok], pk[ok], power[ok], cidx[ok]
        if len(power) > NB_MAX_DET:
            sel = np.argsort(-power)[:NB_MAX_DET]
            first, last, nb, power, cidx = first[sel], last[sel], nb[sel], power[sel], cidx[sel]
        # noise: median per-bin power of a region around the centre, width = band width x (1 + guard)
        hw = 0.5 * nb * (1.0 + self.guard)
        lo = np.maximum(np.rint(cidx - hw).astype(np.intp), 0)
        hi = np.minimum(np.rint(cidx + hw).astype(np.intp), n - 1)
        cnt = hi - lo + 1
        gidx = lo[:, None] + np.arange(int(cnt.max()), dtype=np.intp)[None, :]
        vals = np.where(gidx <= hi[:, None], p_bin[np.minimum(gidx, n - 1)], np.inf)
        vals.sort(axis=1)
        r = np.arange(len(cnt))
        med_p = 0.5 * (vals[r, (cnt - 1) // 2] + vals[r, cnt // 2])
        noise = np.maximum(med_p * nb, 1e-30)
        return (self.fc + (cidx - n // 2) * df, nb * df, 10.0 * np.log10(power), 10.0 * np.log10(power / noise))


class RawLog:
    """Raw detections of every PSD window, newest overwrite oldest. Time is counted from samples."""

    def __init__(self, cap=RAW_LOG_CAP):
        self.cap, self.n = cap, 0                          # n: detections ever appended
        self.t, self.center, self.width = np.zeros(cap), np.zeros(cap), np.zeros(cap)
        self.power, self.snr = np.zeros(cap), np.zeros(cap)
        self.ident = np.full(cap, -1, np.int64)

    def append(self, t, center, width, power, snr, ident):
        m = len(center)
        if m == 0:
            return
        pos = (self.n + np.arange(m)) % self.cap
        self.t[pos], self.center[pos], self.width[pos] = t, center, width
        self.power[pos], self.snr[pos], self.ident[pos] = power, snr, ident
        self.n += m

    def count_since(self, t0):
        m = min(self.n, self.cap)
        return int(np.count_nonzero(self.t[:m] >= t0)) if m else 0


class NbTracker:
    """Registry of narrowband signals. Raw detections of every PSD are collected per signal (the last REG_HISTORY of
    them) and the signal's parameters are their MEDIANS, so one noisy detection moves nothing and the centre does not
    wander. A detection joins the signal whose aggregated centre is nearest, within half the larger width plus two
    bins; otherwise it starts a new signal, which is listed only after REG_MIN_HITS detections (`confirmed`).
    Signals may vanish and come back: present = found in the latest PSD, last_seen = time of the last window."""

    def __init__(self, cap=NB_TRACK_LIMIT):
        self.cap = cap
        self.used = np.zeros(cap, bool)
        self.present = np.zeros(cap, bool)
        self.confirmed = np.zeros(cap, bool)
        self.ident = np.zeros(cap, np.int64)
        self.hits = np.zeros(cap, np.int64)
        self.first_win = np.zeros(cap, np.int64)
        self.center, self.width = np.zeros(cap), np.zeros(cap)          # aggregated (medians of the history)
        self.power, self.snr = np.zeros(cap), np.zeros(cap)
        self.first_seen, self.last_seen = np.zeros(cap), np.zeros(cap)
        self.hist = np.full((4, cap, REG_HISTORY), np.nan)              # centre, width, power, snr of the raw detections
        self.hpos = np.zeros(cap, np.int64)
        self.windows = 0                                                # PSD windows seen
        self.next_id = 1

    def clear(self):
        self.used[:] = False
        self.present[:] = False
        self.confirmed[:] = False
        self.windows = 0
        self.next_id = 1

    @property
    def visible(self):
        return self.used & self.confirmed

    def duty(self, k):
        """Fraction of the windows since the signal first appeared in which it was found."""
        return self.hits[k] / max(1, self.windows - self.first_win[k] + 1)

    def _slot(self):
        free = np.flatnonzero(~self.used)
        if len(free):
            return int(free[0])
        absent = np.flatnonzero(self.used & ~self.present)
        if not len(absent):
            return -1
        # an unconfirmed (maybe a false alarm) signal that was seen longest ago goes first
        return int(absent[np.argmin(self.last_seen[absent] + 1e12 * self.confirmed[absent])])

    def update(self, center, width, power, snr, t, df):
        """One PSD: detections (arrays) found in the window that starts at time t (s). Returns the signal slot of every
        detection (-1: not stored)."""
        self.present[:] = False
        self.windows += 1
        slots = np.full(len(center), -1, np.intp)
        taken = np.zeros(self.cap, bool)
        touched = []
        for j in np.argsort(-snr):                       # strongest first
            c, w = center[j], width[j]
            k = -1
            cand = self.used & ~taken
            if cand.any():
                dist = np.abs(self.center - c)
                near = cand & (dist <= 0.5 * np.maximum(self.width, w) + 2.0 * df)
                if near.any():
                    k = int(np.argmin(np.where(near, dist, np.inf)))
            if k < 0:
                k = self._slot()
                if k < 0:
                    continue
                self.used[k] = True
                self.confirmed[k] = False
                self.ident[k] = self.next_id
                self.next_id += 1
                self.first_seen[k] = t
                self.first_win[k] = self.windows
                self.hits[k] = 0
                self.hist[:, k, :] = np.nan
                self.hpos[k] = 0
                self.center[k], self.width[k] = c, w
            taken[k] = True
            self.present[k] = True
            pos = int(self.hpos[k] % REG_HISTORY)
            self.hist[0, k, pos], self.hist[1, k, pos] = c, w
            self.hist[2, k, pos], self.hist[3, k, pos] = power[j], snr[j]
            self.hpos[k] += 1
            self.last_seen[k] = t
            self.hits[k] += 1
            if self.hits[k] >= REG_MIN_HITS:
                self.confirmed[k] = True
            slots[j] = k
            touched.append(k)
        if touched:
            ks = np.array(touched)
            med = np.nanmedian(self.hist[:, ks, :], axis=2)                # aggregate: medians over the history
            self.center[ks], self.width[ks], self.power[ks], self.snr[ks] = med
        return slots


# ------------------------------------------------------------------ band extraction (digital down-converter)

@dataclass
class RatePlan:
    """Output sample rate = input rate x L / M, with whole numbers L and M: the requested rate is a guide, the nearest
    rate that such a ratio gives is used. Block sizes of the FFT-based extractor follow from it."""
    fs: float
    target: float
    L: int
    M: int
    out_rate: float
    error_ppm: float
    k: int                       # N = M k input samples per FFT, N_out = L k output samples per inverse FFT
    N: int
    N_out: int
    V_half: int                  # output samples cut off at each end of a block
    hop_in: int                  # input samples between two blocks
    hop_out: int                 # valid output samples per block


def _smooth7(n):
    for p in (2, 3, 5, 7):
        while n % p == 0:
            n //= p
    return n == 1


def plan_resampler(fs, target, max_den=CAP_RATIO_MAX_DEN):
    """Nearest L/M (M <= max_den) to target/fs and the block geometry of the extractor. ValueError if impossible."""
    if not 0 < target <= 0.9 * fs:
        raise ValueError(f"The output rate {fmt_hz_short(target)} must be above 0 and at most 90 % of the capture rate "
                         f"{fmt_hz_short(fs)}")
    ratio = (Fraction(float(target)) / Fraction(float(fs))).limit_denominator(max_den)
    if ratio.numerator == 0:
        raise ValueError(f"The output rate {fmt_hz_short(target)} is too low for the capture rate {fmt_hz_short(fs)} "
                         f"(the ratio would need a denominator above {max_den})")
    L, M = ratio.numerator, ratio.denominator
    u = -(-CAP_OVERLAP // L)                                    # V_half = L u >= CAP_OVERLAP output samples
    kmin = 16 * u                                               # a block holds at least 8 overlaps
    k = next((c for c in range(kmin, int(kmin * 1.5) + 1) if _smooth7(c * M)), kmin)
    if M * k > CAP_N_MAX:
        raise ValueError(f"The rate ratio {L}/{M} needs an FFT of {M * k} points (limit {CAP_N_MAX}); pick an output rate "
                         f"that divides the capture rate better")
    out_rate = fs * L / M
    return RatePlan(fs, target, L, M, out_rate, (out_rate / target - 1.0) * 1e6, k, M * k, L * k, L * u,
                    M * (k - 2 * u), L * (k - 2 * u))


class Channelizer:
    """One signal band out of a wide IQ stream: mixing, filtering, decimation and the rate change by L/M in one pass
    (FFT overlap-save: the bins around the centre of every block are cut out and inverse-transformed at the output
    size). Output samples carry the phase of an ideal continuous mixer, so they join across blocks and across gaps:
    the n-th output of a run is the signal at the input time (n_first + n M / L) / fs, a pure tone of the signal
    keeps its phase. Preallocated; one instance per band."""

    def __init__(self, rp, fs, fc, f_c):
        self.rp, self.fs = rp, fs
        N, No = rp.N, rp.N_out
        self.k0 = int(round((f_c - fc) / fs * N))               # bin of the centre (signed)
        self.delta = f_c - (fc + self.k0 * fs / N)              # what is left to the exact centre, Hz (< half a bin)
        idx = np.rint(np.fft.fftfreq(No, d=1.0 / No)).astype(np.int64)
        self.kk = (self.k0 + idx) % N
        fr = np.minimum(np.abs(idx) / (No / 2.0), 1.0)
        self.H = np.where(fr <= CAP_PASS_FRAC, 1.0,
                          0.5 * (1.0 + np.cos(np.pi * (fr - CAP_PASS_FRAC) / (1.0 - CAP_PASS_FRAC)))).astype(np.float32)
        self.scale = np.float32(No / N)                         # a tone keeps its amplitude
        self.work = np.empty(N, np.complex64)
        self.spec = np.empty(N, np.complex64)
        self.Y = np.empty(No, np.complex64)
        self.y = np.empty(No, np.complex64)
        self.out = np.empty(rp.hop_out, np.complex64)
        self.t_rel = np.arange(rp.V_half, No - rp.V_half, dtype=np.float64) * (rp.M / rp.L) / fs

    def first_valid(self, n0):
        """Input sample index that the first valid output of a block starting at input index n0 belongs to."""
        return n0 + self.rp.M * (self.rp.V_half // self.rp.L)

    def process(self, iq, n0):
        """iq: int16 (N, 2) = the input samples n0 ... n0+N-1 of the stream. Returns the hop_out valid output samples."""
        rp = self.rp
        np.copyto(self.work.real, iq[:, 0], casting="unsafe")
        np.copyto(self.work.imag, iq[:, 1], casting="unsafe")
        try:
            np.fft.fft(self.work, out=self.spec)
        except TypeError:                                       # numpy < 2.0
            self.spec[...] = np.fft.fft(self.work)
        np.take(self.spec, self.kk, out=self.Y)
        self.Y *= self.H
        try:
            np.fft.ifft(self.Y, out=self.y)
        except TypeError:
            self.y[...] = np.fft.ifft(self.Y)
        phase = (-2.0 * np.pi * (((self.k0 * int(n0)) % rp.N) / rp.N)
                 - 2.0 * np.pi * self.delta * (int(n0) / self.fs + self.t_rel))
        np.multiply(self.y[rp.V_half:rp.N_out - rp.V_half], np.exp(1j * phase).astype(np.complex64), out=self.out)
        self.out *= self.scale
        return self.out


def to_int16(x, dst):
    """Complex samples to interleaved int16 I/Q (counts as they are, rounded, saturated) in dst (n, 2)."""
    for col, part in ((0, x.real), (1, x.imag)):
        dst[:, col] = np.clip(np.rint(part), -32768, 32767)       # clip BEFORE the cast: a cast would wrap around
    return dst


def wav_header(rate_hz, data_bytes):
    """The 44 bytes of a canonical PCM WAV header: 2 channels (I, Q), 16 bits, for data_bytes bytes of data.
    RIFF size = 36 + data size; the sample rate is a whole number of Hz."""
    rate = int(round(rate_hz))
    return (b"RIFF" + struct.pack("<I", 36 + data_bytes) + b"WAVE" + b"fmt "
            + struct.pack("<IHHIIHH", 16, 1, 2, rate, rate * 4, 4, 16) + b"data" + struct.pack("<I", data_bytes))


def write_json_atomic(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)
    os.replace(tmp, path)


def classify_channels(ch_f, ch_bw, fc, fs):
    """Indices of channels fully inside the capture band, partly inside and outside it."""
    lo, hi = fc - fs / 2, fc + fs / 2
    left, right = ch_f - ch_bw / 2, ch_f + ch_bw / 2
    inside = (left >= lo) & (right <= hi)
    outside = (right <= lo) | (left >= hi)
    partial = ~inside & ~outside
    return np.flatnonzero(inside), int(partial.sum()), int(outside.sum())


# =============================================================================
# DEVICE: streaming capture and processing threads
# =============================================================================

class RingStatus(IntEnum):
    """Why a range of the ring can or cannot be read. `late` tells the two cases of "we did not keep up"."""
    OK = 0                   # the samples are there
    OVERWRITTEN = 1          # too late: the writer has already overwritten the start of the range
    LOST_WHILE_COPYING = 2   # too late: the range was there when the copy began and was overwritten before it ended
    NOT_YET = 3              # too early: the end of the range has not been written yet (wait for it; nothing is lost)
    INVALID = 4              # the request cannot be satisfied: start < 0, n <= 0, or n larger than the ring can hold

    @property
    def late(self):
        return self in (RingStatus.OVERWRITTEN, RingStatus.LOST_WHILE_COPYING)


class SampleRing:
    """Raw IQ samples of the last few seconds in a ring buffer.

    Sample index i (counted from the first sample after the start-up discard; the stream sample of the device is i plus the
    discarded ones) lives at buf[i % size]. The size is
    a whole number of rows, so a row is always one contiguous slice. One thread writes (`advance` after filling
    `buf[pos:pos + n]` at `pos = head % size`), any thread reads. The writer works up to `guard` samples ahead of
    `head`, so data is intact only while head + guard <= start + size: check it AFTER a read."""

    def __init__(self, row_samples, rows, chunk=None):
        self.row_samples, self.rows = row_samples, rows
        self.size = row_samples * rows
        self.chunk = chunk or SYNC_CHUNK     # the writer fills at most this many samples at a time
        self.guard = self.chunk
        self.buf = np.empty((self.size, 2), np.int16)
        self.head = 0                      # samples written since the start (never decreases)
        self.cond = threading.Condition()

    def write_span(self):
        """(position, count) of the next piece the writer may fill without crossing the end of the buffer."""
        pos = self.head % self.size
        return pos, min(self.chunk, self.size - pos)

    def advance(self, n):
        with self.cond:
            self.head += n
            self.cond.notify_all()

    def wait_for(self, end, stop, timeout=0.2):
        """Wait until samples up to index `end` are written; False on timeout or stop."""
        with self.cond:
            self.cond.wait_for(lambda: self.head >= end or stop.is_set(), timeout)
            return self.head >= end

    def check(self, start, n=0):
        """RingStatus of the range start..start+n, without copying (n = 0: only the first sample matters)."""
        head = self.head
        if start < 0 or n < 0 or n > self.size - self.guard:
            return RingStatus.INVALID          # can never be intact, however fast we are
        if start + n > head:
            return RingStatus.NOT_YET
        if head + self.guard > start + self.size:
            return RingStatus.OVERWRITTEN
        return RingStatus.OK

    def intact(self, start, n=0):
        return self.check(start, n) == RingStatus.OK

    def row_view(self, k):
        """Zero-copy view of row k (samples k*row_samples ...). Check intact() after using it."""
        a = (k % self.rows) * self.row_samples
        return self.buf[a:a + self.row_samples]

    def copy_range(self, start, n, out=None):
        """Copy samples start..start+n (any range inside the buffer, wrap handled) into out (n, 2).
        Returns (status, out); out is None unless status is RingStatus.OK. Never waits: for a range that is still
        being recorded (NOT_YET) wait with wait_for(start + n, stop) and ask again."""
        if n <= 0:
            return RingStatus.INVALID, None
        st = self.check(start, n)
        if st != RingStatus.OK:
            return st, None
        if out is None:
            out = np.empty((n, 2), np.int16)
        a = start % self.size
        first = min(n, self.size - a)
        out[:first] = self.buf[a:a + first]
        if first < n:
            out[first:] = self.buf[:n - first]
        if self.check(start, n) != RingStatus.OK:     # the writer lapped us while we copied: the copy is torn
            return RingStatus.LOST_WHILE_COPYING, None
        return RingStatus.OK, out


class Radio:
    """All bladeRF calls of one session; used from the capture thread only."""

    def __init__(self, cfg, plan, stop):
        self.cfg, self.plan, self.stop = cfg, plan, stop
        self.dev = None
        self.ch = None
        self.ch_id = bladerf.CHANNEL_RX(0)
        self.enabled = False
        self.errors = 0                  # read errors that were retried (the stream may have lost samples)
        self.last_error = ""

    def open(self):
        cfg, plan = self.cfg, self.plan
        self.dev = bladerf.BladeRF()
        ch = self.ch = self.dev.Channel(self.ch_id)
        ch.gain_mode = _bladerf.GainMode.Manual
        ch.gain = cfg.gain
        fs = self.dev.set_sample_rate(self.ch_id, cfg.rate)
        bw = self.dev.set_bandwidth(self.ch_id, plan.bandwidth)
        ch.frequency = int(round(plan.center))
        fc = ch.frequency
        self.dev.sync_config(layout=_bladerf.ChannelLayout.RX_X1, fmt=_bladerf.Format.SC16_Q11,
                             num_buffers=SYNC_BUFFERS, buffer_size=SYNC_BUF_SAMPLES,
                             num_transfers=SYNC_TRANSFERS, stream_timeout=SYNC_TIMEOUT_MS)
        self.dev.enable_module(self.ch_id, True)
        self.t_enable = time.time()           # wall clock of the moment the wideband stream was started
        self.enabled = True
        try:                                  # honest gain: read after the stream has started
            gain = ch.gain
        except Exception as e:
            info(f"WARNING: gain cannot be read: {e}")
            gain = None
        return {"fs": float(fs), "bw": float(bw), "fc": float(fc), "gain": gain, "t_enable": self.t_enable}

    def discard(self, scratch, fs):
        """Throw the start-up samples away; returns how many (they are part of the stream's time)."""
        left = max(DISCARD_MIN, int(DISCARD_S * fs))
        done = 0
        while left > 0 and not self.stop.is_set():
            k = min(left, scratch.shape[0], SYNC_CHUNK)
            self.dev.sync_rx(scratch[:k], k, SYNC_TIMEOUT_MS)
            left -= k
            done += k
        return done

    def read(self, buf):
        """Fill buf (n, 2) int16. Returns False if stopped. Raises after repeated read errors."""
        n, pos, errors = buf.shape[0], 0, 0
        while pos < n:
            if self.stop.is_set():
                return False
            k = min(n - pos, SYNC_CHUNK)
            try:
                self.dev.sync_rx(buf[pos:pos + k], k, SYNC_TIMEOUT_MS)
            except Exception as e:
                errors += 1
                self.errors += 1
                self.last_error = f"{type(e).__name__}: {e}"
                if errors >= MAX_CONSECUTIVE_READ_ERRORS:
                    raise RuntimeError(f"stream read failed {errors} times in a row: {e}") from None
                continue
            errors = 0
            pos += k
        return True

    def close(self):
        if self.dev is None:
            return
        if self.enabled:
            try:
                self.dev.enable_module(self.ch_id, False)
            except Exception:
                pass
        try:
            self.dev.close()
        finally:
            self.dev = None


class _CapStats:
    """Capture-thread numbers of one interval. Only the capture thread touches them."""

    def __init__(self, now):
        self.prev_end = now
        self.begin(now)

    def begin(self, now):
        self.t_start, self.samples, self.calls = now, 0, 0
        self.in_read = self.read_max = self.gap_max = 0.0
        self.cpu_start = time.thread_time()

    def call(self, t0, t1, n):
        gap = t0 - self.prev_end                  # time spent outside the read: the USB buffers must bridge it
        if gap > self.gap_max:
            self.gap_max = gap
        d = t1 - t0
        self.in_read += d
        if d > self.read_max:
            self.read_max = d
        self.samples += n
        self.calls += 1
        self.prev_end = t1

    def due(self, now):
        return now - self.t_start >= PERF_PERIOD_S

    def snapshot(self, now, errors, last_error, head):
        snap = {"dt": now - self.t_start, "samples": self.samples, "calls": self.calls, "in_read": self.in_read,
                "read_max": self.read_max, "gap_max": self.gap_max, "errors": errors, "last_error": last_error,
                "cpu": time.thread_time() - self.cpu_start, "head": head}
        self.begin(now)
        return snap


class _ProcStats:
    """Processing-thread numbers of one interval. Only the processing thread touches them."""

    def __init__(self, now, row_dur):
        self.row_dur = row_dur
        self.begin(now)

    def begin(self, now):
        self.t_start, self.rows = now, 0
        self.compute = self.compute_max = 0.0
        self.rt_min = float("inf")
        self.lag_max = 0
        self.cpu_start = time.thread_time()

    def row(self, dur):
        self.rows += 1
        self.compute += dur
        if dur > self.compute_max:
            self.compute_max = dur
        rt = self.row_dur / dur if dur > 0 else float("inf")
        if rt < self.rt_min:
            self.rt_min = rt

    def lag(self, samples):
        if samples > self.lag_max:
            self.lag_max = samples

    def due(self, now):
        return now - self.t_start >= PERF_PERIOD_S

    def snapshot(self, now, backlog):
        snap = {"dt": now - self.t_start, "rows": self.rows, "compute": self.compute, "compute_max": self.compute_max,
                "rt_min": self.rt_min, "lag_max": self.lag_max, "backlog": backlog,
                "cpu": time.thread_time() - self.cpu_start}
        self.begin(now)
        return snap


def safe_name(text):
    out = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text.strip())
    return out or "band"


_FIELD_RE = re.compile(r"\{\{|\}\}|\{([a-z_]+)(?::([^{}]*))?\}|[{}]")
PATTERN_FIELDS = ("name", "time", "utc", "freq", "rate", "id", "n")


def _pattern_tokens(pattern):
    """Yields ('text', str), ('field', key, spec) in order; ValueError for a stray brace."""
    pos = 0
    for m in _FIELD_RE.finditer(pattern):
        if m.start() > pos:
            yield ("text", pattern[pos:m.start()])
        pos = m.end()
        tok = m.group(0)
        if tok == "{{":
            yield ("text", "{")
        elif tok == "}}":
            yield ("text", "}")
        elif m.group(1) is None:
            raise ValueError("A single { or } in the pattern: write {{ or }} for a brace")
        else:
            yield ("field", m.group(1), m.group(2))
    if pos < len(pattern):
        yield ("text", pattern[pos:])


def _field_text(key, spec, ctx):
    """Value of one field of the file name pattern; ValueError says what is wrong."""
    try:
        if key == "name":
            if spec:
                raise ValueError("{name} takes no format")
            return safe_name(ctx["name"])
        if key in ("time", "utc"):
            fmt = spec or "%Y%m%d_%H%M%S"
            tm = time.localtime(ctx["t"]) if key == "time" else time.gmtime(ctx["t"])
            out = time.strftime(fmt, tm)
            if not out:
                raise ValueError("the time format gives nothing")
            return out
        if key == "freq":                               # centre in MHz
            return format(ctx["center_hz"] / 1e6, spec or ".4f")
        if key == "rate":                               # output rate in kHz
            return format(ctx["rate_hz"] / 1e3, spec or "g")
        if key == "id":                                 # signal id; empty when the band is not a registry signal
            return "" if ctx.get("ident") is None else format(int(ctx["ident"]), spec or "d")
        if key == "n":                                  # sequence number, counts up while the name is taken
            return format(int(ctx.get("n", 1)), spec or "d")
    except (ValueError, TypeError, OSError) as e:
        raise ValueError(f"Bad format in {{{key}{':' + spec if spec else ''}}}: {e}") from None
    raise ValueError(f"Unknown field {{{key}}}; the fields are " + ", ".join("{" + f + "}" for f in PATTERN_FIELDS))


_WIN_RESERVED = {"CON", "PRN", "AUX", "NUL"} | {f"COM{i}" for i in range(1, 10)} | {f"LPT{i}" for i in range(1, 10)}


def clean_stem(text):
    """A file name stem that every file system takes: no path separators and none of <>:"|?*, no trailing dot or space."""
    out = "".join("_" if (ch in '<>:"/\\|?*' or ord(ch) < 32) else ch for ch in text).strip().rstrip(". ")
    if not out:
        raise ValueError("The file name pattern gives an empty name")
    out = out[:150].rstrip(". ")
    if out.split(".")[0].rstrip(" ").upper() in _WIN_RESERVED:       # CON, NUL, COM1 ...: devices on Windows, with any extension
        out += "_"
    return out


def render_stem(pattern, ctx):
    """File name stem from the pattern. ctx: name, ident (or None), center_hz, rate_hz, t (epoch seconds), n."""
    if not pattern.strip():
        raise ValueError("The file name pattern is empty")
    out = []
    for tok in _pattern_tokens(pattern):
        out.append(tok[1] if tok[0] == "text" else _field_text(tok[1], tok[2], ctx))
    return clean_stem("".join(out))


def stem_candidates(pattern, ctx, limit=10000):
    """Stems to try, best first. A pattern with {n} counts n up; any other pattern gets _2, _3 ... after the first."""
    uses_n = any(t[0] == "field" and t[1] == "n" for t in _pattern_tokens(pattern))
    for k in range(1, limit + 1):
        stem = render_stem(pattern, dict(ctx, n=k))
        yield stem if (uses_n or k == 1) else f"{stem}_{k}"


def free_stem(directory, pattern, ctx, exists=os.path.exists):
    """The first stem of the pattern for which neither the .wav nor the .json file exists yet."""
    for stem in stem_candidates(pattern, ctx):
        base = os.path.join(directory, stem)
        if not exists(base + ".wav") and not exists(base + ".json"):
            return stem
    raise ValueError("No free file name: change the pattern or the folder")


def iso_utc(t):
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t)) + f".{int((t % 1) * 1000):03d}Z"


class _WavFile:
    """One output file. The 44-byte header is written first with empty sizes; patch_header() writes the real ones."""

    def __init__(self, base, rate_hz, continued_from=None):
        self.base, self.rate_hz, self.continued_from = base, rate_hz, continued_from
        self.fh = open(base + ".wav", "xb")                  # "x": never overwrites a file that is there
        self.fh.write(wav_header(rate_hz, 0))
        self.fh.flush()                                      # the header is on disk from the first moment
        self.data_bytes = 0
        self.frames = 0                                      # I/Q pairs written
        self.segments = []
        self.t_open = time.time()
        self.continued_in = None
        self.closed = False

    @property
    def wav(self):
        return self.base + ".wav"

    @property
    def json(self):
        return self.base + ".json"

    def write(self, frames):
        """frames: int16 (n, 2) = I, Q."""
        self.fh.write(frames.astype("<i2").tobytes())
        self.frames += len(frames)
        self.data_bytes += 4 * len(frames)

    def patch_header(self):
        self.fh.seek(0)
        self.fh.write(wav_header(self.rate_hz, self.data_bytes))
        self.fh.seek(0, 2)
        self.fh.flush()

    def close(self):
        self.patch_header()
        self.fh.close()
        self.closed = True


class BandRecorder:
    """Records one band of the ring buffer into <directory>/<stem>.wav (WAV, int16, I left and Q right) with the sidecar
    <stem>.json that describes the capture. A thread reads blocks of the ring, runs them through the Channelizer and
    appends the result. start / pause / resume / stop. The centre is fixed for the life of the file (another centre means
    another file). A new segment (a line of the sidecar: first frame, stream sample, time in samples) begins after every
    pause and every block that was lost because the recorder was too slow; the phase reference of the output is the same
    across segments. Time is counted in samples from the start of the wideband stream: stream sample = ring sample +
    stream_offset (the start-up samples the device delivered before the ring began). At 4 GiB the file is closed and the next one begins."""

    def __init__(self, ring, fs, fc, wall0, spec, ui_q, context=None, stream_offset=0):
        """spec: dict(center Hz, rate Hz (a guide, see RatePlan), name, directory, pattern, ident).
        context: dict of what the sidecar says about the capture (keys capture, analysis, source)."""
        self.ring, self.fs, self.fc, self.wall0, self.spec, self.ui_q = ring, fs, fc, wall0, dict(spec), ui_q
        self.context = context or {}
        self.stream_offset = int(stream_offset)
        self.plan = plan_resampler(fs, spec["rate"])
        if self.plan.N > ring.size - ring.guard:
            raise ValueError(f"The buffer is too short for this band: one block needs {self.plan.N / fs * 1e3:.0f} ms of "
                             f"signal, the buffer holds {(ring.size - ring.guard) / fs * 1e3:.0f} ms; raise Buffer")
        self.chan = Channelizer(self.plan, fs, fc, spec["center"])
        self.blk = np.empty((self.plan.N, 2), np.int16)
        self.dst = np.empty((self.plan.hop_out, 2), np.int16)
        self.stop_ev, self.paused = threading.Event(), threading.Event()
        self.state = "idle"
        self.error = ""
        self.warning = ""                            # something that did not stop the recording (the sidecar, see _write_meta)
        self.out_samples = 0                         # all files together
        self.bytes = 0
        self.blocks = 0
        self.lost = 0
        self.cur = None                              # the file being written
        self.files_done = []                         # the files that were closed before it (4 GiB limit)
        self.t_work = 0.0
        self.t_start = None
        self.thread = threading.Thread(target=self._run, daemon=True, name="recorder")

    @property
    def base(self):
        return self.cur.base if self.cur else None

    @property
    def segments(self):
        return self.cur.segments if self.cur else []

    @property
    def files(self):
        return (self.cur.wav, self.cur.json) if self.cur else ()

    # ---- control
    def _new_file(self, continued_from=None):
        d = self.spec["directory"]
        os.makedirs(d, exist_ok=True)
        ctx = {"name": self.spec["name"], "ident": self.spec.get("ident"), "center_hz": self.spec["center"],
               "rate_hz": self.plan.out_rate, "t": time.time()}
        for stem in stem_candidates(self.spec.get("pattern", DEFAULT_CAP_PATTERN), ctx):
            base = os.path.join(d, stem)
            if os.path.exists(base + ".json"):
                continue
            try:
                return _WavFile(base, self.plan.out_rate, continued_from)
            except FileExistsError:
                continue
        raise OSError("No free file name: change the pattern or the folder")

    def start(self):
        self.cur = self._new_file()
        self.state = "recording"
        self.t_start = time.perf_counter()
        self._write_meta(self.cur, "recording")
        self.thread.start()

    def pause(self):
        if self.state == "recording":
            self.paused.set()
            self.state = "paused"

    def resume(self):
        if self.state == "paused":
            self.paused.clear()
            self.state = "recording"

    def stop(self, timeout=5.0):
        self.stop_ev.set()
        if self.thread.is_alive():
            self.thread.join(timeout)
        if self.state != "error":
            self.state = "stopped"

    # ---- the sidecar
    def _meta(self, f, state):
        rp = self.plan
        hdr = int(round(rp.out_rate))
        now = time.time()
        meta = {
            "schema": "rssi_waterfall.capture.v1",
            "state": state,                           # recording | paused | closed
            "file": {"wav": os.path.basename(f.wav), "header_bytes": WAV_HEADER_BYTES,
                     "continued_from": f.continued_from, "continued_in": f.continued_in},
            "format": {"container": "WAV (RIFF/WAVE, PCM, canonical 44-byte header)", "sample_format": "int16, little-endian",
                       "channels": 2, "channel_0": "I", "channel_1": "Q", "frame": "one I/Q pair = 4 bytes",
                       "scale": "ADC counts as the radio delivers them (full scale 2048), not normalised",
                       "sample_rate_hz": rp.out_rate, "sample_rate_header_hz": hdr,
                       "sample_rate_header_error_ppm": (hdr / rp.out_rate - 1.0) * 1e6,
                       "read_with": "numpy.fromfile(path, '<i2', offset=44).reshape(-1, 2)"},
            "time": {"started_utc": iso_utc(f.t_open),
                     "started_local": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(f.t_open)),
                     "updated_utc": iso_utc(now),
                     "stream_origin": "t = 0 and stream sample 0: the first sample the device delivers after the wideband "
                                      "stream was started at the capture centre; the start-up samples that are thrown away "
                                      "count (stream sample = ring sample + discarded_samples)",
                     "discarded_samples": self.stream_offset,
                     "segment_time": "stream_time_s = stream_sample / capture sample rate"},
            "band": {"center_hz": float(self.spec["center"]), "name": self.spec["name"],
                     "pattern": self.spec.get("pattern", DEFAULT_CAP_PATTERN),
                     "requested_width_hz": rp.target, "output_rate_hz": rp.out_rate, "rate_ratio_L_M": [rp.L, rp.M],
                     "rate_error_ppm": rp.error_ppm, "passband_fraction": CAP_PASS_FRAC,
                     "extractor": "FFT overlap-save: mixing, filtering and the rate change by L/M in one pass",
                     "block_input_samples": rp.N, "block_output_samples": rp.hop_out},
            "segments": f.segments,
            "result": {"frames": f.frames, "duration_s": f.frames / rp.out_rate, "data_bytes": f.data_bytes,
                       "lost_blocks": self.lost},
        }
        meta.update(self.context)                     # capture, analysis, source
        return meta

    def _write_meta(self, f, state, tries=3):
        """The sidecar is secondary to the recording: if it cannot be replaced (on Windows another program may hold it open)
        it is retried, the next update tries again, and the user is told; the recording goes on."""
        meta = self._meta(f, state)
        err = None
        for _ in range(tries):
            try:
                write_json_atomic(f.json, meta)
                self.warning = ""
                return True
            except OSError as e:
                err = e
                time.sleep(0.05)
        tmp = os.path.basename(f.json) + ".tmp"
        newest = f" The newest version is in {tmp}." if os.path.exists(f.json + ".tmp") else ""
        self.warning = (f"{os.path.basename(f.json)} could not be updated ({type(err).__name__}: {err}); "
                        f"the recording goes on.{newest}")
        return False

    def _segment(self, n_in, why):
        """A new segment of the current file: its first frame is the output sample at ring index n_in."""
        f = self.cur
        stream = int(n_in) + self.stream_offset
        f.segments.append({"frame_start": f.frames, "input_sample": int(n_in), "stream_sample": stream,
                           "stream_time_s": stream / self.fs, "reason": why})
        self._write_meta(f, "recording")

    def _finalize(self, f):
        f.close()                                     # the header gets its real sizes
        self._write_meta(f, "closed", tries=20)

    def _roll(self):
        """The file is full (4 GiB): the next one is opened first (its name goes into the sidecar of the old one), the old
        one is finished, the recording goes on without a gap between the blocks."""
        old = self.cur
        new = self._new_file(continued_from=os.path.basename(old.wav))
        old.continued_in = os.path.basename(new.wav)
        self._finalize(old)
        self.files_done.append(old)
        self.cur = new
        self._write_meta(new, "recording")

    def _stats(self):
        dur = self.out_samples / self.plan.out_rate
        self.ui_q.put(("rec", {"state": self.state, "seconds": dur, "samples": self.out_samples, "bytes": self.bytes,
                               "center": self.spec["center"],
                               "blocks": self.blocks, "lost": self.lost, "error": self.error, "warning": self.warning,
                               "rt": (self.blocks * self.plan.N / self.fs) / self.t_work if self.t_work > 0 else float("nan"),
                               "files": self.files}))

    # ---- the thread
    def _run(self):
        ring, rp = self.ring, self.plan
        n0 = None
        new_segment = None
        t_stats = t_hdr = time.perf_counter()
        try:
            while not self.stop_ev.is_set():
                now = time.perf_counter()
                if now - t_stats >= 0.5:
                    self._stats()
                    t_stats = now
                if now - t_hdr >= HEADER_REFRESH_S:          # a crash then still leaves a file that reads
                    self.cur.patch_header()
                    self._write_meta(self.cur, "recording")
                    t_hdr = now
                if self.paused.is_set():
                    if n0 is not None:
                        self.cur.patch_header()
                        self._write_meta(self.cur, "paused")
                        n0 = None                              # after a pause the next block starts at the then-newest data
                    self.stop_ev.wait(0.05)
                    continue
                if n0 is None:
                    n0 = ring.head
                    new_segment = "start" if (self.out_samples == 0 and not self.cur.segments) else "resume"
                if not ring.wait_for(n0 + rp.N, self.stop_ev):
                    continue
                st, blk = ring.copy_range(n0, rp.N, out=self.blk)
                if st.late:                                    # we were too slow: this block is gone
                    self.lost += 1
                    n0 = ring.head - rp.N
                    new_segment = "lost_block"
                    continue
                if st != RingStatus.OK:
                    raise RuntimeError(f"unexpected ring status {st.name} for the block at {n0}")
                t0 = time.perf_counter()
                out = self.chan.process(blk, n0)
                to_int16(out, self.dst)
                self.t_work += time.perf_counter() - t0
                if self.cur.data_bytes + 4 * rp.hop_out > WAV_MAX_DATA:
                    self._roll()
                    new_segment = new_segment or "continued"
                if new_segment:
                    self._segment(self.chan.first_valid(n0), new_segment)
                    new_segment = None
                self.cur.write(self.dst)
                self.out_samples += rp.hop_out
                self.bytes += rp.hop_out * 4
                self.blocks += 1
                n0 += rp.hop_in
        except Exception as e:
            self.state = "error"
            self.error = f"{type(e).__name__}: {e}"
        finally:
            try:
                if self.cur is not None and not self.cur.closed:
                    self._finalize(self.cur)
            except Exception as e:                             # report it, do not hide a failed final write
                self.error = self.error or f"{type(e).__name__}: {e}"
                self.state = "error"
            self._stats()


class Session:
    """One capture run. The capture thread writes the stream into the ring buffer and never waits for anything
    else; the processing thread reads rows from the ring (Welch PSD) and hands them to the window through ui_q."""

    def __init__(self, cfg, plan, ui_q):
        self.cfg, self.plan, self.ui_q = cfg, plan, ui_q
        self.stop = threading.Event()
        self.paused = threading.Event()
        self.dropped_late = 0            # rows overwritten in the ring before the processing got to them
        self.dropped_during = 0          # rows overwritten while they were being processed
        self.dropped_gui = 0             # rows the window could not take (its queue was full)
        self.rows_done = 0
        self.processed_end = 0           # end of the last row taken by the processing (processed or skipped)
        self.wall0 = None                # wall-clock time of ring index 0 (an estimate, to within one USB transfer)
        self.stream_offset = 0           # start-up samples thrown away: stream sample = ring sample + this
        self.alive = 0
        self.ring = SampleRing(plan.row_samples, plan.ring_rows, ring_chunk(cfg.rate))
        self.rows_free = queue.Queue()
        self.rowbufs = [np.empty(plan.nfft, np.float64) for _ in range(ROW_RING)]
        for i in range(ROW_RING):
            self.rows_free.put(i)
        self.threads = [threading.Thread(target=self._capture, daemon=True, name="capture"),
                        threading.Thread(target=self._process, daemon=True, name="process")]

    @property
    def dropped_rows(self):
        return self.dropped_late + self.dropped_during + self.dropped_gui

    def start(self):
        self.alive = len(self.threads)
        for t in self.threads:
            t.start()

    def release_row(self, i):
        self.rows_free.put(i)

    def shutdown(self, timeout=5.0):
        """Stop both threads; returns True if they finished in time."""
        self.stop.set()
        with self.ring.cond:
            self.ring.cond.notify_all()
        end = time.perf_counter() + timeout
        for t in self.threads:
            t.join(max(0.0, end - time.perf_counter()))
        return not any(t.is_alive() for t in self.threads)

    def _capture(self):
        radio = Radio(self.cfg, self.plan, self.stop)
        ring = self.ring
        try:
            ready_info = radio.open()
            self.stream_offset = radio.discard(ring.buf, ready_info["fs"])   # they go to the (still empty) ring and are overwritten
            ready_info["discarded"] = self.stream_offset
            self.wall0 = time.time()
            self.ui_q.put(("ready", ready_info))
            cs = _CapStats(time.perf_counter())
            while not self.stop.is_set():
                pos, n = ring.write_span()
                t0 = time.perf_counter()
                ok = radio.read(ring.buf[pos:pos + n])
                t1 = time.perf_counter()
                if not ok:
                    break
                ring.advance(n)
                cs.call(t0, t1, n)
                if cs.due(t1):
                    self.ui_q.put(("perf_cap", cs.snapshot(t1, radio.errors, radio.last_error, ring.head)))
        except Exception as e:
            if not self.stop.is_set():
                self.ui_q.put(("error", f"{type(e).__name__}: {e}"))
        finally:
            self.stop.set()
            with ring.cond:
                ring.cond.notify_all()
            radio.close()
            self.ui_q.put(("closed", None))

    def _process(self):
        p, ring = self.plan, self.ring
        R = p.row_samples
        ps = _ProcStats(time.perf_counter(), R / self.cfg.rate)
        try:
            welch = Welch(p.nfft, self.cfg.window, p.n_win)
            k = 0                                          # next row to take
            while not self.stop.is_set():
                now = time.perf_counter()
                if ps.due(now):
                    self.ui_q.put(("perf_proc", ps.snapshot(now, max(0, ring.head - self.processed_end))))
                if not ring.wait_for((k + 1) * R, self.stop):
                    continue
                if self.paused.is_set():                   # the ring keeps recording; nothing is queued up meanwhile
                    k = ring.head // R
                    self.processed_end = k * R
                    self.stop.wait(0.05)
                    continue
                st = ring.check(k * R, R)
                if st.late:                                # fell behind by more than the ring holds: jump to the newest row
                    newest = ring.head // R - 1
                    self.dropped_late += max(0, newest - k)
                    k = newest
                    self.processed_end = k * R
                    continue
                if st != RingStatus.OK:                    # not expected after wait_for: a bug, say so
                    raise RuntimeError(f"row {k}: unexpected ring status {st.name}")
                ps.lag(ring.head - (k + 1) * R)
                try:
                    ri = self.rows_free.get_nowait()
                except queue.Empty:                        # the window is behind: skip this row
                    self.dropped_gui += 1
                    k += 1
                    self.processed_end = k * R
                    continue
                t_a = time.perf_counter()
                welch.reset()
                welch.accumulate(ring.row_view(k), p.hop, p.n_win)
                ps.row(time.perf_counter() - t_a)
                st = ring.check(k * R, R)
                if st.late:                                # overwritten while it was being processed
                    self.dropped_during += 1
                    self.rows_free.put(ri)
                    k += 1
                    self.processed_end = k * R
                    continue
                if st != RingStatus.OK:
                    raise RuntimeError(f"row {k}: unexpected ring status {st.name}")
                welch.power_into(self.rowbufs[ri])
                self.rows_done += 1
                self.ui_q.put(("row", (ri, k * R)))
                k += 1
                self.processed_end = k * R
        except Exception as e:
            if not self.stop.is_set():
                self.ui_q.put(("error", f"{type(e).__name__}: {e}"))
            self.stop.set()


# =============================================================================
# WINDOW
# =============================================================================

# Design tokens. Dark ink surfaces keep the waterfall readable and glare-free; one cool accent for focus and
# the primary action, amber only for measured values and hover.
C = dict(bg="#11151C", panel="#171D26", raised="#1E2632", input="#222B38", line="#2E3948",
         text="#DCE3EE", muted="#8B97A8", accent="#6FB1D1", accent_ink="#0E1A22", amber="#E3AE55",
         error="#E57373", ok="#7FBF9A", band="#EAF1FA", psd="#CFE0EE", cap="#B79CED", sel="#26475A")
S1, S2, S3 = 4, 8, 16          # spacing scale: within a control, within a group, between groups

RSSI_TEXT_W, NB_TEXT_W = 44, 63   # width of the table in characters
REG_RAW_WINDOW_S = 30.0          # raw detections counted in the caption of the list

SPANS_MODES = {"Ranges": "spans", "List": "list"}
SPANS_HINT = {"spans": "[{flo~fhi,step,width},...]", "list": "[{freq,width},...]"}

RESTART_FIELDS = ("center", "rate", "bandwidth", "gain", "ring_dur", "fft_n", "fft_df", "avr_n", "avr_dur",
                  "window", "overlap", "spans_mode", "spans")
FIELD_OF = {"center": "center", "rate": "rate", "bandwidth": "bandwidth", "gain": "gain", "fft": "fft_n",
            "avr": "avr_n", "window": "window", "overlap": "overlap", "spans": "spans", "rows": "rows",
            "ring_dur": "ring_dur", "spans_mode": "spans_mode",
            "nb_max_width": "nb_max_width", "nb_guard": "nb_guard", "nb_thr": "nb_thr"}
NOTE_DEFAULT = "Radio, Spectrum and Channels apply on restart."

TIPS = {
    "center": "Capture center, Hz (k/M/G or 48e6). auto = middle of the channel ranges.",
    "rate": "Sample rate, Hz. The capture band is the rate (0.52 to 61.44 MS/s).",
    "bandwidth": "Analog bandwidth, Hz. auto = the sample rate (the device limit is 56 MHz).",
    "gain": "Manual RX gain, dB. The status bar shows what the device reports.",
    "fft": "FFT size (points per window) and the resolution (bin spacing) are one setting: edit either and the other "
           "follows.\nThe size is a fast FFT size (2^a 3^b 5^c 7^d), so a typed resolution becomes the nearest one that can be reached.",
    "avr": "Windows averaged per row and the time of one row are one setting: edit either and the other follows.\n"
           "The time is a whole number of windows, so a typed time becomes the nearest one that can be reached.",
    "window": "FFT window. rect = no multiplication.",
    "overlap": "Overlap of consecutive windows, 0 to below 1.",
    "spans": "Channels, in the form chosen above. Only channels fully inside the capture band get an RSSI.",
    "spans_mode": "Ranges: [{flo~fhi,step,width},...]  channels from flo to fhi with a step.\n"
                  "List: [{freq,width},...]  every channel by its centre and width.\nOne or the other.",
    "ring": "Seconds of raw samples kept in the ring buffer in memory (not written to disk).\n"
            "The signal search reads its rows from this buffer.",
    "psd": "Which spectrum is drawn and used for the RSSI.",
    "rows": "Waterfall depth in rows. Applies immediately.",
    "range": "Colour range of the waterfall and vertical range of the PSD, dBFS/Hz.",
    "overlay": "Channel bands in RSSI mode, detected signals in NarrowBand Det mode.",
    "thr": "On the PSD: dashed = the level a bin must exceed to count as a signal (the local noise floor + Threshold),\n"
           "dotted = the noise floor. The floor is a block median, so it rises under wide signals.",
    "cap_center": "Centre of the band to record, Hz (k/M/G or 428e6). Dragging the overlay changes it.",
    "cap_rate": "Width of the band = output sample rate, Hz: a guide. The nearest rate that a ratio of whole numbers\n"
                "of the capture rate gives is used, it is shown below. Dragging an edge of the overlay changes it.",
    "cap_name": "Start of the file name.",
    "cap_dir": "Folder of the recordings (created if it does not exist).",
    "cap_pattern": "File name without the extension (.wav and .json are added). Fields: {name} {time} {utc} {freq} {rate} {id} {n}.\n"
                   "A format follows a colon: {time:%Y-%m-%d_%H-%M-%S}  {freq:.3f} (MHz)  {n:03d}.\n"
                   "{time} is local, {utc} UTC; {rate} is in kHz; {id} is the signal's id; {n} counts up while the name is taken.\n"
                   "Any other text stays as it is ({{ and }} give braces). An existing file is never overwritten.",
    "zoom": "Mouse wheel: zoom the frequency axis around the cursor.\nCtrl+wheel on the PSD: zoom its levels.\n"
            "Drag: pan.  Double-click a plot: reset.",
    "nb_max_width": "Widest signal that counts as narrowband, Hz. Applies immediately and clears the tracked signals.",
    "nb_guard": "The noise for the SNR is the median PSD in a region centred on the signal, as wide as\n"
                "signal width x (1 + guard). Applies immediately and clears the tracked signals.",
    "nb_thr": "Detection threshold above the local noise floor, dB. Applies immediately and clears the tracked signals.",
}


class Tooltip:
    """Small hover hint for a widget."""

    def __init__(self, widget, text):
        self.widget, self.text, self.tip, self.job = widget, text, None, None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _e):
        self._hide()
        self.job = self.widget.after(500, self._show)

    def _show(self):
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f"+{x}+{y}")
        text = self.text() if callable(self.text) else self.text
        tk.Label(self.tip, text=text, justify="left", bg=C["raised"], fg=C["text"], bd=1, relief="solid",
                 padx=S2, pady=S1, wraplength=300).pack()

    def _hide(self, _e=None):
        if self.job:
            self.widget.after_cancel(self.job)
            self.job = None
        if self.tip:
            self.tip.destroy()
            self.tip = None


class ScrollFrame(ttk.Frame):
    """Vertically scrollable container; its scrollbar appears only when the content is taller than the view."""

    def __init__(self, parent, width=326):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, bg=C["panel"], highlightthickness=0, bd=0, width=width)
        self.sb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.sb.set)
        self.inner = ttk.Frame(self.canvas)
        self.win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.inner.bind("<Configure>", lambda e: self._resize())
        self.canvas.bind("<Configure>", self._on_canvas)
        self.bind("<Enter>", self._bind_wheel)
        self.bind("<Leave>", self._unbind_wheel)

    def _on_canvas(self, e):
        self.canvas.itemconfigure(self.win, width=e.width)
        self._resize()

    def _resize(self):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        need = self.inner.winfo_reqheight() > self.canvas.winfo_height() > 1
        if need and not self.sb.winfo_ismapped():
            self.sb.pack(side="right", fill="y")
        elif not need and self.sb.winfo_ismapped():
            self.sb.pack_forget()

    def _bind_wheel(self, _e):
        self.bind_all("<MouseWheel>", self._wheel)
        self.bind_all("<Button-4>", self._wheel)
        self.bind_all("<Button-5>", self._wheel)

    def _unbind_wheel(self, _e):
        self.unbind_all("<MouseWheel>")
        self.unbind_all("<Button-4>")
        self.unbind_all("<Button-5>")

    def _wheel(self, e):
        if self.inner.winfo_reqheight() <= self.canvas.winfo_height():
            return
        up = getattr(e, "delta", 0) > 0 or getattr(e, "num", 0) == 4
        self.canvas.yview_scroll(-2 if up else 2, "units")


def pick_font(root, candidates, fallback):
    have = set(tkfont.families(root))
    for name in candidates:
        if name in have:
            return name
    return fallback


PERF_COLS = ("session", "t_s", "stream_s", "rate_ratio", "cap_read_max_ms", "cap_gap_max_ms", "cap_errors", "cap_cpu",
             "rows", "proc_ms_mean", "proc_ms_max", "proc_rt_min", "proc_rt_mean", "proc_busy", "proc_cpu",
             "lag_ms_max", "backlog_ms", "buffer_s", "drop_late", "drop_during", "drop_gui",
             "poll_ms_mean", "poll_ms_max", "draw_n", "draw_ms_mean", "draw_ms_max", "detect_ms_mean", "detect_ms_max",
             "gui_cpu", "ui_q_max", "cpu_total")


def _num(x, spec):
    return "nan" if x is None or not math.isfinite(x) else format(x, spec)


class PerfLog:
    """Collects the per-second numbers of the three threads (capture and processing send snapshots, the window
    measures itself), writes one TSV row per second and keeps totals for the summary. Window thread only.

    rate_ratio   samples received / (sample rate x time): a dip with no surplus afterwards means the stream lost
                 samples; a dip followed by a surplus is a pause that the USB buffers bridged (nothing lost)
    cap_gap_ms   longest pause between two reads; the USB buffers must bridge it (usb_ms), or samples are lost
    proc_rt_*    signal time of a row / time to process it: below 1 the processing cannot keep up
    lag_ms_max   worst distance, at the start of processing a row, between the end of that row and the stream head
    backlog_ms   recorded but not yet processed at the end of the interval
    poll_ms_*    one pass of the window's polling loop without the redraw (draw_* is the redraw)
    *_cpu        CPU seconds per second of that thread (cpu_total: the whole process, in cores)"""

    def __init__(self, fh, session_no, fs, row_dur, usb_ms):
        self.fh, self.no, self.fs, self.row_dur, self.usb_ms = fh, session_no, fs, row_dur, usb_ms
        self.t0 = self.t_last = time.perf_counter()
        self.cpu_proc0, self.cpu_gui0 = time.process_time(), time.thread_time()
        self.cap = self.proc = None
        self.errors_seen = 0
        self.drops_seen = (0, 0, 0)
        self.poll = [0, 0.0, 0.0]        # count, sum, max (s)
        self.draw = [0, 0.0, 0.0]
        self.det = [0, 0.0, 0.0]
        self.q_max = 0
        self.last = {}
        self.tot = {"intervals": 0, "ratio_sum": 0.0, "ratio_min": float("inf"), "ratio_bad": 0, "gap_max": 0.0,
                    "gap_bad": 0, "errors": 0, "rows": 0, "compute": 0.0, "rt_min": float("inf"), "slow": 0,
                    "busy_sum": 0.0, "busy_max": 0.0, "lag_max": 0, "backlog_max": 0, "poll_max": 0.0,
                    "draw_n": 0, "draw_sum": 0.0, "draw_max": 0.0, "cpu_sum": 0.0, "cpu_max": 0.0,
                    "last_error": "", "samples": 0, "cap_dt": 0.0, "read_max": 0.0}

    @staticmethod
    def _add(acc, dt):
        acc[0] += 1
        acc[1] += dt
        if dt > acc[2]:
            acc[2] = dt

    def add_poll(self, dt, qsize):
        self._add(self.poll, dt)
        if qsize > self.q_max:
            self.q_max = qsize

    def add_draw(self, dt):
        self._add(self.draw, dt)

    def add_detect(self, dt):
        self._add(self.det, dt)

    def event(self, name, detail=""):
        if self.fh:
            self.fh.write(f"# event\t{time.perf_counter() - self.t0:.3f}\t{name}\t{detail}\n")
            self.fh.flush()

    def tick(self, now, session):
        """One row as soon as the capture and the processing thread have each sent their snapshot (about once per
        PERF_PERIOD_S; the intervals of the three threads are not aligned to the millisecond). If one of them
        stays silent for 2.5 periods, the row is written without it."""
        both = self.cap is not None and self.proc is not None
        stale = now - self.t_last >= 2.5 * PERF_PERIOD_S and (self.cap is not None or self.proc is not None)
        if not (both or stale):
            return
        if not both:
            self.event("missing_snapshot", "capture" if self.cap is None else "processing")
        dt = now - self.t_last
        cap, proc = self.cap or {}, self.proc or {}
        nan = float("nan")
        ratio = cap["samples"] / (self.fs * cap["dt"]) if cap and cap["dt"] > 0 else nan
        errors = 0
        if cap:
            errors = cap["errors"] - self.errors_seen
            self.errors_seen = cap["errors"]
        if errors > 0:
            self.event("read_error", f"{errors} in the last interval, last: {cap['last_error']}")
        rows = proc.get("rows", 0)
        rt_mean = rows * self.row_dur / proc["compute"] if proc and proc["compute"] > 0 else nan
        rt_min = proc.get("rt_min", nan)
        busy = proc["compute"] / proc["dt"] if proc and proc["dt"] > 0 else nan
        ring = session.ring
        drops = (session.dropped_late, session.dropped_during, session.dropped_gui)
        dd = tuple(a - b for a, b in zip(drops, self.drops_seen))
        self.drops_seen = drops
        cpu_total = (time.process_time() - self.cpu_proc0) / dt
        gui_cpu = (time.thread_time() - self.cpu_gui0) / dt
        ms = lambda v: v * 1e3  # noqa: E731
        mean = lambda a: a[1] / a[0] if a[0] else nan  # noqa: E731
        row = {
            "session": self.no, "t_s": now - self.t0, "stream_s": ring.head / self.fs, "rate_ratio": ratio,
            "cap_read_max_ms": ms(cap.get("read_max", nan)), "cap_gap_max_ms": ms(cap.get("gap_max", nan)),
            "cap_errors": errors, "cap_cpu": cap["cpu"] / cap["dt"] if cap and cap["dt"] > 0 else nan,
            "rows": rows, "proc_ms_mean": ms(proc["compute"] / rows) if rows else nan,
            "proc_ms_max": ms(proc.get("compute_max", nan)), "proc_rt_min": rt_min, "proc_rt_mean": rt_mean,
            "proc_busy": busy, "proc_cpu": proc["cpu"] / proc["dt"] if proc and proc["dt"] > 0 else nan,
            "lag_ms_max": ms(proc.get("lag_max", 0) / self.fs) if proc else nan,
            "backlog_ms": ms(proc.get("backlog", 0) / self.fs) if proc else nan,
            "buffer_s": min(ring.head, ring.size) / self.fs, "drop_late": dd[0], "drop_during": dd[1], "drop_gui": dd[2],
            "poll_ms_mean": ms(mean(self.poll)), "poll_ms_max": ms(self.poll[2]) if self.poll[0] else nan,
            "draw_n": self.draw[0], "draw_ms_mean": ms(mean(self.draw)), "draw_ms_max": ms(self.draw[2]) if self.draw[0] else nan,
            "detect_ms_mean": ms(mean(self.det)), "detect_ms_max": ms(self.det[2]) if self.det[0] else nan,
            "gui_cpu": gui_cpu, "ui_q_max": self.q_max, "cpu_total": cpu_total}
        self.last = row
        if self.fh:
            fm = {"session": "d", "cap_errors": "d", "rows": "d", "drop_late": "d", "drop_during": "d", "drop_gui": "d",
                  "draw_n": "d", "ui_q_max": "d"}
            self.fh.write("\t".join(str(row[c]) if fm.get(c) == "d" else _num(row[c], ".4g") for c in PERF_COLS) + "\n")
            self.fh.flush()
        t = self.tot
        t["intervals"] += 1
        if math.isfinite(ratio):
            t["ratio_sum"] += ratio
            t["ratio_min"] = min(t["ratio_min"], ratio)
            t["ratio_bad"] += ratio < PERF_SLOW_RATIO
            t["samples"] += cap["samples"]
            t["cap_dt"] += cap["dt"]
        if cap:
            t["gap_max"] = max(t["gap_max"], cap["gap_max"])
            t["read_max"] = max(t["read_max"], cap["read_max"])
            t["gap_bad"] += ms(cap["gap_max"]) > self.usb_ms
        t["errors"] += errors
        if cap and cap["last_error"]:
            t["last_error"] = cap["last_error"]
        t["rows"] += rows
        if proc:
            t["compute"] += proc["compute"]
            if rows:
                t["rt_min"] = min(t["rt_min"], rt_min)
                t["slow"] += rt_min < 1.0
            if math.isfinite(busy):
                t["busy_sum"] += busy
                t["busy_max"] = max(t["busy_max"], busy)
            t["lag_max"] = max(t["lag_max"], proc["lag_max"])
            t["backlog_max"] = max(t["backlog_max"], proc["backlog"])
        t["poll_max"] = max(t["poll_max"], self.poll[2])
        t["draw_n"] += self.draw[0]
        t["draw_sum"] += self.draw[1]
        t["draw_max"] = max(t["draw_max"], self.draw[2])
        t["cpu_sum"] += cpu_total
        t["cpu_max"] = max(t["cpu_max"], cpu_total)
        self.poll, self.draw, self.det, self.q_max = [0, 0.0, 0.0], [0, 0.0, 0.0], [0, 0.0, 0.0], 0
        self.cpu_proc0, self.cpu_gui0, self.t_last = time.process_time(), time.thread_time(), now
        self.cap = self.proc = None

    def tooltip_text(self):
        r = self.last
        if not r:
            return "No numbers yet; they appear one second after the start."
        return (f"Capture: {_num(r['rate_ratio'], '.3f')} x the real-time sample count; longest pause between reads "
                f"{_num(r['cap_gap_max_ms'], '.1f')} ms (the USB buffers hold {self.usb_ms:.0f} ms); read errors {r['cap_errors']}\n"
                f"Processing: a row of {self.row_dur * 1e3:.0f} ms is done in {_num(r['proc_ms_mean'], '.1f')} ms "
                f"({_num(r['proc_rt_mean'], '.1f')}x real time, worst {_num(r['proc_rt_min'], '.1f')}x); busy {_num(r['proc_busy'] * 100, '.0f')} %\n"
                f"Backlog: {_num(r['backlog_ms'], '.0f')} ms not yet processed; worst delay at the start of a row "
                f"{_num(r['lag_ms_max'], '.0f')} ms\n"
                f"Window: poll {_num(r['poll_ms_mean'], '.1f')} / {_num(r['poll_ms_max'], '.1f')} ms (mean / max), redraw "
                f"{_num(r['draw_ms_mean'], '.0f')} ms x {r['draw_n']} per s, detection {_num(r['detect_ms_mean'], '.1f')} ms\n"
                f"CPU, in cores: capture {_num(r['cap_cpu'], '.2f')}, processing {_num(r['proc_cpu'], '.2f')}, "
                f"window {_num(r['gui_cpu'], '.2f')}, whole process {_num(r['cpu_total'], '.2f')}")

    def summary_lines(self, session):
        t, n = self.tot, max(self.tot["intervals"], 1)
        wall = time.perf_counter() - self.t0
        drops = (session.dropped_late, session.dropped_during, session.dropped_gui) if session else (0, 0, 0)
        stream = session.ring.head / self.fs if session else float("nan")
        rt_mean = t["rows"] * self.row_dur / t["compute"] if t["compute"] > 0 else float("nan")
        return [
            f"SUMMARY session {self.no}: {wall:.1f} s wall, {_num(stream, '.1f')} s of stream, {t['rows']} rows "
            f"({t['rows'] / max(wall, 1e-9):.1f}/s), {t['intervals']} intervals",
            f"  capture: received/expected samples mean {_num(t['ratio_sum'] / n, '.4f')}, worst interval {_num(t['ratio_min'], '.4f')}; "
            f"longest pause between reads {t['gap_max'] * 1e3:.1f} ms (USB buffers {self.usb_ms:.0f} ms); longest read "
            f"{t['read_max'] * 1e3:.0f} ms (a full piece of {ring_chunk(self.fs) / self.fs * 1e3:.0f} ms of signal normally takes about that long); read errors {t['errors']}"
            + (f" (last: {t['last_error']})" if t['last_error'] else ""),
            f"  processing: {_num(rt_mean, '.1f')}x real time on average, worst interval {_num(t['rt_min'], '.1f')}x; "
            f"busy mean {t['busy_sum'] / n * 100:.0f} %, max {t['busy_max'] * 100:.0f} %; worst delay "
            f"{t['lag_max'] / self.fs * 1e3:.0f} ms; backlog max {t['backlog_max'] / self.fs * 1e3:.0f} ms",
            f"  dropped rows: overwritten before processing {drops[0]}, during {drops[1]}, window {drops[2]}",
            f"  window: poll max {t['poll_max'] * 1e3:.1f} ms, redraw mean {t['draw_sum'] / max(t['draw_n'], 1) * 1e3:.0f} ms, "
            f"max {t['draw_max'] * 1e3:.0f} ms; CPU of the process mean {t['cpu_sum'] / n:.2f} cores, max {t['cpu_max']:.2f}",
            f"  ANOMALIES: intervals with lost samples (ratio < {PERF_SLOW_RATIO}) {t['ratio_bad']}; "
            f"pauses between reads longer than the USB buffers {t['gap_bad']}; intervals with processing slower than real time "
            f"{t['slow']}; read errors {t['errors']}; dropped rows {sum(drops)}"]


class App:
    def __init__(self, root, fields, metrics_fh=None, cap_defaults=None):
        self.root = root
        cd = cap_defaults or {}
        self.cap_center = None           # Hz of the band to record (None until the capture is up)
        self.cap_rate = float(cd.get("rate", CAP_DEFAULT_RATE))   # its width = output rate, a guide, Hz
        self.cap_name = ""
        self.cap_dir0 = cd.get("dir", "captures")
        self.cap_pattern0 = cd.get("pattern", DEFAULT_CAP_PATTERN)
        self.cap_plan = None             # RatePlan of the current centre/width, or None if they are not valid
        self.cap_patches = []            # the band overlay: one rectangle per plot
        self.cap_drag = None
        self._cap_cursor_now = ""
        self.recorder = None
        self._confirming = False         # a confirmation pop-up is open
        self.sections = {}               # fold-able groups of controls by title
        self.fft_primary, self.fft_guide = "n", None     # which of FFT size / resolution was typed last, and the resolution typed
        self.avr_primary, self.dur_guide = "dur", None   # which of windows / duration was typed last, and the duration typed
        self._written = {}               # what the program last wrote into a paired field
        self.sel = None                  # selected list row: ("sig", ident) or ("ch", index)
        self.metrics_fh = metrics_fh
        self.perflog = None
        self.session_no = 0
        self.session = None
        self.ui_q = queue.Queue()
        self.cfg = self.plan = None
        self.view = None             # display model of the running session (None until the device is ready)
        self.state = "Starting"
        self.error_text = ""
        self.paused = False
        self.rows_total = 0
        self.period = 0.0            # seconds per row, smoothed
        self.last_row_t = None
        self.last_draw = self.last_text = self.last_status = 0.0
        self.draw_cost = 0.0         # duration of the last redraw, s
        self.zoom_x = None           # (lo, hi) MHz while the frequency axis is zoomed, else None
        self.zoom_y = None           # (lo, hi) dBFS/Hz while the PSD levels are zoomed, else None
        self._drag = None            # pan state while the left button is held on a plot
        self.need_draw = False
        self.hover_idx = -1
        self.vmin = self.vmax = None
        self.sigint = False
        self.closing = False
        self.spans_invalid = False

        self._fonts()
        self._theme()
        self.vars = {}
        self.entries = {}
        self.applied = {}
        self.mode = tk.StringVar(value="rssi")      # 'rssi' | 'nb': what the overlay and the table show
        self._build_layout()
        self._set_fields(fields)
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.bind("<Control-Return>", lambda e: self.apply_and_restart())
        self.root.bind("<space>", self._on_space)
        self.root.bind("<Escape>", lambda e: self.on_close())
        self.apply_and_restart()
        self.root.after(POLL_MS, self.poll)

    # ------------------------------------------------------------------ look

    def _fonts(self):
        ui = pick_font(self.root, ("Segoe UI", "SF Pro Text", "Noto Sans", "DejaVu Sans"), "TkDefaultFont")
        mono = pick_font(self.root, ("Cascadia Mono", "Consolas", "DejaVu Sans Mono", "Courier New"), "TkFixedFont")
        self.f_ui, self.f_mono = ui, mono
        for name in ("TkDefaultFont", "TkTextFont", "TkMenuFont"):
            tkfont.nametofont(name).configure(family=ui, size=9)

    def _theme(self):
        r = self.root
        r.configure(bg=C["bg"])
        st = ttk.Style(r)
        st.theme_use("clam")
        ui = (self.f_ui, 9)
        st.configure(".", background=C["panel"], foreground=C["text"], fieldbackground=C["input"],
                     bordercolor=C["line"], lightcolor=C["panel"], darkcolor=C["panel"], troughcolor=C["bg"],
                     focuscolor=C["accent"], font=ui)
        st.configure("TFrame", background=C["panel"])
        st.configure("Bar.TFrame", background=C["bg"])
        st.configure("TLabel", background=C["panel"], foreground=C["text"])
        st.configure("Muted.TLabel", foreground=C["muted"])
        st.configure("Bar.TLabel", background=C["bg"], foreground=C["text"])
        st.configure("BarMuted.TLabel", background=C["bg"], foreground=C["muted"])
        st.configure("Section.TLabelframe", background=C["panel"], bordercolor=C["line"], relief="solid")
        st.configure("Section.TLabelframe.Label", background=C["panel"], foreground=C["muted"],
                     font=(self.f_ui, 9, "bold"))
        st.configure("TEntry", fieldbackground=C["input"], foreground=C["text"], insertcolor=C["text"],
                     bordercolor=C["line"], lightcolor=C["input"], darkcolor=C["input"], padding=(S1, 2))
        st.map("TEntry", bordercolor=[("focus", C["accent"])], lightcolor=[("focus", C["accent"])],
               fieldbackground=[("disabled", C["panel"])], foreground=[("disabled", C["muted"])])
        st.configure("Dirty.TEntry", bordercolor=C["accent"], lightcolor=C["accent"], darkcolor=C["accent"])
        st.configure("Invalid.TEntry", bordercolor=C["error"], lightcolor=C["error"], darkcolor=C["error"])
        st.configure("TCombobox", fieldbackground=C["input"], background=C["raised"], foreground=C["text"],
                     arrowcolor=C["muted"], bordercolor=C["line"], lightcolor=C["input"], darkcolor=C["input"],
                     padding=(S1, 2))
        st.map("TCombobox", fieldbackground=[("readonly", C["input"])], foreground=[("readonly", C["text"])],
               bordercolor=[("focus", C["accent"])], lightcolor=[("focus", C["accent"])])
        st.configure("Dirty.TCombobox", bordercolor=C["accent"], lightcolor=C["accent"], darkcolor=C["accent"])
        r.option_add("*TCombobox*Listbox.background", C["input"])
        r.option_add("*TCombobox*Listbox.foreground", C["text"])
        r.option_add("*TCombobox*Listbox.selectBackground", C["accent"])
        r.option_add("*TCombobox*Listbox.selectForeground", C["accent_ink"])
        st.configure("TCheckbutton", background=C["panel"], foreground=C["text"], indicatorbackground=C["input"],
                     indicatorforeground=C["accent"])
        st.map("TCheckbutton", background=[("active", C["panel"])], indicatorbackground=[("selected", C["accent"])])
        st.configure("TButton", background=C["raised"], foreground=C["text"], bordercolor=C["line"],
                     lightcolor=C["raised"], darkcolor=C["raised"], padding=(S2, S1 + 1))
        st.map("TButton", background=[("active", C["line"]), ("disabled", C["panel"])],
               foreground=[("disabled", C["muted"])])
        st.configure("Accent.TButton", background=C["accent"], foreground=C["accent_ink"], bordercolor=C["accent"],
                     lightcolor=C["accent"], darkcolor=C["accent"], font=(self.f_ui, 9, "bold"))
        st.map("Accent.TButton", background=[("active", "#8CC4DE"), ("disabled", C["line"])],
               foreground=[("disabled", C["muted"])])
        st.configure("Vertical.TScrollbar", background=C["raised"], troughcolor=C["panel"], bordercolor=C["panel"],
                     arrowcolor=C["muted"], lightcolor=C["raised"], darkcolor=C["raised"])
        st.map("Vertical.TScrollbar", background=[("active", C["line"]), ("disabled", C["panel"])],
               troughcolor=[("disabled", C["panel"])])
        st.configure("Seg.Toolbutton", background=C["raised"], foreground=C["text"], bordercolor=C["line"],
                     lightcolor=C["raised"], darkcolor=C["raised"], padding=(S2, S1 + 1), anchor="center")
        st.map("Seg.Toolbutton", background=[("selected", C["accent"]), ("active", C["line"])],
               foreground=[("selected", C["accent_ink"])])
        st.configure("Bar.TSeparator", background=C["line"])
        st.configure("SectionHead.TButton", background=C["panel"], foreground=C["muted"], borderwidth=0, relief="flat",
                     anchor="w", padding=(S1, 0), font=(self.f_ui, 9, "bold"), bordercolor=C["panel"],
                     lightcolor=C["panel"], darkcolor=C["panel"], focuscolor=C["accent"])
        st.map("SectionHead.TButton", background=[("active", C["panel"])],
               foreground=[("active", C["accent"]), ("focus", C["text"])])

    # ---------------------------------------------------------------- layout

    def _build_layout(self):
        root = self.root
        root.title("bladeRF waterfall")
        root.geometry("1380x820")
        root.minsize(1120, 660)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        self._build_controls(root)
        self._build_plots(root)
        self._build_rssi(root)
        self._build_status(root)

    def _entry(self, parent, row, label, name, unit, width=10, tip=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, S2), pady=S1)
        var = tk.StringVar()
        e = ttk.Entry(parent, textvariable=var, width=width)
        e.grid(row=row, column=1, sticky="ew", pady=S1)
        ttk.Label(parent, text=unit, style="Muted.TLabel").grid(row=row, column=2, sticky="w", padx=(S1, 0))
        self.vars[name], self.entries[name] = var, e
        var.trace_add("write", lambda *_, n=name: self._on_edit(n))
        e.bind("<FocusOut>", lambda ev, n=name: self._validate_quiet(n), add="+")
        if tip:
            Tooltip(e, tip)
        return e

    def _section(self, parent, title, pack=None, expanded=True):
        """A group of controls that folds and unfolds by a click on its title (a button: Space and Enter work too). The
        title is always there; only the framed body goes. Returns the frame to put the controls into."""
        holder = ttk.Frame(parent)
        holder.pack(**(pack or {"fill": "x", "pady": (0, S2)}))
        head = ttk.Button(holder, style="SectionHead.TButton", cursor="hand2")
        head.pack(fill="x")
        card = ttk.LabelFrame(holder, style="Section.TLabelframe", padding=(S2, S1, S2, S2))
        body = ttk.Frame(card)
        body.pack(fill="x")
        state = {"open": expanded}

        def show():
            head.configure(text=("\u25BE  " if state["open"] else "\u25B8  ") + title)
            if state["open"]:
                card.pack(fill="x")
            else:
                card.pack_forget()

        def toggle():
            state["open"] = not state["open"]
            show()
        head.configure(command=toggle)
        show()
        self.sections[title] = {"frame": holder, "card": card, "body": body, "head": head, "state": state, "toggle": toggle}
        return body

    def _pair_row(self, parent, row, label, a, b, tip):
        """Two linked values in one row: a = (name, unit), b = (name, unit). Edit either, the other follows."""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, S2), pady=S1)
        box = ttk.Frame(parent)
        box.grid(row=row, column=1, columnspan=2, sticky="ew", pady=S1)
        box.columnconfigure(0, weight=1)
        box.columnconfigure(3, weight=1)
        for i, (name, unit) in enumerate((a, b)):
            col = 0 if i == 0 else 3
            var = tk.StringVar()
            e = ttk.Entry(box, textvariable=var, width=8)
            e.grid(row=0, column=col, sticky="ew")
            ttk.Label(box, text=unit, style="Muted.TLabel").grid(row=0, column=col + 1, sticky="w", padx=(S1, 0))
            self.vars[name], self.entries[name] = var, e
            var.trace_add("write", lambda *_, n=name: self._on_edit(n))
            e.bind("<Return>", lambda ev, n=name: self._pair_commit(n, force=True), add="+")
            e.bind("<FocusOut>", lambda ev, n=name: self._pair_commit(n), add="+")
            Tooltip(e, tip)
        ttk.Label(box, text="\u2194", style="Muted.TLabel").grid(row=0, column=2, padx=S1)

    def _put(self, name, text):
        self._written[name] = text
        self.vars[name].set(text)

    def _pair_commit(self, name, force=False):
        """A paired field was edited (or Enter was pressed in it): it becomes the guide of its pair, the other value follows."""
        txt = self.vars[name].get().strip()
        if txt == self._written.get(name) and not force:
            return                                            # leaving a field that was not edited changes nothing
        try:
            if name == "fft_n":
                try:
                    n = int(txt)
                except ValueError:
                    raise ValueError("FFT size must be a whole number") from None
                if n < 16:
                    raise ValueError("FFT size must be at least 16")
                self.fft_primary = "n"
            elif name == "fft_df":
                df = parse_hz(txt)
                if df <= 0:
                    raise ValueError("Resolution must be above 0")
                self.fft_primary, self.fft_guide = "df", df
            elif name == "avr_n":
                try:
                    n = int(txt)
                except ValueError:
                    raise ValueError("Windows must be a whole number") from None
                if n < 1:
                    raise ValueError("Windows must be at least 1")
                self.avr_primary = "n"
            else:
                try:
                    dur = float(txt)
                except ValueError:
                    raise ValueError("The time must be a number of seconds") from None
                if dur <= 0:
                    raise ValueError("The time must be above 0")
                self.avr_primary, self.dur_guide = "dur", dur
        except ValueError as e:
            self._mark_invalid(name, str(e))
            return
        self._clear_invalid()
        self._spectrum_sync()

    def _spectrum_sync(self):
        """Keep the pairs in step: FFT size <-> resolution, then windows <-> row time (which depends on the size, the
        overlap and the rate). The value typed last is the guide and keeps its meaning; the other follows, and the guide's own
        field shows what can be reached."""
        try:
            rate = parse_hz(self.vars["rate"].get())
            overlap = float(self.vars["overlap"].get())
        except ValueError:
            return
        if not (SR_MIN <= rate <= SR_MAX and 0.0 <= overlap < 1.0):
            return
        try:
            if self.fft_primary == "df" and self.fft_guide:
                nfft = nearest_fft_size(rate / self.fft_guide)
            else:
                n = int(self.vars["fft_n"].get().strip())
                nfft = n if is_fft_size(n) else nearest_fft_size(n)
        except ValueError:
            return
        nfft = min(max(nfft, 16), FFT_MAX)
        hop = max(1, int(round(nfft * (1.0 - overlap))))
        try:
            if self.avr_primary == "dur" and self.dur_guide:
                n_win = max(1, int(round((self.dur_guide * rate - nfft) / hop)) + 1)
            else:
                n_win = max(1, int(self.vars["avr_n"].get().strip()))
        except ValueError:
            return
        self._put("fft_n", str(nfft))
        self._put("fft_df", fmt_hz_short(rate / nfft, 5))
        self._put("avr_n", str(n_win))
        self._put("avr_dur", f"{welch_samples(nfft, hop, n_win) / rate:.4g}")

    def _build_controls(self, root):
        side = ttk.Frame(root, padding=(S2, S2, S1, S2))
        side.grid(row=0, column=0, sticky="ns")
        # The primary actions stay at the top; the groups below scroll when the window is short.
        actions = ttk.Frame(side)
        actions.pack(fill="x", pady=(0, S2))
        actions.columnconfigure(0, weight=3)
        actions.columnconfigure(1, weight=2)
        self.btn_apply = ttk.Button(actions, text="Apply and restart", style="Accent.TButton",
                                    command=self.apply_and_restart)
        self.btn_apply.grid(row=0, column=0, sticky="ew", padx=(0, S1))
        Tooltip(self.btn_apply, "Ctrl+Enter")
        self.btn_view_pause = ttk.Button(actions, text="Pause", command=self.toggle_pause)
        self.btn_view_pause.grid(row=0, column=1, sticky="ew")
        Tooltip(self.btn_view_pause, "Space: freeze the display (the capture keeps recording into the ring buffer)")
        self.note = ttk.Label(actions, text=NOTE_DEFAULT, style="Muted.TLabel", wraplength=270, justify="left")
        self.note.grid(row=1, column=0, columnspan=2, sticky="w", pady=(S1 + 2, 0))
        self.msg_kind = None
        scroll = ScrollFrame(side)
        scroll.pack(fill="both", expand=True)
        outer = scroll.inner

        radio = self._section(outer, "Radio")
        radio.columnconfigure(1, weight=1)
        self._entry(radio, 0, "Center", "center", "Hz", tip=TIPS["center"])
        self._entry(radio, 1, "Rate", "rate", "S/s", tip=TIPS["rate"])
        self._entry(radio, 2, "Bandwidth", "bandwidth", "Hz", tip=TIPS["bandwidth"])
        self._entry(radio, 3, "Gain", "gain", "dB", tip=TIPS["gain"])
        self._entry(radio, 4, "Buffer", "ring_dur", "s", tip=TIPS["ring"])

        spec = self._section(outer, "Spectrum")
        spec.columnconfigure(1, weight=1)
        self._pair_row(spec, 0, "FFT", ("fft_n", "pts"), ("fft_df", "Hz"), TIPS["fft"])
        self._pair_row(spec, 1, "Average", ("avr_n", "win"), ("avr_dur", "s"), TIPS["avr"])
        ttk.Label(spec, text="Window").grid(row=2, column=0, sticky="w", padx=(0, S2), pady=S1)
        wvar = tk.StringVar()
        wc = ttk.Combobox(spec, textvariable=wvar, values=list(WINDOWS), state="readonly", width=10)
        wc.grid(row=2, column=1, sticky="ew", pady=S1)
        self.vars["window"], self.entries["window"] = wvar, wc
        wvar.trace_add("write", lambda *_: self._on_edit("window"))
        Tooltip(wc, TIPS["window"])
        self._entry(spec, 3, "Overlap", "overlap", "", tip=TIPS["overlap"])
        for name in ("rate", "overlap"):                      # the pairs follow the rate and the overlap
            self.entries[name].bind("<Return>", lambda e: self._spectrum_sync(), add="+")
            self.entries[name].bind("<FocusOut>", lambda e: self._spectrum_sync(), add="+")

        chan = self._section(outer, "Channels")
        hdr = ttk.Frame(chan)
        hdr.pack(fill="x")
        mvar = tk.StringVar()
        mc = ttk.Combobox(hdr, textvariable=mvar, values=list(SPANS_MODES), state="readonly", width=8)
        mc.pack(side="left")
        self.spans_hint = ttk.Label(hdr, text="", style="Muted.TLabel")
        self.spans_hint.pack(side="left", padx=(S2, 0))
        self.vars["spans_mode"], self.entries["spans_mode"] = mvar, mc

        def on_spans_mode(*_):
            self.spans_hint.configure(text=SPANS_HINT[SPANS_MODES.get(mvar.get(), "spans")])
            self._on_edit("spans_mode")
        mvar.trace_add("write", on_spans_mode)
        Tooltip(mc, TIPS["spans_mode"])
        self.spans_text = tk.Text(chan, height=3, width=30, wrap="word", bg=C["input"], fg=C["text"],
                                  insertbackground=C["text"], relief="flat", highlightthickness=1,
                                  highlightbackground=C["line"], highlightcolor=C["accent"], padx=S1, pady=S1,
                                  font=(self.f_mono, 9), undo=True)
        self.spans_text.pack(fill="x", pady=(S1, S1))
        self.spans_text.bind("<<Modified>>", self._on_spans_modified)
        self.spans_text.bind("<FocusOut>", lambda ev: self._validate_quiet("spans"), add="+")
        self.spans_text.bind("<Tab>", lambda ev: (self.spans_text.tk_focusNext().focus(), "break")[1])
        Tooltip(self.spans_text, TIPS["spans"])
        self.vars["spans"] = tk.StringVar()
        self.entries["spans"] = self.spans_text
        line = tkfont.nametofont("TkDefaultFont").metrics("linespace")
        box = ttk.Frame(chan, height=3 * line + 2)           # fixed height: the layout does not jump
        box.pack(fill="x")
        box.pack_propagate(False)
        self.chan_summary = ttk.Label(box, text="", style="Muted.TLabel", wraplength=262, justify="left")
        self.chan_summary.pack(anchor="nw")

        nbf = self._section(outer, "NarrowBand Det")
        nbf.columnconfigure(1, weight=1)
        self._entry(nbf, 0, "Max width", "nb_max_width", "Hz", tip=TIPS["nb_max_width"])
        self._entry(nbf, 1, "Guard", "nb_guard", "x width", tip=TIPS["nb_guard"])
        self._entry(nbf, 2, "Threshold", "nb_thr", "dB", tip=TIPS["nb_thr"])
        for name in ("nb_max_width", "nb_guard", "nb_thr"):
            self.entries[name].bind("<Return>", lambda e: self._nb_changed(), add="+")
            self.entries[name].bind("<FocusOut>", lambda e: self._nb_changed(), add="+")

        disp = self._section(outer, "Display")
        disp.columnconfigure(1, weight=1)
        ttk.Label(disp, text="PSD shows").grid(row=0, column=0, sticky="w", padx=(0, S2), pady=S1)
        self.psd_source = tk.StringVar(value=PSD_SOURCES[0])
        pc = ttk.Combobox(disp, textvariable=self.psd_source, values=list(PSD_SOURCES), state="readonly", width=16)
        pc.grid(row=0, column=1, columnspan=2, sticky="ew", pady=S1)
        pc.bind("<<ComboboxSelected>>", lambda e: self._display_changed(refocus=True))
        Tooltip(pc, TIPS["psd"])
        self._entry(disp, 1, "Rows", "rows", "", tip=TIPS["rows"])
        self.vars["rows"].trace_add("write", lambda *_: None)
        self.entries["rows"].bind("<Return>", lambda e: self._rows_changed(), add="+")
        self.entries["rows"].bind("<FocusOut>", lambda e: self._rows_changed(), add="+")
        ttk.Label(disp, text="Colormap").grid(row=2, column=0, sticky="w", padx=(0, S2), pady=S1)
        self.cmap = tk.StringVar(value=CMAPS[0])
        cc = ttk.Combobox(disp, textvariable=self.cmap, values=list(CMAPS), state="readonly", width=10)
        cc.grid(row=2, column=1, columnspan=2, sticky="ew", pady=S1)
        cc.bind("<<ComboboxSelected>>", lambda e: self._display_changed(refocus=True))
        self.auto_range = tk.BooleanVar(value=False)         # a fixed range is the default
        self.show_thr = tk.BooleanVar(value=True)            # the detection threshold on the PSD
        ar = ttk.Checkbutton(disp, text="Auto range", variable=self.auto_range,
                             command=lambda: (self._range_mode(), self._focus_plots()))
        ar.grid(row=3, column=0, sticky="w", pady=S1)
        Tooltip(ar, TIPS["range"])
        rng = ttk.Frame(disp)
        rng.grid(row=3, column=1, columnspan=2, sticky="ew", pady=S1)
        self.rmin, self.rmax = tk.StringVar(value="-120"), tk.StringVar(value="-60")
        self.e_rmin = ttk.Entry(rng, textvariable=self.rmin, width=6)
        self.e_rmax = ttk.Entry(rng, textvariable=self.rmax, width=6)
        self.e_rmin.pack(side="left")
        ttk.Label(rng, text=" to ", style="Muted.TLabel").pack(side="left")
        self.e_rmax.pack(side="left")
        ttk.Label(rng, text=" dB", style="Muted.TLabel").pack(side="left")
        for e in (self.e_rmin, self.e_rmax):
            e.bind("<Return>", lambda ev: self._display_changed(), add="+")
            e.bind("<FocusOut>", lambda ev: self._display_changed(), add="+")
        self.overlay = tk.BooleanVar(value=True)
        ov = ttk.Checkbutton(disp, text="Overlay", variable=self.overlay,
                             command=lambda: self._display_changed(refocus=True))
        ov.grid(row=4, column=0, columnspan=3, sticky="w", pady=S1)
        Tooltip(ov, TIPS["overlay"])
        th = ttk.Checkbutton(disp, text="Detection threshold", variable=self.show_thr,
                             command=lambda: self._display_changed(refocus=True))
        th.grid(row=5, column=0, columnspan=3, sticky="w", pady=S1)
        Tooltip(th, TIPS["thr"])
        rz = ttk.Button(disp, text="Reset zoom", command=self.reset_zoom)
        rz.grid(row=6, column=0, columnspan=3, sticky="ew", pady=(S2, S1))
        Tooltip(rz, TIPS["zoom"])
        self._range_mode()

    def _build_plots(self, root):
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure
        matplotlib.rcParams["font.family"] = [self.f_ui, "DejaVu Sans"]

        frame = ttk.Frame(root, style="Bar.TFrame", padding=(S1, S2, S1, S2))
        frame.grid(row=0, column=1, sticky="nsew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        self.fig = Figure(figsize=(8, 6), dpi=100, facecolor=C["bg"])
        gs = self.fig.add_gridspec(2, 2, width_ratios=[1, 0.022], height_ratios=[3, 2], hspace=0.07, wspace=0.025,
                                   left=0.095, right=0.875, top=0.955, bottom=0.09)
        self.ax_wf = self.fig.add_subplot(gs[0, 0])
        self.cax = self.fig.add_subplot(gs[0, 1])
        self.ax_psd = self.fig.add_subplot(gs[1, 0], sharex=self.ax_wf)
        self.spacer = self.fig.add_subplot(gs[1, 1])     # same width as the colorbar: the x axes line up
        self.spacer.axis("off")
        for ax in (self.ax_wf, self.ax_psd, self.cax):
            ax.set_facecolor(C["panel"])
            for sp in ax.spines.values():
                sp.set_color(C["line"])
            ax.tick_params(colors=C["muted"], labelsize=8, length=3)
        self.ax_wf.tick_params(labelbottom=False)
        self.ax_wf.set_ylabel("Age, s", color=C["muted"], fontsize=9)
        self.ax_psd.set_ylabel("PSD, dBFS/Hz", color=C["muted"], fontsize=9)
        self.ax_psd.set_xlabel("Frequency, MHz", color=C["muted"], fontsize=9)
        self.ax_psd.grid(True, color=C["line"], alpha=0.7, lw=0.6)
        self.ax_psd.xaxis.get_major_formatter().set_useOffset(False)
        self.gs = gs
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.bind("<Configure>", self._on_canvas_resize, add="+")
        self.canvas_widget.configure(bg=C["bg"], highlightthickness=0)
        self.canvas_widget.grid(row=0, column=0, sticky="nsew")
        self.empty_text = self.ax_wf.text(0.5, 0.5, "", transform=self.ax_wf.transAxes, ha="center", va="center",
                                          color=C["muted"], fontsize=11)
        self.tip = tk.Label(self.canvas_widget, text="", justify="left", bg=C["raised"], fg=C["text"], bd=1,
                            relief="solid", padx=S2, pady=S1, font=(self.f_ui, 9))
        self.canvas_widget.bind("<Button-1>", lambda e: self.canvas_widget.focus_set(), add="+")
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("button_release_event", self.on_release)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas.mpl_connect("figure_leave_event", self.on_leave)
        self.im = self.psd_line = self.levels = self.cbar = None
        self.thr_line = self.floor_line = self.thr_text = None
        self.band_pcs = []
        self.hl_patches = []
        self.nb_lcs = []

    def _on_canvas_resize(self, e):
        """Margins in pixels (room for the tick labels, the axis titles and the colour bar), not fractions of a width."""
        w, h = max(e.width, 240), max(e.height, 240)
        self.gs.update(left=min(0.3, 84.0 / w), right=1.0 - 92.0 / w, top=1.0 - 36.0 / h, bottom=48.0 / h)
        for ax in self.fig.axes:                 # GridSpec.update moves only the axes of pyplot figures: do it here
            ss = ax.get_subplotspec()
            if ss is not None:
                ax.set_position(ss.get_position(self.fig))

    def _build_rssi(self, root):
        frame = ttk.Frame(root, padding=(S1, S2, S2, S2))
        frame.grid(row=0, column=2, sticky="ns")
        self._build_capture(frame)
        box = ttk.LabelFrame(frame, text="RSSI per channel", style="Section.TLabelframe", padding=(S2, S1, S2, S2))
        box.pack(fill="both", expand=True)
        self.rssi_box = box
        seg = ttk.Frame(box)
        seg.pack(fill="x", pady=(0, S2))
        for text, val in (("RSSI", "rssi"), ("NarrowBand Det", "nb")):
            ttk.Radiobutton(seg, text=text, value=val, variable=self.mode, style="Seg.Toolbutton",
                            command=self._mode_changed).pack(side="left", fill="x", expand=True)
        head = ttk.Frame(box)
        head.pack(fill="x")
        self.rssi_caption = ttk.Label(head, text="dBFS from the displayed PSD", style="Muted.TLabel")
        self.rssi_caption.pack(side="left")
        ttk.Button(head, text="Copy", command=self.copy_rssi).pack(side="right")
        body = ttk.Frame(box)
        body.pack(fill="both", expand=True, pady=(S2, 0))
        # the column titles sit outside the scrolled text, so they stay in view while the rows scroll
        self.rssi_head = tk.Text(body, width=RSSI_TEXT_W, height=1, wrap="none", bg=C["input"], fg=C["muted"], relief="flat",
                                 highlightthickness=1, highlightbackground=C["line"], highlightcolor=C["line"],
                                 font=(self.f_mono, 9), padx=S2, pady=S1, state="disabled", cursor="arrow")
        self.rssi_head.pack(side="top", fill="x")
        rows = ttk.Frame(body)
        rows.pack(side="top", fill="both", expand=True)
        self.rssi_text = tk.Text(rows, width=RSSI_TEXT_W, wrap="none", bg=C["input"], fg=C["text"], relief="flat",
                                 highlightthickness=1, highlightbackground=C["line"], highlightcolor=C["line"],
                                 font=(self.f_mono, 9), padx=S2, pady=S1, state="disabled", cursor="arrow")
        sb = ttk.Scrollbar(rows, orient="vertical", command=self.rssi_text.yview)
        self.rssi_text.configure(yscrollcommand=sb.set)
        self.rssi_text.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.rssi_text.tag_configure("head", foreground=C["muted"])
        self.rssi_text.tag_configure("absent", foreground=C["muted"])
        self.rssi_text.tag_configure("hover", background=C["raised"], foreground=C["amber"])
        self.rssi_text.tag_configure("empty", foreground=C["muted"])
        self.rssi_text.tag_configure("selected", background=C["sel"], foreground=C["text"])
        self.rssi_text.tag_raise("hover")
        self.rssi_text.bind("<Button-1>", self._on_table_click)

    # ------------------------------------------------------------------- capture panel

    def _build_capture(self, parent):
        cap = self._section(parent, "Capture", pack={"side": "bottom", "fill": "x", "pady": (S2, 0)})
        cap.columnconfigure(1, weight=1)
        self.cap_src = ttk.Label(cap, text="", style="Muted.TLabel", wraplength=430, justify="left")
        self.cap_src.grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, S1))
        self.cap_vars, self.cap_entries = {}, {}
        for i, (key, label, unit) in enumerate((("center", "Center", "Hz"), ("rate", "Width", "Hz"), ("name", "Name", "")),
                                               start=1):
            ttk.Label(cap, text=label).grid(row=i, column=0, sticky="w", padx=(0, S2), pady=S1)
            var = tk.StringVar()
            e = ttk.Entry(cap, textvariable=var, width=18)
            e.grid(row=i, column=1, sticky="ew", pady=S1)
            ttk.Label(cap, text=unit, style="Muted.TLabel").grid(row=i, column=2, sticky="w", padx=(S1, 0))
            e.bind("<Return>", lambda ev: self._cap_fields_changed(), add="+")
            e.bind("<FocusOut>", lambda ev: self._cap_fields_changed(), add="+")
            Tooltip(e, TIPS["cap_" + key])
            self.cap_vars[key], self.cap_entries[key] = var, e
        ttk.Label(cap, text="Pattern").grid(row=4, column=0, sticky="w", padx=(0, S2), pady=S1)
        self.cap_pattern = tk.StringVar(value=self.cap_pattern0)
        pe = self.cap_pattern_entry = ttk.Combobox(cap, textvariable=self.cap_pattern, values=list(CAP_PATTERN_PRESETS), width=18)
        pe.grid(row=4, column=1, columnspan=2, sticky="ew", pady=S1)
        Tooltip(pe, TIPS["cap_pattern"])
        self.cap_pattern.trace_add("write", lambda *_: self._cap_pattern_changed())
        self.cap_preview = ttk.Label(cap, text="", style="Muted.TLabel", wraplength=430, justify="left")
        self.cap_preview.grid(row=5, column=0, columnspan=3, sticky="w")
        self.cap_info = ttk.Label(cap, text="", style="Muted.TLabel", wraplength=430, justify="left")
        self.cap_info.grid(row=6, column=0, columnspan=3, sticky="w")
        ttk.Label(cap, text="Folder").grid(row=7, column=0, sticky="w", padx=(0, S2), pady=S1)
        self.cap_dir = tk.StringVar(value=self.cap_dir0)
        de = self.cap_dir_entry = ttk.Entry(cap, textvariable=self.cap_dir, width=18)
        de.grid(row=7, column=1, sticky="ew", pady=S1)
        Tooltip(de, TIPS["cap_dir"])
        self.cap_dir.trace_add("write", lambda *_: (self._cap_preview(), self._rec_ui()))
        self.cap_browse_btn = ttk.Button(cap, text="Browse...", command=self._cap_browse)
        self.cap_browse_btn.grid(row=7, column=2, padx=(S1, 0))
        bar = ttk.Frame(cap)
        bar.grid(row=8, column=0, columnspan=3, sticky="ew", pady=(S2, 0))
        self.btn_rec = ttk.Button(bar, text="Start", style="Accent.TButton", command=self._rec_start, width=13)
        self.btn_rec.pack(side="left")
        Tooltip(self.btn_rec, "Record the band from the ring buffer. After a pause: continue in the same file, or start a new "
                              "one if the centre, the width, the name, the pattern or the folder was changed.")
        self.btn_pause = ttk.Button(bar, text="Pause", command=self._rec_pause, state="disabled")
        self.btn_pause.pack(side="left", padx=(S1, 0))
        Tooltip(self.btn_pause, "Stop writing but keep the file open.")
        self.btn_stop = ttk.Button(bar, text="Stop", command=self._rec_stop_clicked, state="disabled")
        self.btn_stop.pack(side="left", padx=(S1, 0))
        Tooltip(self.btn_stop, "Close the file.")
        self.rec_lbl = ttk.Label(cap, text="", style="Muted.TLabel", wraplength=430, justify="left")
        self.rec_lbl.grid(row=9, column=0, columnspan=3, sticky="w", pady=(S1, 0))
        self.cap_msg = ttk.Label(cap, text="", wraplength=430, justify="left", foreground=C["error"])
        self.cap_msg.grid(row=10, column=0, columnspan=3, sticky="w")
        self.cap_src.configure(text="Shift+click on a plot places the band; or pick a signal or a channel in the list, or type a centre.")

    def _cap_say(self, text="", kind="error"):
        self.cap_msg.configure(text=text, foreground=C["amber"] if kind == "warn" else C["error"])

    def _cap_ctx(self):
        """What the file name pattern is given for the band as it is now (the time is now)."""
        v = self.view
        center = self.cap_center if self.cap_center is not None else (v.fc if v is not None else 0.0)
        rate = self.cap_plan.out_rate if self.cap_plan is not None else self.cap_rate
        ident = self.sel[1] if self.sel is not None and self.sel[0] == "sig" else None
        return {"name": self.cap_name or "band", "ident": ident, "center_hz": center, "rate_hz": rate, "t": time.time(), "n": 1}

    def _cap_preview(self):
        """The file name the next recording would get, or what is wrong with the pattern."""
        if not hasattr(self, "cap_preview"):
            return None
        try:
            stem = free_stem(self.cap_dir.get().strip() or ".", self.cap_pattern.get(), self._cap_ctx())
        except ValueError as e:
            self.cap_preview.configure(text=str(e), foreground=C["error"])
            return None
        self.cap_preview.configure(text=f"\u2192 {stem}.wav  (+ .json)", foreground=C["muted"])
        return stem

    def _cap_pattern_changed(self):
        self._cap_preview()
        self._rec_ui()

    def _cap_recording(self):
        return self.recorder is not None and self.recorder.state == "recording"

    def _cap_has_file(self):
        return self.recorder is not None and self.recorder.state in ("recording", "paused")

    def _cap_check(self, c, r):
        """(ok, message) for centre c and width r (Hz); on success self.cap_plan is the plan of this width."""
        v = self.view
        if v is None:
            return False, "The device is not running."
        try:
            rp = plan_resampler(v.fs, r)
        except ValueError as e:
            return False, str(e)
        if c - rp.out_rate / 2 < v.left * 1e6 or c + rp.out_rate / 2 > v.right * 1e6:
            return False, (f"The band {fmt_hz_short(c - rp.out_rate / 2)}Hz .. {fmt_hz_short(c + rp.out_rate / 2)}Hz reaches "
                           f"outside the capture band; move it or make it narrower.")
        self.cap_plan = rp
        return True, ""

    def _cap_set(self, center=None, rate=None, name=None, source=None):
        c = self.cap_center if center is None else float(center)
        r = self.cap_rate if rate is None else float(rate)
        ok, msg = self._cap_check(c, r)
        if not ok:
            self._cap_say(msg)
            return False
        self.cap_center, self.cap_rate = c, r
        if name is not None:
            self.cap_name = name
        if source is not None:
            self.cap_src.configure(text=source)
        self._cap_sync()
        self._cap_say("")
        self._cap_overlay()
        self._cap_preview()
        self._rec_ui()
        return True

    def _cap_sync(self):
        """The fields and the info line follow the band."""
        self.cap_vars["rate"].set(fmt_hz_short(self.cap_rate))
        self.cap_vars["name"].set(self.cap_name)
        if self.cap_center is None:
            self.cap_vars["center"].set("")
            self.cap_info.configure(text="")
            return
        self.cap_vars["center"].set(fmt_hz_short(self.cap_center))
        rp, v = self.cap_plan, self.view
        if rp is not None and v is not None:
            self.cap_info.configure(
                text=f"Output rate {rp.out_rate:.3f} Hz = {v.fs * 1e-6:g} MS/s x {rp.L}/{rp.M}  ({rp.error_ppm:+.1f} ppm from "
                     f"the {fmt_hz_short(rp.target)}Hz asked); one block {rp.N / v.fs * 1e3:.0f} ms")

    def _cap_fields_changed(self):
        if self.view is None:
            return
        c_txt = self.cap_vars["center"].get().strip()
        try:
            r = parse_hz(self.cap_vars["rate"].get())
            c = parse_hz(c_txt) if c_txt else None
        except ValueError as e:
            self._cap_say(str(e))
            return
        name = self.cap_vars["name"].get().strip()
        if self._cap_recording():                    # the band belongs to the file: another centre asks for a new file
            if self._confirming:
                return
            if c is not None and c != self.recorder.spec["center"]:
                self._rec_center_request(center=c, source="Manual band")
            else:
                self._cap_sync()
            return
        if c is None:                                # no band yet: the width and the name are remembered for the next one
            self.cap_rate, self.cap_name = r, name
            self._cap_preview()
            return
        if (c, r, name) == (self.cap_center, self.cap_rate, self.cap_name):
            return
        if self._cap_set(center=c, rate=r, name=name, source="Manual band"):
            self.sel = None
            self.last_text = 0.0

    def _cap_overlay(self):
        v = self.view
        if v is None or self.cap_center is None or not self.cap_patches:
            return
        rate = self.cap_plan.out_rate if self.cap_plan is not None else self.cap_rate
        w = max(rate * 1e-6, 0.003 * self._view_span())               # never thinner than a hairline
        for r in self.cap_patches:
            r.set_x(self.cap_center * 1e-6 - w / 2)
            r.set_width(w)
            r.set_visible(True)
        self.need_draw = True
        self.last_draw = 0.0

    def _cap_browse(self):
        from tkinter import filedialog
        d = filedialog.askdirectory(initialdir=self.cap_dir.get() if os.path.isdir(self.cap_dir.get()) else os.getcwd(),
                                    title="Folder of the recordings")
        if d:
            self.cap_dir.set(d)

    # ---- the band on the plots
    def _cap_half(self):
        return (self.cap_plan.out_rate if self.cap_plan is not None else self.cap_rate) / 2.0

    def _cap_hit(self, event, edges=True):
        """What of the band is under the cursor: 'left' / 'right' edge (only if edges), 'move' (inside) or None."""
        v = self.view
        if v is None or self.cap_center is None or event.inaxes not in (self.ax_wf, self.ax_psd) or event.xdata is None:
            return None
        half, c, x = self._cap_half() * 1e-6, self.cap_center * 1e-6, event.xdata
        lo, hi = c - half, c + half
        span = self._view_span()
        tol = 6.0 * span / max(self.ax_psd.bbox.width, 1.0)           # 6 pixels
        if hi - lo < 2 * tol:                                         # a hairline: grab it anywhere near
            return "move" if lo - tol <= x <= hi + tol else None
        if edges:
            if abs(x - lo) <= tol:
                return "left"
            if abs(x - hi) <= tol:
                return "right"
        return "move" if (lo < x < hi and hi - lo <= 0.6 * span) else None   # a band that fills the view: pan instead

    def _cap_spawn(self, event):
        """Shift + left click: the band to record is placed at the click (a new one, or the existing one is moved there),
        with the width it has. During a recording this stops the capture, after a question."""
        v = self.view
        if v is None or event.xdata is None or self._confirming:
            return
        try:
            rp = plan_resampler(v.fs, self.cap_rate)
        except ValueError as e:
            self._cap_say(str(e))
            return
        half = rp.out_rate / 2.0
        c = float(round(min(max(event.xdata * 1e6, v.left * 1e6 + half), v.right * 1e6 - half)))   # inside the spectrum
        if self._cap_recording():
            r = self.recorder
            self._confirming = True
            try:
                ok = self._confirm("Stop the capture?",
                                   f"A recording is in progress.\n\nA new ROI stops the current capture: the file "
                                   f"{os.path.basename(r.files[0])} will be closed.\n\nStop the capture and place the ROI "
                                   f"at {fmt_freq(c, 4).strip()} MHz?")
            finally:
                self._confirming = False
            if not ok or self.recorder is not r or r.state != "recording":
                return
            self._rec_stop()
        if self._cap_set(center=c, source="Manual band"):
            self.sel = None
            self.last_text = 0.0

    def _cap_cursor(self, hit):
        cur = {"left": "sb_h_double_arrow", "right": "sb_h_double_arrow", "move": "fleur"}.get(hit, "")
        if cur != self._cap_cursor_now:
            self._cap_cursor_now = cur
            self.canvas_widget.configure(cursor=cur)

    def _cap_drag_move(self, event):
        v, d = self.view, self.cap_drag
        x_hz = self.ax_psd.transData.inverted().transform((event.x, event.y))[0] * 1e6
        c0, h0 = d["center"], d["half"]
        lo_lim, hi_lim = v.left * 1e6, v.right * 1e6
        min_w = max(500.0, 4.0 * v.df)
        if d["mode"] == "move":
            c = min(max(c0 + (x_hz - d["x0"] * 1e6), lo_lim + h0), hi_lim - h0)
            self._cap_set(center=round(c), source="Manual band")
        else:
            lo, hi = c0 - h0, c0 + h0
            if d["mode"] == "left":
                lo = min(max(x_hz, lo_lim), hi - min_w)
            else:
                hi = max(min(x_hz, hi_lim), lo + min_w)
            self._cap_set(center=round((lo + hi) / 2), rate=round(hi - lo), source="Manual band")
        self._show_tip(event, self._cap_tip_text())

    # ---- selection from the lists
    def _on_table_click(self, event):
        v = self.view
        if v is None:
            return
        row = int(self.rssi_text.index(f"@{event.x},{event.y}").split(".")[0]) - 1      # line 1 is the first row
        if self.mode.get() == "nb":
            if 0 <= row < len(v.nb_order):
                self._select_signal(int(v.nb_order[row]))
        elif 0 <= row < len(v.idx):
            self._select_channel(row)

    def _select_signal(self, slot):
        tr = self.view.tracker
        ident = int(tr.ident[slot])
        src = (f"Signal {ident}: {fmt_freq(tr.center[slot], 4).strip()} MHz, width {tr.width[slot] * 1e-3:.2f} kHz, "
               f"SNR {tr.snr[slot]:.1f} dB, seen {int(tr.hits[slot])} times")
        c = round(float(tr.center[slot]))
        if self._cap_recording():                         # another centre during a recording: confirm, new file
            self._rec_center_request(center=c, name=f"sig{ident}", source=src, sel=("sig", ident))
        elif self._cap_set(center=c, name=f"sig{ident}", source=src):
            self.sel = ("sig", ident)
            self.last_text = 0.0

    def _select_channel(self, k):
        v = self.view
        n = int(v.idx[k]) + 1
        src = f"Channel {n}: {fmt_freq(v.f_in[k], 3).strip()} MHz, width {fmt_freq(v.bw_in[k], 3).strip()} MHz"
        c = round(float(v.f_in[k]))
        if self._cap_recording():
            self._rec_center_request(center=c, rate=float(v.bw_in[k]), name=f"ch{n}", source=src, sel=("ch", k))
        elif self._cap_set(center=c, rate=float(v.bw_in[k]), name=f"ch{n}", source=src):
            self.sel = ("ch", k)
            self.last_text = 0.0

    # ---- recording
    def _cap_spec(self):
        """The band to record as the recorder takes it, or ValueError with the reason."""
        if self.cap_center is None or self.cap_plan is None:
            raise ValueError("There is no valid band: " + (self.cap_msg.cget("text") or "set a centre and a width."))
        d = self.cap_dir.get().strip()
        if not d:
            raise ValueError("Set the folder of the recordings.")
        pattern = self.cap_pattern.get()
        render_stem(pattern, self._cap_ctx())                       # ValueError with the reason if the pattern is wrong
        ident = self.sel[1] if self.sel is not None and self.sel[0] == "sig" else None
        return {"center": self.cap_center, "rate": self.cap_rate, "name": self.cap_name or "band", "directory": d,
                "pattern": pattern, "ident": ident}

    def _cap_context(self):
        """What the sidecar says about the capture: the radio, the analysis, and what the band was taken from."""
        v, cfg, plan, s = self.view, self.cfg, self.plan, self.session
        source = {"kind": "manual"}
        if self.sel is not None and self.sel[0] == "sig":
            tr = v.tracker
            source = {"kind": "signal", "id": int(self.sel[1])}
            hit = np.flatnonzero(tr.used & (tr.ident == self.sel[1]))
            if len(hit):
                k = int(hit[0])
                source["registry"] = {"center_hz": float(tr.center[k]), "width_hz": float(tr.width[k]),
                                      "power_dbfs": float(tr.power[k]), "snr_db": float(tr.snr[k]),
                                      "detections": int(tr.hits[k]), "present": bool(tr.present[k]),
                                      "first_seen_s": float(tr.first_seen[k]), "last_seen_s": float(tr.last_seen[k]),
                                      "duty": float(tr.duty(k))}
        elif self.sel is not None and self.sel[0] == "ch":
            k = int(self.sel[1])
            source = {"kind": "channel", "number": int(v.idx[k]) + 1, "center_hz": float(v.f_in[k]),
                      "width_hz": float(v.bw_in[k])}
        return {
            "capture": {"device": "bladeRF 2.0 (libbladeRF)", "rx_channel": 0, "center_hz": float(v.fc),
                        "sample_rate_hz": float(v.fs), "analog_bandwidth_hz": float(v.bw), "gain_mode": "manual",
                        "gain_db_requested": cfg.gain, "gain_db_reported": v.gain, "stream_format": "SC16_Q11",
                        "stream_start_utc": iso_utc(v.t_enable) if v.t_enable else None,
                        "stream_origin_samples_discarded": v.discard},
            "analysis": {"fft_size": plan.nfft, "window": cfg.window, "overlap": cfg.overlap,
                         "windows_per_row": plan.n_win, "row_ms": plan.row_samples / v.fs * 1e3,
                         "ring_seconds": s.ring.size / v.fs,
                         "narrowband": {"max_width_hz": v.nb_params[0], "guard": v.nb_params[1],
                                        "threshold_db": v.nb_params[2]}},
            "source": source}

    def _same_file(self, spec, r):
        """Would a recording of `spec` continue the file of recorder r? (Only if the whole band is the same.)"""
        keys = ("center", "rate", "name", "directory", "pattern")
        return tuple(spec[k] for k in keys) == tuple(r.spec[k] for k in keys)

    def _rec_start(self):
        """Start, Resume, or Start new file (after a pause with the width, the name or the folder changed)."""
        self._focus_plots()
        s, v, r = self.session, self.view, self.recorder
        if s is None or v is None:
            self._cap_say("The device is not running.")
            return
        if self.cap_center is None and not self.cap_vars["center"].get().strip():
            self._cap_say("Pick a signal or a channel in the list, drag a band, or type a centre first.")
            return
        try:
            try:
                typed = (parse_hz(self.cap_vars["center"].get()), parse_hz(self.cap_vars["rate"].get()))
            except ValueError as e:
                raise ValueError(f"The fields do not hold a valid band: {e}") from None
            if typed != (self.cap_center, self.cap_rate):
                self._cap_fields_changed()
                if typed != (self.cap_center, self.cap_rate):
                    raise ValueError("The fields hold a band that is not valid: " + (self.cap_msg.cget("text") or "fix them"))
            spec = self._cap_spec()
            if r is not None and r.state == "paused" and self._same_file(spec, r):
                r.resume()
                self._rec_ui()
                return
        except ValueError as e:
            self._cap_say(str(e))
            return
        self._rec_switch(spec)                            # a new file; the paused one (if any) is closed once this one runs

    def _rec_switch(self, spec):
        """Open a new file for spec and close the one that was being written (if any). True if the new recording runs."""
        s, v, old = self.session, self.view, self.recorder
        try:
            rec = BandRecorder(s.ring, v.fs, v.fc, s.wall0, spec, self.ui_q, self._cap_context(), v.discard)
            rec.start()                                   # may fail (folder): the earlier recording is left as it is then
        except (ValueError, OSError) as e:
            self._cap_say(str(e))
            return False
        if old is not None:
            self._rec_stop()
        self.recorder = rec
        self._cap_say("")
        self._rec_ui()
        return True

    def _confirm(self, title, message):
        """A modal question; True = go on. (A method of its own, so that it can be replaced in tests.)"""
        from tkinter import messagebox
        return messagebox.askokcancel(title, message, icon="question", parent=self.root)

    def _cap_restore(self):
        """Back to the band that is being recorded (the request was cancelled or could not be carried out)."""
        r = self.recorder
        if r is None:
            return
        self.cap_center, self.cap_rate, self.cap_name = r.spec["center"], r.spec["rate"], r.spec["name"]
        self._cap_check(self.cap_center, self.cap_rate)
        self._cap_sync()
        self._cap_overlay()
        self._cap_preview()

    def _rec_center_request(self, center, rate=None, name=None, source=None, sel=None):
        """Another centre was asked for while a file is being written: confirm in a pop-up, then close this file and open a
        new one for the new band. The recording goes on while the question is open."""
        r = self.recorder
        if r is None or self._confirming:
            return False
        old_c, new_c = r.spec["center"], float(center)
        if new_c == old_c and (rate is None or float(rate) == r.spec["rate"]):
            self._cap_restore()
            return False
        what = f"from {fmt_freq(old_c, 4).strip()} MHz to {fmt_freq(new_c, 4).strip()} MHz"
        self._confirming = True
        try:
            ok = self._confirm("Change the centre?", f"A recording is in progress.\n\nChange the centre {what}?\n\n"
                                                      f"The file {os.path.basename(r.files[0])} will be closed and a new one "
                                                      f"opened for the new band.")
        finally:
            self._confirming = False
        if not ok or self.recorder is not r or r.state != "recording":      # cancelled, or the recording ended meanwhile
            self._cap_restore()
            return False
        if not self._cap_set(center=new_c, rate=rate, name=name, source=source):
            self._cap_restore()
            return False
        if sel is not None:
            self.sel = sel
            self.last_text = 0.0
        try:
            spec = self._cap_spec()
        except ValueError as e:
            self._cap_say(str(e))
            self._cap_restore()
            return False
        if not self._rec_switch(spec):
            self._cap_restore()
            return False
        return True

    def _rec_pause(self):
        self._focus_plots()
        if self.recorder is not None:
            self.recorder.pause()
        self._rec_ui()

    def _rec_stop_clicked(self):
        self._focus_plots()
        self._rec_stop()

    def _rec_stop(self):
        """Close the file (Stop)."""
        r = self.recorder
        self.recorder = None
        if r is not None:
            r.stop()
            if r.files:
                more = f" and {len(r.files_done)} earlier" if r.files_done else ""
                self.rec_lbl.configure(text=f"Saved {os.path.basename(r.files[0])} (+ .json){more}  "
                                            f"({r.out_samples / r.plan.out_rate:.1f} s, {r.bytes / 1e6:.2f} MB, "
                                            f"{len(r.segments)} segment{'s' if len(r.segments) != 1 else ''}, lost blocks {r.lost})")
        self._rec_ui()

    def _set_cap_entries(self):
        rec = self._cap_recording()
        for key, e in self.cap_entries.items():            # the centre can be moved while recording, the rest is the file's
            e.configure(state="normal" if key == "center" or not rec else "disabled")
        self.cap_dir_entry.configure(state="disabled" if rec else "normal")
        self.cap_pattern_entry.configure(state="disabled" if rec else "normal")
        self.cap_browse_btn.configure(state="disabled" if rec else "normal")

    def _rec_ui(self):
        if not hasattr(self, "btn_stop"):
            return
        r = self.recorder
        st = r.state if r is not None else "idle"
        if st == "paused":
            try:
                new_file = not self._same_file(self._cap_spec(), r)
            except ValueError:
                new_file = False
            first = "Start new file" if new_file else "Resume"
        else:
            first = "Start"
        self._cap_preview()                              # the file just created (or closed) changes which name is free next
        self.btn_rec.configure(text=first, state="disabled" if st == "recording" else "normal")
        self.btn_pause.configure(state="normal" if st == "recording" else "disabled")
        self.btn_stop.configure(state="normal" if st in ("recording", "paused") else "disabled")
        self._set_cap_entries()

    def _rec_update(self, st):
        r = self.recorder
        if r is None:
            return
        rt = "" if st["rt"] != st["rt"] else f", x{st['rt']:.1f} real time"
        self.rec_lbl.configure(text=f"{st['state']} at {fmt_freq(st['center'], 4).strip()} MHz: {st['seconds']:.1f} s, "
                                    f"{st['samples'] / 1e3:.0f}k samples, {st['bytes'] / 1e6:.2f} MB, "
                                    f"lost blocks {st['lost']}{rt}")
        if st.get("warning"):
            self._cap_say(st["warning"], "warn")
        if st["state"] == "error":
            self._cap_say(f"Recording stopped: {st['error']}")
            self.recorder = None
            self._rec_ui()

    def _build_status(self, root):
        bar = ttk.Frame(root, style="Bar.TFrame", padding=(S2, S1 + 1))
        bar.grid(row=1, column=0, columnspan=3, sticky="ew")
        self.state_dot = tk.Label(bar, text="\u25CF", bg=C["bg"], fg=C["muted"], font=(self.f_ui, 9))
        self.state_dot.grid(row=0, column=0)
        self.state_lbl = ttk.Label(bar, text="Starting", style="Bar.TLabel", font=(self.f_ui, 9, "bold"))
        self.state_lbl.grid(row=0, column=1, padx=(S1, S2))
        self.stat_lbls, self.stat_seps = {}, {}
        col = 2
        for key in ("center", "rate", "res", "gain", "rows", "ring", "dropped", "chan", "rec", "zoom"):
            sep = ttk.Separator(bar, orient="vertical", style="Bar.TSeparator")
            sep.grid(row=0, column=col, sticky="ns", padx=S2)
            lbl = ttk.Label(bar, text="", style="BarMuted.TLabel")
            lbl.grid(row=0, column=col + 1)
            self.stat_seps[key], self.stat_lbls[key] = sep, lbl
            sep.grid_remove()
            lbl.grid_remove()
            col += 2
        Tooltip(self.stat_lbls["dropped"], self._dropped_text)
        Tooltip(self.stat_lbls["rows"], self._perf_text)
        bar.columnconfigure(col, weight=1)
        self.hover_lbl = ttk.Label(bar, text="", style="BarMuted.TLabel", font=(self.f_mono, 9))
        self.hover_lbl.grid(row=0, column=col + 1, sticky="e")

    def _perf_text(self):
        return self.perflog.tooltip_text() if self.perflog else "No numbers yet."

    def _finish_perf(self, session):
        """Summary of a finished run: stderr and the metrics file."""
        pl = self.perflog
        self.perflog = None
        if pl is None or pl.tot["intervals"] == 0:
            return
        lines = pl.summary_lines(session)
        for line in lines:
            info(line)
        if self.metrics_fh:
            for line in lines:
                self.metrics_fh.write("# " + line + "\n")
            self.metrics_fh.flush()

    def _dropped_text(self):
        s = self.session
        if s is None:
            return "No rows skipped."
        return ("Rows skipped:\n"
                f"  overwritten in the ring before they were processed: {s.dropped_late}\n"
                f"  overwritten while they were being processed: {s.dropped_during}\n"
                f"  the window could not take them: {s.dropped_gui}")

    def _stat(self, key, text, color=None):
        """Show a status field, or hide it (with its separator) when there is nothing to say."""
        lbl, sep = self.stat_lbls[key], self.stat_seps[key]
        if text:
            lbl.configure(text=text, foreground=color or C["muted"])
            lbl.grid()
            sep.grid()
        else:
            lbl.grid_remove()
            sep.grid_remove()

    # --------------------------------------------------------------- fields

    def _set_fields(self, f):
        """Fill the controls from validate()-style strings."""
        for name in ("center", "rate", "bandwidth", "gain", "overlap", "rows", "nb_max_width", "nb_guard", "nb_thr",
                     "ring_dur"):
            self.vars[name].set(f[name])
        self.vars["spans_mode"].set({v: k for k, v in SPANS_MODES.items()}[f["spans_mode"]])
        self.vars["window"].set(f["window"])
        if f["fft_mode"] == "df":                            # what the command line gave is the guide of its pair
            self.fft_primary, self.fft_guide = "df", parse_hz(f["fft_value"])
        else:
            self.fft_primary = "n"
            self.vars["fft_n"].set(f["fft_value"])
        if f["avr_mode"] == "dur":
            self.avr_primary, self.dur_guide = "dur", float(f["avr_value"])
        else:
            self.avr_primary = "n"
            self.vars["avr_n"].set(f["avr_value"])
        self._spectrum_sync()
        self._setting_spans = True
        self.spans_text.delete("1.0", "end")
        self.spans_text.insert("1.0", f["spans"])
        self.spans_text.edit_modified(False)
        self._setting_spans = False

    def _fields(self):
        """Current control values as validate()-style strings. The paired fields are always consistent: the size and the
        number of windows are passed, the resolution and the time follow from them."""
        d = {n: self.vars[n].get() for n in ("center", "rate", "bandwidth", "gain", "overlap", "rows",
                                             "nb_max_width", "nb_guard", "nb_thr", "ring_dur",
                                             "fft_n", "fft_df", "avr_n", "avr_dur")}
        d["spans_mode"] = SPANS_MODES[self.vars["spans_mode"].get()]
        d["fft_mode"], d["fft_value"] = "fft_n", d["fft_n"]
        d["avr_mode"], d["avr_value"] = "n", d["avr_n"]
        d["window"] = self.vars["window"].get()
        d["spans"] = self.spans_text.get("1.0", "end").strip()
        return d

    def _on_spans_modified(self, _e):
        if getattr(self, "_setting_spans", False) or not self.spans_text.edit_modified():
            return
        self.spans_text.edit_modified(False)
        self._on_edit("spans")

    def _show_message(self, text="", kind="invalid"):
        """One place for feedback, next to the primary action. kind: invalid | error (red), warn (amber)."""
        if not text:
            self.msg_kind = None
            self.note.configure(text=NOTE_DEFAULT, foreground=C["muted"])
            return
        self.msg_kind = kind
        self.note.configure(text=text, foreground=C["amber"] if kind == "warn" else C["error"])

    def _on_edit(self, name=None):
        """The user changed a control: drop its error mark, then mark the controls that differ from the run."""
        if name:
            w = self.entries.get(name)
            if name in ("spans", "spans_mode"):
                self.spans_invalid = False
                if self.msg_kind == "invalid":
                    self._show_message("")
            elif w is not None and not isinstance(w, ttk.Combobox) and str(w.cget("style")) == "Invalid.TEntry":
                w.configure(style="TEntry")
                if self.msg_kind == "invalid":
                    self._show_message("")
        if not self.applied:
            return
        cur = self._fields()
        for key in RESTART_FIELDS:
            w = self.entries.get(key)
            changed = cur[key] != self.applied.get(key)
            if key == "spans":
                col = C["error"] if self.spans_invalid else (C["accent"] if changed else C["line"])
                self.spans_text.configure(highlightbackground=col)
            elif isinstance(w, ttk.Combobox):
                w.configure(style="Dirty.TCombobox" if changed else "TCombobox")
            elif w is not None and str(w.cget("style")) != "Invalid.TEntry":
                w.configure(style="Dirty.TEntry" if changed else "TEntry")

    def _mark_invalid(self, field, message):
        name = FIELD_OF.get(field, field)
        w = self.entries.get(name)
        if name == "spans":
            self.spans_invalid = True
            self.spans_text.configure(highlightbackground=C["error"])
        elif w is not None and not isinstance(w, ttk.Combobox):
            w.configure(style="Invalid.TEntry")
        self._show_message(message, "invalid")

    def _clear_invalid(self):
        if self.msg_kind == "invalid":
            self._show_message("")
        self.spans_invalid = False
        for name in RESTART_FIELDS + ("rows", "nb_max_width", "nb_guard", "nb_thr"):
            w = self.entries.get(name)
            if name == "spans":
                continue
            if w is not None and not isinstance(w, ttk.Combobox) and str(w.cget("style")) == "Invalid.TEntry":
                w.configure(style="TEntry")
        self._on_edit()

    def _validate_quiet(self, name):
        """Check the form when a control loses focus; only the error of that control is shown."""
        try:
            validate(self._fields())
        except ConfigError as e:
            if FIELD_OF.get(e.field) == name:
                self._mark_invalid(e.field, str(e))
            return
        self._clear_invalid()

    # ----------------------------------------------------------- display opts

    def _range_mode(self):
        state = "disabled" if self.auto_range.get() else "normal"
        self.e_rmin.configure(state=state)
        self.e_rmax.configure(state=state)
        self._display_changed()

    def _display_changed(self, refocus=False):
        if refocus:
            self._focus_plots()
        if not self.auto_range.get():
            try:
                lo, hi = float(self.rmin.get()), float(self.rmax.get())
                if hi > lo:
                    self.vmin, self.vmax = lo, hi
            except ValueError:
                pass
        self.need_draw = True
        self.last_text = 0.0
        self.last_draw = 0.0

    def _rows_changed(self):
        if self.view is None:
            return
        try:
            rows = int(self.vars["rows"].get().strip())
        except ValueError:
            self._mark_invalid("rows", "Rows must be a whole number")
            return
        if not ROWS_MIN <= rows <= ROWS_MAX:
            self._mark_invalid("rows", f"Rows must be {ROWS_MIN}..{ROWS_MAX}")
            return
        self._clear_invalid()
        if rows != self.view.wf.rows:
            self.view.wf = self.view.wf.resized(rows)
            self.im.set_data(self.view.wf.view())
            self._fit_waterfall_extent()
            self.need_draw = True

    # ---------------------------------------------------------- run control

    def apply_and_restart(self):
        for name in ("fft_n", "fft_df", "avr_n", "avr_dur"):    # a value typed but not yet left counts
            self._pair_commit(name)
        if any(str(self.entries[n].cget("style")) == "Invalid.TEntry" for n in ("fft_n", "fft_df", "avr_n", "avr_dur")):
            return
        self._spectrum_sync()
        try:
            cfg, plan = validate(self._fields())
        except ConfigError as e:
            self._mark_invalid(e.field, str(e))
            return
        self._clear_invalid()
        self.applied = {n: v for n, v in self._fields().items() if n in RESTART_FIELDS}
        self._on_edit()
        self._stop_session(wait=True)
        self.cfg, self.plan = cfg, plan
        for note in plan.notes:
            info(note)
        self.view = None
        self.rows_total = 0
        self.period = 0.0
        self.last_row_t = None
        self.hover_idx = -1
        self.paused = False
        self.btn_view_pause.configure(text="Pause")
        self._reset_artists()
        self._set_state("Starting", C["muted"])
        self.empty_text.set_text("Waiting for the device")
        self._update_status(force=True)
        try:
            self.session = Session(cfg, plan, self.ui_q)
        except MemoryError:
            self.session = None
            self._set_state("Error", C["error"])
            self._show_message("Not enough memory for the buffer; shorten Buffer or lower the rate", "error")
            self.empty_text.set_text("No data. Shorten Buffer or lower the rate, then press Apply and restart.")
            return
        self.session.start()
        self.need_draw = True
        self._focus_plots()

    def _stop_session(self, wait):
        s = self.session
        if s is None:
            return
        self.session = None
        self._rec_stop()
        self._finish_perf(s)
        if not s.shutdown(5.0):
            self._show_message("The device did not stop in 5 s; if it fails to open, unplug and reconnect it", "warn")
        while True:                       # drop messages of the old session
            try:
                self.ui_q.get_nowait()
            except queue.Empty:
                break

    def toggle_pause(self):
        s = self.session
        if s is None or self.view is None:
            return
        self.paused = not self.paused
        if self.paused:
            s.paused.set()
        else:
            s.paused.clear()
        self.btn_view_pause.configure(text="Resume" if self.paused else "Pause")
        self._set_state("Paused" if self.paused else "Running", C["amber"] if self.paused else C["ok"])
        self._focus_plots()

    def _focus_plots(self):
        """Keyboard focus to the plots, so that Space keeps working after a click on a control."""
        self.canvas_widget.focus_set()

    def _on_space(self, event):
        """Space pauses, unless the focused widget uses Space itself (typing, buttons, check boxes)."""
        try:
            w = self.root.focus_get()
        except KeyError:                       # focus is in a pop-down list
            return
        if w is not None and w.winfo_class() in TEXT_INPUT_CLASSES:
            return
        self.toggle_pause()

    def _set_state(self, text, color):
        self.state = text
        self.state_lbl.configure(text=text)
        self.state_dot.configure(fg=color)

    def on_close(self):
        if self.closing:
            return
        self.closing = True
        self._stop_session(wait=True)
        if self.metrics_fh:
            self.metrics_fh.close()
            self.metrics_fh = None
        self.root.destroy()

    # --------------------------------------------------------------- polling

    def poll(self):
        if self.closing:
            return
        try:
            if self.sigint:
                self.on_close()
                return
            self._poll()
        except Exception as e:            # keep the window alive and say what went wrong
            import traceback
            traceback.print_exc()
            self._show_message(f"Window error: {type(e).__name__}: {e}", "error")
        if not self.closing:
            self.root.after(POLL_MS, self.poll)

    def _poll(self):
        s = self.session
        t_poll = time.perf_counter()
        qsize = self.ui_q.qsize()
        self._draw_in_poll = 0.0
        rows_seen = 0
        try:
            while True:
                kind, payload = self.ui_q.get_nowait()
                if kind == "ready":
                    self.on_ready(payload)
                elif kind == "row":
                    rows_seen += 1
                    ri, start = payload
                    if s is not None and self.view is not None:
                        self.on_row(s, ri, start)
                    elif s is not None:
                        s.release_row(ri)
                elif kind == "perf_cap" and self.perflog is not None:
                    self.perflog.cap = payload
                elif kind == "perf_proc" and self.perflog is not None:
                    self.perflog.proc = payload
                elif kind == "rec":
                    self._rec_update(payload)
                elif kind == "error":
                    self.on_error(payload)
        except queue.Empty:
            pass
        now = time.perf_counter()
        if (self.need_draw or (self.view is not None and rows_seen)) and \
                now - self.last_draw >= max(DRAW_MIN_S, 2.0 * self.draw_cost):
            if self.view is not None:
                self.refresh_plots()
            self.last_draw = now
            t_draw = time.perf_counter()
            self.canvas.draw()
            self.draw_cost = time.perf_counter() - t_draw
            self._draw_in_poll = self.draw_cost
            if self.perflog is not None:
                self.perflog.add_draw(self.draw_cost)
            self.need_draw = False
        if self.view is not None and now - self.last_text >= TEXT_MIN_S:
            self.refresh_text()
            self.last_text = now
        if now - self.last_status >= STATUS_MIN_S:
            self._update_status()
            self.last_status = now
        if self.perflog is not None and s is not None:
            self.perflog.add_poll(time.perf_counter() - t_poll - self._draw_in_poll, qsize)
            self.perflog.tick(time.perf_counter(), s)

    def on_error(self, text):
        self._set_state("Error", C["error"])
        self._show_message(text, "error")
        self.empty_text.set_text("No data. Check the device, then press Apply and restart.")
        self._rssi_message("No data.\nCheck the device, then press\nApply and restart.")
        self.need_draw = True
        if self.session is not None:
            self.session.stop.set()

    # ------------------------------------------------------ display model

    def _reset_artists(self):
        for art in (self.im, self.psd_line, self.levels, self.thr_line, self.floor_line, self.thr_text, *self.band_pcs,
                    *self.hl_patches, *self.nb_lcs, *self.cap_patches):
            if art is not None:
                try:
                    art.remove()
                except (ValueError, NotImplementedError, AttributeError):
                    pass
        self.im = self.psd_line = self.levels = None
        self.thr_line = self.floor_line = self.thr_text = None
        self.band_pcs, self.hl_patches, self.nb_lcs, self.cap_patches = [], [], [], []
        self.cax.clear()
        self.cax.set_facecolor(C["panel"])
        self.cax.set_visible(False)
        self._rssi_message("Waiting for the device")
        self.chan_summary.configure(text="")
        for key in self.stat_lbls:
            self._stat(key, "")
        self.hover_lbl.configure(text="")

    def on_ready(self, ready):
        from matplotlib.collections import LineCollection, PolyCollection
        from matplotlib.colors import to_rgba
        from matplotlib.patches import Rectangle

        plan, cfg = self.plan, self.cfg
        self.zoom_x = self.zoom_y = None
        self._drag = None
        fs, fc = ready["fs"], ready["fc"]
        self.session_no += 1
        self.perflog = PerfLog(self.metrics_fh, self.session_no, fs, plan.row_samples / fs,
                               SYNC_BUFFERS * SYNC_BUF_SAMPLES / fs * 1e3)
        self.perflog.event("session_start", f"center {fc:.0f} Hz, rate {fs:.0f}, bandwidth {ready['bw']:.0f}, gain "
                                            f"{ready['gain']}, FFT {plan.nfft}, row {plan.row_samples / fs * 1e3:.1f} ms, "
                                            f"ring {plan.ring_rows} rows")
        n = plan.nfft
        df = fs / n
        left = fc - (n // 2) * df - df / 2
        right = left + n * df
        fgrid = (fc + (np.arange(n) - n // 2) * df) * 1e-6
        idx, partial, outside = classify_channels(plan.ch_f, plan.ch_bw, fc, fs)
        f_in, bw_in = plan.ch_f[idx], plan.ch_bw[idx]
        v = _View()
        v.fs, v.fc, v.df, v.n, v.bw, v.gain = fs, fc, df, n, ready["bw"], ready["gain"]
        v.row_dur = plan.row_samples / fs      # seconds of signal in one waterfall row
        v.discard = int(ready.get("discarded", 0))   # start-up samples before ring sample 0: time counts from the stream start
        v.t_enable = ready.get("t_enable")
        v.left, v.right = left * 1e-6, right * 1e-6
        v.fgrid = fgrid
        v.idx, v.f_in, v.bw_in = idx, f_in, bw_in
        v.n_partial, v.n_outside, v.n_total = partial, outside, len(plan.ch_f)
        v.integ = BandIntegrator(n, fs, fc, f_in, bw_in)
        v.nb = NbDetector(n, fs, fc, cfg.window)
        v.nb.set_params(cfg.nb_max_width, cfg.nb_guard, cfg.nb_thr)
        v.nb_params = (cfg.nb_max_width, cfg.nb_guard, cfg.nb_thr)
        v.tracker = NbTracker()
        v.t_win = 0.0                      # start of the displayed PSD window, s from samples
        v.nb_order = np.empty(0, np.intp)  # track slots in table order
        v.rawlog = RawLog()
        v.log_bw = 10.0 * np.log10(np.maximum(bw_in, 1.0))
        v.latest = np.zeros(n)
        v.mean_dens = np.zeros(n)
        v.mean_p = np.zeros(n)
        v.dens_db = np.empty(n)
        v.tmp_row = np.empty(n)
        v.rssi = np.full(len(idx), np.nan)
        v.wf = WaterfallBuffer(self._rows_value(cfg.rows), n)
        v.have_row = False
        self.view = v

        self.cax.set_visible(True)
        cmap = self._cmap(self.cmap.get())
        self.im = self.ax_wf.imshow(v.wf.view(), aspect="auto", interpolation="nearest", cmap=cmap,
                                    extent=(v.left, v.right, 1.0, 0.0), vmin=-120, vmax=-60, zorder=1)
        self.cbar = self.fig.colorbar(self.im, cax=self.cax)
        self.cax.set_title("dBFS/Hz", loc="left", fontsize=8, color=C["muted"], pad=4)
        self.cbar.ax.tick_params(colors=C["muted"], labelsize=8, length=3)
        self.cbar.outline.set_edgecolor(C["line"])
        (self.psd_line,) = self.ax_psd.plot(fgrid, np.full(n, np.nan), lw=1.0, color=C["psd"], zorder=2)
        # what the narrowband detector uses: the noise floor (dotted) and floor + threshold (dashed)
        (self.floor_line,) = self.ax_psd.plot(fgrid, np.full(n, np.nan), lw=0.9, ls=":", color=C["muted"], zorder=2)
        (self.thr_line,) = self.ax_psd.plot(fgrid, np.full(n, np.nan), lw=1.1, ls=(0, (5, 3)), color=C["ok"], zorder=3)
        self.thr_text = self.ax_psd.text(0.995, 0.97, "", transform=self.ax_psd.transAxes, ha="right", va="top",
                                         fontsize=8, color=C["ok"], zorder=6)
        self.ax_wf.set_xlim(v.left, v.right)
        self.ax_psd.set_ylim(-120, -60)

        m = len(idx)
        if m:
            x0, x1 = (f_in - bw_in / 2) * 1e-6, (f_in + bw_in / 2) * 1e-6
            verts = np.empty((m, 4, 2))
            verts[:, :, 1] = np.array([0, 1, 1, 0])
            verts[:, 0, 0] = verts[:, 1, 0] = x0
            verts[:, 2, 0] = verts[:, 3, 0] = x1
            face = np.array([to_rgba(C["band"], 0.10 if i % 2 else 0.17) for i in range(m)])
            edge = to_rgba(C["band"], 0.55)
            for ax in (self.ax_wf, self.ax_psd):
                pc = PolyCollection(verts, transform=ax.get_xaxis_transform(), facecolors=face, edgecolors=edge,
                                    linewidths=0.8, zorder=3)
                ax.add_artist(pc)
                self.band_pcs.append(pc)
            segs = np.zeros((m, 2, 2))
            segs[:, 0, 0], segs[:, 1, 0] = x0, x1
            segs[:, :, 1] = -200.0
            self.levels = LineCollection(segs, colors=C["amber"], linewidths=2.4, zorder=4, capstyle="butt")
            self.ax_psd.add_collection(self.levels, autolim=False)
            v.segs = segs
        v.nb_segs = np.full((NB_TRACK_LIMIT, 2, 2), np.nan)     # one vertical marker per track slot
        v.nb_segs[:, 0, 1] = 0.0
        v.nb_segs[:, 1, 1] = 1.0
        for ax in (self.ax_wf, self.ax_psd):
            lc = LineCollection(v.nb_segs, colors=C["accent"], linewidths=1.3, zorder=4,
                                transform=ax.get_xaxis_transform())
            ax.add_collection(lc, autolim=False)
            self.nb_lcs.append(lc)
        for ax in (self.ax_wf, self.ax_psd):
            r = Rectangle((0, 0), 0, 1, transform=ax.get_xaxis_transform(), facecolor=C["amber"], alpha=0.22,
                          edgecolor=C["amber"], linewidth=1.2, zorder=5, visible=False)
            ax.add_artist(r)
            self.hl_patches.append(r)

        for ax in (self.ax_wf, self.ax_psd):                    # the band to record: one object on both plots
            r = Rectangle((0, 0), 0, 1, transform=ax.get_xaxis_transform(), facecolor=C["cap"], alpha=0.25,
                          edgecolor=C["cap"], linewidth=1.5, zorder=7, visible=False)
            ax.add_artist(r)
            self.cap_patches.append(r)
        if self.cap_center is not None and not self._cap_check(self.cap_center, self.cap_rate)[0]:
            self.cap_center = None                              # a band set earlier stays only if it still fits
            self.cap_plan = None
        self._cap_sync()
        self._cap_overlay()
        self._fit_waterfall_extent()
        self.empty_text.set_text("")
        self._build_rssi_head()
        self._set_state("Paused" if self.paused else "Running", C["amber"] if self.paused else C["ok"])
        self._apply_overlay_visibility()
        gain_txt = "?" if ready["gain"] is None else f"{ready['gain']} dB"
        if ready["gain"] is not None and ready["gain"] != cfg.gain:
            self._show_message(f"Gain: the device reports {ready['gain']} dB, requested {cfg.gain} dB", "warn")
        self.stat_gain_warn = ready["gain"] is not None and ready["gain"] != cfg.gain
        info(f"capture: center {fmt_freq(fc, 3).strip()} MHz, rate {fs * 1e-6:.3f} MS/s, bandwidth "
             f"{ready['bw'] * 1e-6:.3f} MHz, gain {gain_txt}, {m} of {v.n_total} channels inside")
        self._update_chan_summary()
        self._update_status(force=True)
        self.need_draw = True

    @staticmethod
    def _rows_value(default):
        return default

    def _fit_waterfall_extent(self):
        v = self.view
        if v is None or self.im is None:
            return
        height = v.wf.rows * v.row_dur          # exact row duration: the Age axis stays fixed between restarts
        self.im.set_extent((v.left, v.right, height, 0.0))
        self.ax_wf.set_ylim(height, 0.0)
        self.ax_wf.set_xlim(*(self.zoom_x or (v.left, v.right)))

    def _apply_overlay_visibility(self):
        on = self.overlay.get()
        nb = self.mode.get() == "nb"
        for art in (*self.band_pcs, self.levels):
            if art is not None:
                art.set_visible(on and not nb)
        for lc in self.nb_lcs:
            lc.set_visible(on and nb)

    def _update_chan_summary(self):
        v = self.view
        if v is None:
            return
        if v.n_total == 0:
            text = "No channel ranges set. Add some to get an RSSI."
        else:
            text = f"{len(v.idx)} of {v.n_total} channels are inside the capture band."
            if v.n_partial:
                text += f" {v.n_partial} partly outside are skipped."
            if v.n_outside:
                text += f" {v.n_outside} outside."
        self.chan_summary.configure(text=text)

    def _build_rssi_head(self):
        v = self.view
        if self.mode.get() == "nb" and v is not None:
            tr = v.tracker
            self.rssi_caption.configure(text=f"window start t = {v.t_win:.3f} s,  {int((tr.present & tr.confirmed).sum())} "
                                             f"present of {int(tr.visible.sum())} listed, "
                                             f"raw {v.rawlog.count_since(v.t_win - REG_RAW_WINDOW_S)}/{REG_RAW_WINDOW_S:g} s")
        else:
            self.rssi_caption.configure(text=f"dBFS, PSD: {self.psd_source.get().lower()}")

    # ------------------------------------------------------------- data flow

    def on_row(self, s, ri, start):
        v = self.view
        now = time.perf_counter()
        np.copyto(v.latest, s.rowbufs[ri])
        s.release_row(ri)
        np.divide(v.latest, v.df, out=v.tmp_row)
        np.maximum(v.tmp_row, 1e-30, out=v.tmp_row)
        np.log10(v.tmp_row, out=v.tmp_row)
        v.tmp_row *= 10.0
        v.wf.push(v.tmp_row, start)
        v.have_row = True
        self.rows_total += 1
        if self.last_row_t is not None:
            dt = now - self.last_row_t
            self.period = dt if self.period == 0.0 else 0.9 * self.period + 0.1 * dt
        self.last_row_t = now

    def _displayed_power(self):
        """Per-bin power of the representative PSD (latest row or mean of the waterfall rows)."""
        v = self.view
        if self.psd_source.get() == PSD_SOURCES[0]:
            return v.latest
        v.wf.mean_linear_into(v.mean_dens)
        np.multiply(v.mean_dens, v.df, out=v.mean_p)
        return v.mean_p

    def refresh_plots(self):
        v = self.view
        if not v.have_row:
            return
        p = self._displayed_power()
        np.divide(p, v.df, out=v.dens_db)
        np.maximum(v.dens_db, 1e-30, out=v.dens_db)
        np.log10(v.dens_db, out=v.dens_db)
        v.dens_db *= 10.0
        self.psd_line.set_ydata(v.dens_db)
        wf = v.wf.view()
        self.im.set_data(wf)
        if self.auto_range.get():
            vals = wf[:v.wf.count]
            lo, hi = np.percentile(vals, (1.0, 99.7))
            hi = max(hi, lo + 3.0)
            if self.vmin is None or not np.isfinite(self.vmin):
                self.vmin, self.vmax = lo, hi
            else:
                self.vmin += 0.3 * (lo - self.vmin)
                self.vmax += 0.3 * (hi - self.vmax)
        if self.vmin is not None:
            self.im.set_clim(self.vmin, self.vmax)
            peak = float(np.max(v.dens_db))
            if self.zoom_y is None:
                self.ax_psd.set_ylim(self.vmin - 3.0, max(self.vmax, peak) + 3.0)
            else:
                self.ax_psd.set_ylim(*self.zoom_y)
        self.im.set_cmap(self._cmap(self.cmap.get()))
        self._detect(p)
        self._update_threshold()
        if self.mode.get() == "rssi" and len(v.idx):
            v.rssi[:] = v.integ.db(p)
            lvl = v.rssi - v.log_bw
            v.segs[:, 0, 1] = lvl
            v.segs[:, 1, 1] = lvl
            self.levels.set_segments(v.segs)
        self._apply_overlay_visibility()
        self._update_highlight()

    def _update_threshold(self):
        """The detection threshold on the PSD: floor + Threshold (dashed), the floor (dotted). Shown only when the detector
        works (a bin must fit into Max width) and the check box is on."""
        v = self.view
        on = bool(self.show_thr.get() and v.nb.possible)
        for art in (self.thr_line, self.floor_line, self.thr_text):
            art.set_visible(on)
        if on:
            self.thr_line.set_ydata(v.nb.thr_db)
            self.floor_line.set_ydata(v.nb.floor)
            self.thr_text.set_text(f"threshold = floor + {v.nb.thr:g} dB")

    def _cmap(self, name):
        import matplotlib
        if getattr(self, "_cmap_name", None) != name:
            cm = matplotlib.colormaps[name].copy()
            cm.set_bad(cm(0.0))              # rows that have no data yet look like the noise floor
            self._cmap_obj, self._cmap_name = cm, name
        return self._cmap_obj

    def refresh_text(self):
        v = self.view
        t = self.rssi_text
        top, _ = t.yview()
        t.configure(state="normal")
        t.delete("1.0", "end")
        self._set_head("")                      # set again below when there are rows
        self._build_rssi_head()
        if self.mode.get() == "nb":
            self._fill_nb(t)
        elif v.n_total == 0:
            t.insert("end", "No channel ranges set.\nEnter ranges in the Channels group,\nthen press Apply and restart.",
                     "empty")
        elif not len(v.idx):
            t.insert("end", "No channel lies inside the capture band.\nChange the center or the rate.", "empty")
        else:
            self._set_head(f"{'Ch':>4} {'Center, MHz':>12} {'Width, MHz':>11} {'RSSI, dBFS':>11}")
            lines = []
            for k in range(len(v.idx)):
                val = v.rssi[k]
                txt = "      n/a" if not np.isfinite(val) or not v.have_row else f"{val:11.2f}"
                lines.append(f"{int(v.idx[k]) + 1:>4} {fmt_freq(v.f_in[k], 3):>12} {fmt_freq(v.bw_in[k], 3):>11} {txt}\n")
            t.insert("end", "".join(lines))
            if self.sel is not None and self.sel[0] == "ch":
                t.tag_add("selected", f"{self.sel[1] + 1}.0", f"{self.sel[1] + 1}.end")
            if self.hover_idx >= 0:
                line = self.hover_idx + 1
                t.tag_add("hover", f"{line}.0", f"{line}.end")
        t.yview_moveto(top)
        t.configure(state="disabled")

    def _fill_nb(self, t):
        v = self.view
        tr = v.tracker
        if not v.nb.possible:
            t.insert("end", f"Resolution {v.df / 1e3:.3f} kHz is coarser than\nMax width {fmt_hz_short(v.nb.max_width)}Hz: nothing "
                            f"can be detected.\nUse a finer resolution or a larger Max width.", "empty")
            v.nb_order = np.empty(0, np.intp)
            return
        vis = np.flatnonzero(tr.visible)
        if not len(vis):
            t.insert("end", f"No narrowband signal listed yet (a signal is listed after {REG_MIN_HITS} detections).\n"
                            f"Threshold: {v.nb.thr:g} dB above the local noise floor.", "empty")
            v.nb_order = np.empty(0, np.intp)
            return
        order = vis[np.argsort(tr.center[vis])]
        v.nb_order = order
        self._set_head(f"{'Id':>4} {'Center, MHz':>11} {'W, kHz':>7} {'P, dBFS':>8} {'SNR, dB':>7} {'Pres':>4} "
                       f"{'Seen, s':>9} {'N':>5}")
        for pos, k in enumerate(order):
            line = (f"{int(tr.ident[k]):>4} {fmt_freq(tr.center[k], 4):>11} {tr.width[k] * 1e-3:7.2f} "
                    f"{tr.power[k]:8.2f} {tr.snr[k]:7.1f} {'yes' if tr.present[k] else 'no':>4} "
                    f"{tr.last_seen[k]:9.3f} {int(tr.hits[k]):>5}\n")
            t.insert("end", line, () if tr.present[k] else "absent")
        if self.sel is not None and self.sel[0] == "sig":
            pos = np.flatnonzero(tr.ident[order] == self.sel[1])
            if len(pos):
                t.tag_add("selected", f"{int(pos[0]) + 1}.0", f"{int(pos[0]) + 1}.end")
        if self.hover_idx >= 0:
            pos = np.flatnonzero(order == self.hover_idx)
            if len(pos):
                line = int(pos[0]) + 1
                t.tag_add("hover", f"{line}.0", f"{line}.end")

    def _detect(self, p):
        """Narrowband detection on the displayed PSD; the time of the window comes from sample counts."""
        v = self.view
        j = 0 if self.psd_source.get() == PSD_SOURCES[0] else max(v.wf.count - 1, 0)
        v.t_win = (int(v.wf.starts()[j]) + v.discard) / v.fs      # samples since the stream was started
        t_det = time.perf_counter()
        center, width, power, snr = v.nb.detect(p, v.dens_db)
        tr = v.tracker
        slots = tr.update(center, width, power, snr, v.t_win, v.df)
        v.rawlog.append(v.t_win, center, width, power, snr, np.where(slots >= 0, tr.ident[np.maximum(slots, 0)], -1))
        if self.perflog is not None:
            self.perflog.add_detect(time.perf_counter() - t_det)
        x = v.nb_segs[:, :, 0]
        x.fill(np.nan)
        pres = np.flatnonzero(tr.present & tr.confirmed)
        if len(pres):
            x[pres, 0] = x[pres, 1] = tr.center[pres] * 1e-6
        for lc in self.nb_lcs:
            lc.set_segments(v.nb_segs)

    def _mode_changed(self):
        nb = self.mode.get() == "nb"
        self.rssi_box.configure(text="Narrowband signals" if nb else "RSSI per channel")
        self.rssi_text.configure(width=NB_TEXT_W if nb else RSSI_TEXT_W)
        self.rssi_head.configure(width=NB_TEXT_W if nb else RSSI_TEXT_W)
        self.hover_idx = -1
        self.tip.place_forget()
        self._focus_plots()
        self._apply_overlay_visibility()
        self._update_highlight()
        self.last_text = 0.0
        self.last_draw = 0.0
        self.need_draw = True

    def _nb_changed(self):
        """Max width, guard or threshold changed: apply at once; the tracked signals start over."""
        try:
            params = parse_nb(self._fields())
        except ConfigError as e:
            self._mark_invalid(e.field, str(e))
            return
        self._clear_invalid()
        v = self.view
        if v is None or params == v.nb_params:
            return
        v.nb_params = params
        v.nb.set_params(*params)
        v.tracker.clear()
        self.last_text = 0.0
        self.last_draw = 0.0
        self.need_draw = True

    def _set_head(self, text):
        """The fixed line with the column titles above the rows."""
        h = self.rssi_head
        h.configure(state="normal")
        h.delete("1.0", "end")
        h.insert("end", text)
        h.configure(state="disabled")

    def _rssi_message(self, text):
        self._set_head("")
        t = self.rssi_text
        t.configure(state="normal")
        t.delete("1.0", "end")
        t.insert("end", text, "empty")
        t.configure(state="disabled")

    def copy_rssi(self):
        head = self.rssi_head.get("1.0", "end").strip("\n")
        body = self.rssi_text.get("1.0", "end").strip()
        text = f"{head}\n{body}" if head.strip() else body
        self.root.clipboard_clear()
        self.root.clipboard_append(text)

    # ---------------------------------------------------------------- status

    def _update_status(self, force=False):
        s = self.session
        if self.cfg is None:
            return
        v = self.view
        if v is None:
            self._stat("center", f"Center {fmt_freq(self.plan.center, 3).strip()} MHz")
            self._stat("rate", f"Rate {self.cfg.rate * 1e-6:.3f} MS/s")
            return
        self._stat("center", f"Center {fmt_freq(v.fc, 3).strip()} MHz")
        self._stat("rate", f"Rate {v.fs * 1e-6:.3f} MS/s")
        self._stat("res", f"Resolution {v.df / 1e3:.3f} kHz")
        g = "?" if v.gain is None else f"{v.gain} dB"
        self._stat("gain", f"Gain {g}", C["amber"] if getattr(self, "stat_gain_warn", False) else None)
        rate = 1.0 / self.period if self.period > 0 else 0.0
        self._stat("rows", f"Rows {rate:.1f}/s")
        if s is not None:
            r = s.ring
            lag = max(0, r.head - s.processed_end - r.row_samples) / v.fs      # behind by more than the row in progress
            text = f"Buffer {min(r.head, r.size) / v.fs:.1f}/{r.size / v.fs:.1f} s"
            self._stat("ring", text + (f"  lag {lag:.1f} s" if lag >= 0.2 else ""), C["amber"] if lag >= 0.2 else None)
        dropped = 0 if s is None else s.dropped_rows
        self._stat("dropped", f"Dropped {dropped}", C["amber"] if dropped else None)
        if self.mode.get() == "nb":
            self._stat("chan", f"Narrowband {int((v.tracker.present & v.tracker.confirmed).sum())} present")
        else:
            self._stat("chan", f"Channels {len(v.idx)}")
        rec = self.recorder
        if rec is not None and rec.state in ("recording", "paused"):
            self._stat("rec", f"{'Rec' if rec.state == 'recording' else 'Rec paused'} "
                              f"{rec.out_samples / rec.plan.out_rate:.1f} s", C["error"] if rec.state == "recording" else C["amber"])
        else:
            self._stat("rec", "")
        parts = []
        lo, hi = self.ax_wf.get_xlim()
        if (v.right - v.left) / max(hi - lo, 1e-12) > 1.001:
            parts.append(f"Zoom x{(v.right - v.left) / (hi - lo):.1f}")
        if self.zoom_y is not None:
            parts.append("PSD levels")
        self._stat("zoom", "  ".join(parts), C["accent"] if parts else None)

    # ----------------------------------------------------------------- hover

    def _channel_at(self, x_mhz):
        v = self.view
        if v is None or not len(v.idx):
            return -1
        half = v.bw_in * 0.5e-6
        d = np.abs(x_mhz - v.f_in * 1e-6)
        ok = np.flatnonzero(d <= half)
        if not len(ok):
            return -1
        return int(ok[np.argmin(d[ok])])

    def _nb_at(self, x_mhz):
        """Track slot of the present narrowband signal under the cursor, or -1."""
        v = self.view
        tr = v.tracker
        pres = np.flatnonzero(tr.present & tr.confirmed)
        if v is None or not len(pres):
            return -1
        tol = np.maximum(0.5e-6 * tr.width[pres], 0.004 * self._view_span())
        d = np.abs(x_mhz - tr.center[pres] * 1e-6)
        ok = np.flatnonzero(d <= tol)
        if not len(ok):
            return -1
        return int(pres[ok[np.argmin(d[ok])]])

    def _hover_target(self, x_mhz):
        if not self.overlay.get():
            return -1
        return self._nb_at(x_mhz) if self.mode.get() == "nb" else self._channel_at(x_mhz)

    def _view_span(self):
        """Visible frequency span, MHz."""
        lo, hi = self.ax_wf.get_xlim()
        return hi - lo

    def _set_xlim(self, lo, hi):
        """Zoom/pan the frequency axis of both plots (they share it), kept inside the capture band."""
        v = self.view
        full = v.right - v.left
        span = min(max(hi - lo, ZOOM_MIN_BINS * v.df * 1e-6), full)
        lo = min(max(lo, v.left), v.right - span)
        if span >= full * (1.0 - 1e-9):
            self.zoom_x, lo, span = None, v.left, full
        else:
            self.zoom_x = (lo, lo + span)
        self.ax_wf.set_xlim(lo, lo + span)
        self._after_zoom()

    def _after_zoom(self):
        self.need_draw = True
        self.last_draw = 0.0
        self._update_status()

    def reset_zoom(self):
        v = self.view
        self.zoom_x = self.zoom_y = None
        self._drag = None
        if v is not None:
            self.ax_wf.set_xlim(v.left, v.right)
            self._after_zoom()
        self._focus_plots()

    def on_scroll(self, event):
        v = self.view
        if v is None or event.inaxes not in (self.ax_wf, self.ax_psd) or event.xdata is None:
            return
        zoom_in = event.button == "up" or getattr(event, "step", 0) > 0
        f = 1.0 / ZOOM_STEP if zoom_in else ZOOM_STEP
        gui = getattr(event, "guiEvent", None)
        ctrl = bool(getattr(gui, "state", 0) & 0x4)
        if ctrl and event.inaxes is self.ax_psd and event.ydata is not None:
            lo, hi = self.ax_psd.get_ylim()
            y = event.ydata
            self.zoom_y = (y - (y - lo) * f, y + (hi - y) * f)
            self.ax_psd.set_ylim(*self.zoom_y)
            self._after_zoom()
        else:
            lo, hi = self.ax_wf.get_xlim()
            x = event.xdata                       # the point under the cursor stays where it is
            self._set_xlim(x - (x - lo) * f, x + (hi - x) * f)

    def on_press(self, event):
        if self.view is None or event.button != 1 or event.inaxes not in (self.ax_wf, self.ax_psd):
            return
        gui = getattr(event, "guiEvent", None)
        if gui is not None and (getattr(gui, "state", 0) & 0x1):       # Shift + left click
            self._cap_spawn(event)
            return
        if getattr(event, "dblclick", False):
            self.reset_zoom()
            return
        hit = self._cap_hit(event, edges=not self._cap_recording())    # while recording the width is the file's: no edges
        if hit:                                                    # the band to record is grabbed, not the plot
            self.cap_drag = {"mode": hit, "center": self.cap_center, "half": self._cap_half(), "x0": event.xdata}
            self.tip.place_forget()
            return
        self._drag = {"px": event.x, "py": event.y, "xlim": self.ax_wf.get_xlim(),
                      "ylim": self.ax_psd.get_ylim(), "ax": event.inaxes}
        self.tip.place_forget()

    def on_release(self, _event):
        self._drag = None
        if self.cap_drag is not None:
            self.cap_drag = None
            if self._cap_recording() and self.cap_center != self.recorder.spec["center"]:
                self._rec_center_request(center=self.cap_center, source="Manual band")   # the band was dragged: ask

    def _pan(self, event):
        d = self._drag
        if event.x is None:
            return
        x0, x1 = d["xlim"]
        dx = (event.x - d["px"]) / max(self.ax_psd.bbox.width, 1.0) * (x1 - x0)
        self._set_xlim(x0 - dx, x1 - dx)
        if d["ax"] is self.ax_psd and self.zoom_y is not None:
            y0, y1 = d["ylim"]
            dy = (event.y - d["py"]) / max(self.ax_psd.bbox.height, 1.0) * (y1 - y0)
            self.zoom_y = (y0 - dy, y1 - dy)
            self.ax_psd.set_ylim(*self.zoom_y)
            self._after_zoom()

    def _show_tip(self, event, text):
        """The label next to the cursor; it flips to the other side near the right and bottom edge of the plots."""
        gui = getattr(event, "guiEvent", None)
        ch = self.canvas_widget.winfo_height()
        gx, gy = (gui.x, gui.y) if gui is not None else (event.x, ch - event.y)    # Tk's y runs downwards
        self.tip.configure(text=text)
        self.tip.update_idletasks()
        w, h = self.tip.winfo_reqwidth(), self.tip.winfo_reqheight()
        cw = self.canvas_widget.winfo_width()
        x, y = gx + 14, gy + 14
        if x + w > cw:
            x = max(0, gx - w - 14)
        if y + h > ch:
            y = max(0, gy - h - 14)
        self.tip.place(x=x, y=y)
        self.tip.lift()

    def _cap_tip_text(self):
        c, half = self.cap_center, self._cap_half()
        return (f"Capture band\nCenter {fmt_freq(c, 4).strip()} MHz\nLow     {fmt_freq(c - half, 4).strip()} MHz\n"
                f"High    {fmt_freq(c + half, 4).strip()} MHz\nWidth {2 * half / 1e3:.3f} kHz")

    def on_motion(self, event):
        v = self.view
        if self.cap_drag is not None and v is not None:
            self._cap_drag_move(event)
            return
        if self._drag is not None and v is not None:
            self._pan(event)
            return
        hit = self._cap_hit(event, edges=not self._cap_recording()) if v is not None else None
        if v is not None:
            self._cap_cursor(hit)
        if hit:                                      # over the band to record: its frequencies, and nothing else
            if self.hover_idx != -1:
                self.hover_idx = -1
                self._update_highlight()
                self.need_draw = True
                self.last_draw = 0.0
                self.last_text = 0.0
            self._show_tip(event, self._cap_tip_text())
            self.hover_lbl.configure(text=self._cursor_text(event))
            return
        if v is None or event.inaxes not in (self.ax_wf, self.ax_psd) or event.xdata is None:
            self._clear_hover()
            return
        k = self._hover_target(event.xdata)
        if k != self.hover_idx:
            self.hover_idx = k
            self._update_highlight()
            self.need_draw = True
            self.last_draw = 0.0
            self.last_text = 0.0
        if k >= 0 and v.have_row:
            self._show_tip(event, self._tip_text(k))
        else:
            self.tip.place_forget()
        self.hover_lbl.configure(text=self._cursor_text(event))

    def on_leave(self, _event=None):
        self._clear_hover()

    def _clear_hover(self):
        self.tip.place_forget()
        if self.hover_idx != -1:
            self.hover_idx = -1
            self._update_highlight()
            self.need_draw = True
            self.last_draw = 0.0
            self.last_text = 0.0
        self.hover_lbl.configure(text="")

    def _tip_text(self, k):
        v = self.view
        if self.mode.get() == "nb":
            tr = v.tracker
            return (f"Signal {int(tr.ident[k])}\nCenter {fmt_freq(tr.center[k], 4).strip()} MHz\n"
                    f"Width {tr.width[k] * 1e-3:.2f} kHz\nPower {tr.power[k]:.2f} dBFS\nSNR {tr.snr[k]:.1f} dB\n"
                    f"Present: {'yes' if tr.present[k] else 'no'}\nLast seen: {tr.last_seen[k]:.3f} s")
        val = v.rssi[k]
        rssi = "n/a" if not np.isfinite(val) else f"{val:.2f} dBFS"
        return (f"Channel {int(v.idx[k]) + 1}\nCenter {fmt_freq(v.f_in[k], 3).strip()} MHz\n"
                f"Width {fmt_freq(v.bw_in[k], 3).strip()} MHz\nRSSI {rssi}")

    def _update_highlight(self):
        v = self.view
        k = self.hover_idx
        x0 = w = None
        if k >= 0 and v is not None:
            if self.mode.get() == "nb":
                tr = v.tracker
                if tr.visible[k] and tr.present[k]:
                    w = max(tr.width[k] * 1e-6, 0.004 * self._view_span())
                    x0 = tr.center[k] * 1e-6 - w / 2
            elif k < len(v.idx):
                x0 = (v.f_in[k] - v.bw_in[k] / 2) * 1e-6
                w = v.bw_in[k] * 1e-6
        for r in self.hl_patches:
            if x0 is None:
                r.set_visible(False)
            else:
                r.set_x(x0)
                r.set_width(w)
                r.set_visible(True)

    def _cursor_text(self, event):
        v = self.view
        x = event.xdata
        f = f"{x:10.3f} MHz"
        if event.inaxes is self.ax_psd and v.have_row:
            j = int(round((x * 1e6 - (v.fc - (v.n // 2) * v.df)) / v.df))
            if 0 <= j < v.n:
                return f"{f}  {v.dens_db[j]:8.2f} dBFS/Hz"
        elif event.inaxes is self.ax_wf and v.have_row and event.ydata is not None:
            j = int(round((x * 1e6 - (v.fc - (v.n // 2) * v.df)) / v.df))
            r = int(event.ydata / v.row_dur)
            if 0 <= j < v.n and 0 <= r < v.wf.count:
                return f"{f}  {event.ydata:6.2f} s  {v.wf.view()[r, j]:8.2f} dBFS/Hz"
        return f


class _View:
    """Display model of the running session (filled by App.on_ready)."""


# ---------------------------------------------------------------------- entry

def main():
    args = parse_args()
    try:
        validate(args_to_fields(args))
    except ConfigError as e:
        info(f"ERROR: {e}")
        return 2
    try:
        cap_rate = parse_hz(args.cap_rate)
        if cap_rate <= 0:
            raise ValueError("The capture width must be above 0")
    except ValueError as e:
        info(f"ERROR: --cap-rate: {e}")
        return 2
    try:
        render_stem(args.cap_pattern, {"name": "x", "ident": 1, "center_hz": 1e9, "rate_hz": 48e3, "t": time.time(), "n": 1})
    except ValueError as e:
        info(f"ERROR: --cap-pattern: {e}")
        return 2
    metrics_fh = None
    if args.metrics:
        try:
            metrics_fh = open(args.metrics, "w", encoding="utf-8")
        except OSError as e:
            info(f"ERROR: cannot write the metrics file: {e}")
            return 2
        metrics_fh.write("# rssi_waterfall performance metrics; one row per second per session, '# event' and "
                         "'# summary' lines in between\n" + "\t".join(PERF_COLS) + "\n")
    root = tk.Tk()
    app = App(root, args_to_fields(args), metrics_fh, {"rate": cap_rate, "dir": args.cap_dir, "pattern": args.cap_pattern})

    def on_sigint(_s, _f):
        app.sigint = True
    signal.signal(signal.SIGINT, on_sigint)
    root.mainloop()
    return 0


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