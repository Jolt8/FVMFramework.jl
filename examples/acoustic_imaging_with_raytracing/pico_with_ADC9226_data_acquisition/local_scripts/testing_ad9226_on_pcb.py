"""
Triggered AD9226 diagnostic for the custom tomography PCB.

This is deliberately a signal-integrity test, not a Time-of-Flight detector.
It answers four progressively stronger questions:

1. Are the 12 ADC data bits changing?
2. Is a response repeatably locked to the transmit pulse?
3. Is there energy in the expected ultrasonic frequency band?
4. Is that response usable, or is the ADC/amplifier clipping?

The script captures several waveforms through PicoADCInterface, prints a
plain-language verdict, and saves a plot, JSON report, and compressed raw data.

Examples:
    python testing_ad9226_on_pcb.py
    python testing_ad9226_on_pcb.py --gain-code 64 --captures 12
    python testing_ad9226_on_pcb.py --band-low-mhz 1.5 --band-high-mhz 2.5
    python testing_ad9226_on_pcb.py --port COM5 --tx 1 --rx 2
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime

import numpy as np

try:
    from scipy.signal import butter, hilbert, sosfiltfilt, welch
except ImportError:
    print("[X] scipy is required. Install it with: pip install scipy")
    raise

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface


ADC_MIN_CODE = 0
ADC_MAX_CODE = 4095
ADC_MIDSCALE = 2048.0


def bounded_int(minimum, maximum, label):
    """Build an argparse integer validator."""
    def parse(value):
        parsed = int(value)
        if not minimum <= parsed <= maximum:
            raise argparse.ArgumentTypeError(
                "{} must be from {} to {}".format(label, minimum, maximum)
            )
        return parsed

    return parse


def positive_float(value):
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return parsed


def volts_to_codes(volts, vref):
    """Undo PicoADCInterface's 12-bit code-to-voltage conversion exactly."""
    codes = np.rint((volts / vref) * ADC_MIDSCALE + ADC_MIDSCALE)
    return np.clip(codes, ADC_MIN_CODE, ADC_MAX_CODE).astype(np.uint16)


def rms(values):
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(values * values)))


def make_time_mask(times_us, start_us, end_us):
    return (times_us >= start_us) & (times_us <= end_us)


def bandpass_captures(captures, sample_rate_hz, low_hz, high_hz):
    """Zero-phase bandpass each capture, preserving trigger alignment."""
    nyquist_hz = sample_rate_hz / 2.0
    if low_hz <= 0 or high_hz >= nyquist_hz or low_hz >= high_hz:
        raise ValueError(
            "band must satisfy 0 < low < high < {:.3f} MHz".format(
                nyquist_hz / 1e6
            )
        )

    sos = butter(
        4,
        [low_hz / nyquist_hz, high_hz / nyquist_hz],
        btype="bandpass",
        output="sos",
    )
    return sosfiltfilt(sos, captures, axis=1)


def bit_activity(codes):
    """Return the fraction of samples in which each ADC result bit is one."""
    flattened = codes.reshape(-1).astype(np.uint16)
    return np.array(
        [np.mean((flattened >> bit) & 1) for bit in range(12)],
        dtype=np.float64,
    )


def pairwise_correlation(captures):
    """Median correlation between independently triggered captures."""
    if captures.shape[0] < 2 or captures.shape[1] < 2:
        return None

    row_std = np.std(captures, axis=1)
    usable = captures[row_std > 1e-12]
    if usable.shape[0] < 2:
        return None

    matrix = np.corrcoef(usable)
    values = matrix[np.triu_indices(matrix.shape[0], k=1)]
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else None


def analyse_captures(
    captures,
    sample_rate_hz,
    dt_us,
    vref,
    signal_start_us,
    signal_end_us,
    band_low_hz,
    band_high_hz,
):
    """Calculate digital, clipping, timing, and frequency diagnostics."""
    capture_count, sample_count = captures.shape
    times_us = np.arange(sample_count, dtype=np.float64) * dt_us
    duration_us = times_us[-1] if sample_count else 0.0

    if signal_start_us >= duration_us:
        raise ValueError(
            "signal window starts after the {:.1f} us capture".format(duration_us)
        )
    signal_end_us = min(signal_end_us, duration_us)
    if signal_end_us <= signal_start_us:
        raise ValueError("signal window end must be after its start")

    # Use the quiet tail as the noise reference. This avoids assuming that the
    # firmware provides an external-scope-style pre-trigger region.
    noise_start_us = max(signal_end_us + 50.0, duration_us * 0.60)
    noise_end_us = duration_us * 0.95
    if noise_start_us >= noise_end_us:
        noise_start_us = max(signal_end_us, duration_us * 0.75)
        noise_end_us = duration_us

    signal_mask = make_time_mask(times_us, signal_start_us, signal_end_us)
    noise_mask = make_time_mask(times_us, noise_start_us, noise_end_us)
    if np.count_nonzero(signal_mask) < 16 or np.count_nonzero(noise_mask) < 16:
        raise ValueError("capture is too short for the requested analysis windows")

    codes = volts_to_codes(captures, vref)
    flattened_codes = codes.reshape(-1)
    unique_codes = int(np.unique(flattened_codes).size)
    minimum_code = int(np.min(flattened_codes))
    maximum_code = int(np.max(flattened_codes))
    code_span = maximum_code - minimum_code

    exact_low_fraction = float(np.mean(flattened_codes == ADC_MIN_CODE))
    exact_high_fraction = float(np.mean(flattened_codes == ADC_MAX_CODE))
    near_rail_codes = 41  # Approximately 1% of the 12-bit range.
    near_low_fraction = float(np.mean(flattened_codes <= near_rail_codes))
    near_high_fraction = float(
        np.mean(flattened_codes >= ADC_MAX_CODE - near_rail_codes)
    )
    near_rail_fraction = near_low_fraction + near_high_fraction

    transitions = np.diff(codes.astype(np.int32), axis=1)
    changing_fraction = float(np.mean(transitions != 0))
    largest_step_codes = int(np.max(np.abs(transitions))) if transitions.size else 0
    ones_fraction = bit_activity(codes)

    filtered = bandpass_captures(
        captures, sample_rate_hz, band_low_hz, band_high_hz
    )
    envelopes = np.abs(hilbert(filtered, axis=1))
    mean_trace = np.mean(captures, axis=0)
    mean_filtered = np.mean(filtered, axis=0)
    mean_envelope = np.mean(envelopes, axis=0)

    signal_band_rms = rms(filtered[:, signal_mask])
    noise_band_rms = rms(filtered[:, noise_mask])
    band_rms_ratio = signal_band_rms / max(noise_band_rms, 1e-15)

    signal_indices = np.flatnonzero(signal_mask)
    local_peak = int(np.argmax(mean_envelope[signal_mask]))
    peak_index = int(signal_indices[local_peak])
    peak_time_us = float(times_us[peak_index])
    peak_envelope_v = float(mean_envelope[peak_index])
    noise_envelope_median = float(np.median(envelopes[:, noise_mask]))
    peak_to_noise = peak_envelope_v / max(noise_envelope_median, 1e-15)

    signal_filtered = filtered[:, signal_mask]
    coherent_rms = rms(np.mean(signal_filtered, axis=0))
    total_rms = rms(signal_filtered)
    coherent_fraction = coherent_rms / max(total_rms, 1e-15)
    median_correlation = pairwise_correlation(signal_filtered)

    # Average Welch spectra across captures. Detrending suppresses the large,
    # slow recovery transient visible in the oscilloscope screenshot.
    signal_length = int(np.count_nonzero(signal_mask))
    nperseg = min(2048, signal_length)
    frequencies_hz, psd = welch(
        captures[:, signal_mask],
        fs=sample_rate_hz,
        window="hann",
        nperseg=nperseg,
        detrend="linear",
        axis=1,
        scaling="density",
    )
    mean_psd = np.mean(psd, axis=0)
    non_dc = frequencies_hz >= 50_000.0
    band_bins = (frequencies_hz >= band_low_hz) & (
        frequencies_hz <= band_high_hz
    )
    if np.any(non_dc):
        peak_spectrum_index = int(np.argmax(mean_psd[non_dc]))
        peak_frequency_hz = float(frequencies_hz[non_dc][peak_spectrum_index])
    else:
        peak_frequency_hz = 0.0

    total_power = float(np.trapezoid(mean_psd[non_dc], frequencies_hz[non_dc]))
    band_power = float(np.trapezoid(mean_psd[band_bins], frequencies_hz[band_bins]))
    ultrasonic_power_fraction = band_power / max(total_power, 1e-30)

    adc_bus_active = (
        unique_codes >= 16 and code_span >= 16 and changing_fraction >= 0.001
    )
    severe_clipping = near_rail_fraction >= 0.10 or (
        exact_low_fraction + exact_high_fraction
    ) >= 0.01
    band_response = band_rms_ratio >= 3.0 and peak_to_noise >= 5.0
    repeatable_response = coherent_fraction >= 0.50 and (
        median_correlation is None or median_correlation >= 0.20
    )

    if not adc_bus_active:
        verdict = (
            "NO MEANINGFUL ADC ACTIVITY: the captured codes are nearly static. "
            "Check the ADC clock, data bus, output-enable state, and analog input."
        )
    elif severe_clipping and band_response:
        verdict = (
            "ADC IS RESPONDING, BUT IT IS OVERDRIVEN: a trigger-related ultrasonic-"
            "band response exists, but rail clipping makes it unsuitable for ToF. "
            "Reduce the MCP4161 gain and/or analog input amplitude."
        )
    elif severe_clipping:
        verdict = (
            "ADC DATA IS DEFINITELY CHANGING, BUT IT IS HEAVILY CLIPPED: this proves "
            "activity, not a clean ultrasonic echo. Reduce gain before tuning ToF."
        )
    elif band_response and repeatable_response:
        verdict = (
            "USABLE TRIGGER-LOCKED RESPONSE DETECTED: the ADC sees repeatable energy "
            "inside the requested ultrasonic band without severe clipping."
        )
    elif band_response:
        verdict = (
            "ULTRASONIC-BAND ENERGY DETECTED, BUT REPEATABILITY IS WEAK: inspect the "
            "overlay and check triggering, mux selection, grounding, and noise."
        )
    else:
        verdict = (
            "ADC DATA IS CHANGING, BUT NO CLEAR ULTRASONIC-BAND RESPONSE WAS FOUND "
            "in the selected time/frequency window. Inspect the raw plot and adjust "
            "the band or signal window if necessary."
        )

    per_capture = []
    for index in range(capture_count):
        capture_codes = codes[index]
        per_capture.append(
            {
                "capture": index + 1,
                "min_code": int(np.min(capture_codes)),
                "max_code": int(np.max(capture_codes)),
                "mean_v": float(np.mean(captures[index])),
                "std_v": float(np.std(captures[index])),
                "near_rail_percent": float(
                    100.0
                    * np.mean(
                        (capture_codes <= near_rail_codes)
                        | (capture_codes >= ADC_MAX_CODE - near_rail_codes)
                    )
                ),
            }
        )

    metrics = {
        "capture_count": capture_count,
        "sample_count": sample_count,
        "sample_rate_hz": int(sample_rate_hz),
        "sample_interval_us": float(dt_us),
        "capture_duration_us": float(duration_us),
        "signal_window_us": [float(signal_start_us), float(signal_end_us)],
        "noise_window_us": [float(noise_start_us), float(noise_end_us)],
        "band_hz": [float(band_low_hz), float(band_high_hz)],
        "minimum_code": minimum_code,
        "maximum_code": maximum_code,
        "code_span": code_span,
        "unique_codes": unique_codes,
        "changing_sample_percent": changing_fraction * 100.0,
        "largest_adjacent_step_codes": largest_step_codes,
        "exact_low_percent": exact_low_fraction * 100.0,
        "exact_high_percent": exact_high_fraction * 100.0,
        "near_low_percent": near_low_fraction * 100.0,
        "near_high_percent": near_high_fraction * 100.0,
        "near_rail_percent": near_rail_fraction * 100.0,
        "bit_one_percent_lsb_to_msb": (ones_fraction * 100.0).tolist(),
        "signal_band_rms_v": signal_band_rms,
        "noise_band_rms_v": noise_band_rms,
        "band_rms_ratio": band_rms_ratio,
        "peak_envelope_v": peak_envelope_v,
        "peak_time_us": peak_time_us,
        "peak_to_noise_ratio": peak_to_noise,
        "coherent_fraction": coherent_fraction,
        "median_pairwise_correlation": median_correlation,
        "peak_spectrum_frequency_hz": peak_frequency_hz,
        "ultrasonic_power_fraction": ultrasonic_power_fraction,
        "adc_bus_active": bool(adc_bus_active),
        "severe_clipping": bool(severe_clipping),
        "band_response": bool(band_response),
        "repeatable_response": bool(repeatable_response),
        "verdict": verdict,
        "per_capture": per_capture,
    }

    arrays = {
        "times_us": times_us,
        "codes": codes,
        "mean_trace": mean_trace,
        "filtered": filtered,
        "mean_filtered": mean_filtered,
        "envelopes": envelopes,
        "mean_envelope": mean_envelope,
        "frequencies_hz": frequencies_hz,
        "mean_psd": mean_psd,
        "signal_mask": signal_mask,
        "noise_mask": noise_mask,
    }
    return metrics, arrays


def print_report(metrics, gain_code, tx, rx):
    correlation = metrics["median_pairwise_correlation"]
    correlation_text = "n/a" if correlation is None else "{:.3f}".format(correlation)

    print("\n" + "=" * 72)
    print("AD9226 PCB DIAGNOSTIC REPORT")
    print("=" * 72)
    print("Path                    : TX {} -> RX {}".format(tx, rx))
    print("MCP4161 wiper           : {}".format(gain_code))
    print(
        "Capture                 : {} x {} samples at {:.3f} MSPS".format(
            metrics["capture_count"],
            metrics["sample_count"],
            metrics["sample_rate_hz"] / 1e6,
        )
    )
    print(
        "ADC code range          : {} to {} (span {}, {} unique)".format(
            metrics["minimum_code"],
            metrics["maximum_code"],
            metrics["code_span"],
            metrics["unique_codes"],
        )
    )
    print(
        "Samples that change     : {:.2f}%".format(
            metrics["changing_sample_percent"]
        )
    )
    print(
        "Exact rails             : low {:.3f}% | high {:.3f}%".format(
            metrics["exact_low_percent"], metrics["exact_high_percent"]
        )
    )
    print(
        "Within 1% of rails      : low {:.2f}% | high {:.2f}% | total {:.2f}%".format(
            metrics["near_low_percent"],
            metrics["near_high_percent"],
            metrics["near_rail_percent"],
        )
    )
    print(
        "Band RMS signal/noise   : {:.6f} / {:.6f} V ({:.2f}x)".format(
            metrics["signal_band_rms_v"],
            metrics["noise_band_rms_v"],
            metrics["band_rms_ratio"],
        )
    )
    print(
        "Band-envelope peak      : {:.6f} V at {:.3f} us ({:.2f}x noise)".format(
            metrics["peak_envelope_v"],
            metrics["peak_time_us"],
            metrics["peak_to_noise_ratio"],
        )
    )
    print(
        "Trigger repeatability   : coherent {:.1f}% | median correlation {}".format(
            100.0 * metrics["coherent_fraction"], correlation_text
        )
    )
    print(
        "Strongest AC frequency  : {:.4f} MHz".format(
            metrics["peak_spectrum_frequency_hz"] / 1e6
        )
    )
    print(
        "Power inside test band  : {:.2f}%".format(
            100.0 * metrics["ultrasonic_power_fraction"]
        )
    )
    print("\nVERDICT:")
    print("  " + metrics["verdict"])
    print("=" * 72)


def save_plot(captures, metrics, arrays, output_path, gain_code, tx, rx, vref):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    times_us = arrays["times_us"]
    signal_mask = arrays["signal_mask"]
    noise_mask = arrays["noise_mask"]
    display_end_us = min(metrics["signal_window_us"][1] + 50.0, times_us[-1])
    display_mask = times_us <= display_end_us

    fig, axes = plt.subplots(3, 2, figsize=(16, 13))
    fig.suptitle(
        "AD9226 PCB Diagnostic - TX {} -> RX {}, MCP4161 wiper {}".format(
            tx, rx, gain_code
        ),
        fontsize=15,
        fontweight="bold",
    )

    ax = axes[0, 0]
    for capture in captures:
        ax.plot(times_us[display_mask], capture[display_mask], alpha=0.25, lw=0.6)
    ax.plot(
        times_us[display_mask],
        arrays["mean_trace"][display_mask],
        color="black",
        lw=1.2,
        label="Mean",
    )
    ax.axhline(vref, color="red", ls="--", alpha=0.6)
    ax.axhline(-vref, color="red", ls="--", alpha=0.6, label="ADC rails")
    ax.set_title("Triggered raw captures (zoom)")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Voltage (V)")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(times_us, arrays["mean_trace"], color="#1665A7", lw=0.8)
    ax.axvspan(*metrics["signal_window_us"], color="orange", alpha=0.18, label="Signal")
    ax.axvspan(*metrics["noise_window_us"], color="green", alpha=0.14, label="Noise")
    ax.set_title("Mean raw capture (complete window)")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Voltage (V)")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1, 0]
    for filtered_capture in arrays["filtered"]:
        ax.plot(
            times_us[display_mask],
            filtered_capture[display_mask],
            alpha=0.20,
            lw=0.55,
        )
    ax.plot(
        times_us[display_mask],
        arrays["mean_filtered"][display_mask],
        color="black",
        lw=1.0,
        label="Mean bandpassed",
    )
    ax.set_title(
        "Bandpassed {:.2f}-{:.2f} MHz".format(
            metrics["band_hz"][0] / 1e6, metrics["band_hz"][1] / 1e6
        )
    )
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Voltage (V)")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1, 1]
    ax.plot(times_us, arrays["mean_envelope"], color="#D34A24", lw=1.0)
    ax.axvspan(*metrics["signal_window_us"], color="orange", alpha=0.18)
    ax.axvspan(*metrics["noise_window_us"], color="green", alpha=0.14)
    ax.axvline(
        metrics["peak_time_us"],
        color="purple",
        ls="--",
        label="Peak {:.2f} us".format(metrics["peak_time_us"]),
    )
    ax.set_title("Bandpassed Hilbert envelope")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Envelope (V)")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[2, 0]
    positive = arrays["frequencies_hz"] > 0
    ax.semilogy(
        arrays["frequencies_hz"][positive] / 1e6,
        np.maximum(arrays["mean_psd"][positive], 1e-30),
        color="#4C8B3B",
        lw=1.0,
    )
    ax.axvspan(
        metrics["band_hz"][0] / 1e6,
        metrics["band_hz"][1] / 1e6,
        color="orange",
        alpha=0.20,
        label="Test band",
    )
    ax.axvline(
        metrics["peak_spectrum_frequency_hz"] / 1e6,
        color="purple",
        ls="--",
        label="Peak {:.3f} MHz".format(
            metrics["peak_spectrum_frequency_hz"] / 1e6
        ),
    )
    ax.set_xlim(0, metrics["sample_rate_hz"] / 2e6)
    ax.set_title("Mean power spectrum in signal window")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("PSD (V^2/Hz)")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[2, 1]
    ax.hist(
        arrays["codes"].reshape(-1),
        bins=128,
        range=(ADC_MIN_CODE, ADC_MAX_CODE),
        color="#3477B8",
        alpha=0.85,
    )
    ax.axvline(41, color="red", ls="--")
    ax.axvline(ADC_MAX_CODE - 41, color="red", ls="--", label="1% rail zones")
    ax.set_title(
        "ADC code histogram - {:.2f}% near rails".format(
            metrics["near_rail_percent"]
        )
    )
    ax.set_xlabel("12-bit ADC code")
    ax.set_ylabel("Count")
    ax.grid(alpha=0.25)
    ax.legend()

    verdict_color = "darkred" if metrics["severe_clipping"] else "darkgreen"
    fig.text(
        0.5,
        0.012,
        metrics["verdict"],
        ha="center",
        va="bottom",
        fontsize=10,
        color=verdict_color,
        wrap=True,
    )
    plt.tight_layout(rect=(0, 0.045, 1, 0.96))
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def acquire(args):
    adc = PicoADCInterface(port=args.port, vref=args.vref)
    if not adc.connect():
        raise RuntimeError("could not connect to the Pico")

    captures = []
    actual_rate = None
    actual_dt_us = None
    active_gain_code = None
    try:
        status = adc.get_status()
        if status:
            print("[+] Pico status: {}".format(status))
            active_gain_code = adc.last_wiper_code

        if not adc.configure(args.sample_rate, args.samples):
            raise RuntimeError("Pico rejected the ADC capture configuration")
        if args.gain_code is None:
            print("[+] Leaving the MCP4161 wiper unchanged")
        else:
            if not adc.set_gain(args.gain_code):
                raise RuntimeError("Pico rejected the MCP4161 wiper setting")
            active_gain_code = args.gain_code

        for warmup_index in range(args.warmup_captures):
            print(
                "[*] Warm-up capture {}/{}...".format(
                    warmup_index + 1, args.warmup_captures
                )
            )
            warmup_volts, _, _ = adc.trigger_and_capture(tx=args.tx, rx=args.rx)
            if warmup_volts is None:
                raise RuntimeError(
                    "warm-up capture {} failed".format(warmup_index + 1)
                )
            time.sleep(args.delay)

        for capture_index in range(args.captures):
            print(
                "[*] Diagnostic capture {}/{} (TX {} -> RX {})...".format(
                    capture_index + 1, args.captures, args.tx, args.rx
                ),
                end="",
                flush=True,
            )
            volts, sample_rate_hz, dt_us = adc.trigger_and_capture(
                tx=args.tx, rx=args.rx
            )
            if volts is None:
                print(" failed")
                raise RuntimeError("capture {} failed".format(capture_index + 1))
            print(" received")

            if actual_rate is None:
                actual_rate = sample_rate_hz
                actual_dt_us = dt_us
            elif sample_rate_hz != actual_rate or not np.isclose(dt_us, actual_dt_us):
                raise RuntimeError("sample timing changed between captures")

            captures.append(np.asarray(volts, dtype=np.float64))
            time.sleep(args.delay)
    finally:
        adc.close()

    return np.stack(captures), actual_rate, actual_dt_us, active_gain_code


def build_parser():
    parser = argparse.ArgumentParser(
        description="Determine whether the PCB-mounted AD9226 sees a usable signal."
    )
    parser.add_argument("--port", default=None, help="Pico serial port (auto-detect by default)")
    parser.add_argument("--tx", type=bounded_int(0, 15, "TX channel"), default=1)
    parser.add_argument("--rx", type=bounded_int(0, 15, "RX channel"), default=2)
    parser.add_argument(
        "--gain-code",
        type=bounded_int(0, 256, "gain code"),
        default=None,
        metavar="0..256",
        help="set the MCP4161 volatile wiper (unchanged by default)",
    )
    parser.add_argument("--captures", type=bounded_int(2, 100, "captures"), default=8)
    parser.add_argument(
        "--warmup-captures",
        type=bounded_int(0, 20, "warm-up captures"),
        default=0,
        help="discard this many captures before testing (default: 0)",
    )
    parser.add_argument("--sample-rate", type=int, default=10_000_000)
    parser.add_argument("--samples", type=bounded_int(100, 50_000, "samples"), default=20_000)
    parser.add_argument("--vref", type=positive_float, default=1.0)
    parser.add_argument("--signal-start-us", type=float, default=2.0)
    parser.add_argument("--signal-end-us", type=positive_float, default=300.0)
    parser.add_argument("--band-low-mhz", type=positive_float, default=1.0)
    parser.add_argument("--band-high-mhz", type=positive_float, default=3.0)
    parser.add_argument("--delay", type=float, default=0.02)
    parser.add_argument(
        "--output-dir",
        default=None,
        help="output directory (default: local_scripts/diagnostic_results)",
    )
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--no-raw-save", action="store_true")
    return parser


def main():
    args = build_parser().parse_args()
    if not 100_000 <= args.sample_rate <= 65_000_000:
        raise SystemExit("--sample-rate must be from 100000 to 65000000")
    if args.signal_start_us < 0:
        raise SystemExit("--signal-start-us cannot be negative")
    if args.signal_end_us <= args.signal_start_us:
        raise SystemExit("--signal-end-us must be after --signal-start-us")
    if args.band_high_mhz <= args.band_low_mhz:
        raise SystemExit("--band-high-mhz must exceed --band-low-mhz")
    if args.delay < 0:
        raise SystemExit("--delay cannot be negative")

    print("=" * 72)
    print("AD9226 ON-PCB TRIGGERED SIGNAL DIAGNOSTIC")
    print("=" * 72)
    print(
        "Capturing TX {} -> RX {}, requested wiper {}, {:.2f}-{:.2f} MHz test band".format(
            args.tx,
            args.rx,
            args.gain_code if args.gain_code is not None else "unchanged",
            args.band_low_mhz,
            args.band_high_mhz,
        )
    )

    try:
        captures, sample_rate_hz, dt_us, active_gain_code = acquire(args)
        metrics, arrays = analyse_captures(
            captures=captures,
            sample_rate_hz=sample_rate_hz,
            dt_us=dt_us,
            vref=args.vref,
            signal_start_us=args.signal_start_us,
            signal_end_us=args.signal_end_us,
            band_low_hz=args.band_low_mhz * 1e6,
            band_high_hz=args.band_high_mhz * 1e6,
        )
    except (RuntimeError, ValueError) as error:
        print("[X] Diagnostic failed: {}".format(error))
        return 1

    metrics.update(
        {
            "tx_channel": args.tx,
            "rx_channel": args.rx,
            "mcp4161_wiper_code": active_gain_code,
            "vref": args.vref,
        }
    )
    print_report(metrics, active_gain_code, args.tx, args.rx)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = args.output_dir or os.path.join(script_dir, "diagnostic_results")
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
    stem = "ad9226_pcb_diagnostic_{}".format(timestamp)

    report_path = os.path.join(output_dir, stem + ".json")
    with open(report_path, "w", encoding="utf-8") as report_file:
        json.dump(metrics, report_file, indent=2)
    print("[+] JSON report saved: {}".format(report_path))

    if not args.no_raw_save:
        raw_path = os.path.join(output_dir, stem + ".npz")
        np.savez_compressed(
            raw_path,
            volts=captures,
            adc_codes=arrays["codes"],
            times_us=arrays["times_us"],
            filtered_volts=arrays["filtered"],
            envelope_volts=arrays["envelopes"],
            sample_rate_hz=sample_rate_hz,
            mcp4161_wiper_code=(
                active_gain_code if active_gain_code is not None else -1
            ),
            tx_channel=args.tx,
            rx_channel=args.rx,
        )
        print("[+] Raw capture bundle saved: {}".format(raw_path))

    if not args.no_plot:
        plot_path = os.path.join(output_dir, stem + ".png")
        try:
            save_plot(
                captures,
                metrics,
                arrays,
                plot_path,
                active_gain_code,
                args.tx,
                args.rx,
                args.vref,
            )
            print("[+] Diagnostic plot saved: {}".format(plot_path))
        except ImportError:
            print("[!] matplotlib is unavailable; JSON and raw data were still saved")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
