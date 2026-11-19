"""Native-resolution, model-free DNA gel analysis."""

import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d, grey_opening
from scipy.signal import find_peaks, peak_widths

LADDERS = {
    "50bp": [1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50],
    "100bp": [1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
    "1kb": [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 500],
}


@dataclass
class GelResult:
    bands: pd.DataFrame
    lane_centers: list[int]
    warnings: list[str] = field(default_factory=list)
    output_dir: Path | None = None


def load_grayscale(path: str | Path, crop: tuple[int, int, int, int] | None = None) -> np.ndarray:
    image = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise ValueError(f"Could not read gel image: {path}")
    if crop is not None:
        x, y, width, height = crop
        if min(x, y) < 0 or width <= 0 or height <= 0 or x + width > image.shape[1] or y + height > image.shape[0]:
            raise ValueError(f"Crop {crop} lies outside image dimensions {image.shape[::-1]}")
        image = image[y:y + height, x:x + width]
    if min(image.shape) < 20:
        raise ValueError("Gel image or crop is too small")
    return image


def detect_lanes(image: np.ndarray, n_lanes: int, prominence: float = 2.5) -> tuple[list[int], np.ndarray, list[str]]:
    """Detect lane centers from high-frequency column signal at full resolution."""
    height, width = image.shape
    normalized = image.astype(np.float32) / 255.0
    background = cv2.GaussianBlur(normalized, (0, 0), max(3, width / 30), max(3, height / 30))
    residual = normalized - background
    column_signal = np.mean(np.maximum(residual, 0), axis=0)
    column_signal = gaussian_filter1d(column_signal, max(1, width / (n_lanes * 12)))
    noise = np.median(np.abs(column_signal - np.median(column_signal))) / 0.6745 + 1e-8
    peaks, properties = find_peaks(
        column_signal, distance=max(3, int(width / (n_lanes * 1.5))),
        prominence=max(prominence * noise, 0.002),
    )
    warnings = []
    if len(peaks) < n_lanes:
        if n_lanes == 1:
            return [int(np.argmax(column_signal))], residual, []
        # Gel lanes are commonly loaded on a regular comb. Fit offset and
        # spacing to the observed signal rather than assuming image midpoints.
        base = width / n_lanes
        best_score, best_centers = -np.inf, None
        max_spacing = min(base * 1.2, (width - 1) / (n_lanes - 1))
        for spacing in np.linspace(base * 0.8, max_spacing, 81):
            max_start = min(base * 1.5, width - 1 - spacing * (n_lanes - 1))
            min_start = base * 0.4
            if max_start < min_start:
                continue
            for start in np.linspace(min_start, max_start, 61):
                candidate = start + np.arange(n_lanes) * spacing
                score = float(np.interp(candidate, np.arange(width), column_signal).mean())
                if score > best_score:
                    best_score, best_centers = score, candidate
        centers = np.rint(best_centers).astype(int).tolist()
        warnings.append(f"Lane peaks were ambiguous; fitted a regular {n_lanes}-lane grid to the image signal")
    else:
        selected = np.argsort(properties["prominences"])[-n_lanes:]
        centers = []
        radius = max(3, int(width / (n_lanes * 5)))
        for peak in sorted(peaks[selected]):
            local_x = np.arange(max(0, peak - radius), min(width, peak + radius + 1))
            weights = np.maximum(column_signal[local_x] - column_signal[local_x].min(), 0)
            centers.append(round(np.average(local_x, weights=weights)) if weights.sum() else int(peak))
    return centers, residual, warnings


def _subpixel_peak(trace: np.ndarray, index: int) -> float:
    if 0 < index < len(trace) - 1:
        left, center, right = trace[index - 1:index + 2]
        denominator = left - 2 * center + right
        if denominator < 0:
            return float(index + 0.5 * (left - right) / denominator)
    return float(index)


def detect_bands(image: np.ndarray, residual: np.ndarray, centers: Sequence[int], prominence_scale: float = 2.5, min_band_distance: int | None = None) -> list[pd.DataFrame]:
    """Quantify bands from background-corrected, per-lane 1-D traces."""
    height, width = residual.shape
    normalized = image.astype(np.float32) / 255.0
    all_lanes = []
    for lane_index, center in enumerate(centers):
        left_neighbor = centers[lane_index - 1] if lane_index else 2 * center - centers[lane_index + 1]
        right_neighbor = centers[lane_index + 1] if lane_index + 1 < len(centers) else 2 * center - centers[lane_index - 1]
        half_width = max(2, int(min(center - left_neighbor, right_neighbor - center) * 0.28))
        left, right = max(0, center - half_width), min(width, center + half_width + 1)
        trace = np.mean(normalized[:, left:right], axis=1)
        # Morphological opening follows the background without the negative
        # Gaussian halos that can become false peaks between strong bands.
        baseline = grey_opening(trace, size=max(5, height // 20))
        corrected = gaussian_filter1d(trace - baseline, 1.0)
        noise = np.median(np.abs(corrected - np.median(corrected))) / 0.6745 + 1e-8
        dynamic = np.percentile(corrected, 99) - np.percentile(corrected, 10)
        min_prominence = max(prominence_scale * noise, 0.05 * dynamic, 1e-5)
        distance = max(2, height // 300) if min_band_distance is None else min_band_distance
        if distance < 1:
            raise ValueError("min_band_distance must be at least one pixel")
        peaks, properties = find_peaks(corrected, prominence=min_prominence, distance=distance)
        widths = peak_widths(corrected, peaks, rel_height=0.5)[0] if len(peaks) else []
        rows = []
        for band_index, peak in enumerate(peaks):
            radius = max(2, int(np.ceil(widths[band_index])))
            local = corrected[max(0, peak - radius):min(height, peak + radius + 1)]
            rows.append({
                "lane": lane_index + 1, "pixel_x": center,
                "pixel_y": _subpixel_peak(corrected, int(peak)),
                "integrated_intensity": float(np.sum(np.maximum(local, 0))),
                "prominence": float(properties["prominences"][band_index]),
            })
        lane = pd.DataFrame(rows)
        if not lane.empty:
            lane["relative_intensity"] = lane["integrated_intensity"] / lane["integrated_intensity"].sum()
        all_lanes.append(lane)
    return all_lanes


def calibrate_ladder(pixel_y: Sequence[float], sizes_bp: Sequence[int]) -> PchipInterpolator:
    """Monotonic log-size interpolation; a count/order mismatch is an error."""
    if len(pixel_y) != len(sizes_bp):
        raise ValueError(
            f"Detected {len(pixel_y)} ladder bands but {len(sizes_bp)} sizes were supplied. "
            "Inspect the QC image and supply --ladder-y explicitly or adjust detection."
        )
    if len(sizes_bp) < 2 or min(sizes_bp) <= 0:
        raise ValueError("At least two positive ladder sizes are required")
    order = np.argsort(pixel_y)
    y = np.asarray(pixel_y, dtype=float)[order]
    sizes = np.asarray(sizes_bp, dtype=float)[order]
    if np.any(np.diff(y) <= 0) or np.any(np.diff(sizes) >= 0):
        raise ValueError("Ladder sizes must strictly decrease from the top to the bottom of the gel")
    return PchipInterpolator(y, np.log10(sizes), extrapolate=False)



def _draw_qc(path: Path, image: np.ndarray, bands: pd.DataFrame, centers: Sequence[int], ladder_lanes: set[int]) -> None:
    canvas = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)
    for lane, center in enumerate(centers, 1):
        cv2.putText(canvas, str(lane), (max(0, center - 8), 18), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 200, 0), 1, cv2.LINE_AA)
    dense = len(centers) > 8
    for band in bands.itertuples():
        is_ladder = int(band.lane) in ladder_lanes
        color = (255, 180, 0) if is_ladder else (80, 240, 80)
        y = round(band.pixel_y)
        cv2.line(canvas, (int(band.pixel_x) - 8, y), (int(band.pixel_x) + 8, y), color, 1)
        if not dense or is_ladder:
            label = f"{band.fragment_bp:.0f} bp" if pd.notna(band.fragment_bp) else f"y={band.pixel_y:.1f}"
            cv2.putText(canvas, label, (int(band.pixel_x) + 10, y + 3), cv2.FONT_HERSHEY_SIMPLEX, 0.34, color, 1, cv2.LINE_AA)
    cv2.imwrite(str(path), canvas)


def analyze_gel(
    image_path: str | Path,
    enzyme: str,
    output_dir: str | Path,
    n_lanes: int,
    ladder_lane: int | None = None,
    secondary_ladder_lane: int | None = None,
    ladder_sizes: Sequence[int] | None = None,
    ladder_y: Sequence[float] | None = None,
    secondary_ladder_y: Sequence[float] | None = None,
    sample_names: Sequence[str] | None = None,
    crop: tuple[int, int, int, int] | None = None,
    lane_centers: Sequence[int] | None = None,
    polarity: str = "bright",
    prominence: float = 2.5,
    lane_roles: dict[int, str] | None = None,
    sample_types: dict[int, str] | None = None,
    min_band_distance: int | None = None,
    experiment_metadata: dict | None = None,
) -> GelResult:
    """Detect and quantify all lanes, optionally converting migration to bp."""
    if n_lanes < 1:
        raise ValueError("n_lanes must be positive")
    original = load_grayscale(image_path)
    image = load_grayscale(image_path, crop)
    if polarity not in {"bright", "dark"}:
        raise ValueError("polarity must be 'bright' or 'dark'")
    if polarity == "dark":
        image = 255 - image
    if lane_centers is None:
        centers, residual, warnings = detect_lanes(image, n_lanes, prominence)
    else:
        if len(lane_centers) != n_lanes:
            raise ValueError("The number of manual lane centers must equal n_lanes")
        centers = sorted(int(value) for value in lane_centers)
        normalized = image.astype(np.float32) / 255.0
        background = cv2.GaussianBlur(normalized, (0, 0), max(3, image.shape[1] / 30), max(3, image.shape[0] / 30))
        residual, warnings = normalized - background, []
    lanes = detect_bands(image, residual, centers, prominence, min_band_distance)

    calibrations: list[PchipInterpolator | None] = [None] * n_lanes
    ladder_lanes = {lane for lane in (ladder_lane, secondary_ladder_lane) if lane is not None}
    if ladder_sizes is not None:
        if ladder_lane is None or not 1 <= ladder_lane <= n_lanes:
            raise ValueError("ladder_lane must be supplied using one-based lane numbering")
        if secondary_ladder_lane is not None and (
            not 1 <= secondary_ladder_lane <= n_lanes or secondary_ladder_lane == ladder_lane
        ):
            raise ValueError("secondary_ladder_lane must be a different, valid one-based lane")
        primary_y = ladder_y if ladder_y is not None else (
            lanes[ladder_lane - 1]["pixel_y"].tolist() if not lanes[ladder_lane - 1].empty else []
        )
        calibrate_ladder(primary_y, ladder_sizes)  # Validate before processing samples.
        secondary_y = None
        if secondary_ladder_lane is not None:
            secondary_y = secondary_ladder_y if secondary_ladder_y is not None else (
                lanes[secondary_ladder_lane - 1]["pixel_y"].tolist() if not lanes[secondary_ladder_lane - 1].empty else []
            )
            calibrate_ladder(secondary_y, ladder_sizes)
        for index in range(n_lanes):
            reference_y = np.asarray(primary_y, dtype=float)
            if secondary_y is not None:
                fraction = (index + 1 - ladder_lane) / (secondary_ladder_lane - ladder_lane)
                reference_y = reference_y + fraction * (np.asarray(secondary_y) - reference_y)
            calibrations[index] = calibrate_ladder(reference_y, ladder_sizes)
    elif ladder_y is not None or secondary_ladder_y is not None:
        raise ValueError("Ladder positions were supplied without ladder sizes")

    sample_indices = [index for index in range(n_lanes) if index + 1 not in ladder_lanes]
    if sample_names is not None and len(sample_names) != len(sample_indices):
        raise ValueError(f"Expected {len(sample_indices)} sample names (excluding ladders), got {len(sample_names)}")
    names = {index: (sample_names[position] if sample_names else f"sample_{position + 1}") for position, index in enumerate(sample_indices)}
    if len(set(names.values())) != len(names):
        raise ValueError("Sample names must be unique")

    frames = []
    for index, lane in enumerate(lanes):
        if lane.empty:
            warnings.append(f"No bands detected in lane {index + 1}")
            continue
        lane = lane.copy()
        lane["sample"] = "ladder" if index + 1 in ladder_lanes else names[index]
        lane["enzyme"] = enzyme
        default_role = "ladder" if index + 1 in ladder_lanes else "sample"
        lane["role"] = (lane_roles or {}).get(index + 1, default_role)
        lane["sample_type"] = (sample_types or {}).get(index + 1, "unknown")
        calibration = calibrations[index]
        if calibration is None:
            lane["fragment_bp"] = np.nan
        else:
            raw_y = lane["pixel_y"].to_numpy()
            lower, upper = float(calibration.x[0]), float(calibration.x[-1])
            # A manually supplied coordinate can differ by a fraction of a pixel
            # from the detected peak; only genuinely out-of-range bands are lost.
            outside = (raw_y < lower - 1) | (raw_y > upper + 1)
            log_sizes = calibration(np.clip(raw_y, lower, upper))
            lane["fragment_bp"] = np.power(10, log_sizes)
            if outside.any():
                warnings.append(f"Excluded {int(outside.sum())} band(s) outside the ladder range in lane {index + 1}")
            lane = lane.loc[~outside]
        frames.append(lane)
    columns = ["sample", "enzyme", "lane", "role", "sample_type", "pixel_x", "pixel_y", "fragment_bp", "integrated_intensity", "relative_intensity", "prominence"]
    bands = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
    bands = bands[columns]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bands.to_csv(output_dir / "bands.csv", index=False, float_format="%.4f")
    bands.dropna(subset=["fragment_bp"])[["sample", "enzyme", "fragment_bp"]].to_csv(
        output_dir / "profiles.csv", index=False, float_format="%.3f"
    )
    _draw_qc(output_dir / "qc_overlay.png", image, bands, centers, ladder_lanes)
    metadata = {
        "image": str(image_path), "original_shape": list(original.shape), "crop": crop,
        "polarity": polarity, "lane_centers": centers, "ladder_lanes": sorted(ladder_lanes),
        "ladder_sizes": list(ladder_sizes) if ladder_sizes is not None else None,
        "ladder_y": list(ladder_y) if ladder_y is not None else None,
        "secondary_ladder_y": list(secondary_ladder_y) if secondary_ladder_y is not None else None,
        "lane_roles": lane_roles or {}, "sample_types": sample_types or {},
        "min_band_distance": min_band_distance, "experiment": experiment_metadata or {},
        "calibration": (
            "Per-lane PCHIP of log10(bp), linearly corrected between two ladders"
            if secondary_ladder_lane is not None else
            "PCHIP interpolation of log10(bp)" if ladder_sizes is not None else None
        ),
        "warnings": warnings,
    }
    (output_dir / "analysis.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return GelResult(bands=bands, lane_centers=centers, warnings=warnings, output_dir=output_dir)
