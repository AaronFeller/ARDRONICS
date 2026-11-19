from pathlib import Path

import cv2
import numpy as np
import pytest

LADDER_SIZES = [1500, 1000, 700, 500, 300, 100]
SAMPLE_SIZES = [[1200, 700, 250], [1000, 500, 200], [1200, 500]]


def migration_y(size):
    return 1100 - 250 * np.log10(size)


@pytest.fixture
def synthetic_gel(tmp_path: Path):
    height, width = 720, 900
    centers = [120, 340, 560, 780]
    yy, xx = np.mgrid[:height, :width]
    image = np.full((height, width), 12.0)
    image += np.linspace(0, 8, height)[:, None]
    for lane, center in enumerate(centers):
        sizes = LADDER_SIZES if lane == 0 else SAMPLE_SIZES[lane - 1]
        for band, size in enumerate(sizes):
            amplitude = 220 - 15 * band
            image += amplitude * np.exp(
                -0.5 * ((xx - center) / 28) ** 2 - 0.5 * ((yy - migration_y(size)) / 3.0) ** 2
            )
    rng = np.random.default_rng(42)
    image += rng.normal(0, 0.7, image.shape)
    path = tmp_path / "synthetic.png"
    cv2.imwrite(str(path), np.clip(image, 0, 255).astype(np.uint8))
    return path
