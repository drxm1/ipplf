"""Generate random initial-potential fixtures for RandomPotentialHistoryFixture.h.

Each file stores A_nm1 followed by A_n: two potential time levels, each with
scalar potential in component zero and vector potential in components 1-3.
Values are ordered by level, component, x, y, z (z contiguous), as little-endian
binary64. Each case uses a fresh seeded NumPy generator with a standard normal
distribution. The data is independent of boundary conditions; the C++ loader
handles optional boundary-layer masking and the solver prepares the ghosts.
"""

from pathlib import Path
import sys

import numpy as np


def main():
    """Generate the medium comparison and the two long-run box fixtures."""
    destination = Path(sys.argv[1])
    destination.mkdir(parents=True, exist_ok=True)
    for name, seed, shape in (
        ("medium", 271828, (16, 16, 16)),
        ("stability-open", 81413, (16, 16, 16)),
        ("stability-mixed", 81413, (16, 8, 16)),
    ):
        values = np.random.default_rng(seed).standard_normal((2, 4, *shape))
        values.astype("<f8").tofile(destination / f"{name}.bin")


if __name__ == "__main__":
    main()
