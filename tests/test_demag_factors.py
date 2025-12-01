import math

import pytest

from cmtj.utils.demag import rectangular_prism_demag_factors


def test_cube_has_uniform_demag():
    dx, dy, dz = rectangular_prism_demag_factors(1.0, 1.0, 1.0)
    expected = pytest.approx(1.0 / 3.0, rel=1e-9, abs=1e-12)
    assert dx == expected
    assert dy == expected
    assert dz == expected


@pytest.mark.parametrize(
    "width,length,height",
    [
        (1.0, 1.0, 0.01),
        (4.0, 2.0, 0.25),
        (10.0, 10.0, 1.0),
    ],
)
def test_demag_factors_sum_to_one(width, length, height):
    dx, dy, dz = rectangular_prism_demag_factors(width, length, height)
    assert math.isclose(dx + dy + dz, 1.0, rel_tol=1e-10, abs_tol=1e-12)


def test_thin_film_limit_matches_reference_value():
    dx, dy, dz = rectangular_prism_demag_factors(1.0, 1.0, 0.01)
    # Reference values obtained from Aharoni (1998) formula.
    assert dx == pytest.approx(0.01698020892105028, rel=1e-12, abs=1e-12)
    assert dy == pytest.approx(0.01698020892105028, rel=1e-12, abs=1e-12)
    assert dz == pytest.approx(0.9660395821578994, rel=1e-12, abs=1e-12)
