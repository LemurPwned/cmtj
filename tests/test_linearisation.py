from pytest import fixture
import pytest
from cmtj.utils import TtoAm
from cmtj.models.general_sb import LayerSB, Solver
from cmtj.utils.general import VectorObj
import numpy as np


@fixture
def _sample_system_spherical():
    th = 1e-9
    Ms = 0.73 * TtoAm
    Ks = 300e3
    J1 = -180e-6
    l1 = LayerSB(
        _id=0,
        Ms=Ms,
        Kv=VectorObj(0, 0, 0),
        Ks=Ks,
        thickness=th,
        Ndemag=VectorObj(0, 0, 1),
    )

    l2 = LayerSB(
        _id=1,
        Ms=Ms,
        Kv=VectorObj(0, 0, 0),
        Ks=Ks,
        thickness=th,
        Ndemag=VectorObj(0, 0, 1),
    )

    return l1, l2, J1


# Expected positive frequencies for different configurations
# Doing non zero h for better stability
NON_ZERO_H_TRUE_P = [7.12, 10.24]
NON_ZERO_H_TRUE_AP = [13.05, 16.57]


@pytest.mark.parametrize(
    "prefer_numerical_roots,use_LU_decomposition",
    [
        (False, False),  # Analytical + no LU
        (True, False),  # Numerical + no LU (no AP frequencies)
        (True, True),  # Numerical + LU (no AP frequencies)
        (False, True),  # Analytical + LU
    ],
)
def test_linearised_frequencies_positive_only(
    _sample_system_spherical,
    prefer_numerical_roots,
    use_LU_decomposition,
):
    """Test linearised frequency calculation with different solver configurations, checking only positive frequencies."""
    l1, l2, J1 = _sample_system_spherical
    Hz = 50e3
    solver = Solver(
        layers=[l1, l2],
        J1=[J1],
        J2=[0],
        H=VectorObj(0, 0, Hz),
        prefer_numerical_roots=prefer_numerical_roots,
        use_LU_decomposition=use_LU_decomposition,
    )

    freqsP, freqsAP = solver.solve_linearised_frequencies(
        H=None, linearisation_axis="z"
    )

    # Filter for positive frequencies only
    pos_freqsAP = sorted([f for f in freqsAP if f > 0])
    pos_freqsP = sorted([f for f in freqsP if f > 0])

    # Check P frequencies (should be consistent across all methods)
    assert (
        len(pos_freqsP) > 0
    ), f"No positive P frequencies found for config: numerical={prefer_numerical_roots}, LU={use_LU_decomposition}"
    assert np.allclose(
        np.around(pos_freqsP, 2), NON_ZERO_H_TRUE_P,
        rtol=1e-2,
        atol=1e-2,
    ), f"P frequencies mismatch for config: numerical={prefer_numerical_roots}, LU={use_LU_decomposition}. Got: {np.around(pos_freqsP, 2)}, Expected: {NON_ZERO_H_TRUE_P}"

    # Analytical methods should find AP frequencies
    assert (
        len(pos_freqsAP) > 0
    ), f"No positive AP frequencies found for analytical config: numerical={prefer_numerical_roots}, LU={use_LU_decomposition}"
    assert np.allclose( 
        np.around(pos_freqsAP, 2), NON_ZERO_H_TRUE_AP,
        rtol=1e-2,
        atol=1e-2,
    ), f"AP frequencies mismatch for config: numerical={prefer_numerical_roots}, LU={use_LU_decomposition}. Got: {np.around(pos_freqsAP, 2)}, Expected: {NON_ZERO_H_TRUE_AP}"
