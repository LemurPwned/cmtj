import pytest

from cmtj import Layer, Junction, CVector


def test_throws_error_if_layers_have_duplicate_ids():
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    alpha = 0.005
    surface = 0

    layers = [
        Layer(
            "free",
            mag=CVector(0, 0, 1),
            anis=CVector(0, 0, 1),
            Ms=1.65,
            thickness=3e-9,
            cellSurface=surface,
            demagTensor=demag,
            damping=alpha,
        ),
        Layer(
            "free",
            mag=CVector(0, 0, 1),
            anis=CVector(0, 1.0, 0),
            Ms=1.65,
            thickness=3e-9,
            cellSurface=surface,
            demagTensor=demag,
            damping=alpha,
        ),
    ]
    with pytest.raises(ValueError, match="Layers must have unique ids!"):
        Junction(layers)


def test_throws_if_no_layers_are_passed():
    with pytest.raises(ValueError, match="Passed a zero length Layer vector!"):
        Junction([])
