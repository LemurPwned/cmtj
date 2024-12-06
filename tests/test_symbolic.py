import numpy as np
import sympy as sym
from cmtj.models.general_sb import LayerSB, Solver, VectorObj
import sympy as sym


def test_layer_energy():
    # create a test layer
    layer = LayerSB(
        _id=1, thickness=1.0, Kv=VectorObj(theta=0, phi=0.0, mag=1.0), Ks=1.0, Ms=1.0
    )

    # create test values for the input parameters
    H = sym.ImmutableMatrix([0, 0, 1])
    J1top = 0.0
    J1bottom = 0.0
    J2top = 0.0
    J2bottom = 0.0
    top_layer = None
    down_layer = None

    # calculate the energy of the layer
    energy = layer.total_symbolic_layer_energy(
        H, J1top, J1bottom, J2top, J2bottom, top_layer, down_layer
    )

    # check that the energy is a sympy expression
    assert isinstance(energy, sym.Expr)


def test_solver_init():
    # create test layers
    layer1 = LayerSB(
        _id=0, thickness=1.0, Kv=VectorObj(theta=00, phi=0.0, mag=1.0), Ks=1.0, Ms=1.0
    )
    layer2 = LayerSB(
        _id=1, thickness=1.0, Kv=VectorObj(theta=0, phi=0.0, mag=1.0), Ks=1.0, Ms=1.0
    )
    layers = [layer1, layer2]

    # create test values for J1 and J2
    J1 = [0.0]
    J2 = [0.0]

    # create a solver object
    solver = Solver(layers=layers, J1=J1, J2=J2)

    # check that the solver object was created successfully
    assert isinstance(solver, Solver)
