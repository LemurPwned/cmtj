import numpy as np
import pytest
from cmtj import (
    CVector,
    Layer,
    Junction,
    reservoir,
)


def create_test_junction(pos=(0, 0, 0), Ms=1e6):
    """Helper function to create a test junction"""
    mag = CVector(0, 0, 1)
    anis = CVector(0, 0, 1)
    demag = [CVector(0, 0, 0)] * 3

    layer = Layer(
        "free",
        mag,
        anis,
        Ms=Ms,
        thickness=2e-9,
        cellSurface=np.pi * (20e-9) ** 2,
        demagTensor=demag,
    )
    return Junction([layer]), CVector(*pos)


def test_group_interaction_initialization():
    """Test basic initialization of GroupInteraction"""
    j1, pos1 = create_test_junction((0, 0, 0))
    j2, pos2 = create_test_junction((100e-9, 0, 0))

    # Should initialize successfully
    group = reservoir.GroupInteraction([pos1, pos2], [j1, j2])

    # Should fail with mismatched sizes
    with pytest.raises(RuntimeError):
        reservoir.GroupInteraction([pos1], [j1, j2])

    # Should fail with empty lists
    with pytest.raises(RuntimeError):
        reservoir.GroupInteraction([], [])

    # Should fail with duplicate positions
    with pytest.raises(RuntimeError):
        reservoir.GroupInteraction([pos1, pos1], [j1, j2])


def test_dipole_interactions():
    """Test different dipole interaction functions"""
    j1, pos1 = create_test_junction((0, 0, 0))
    j2, pos2 = create_test_junction((100e-9, 0, 0))

    # Test null interaction
    h_null = reservoir.nullDipoleInteraction(
        pos1, pos2, j1.getLayer("free"), j2.getLayer("free")
    )
    assert h_null.x == 0 and h_null.y == 0 and h_null.z == 0

    # Test regular dipole interaction
    h_dipole = reservoir.computeDipoleInteraction(
        pos1, pos2, j1.getLayer("free"), j2.getLayer("free")
    )
    assert isinstance(h_dipole, CVector)

    # Test Noumra dipole interaction
    h_noumra = reservoir.computeDipoleInteractionNoumra(
        pos1, pos2, j1.getLayer("free"), j2.getLayer("free")
    )
    assert isinstance(h_noumra, CVector)


def test_group_simulation():
    """Test running a simulation with group interaction"""
    j1, pos1 = create_test_junction((0, 0, 0))
    j2, pos2 = create_test_junction((100e-9, 0, 0))

    group = reservoir.GroupInteraction([pos1, pos2], [j1, j2])
    # Then test invalid indices
    with pytest.raises(
        (RuntimeError, IndexError, TypeError)
    ):  # Accept either exception type
        group.getLog(-1)  # Test negative index

    with pytest.raises((RuntimeError, IndexError)):  # Accept either exception type
        group.getLog(2)  # Test out of bounds index
    # Run a short simulation
    total_time = 1e-9
    time_step = 1e-12
    group.runSimulation(total_time, time_step)

    # Check logs exist
    log1 = group.getLog(0)
    log2 = group.getLog(1)
    assert len(log1["time"]) > 0
    assert len(log2["time"]) > 0

    # Clear logs
    group.clearLogs()
    log1 = group.getLog(0)
    log2 = group.getLog(1)
    with pytest.raises(KeyError):
        log1["time"]
    with pytest.raises(KeyError):
        log2["time"]


def test_interaction_functions():
    """Test setting different interaction functions"""
    j1, pos1 = create_test_junction((0, 0, 0))
    j2, pos2 = create_test_junction((100e-9, 0, 0))

    group = reservoir.GroupInteraction([pos1, pos2], [j1, j2])

    # Test setting null interaction
    group.setInteractionFunction(reservoir.nullDipoleInteraction)

    # Test setting Noumra interaction
    group.setInteractionFunction(reservoir.computeDipoleInteractionNoumra)

    # Test setting regular dipole interaction
    group.setInteractionFunction(reservoir.computeDipoleInteraction)


def test_invalid_log_access():
    """Test accessing invalid log indices"""
    j1, pos1 = create_test_junction((0, 0, 0))
    j2, pos2 = create_test_junction((100e-9, 0, 0))

    group = reservoir.GroupInteraction([pos1, pos2], [j1, j2])

    # Should raise error for invalid index
    with pytest.raises(RuntimeError):
        group.getLog(2)
