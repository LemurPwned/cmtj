import pytest
import cmtj


def test_reference_layers():
    """Test setting and getting reference layers"""

    # Create a basic SOT layer
    mag = cmtj.CVector(0, 0, 1)
    anis = cmtj.CVector(0, 0, 1)
    demag = [cmtj.CVector(0, 0, 0), cmtj.CVector(0, 0, 0), cmtj.CVector(0, 0, 1)]

    layer = cmtj.Layer(
        id="free",
        mag=mag,
        anis=anis,
        Ms=1,
        thickness=2e-9,
        cellSurface=1e-18,
        demagTensor=demag,
        damping=0.01,
    )

    # Test setting reference layer magnetization
    ref_mag = cmtj.CVector(1, 0, 0)
    layer.setReferenceLayer(ref_mag)
    assert layer.getReferenceLayer() == ref_mag

    # Test setting secondary reference layer
    sec_ref_mag = cmtj.CVector(0, 1, 0)
    layer.setSecondaryReferenceLayer(sec_ref_mag)
    assert layer.getSecondaryReferenceLayer() == sec_ref_mag


def test_secondary_torques():
    """Test setting secondary torque drivers"""

    # Create a basic SOT layer
    mag = cmtj.CVector(0, 0, 1)
    anis = cmtj.CVector(0, 0, 1)
    demag = [cmtj.CVector(0, 0, 0)] * 3

    layer = cmtj.Layer(
        id="free",
        mag=mag,
        anis=anis,
        Ms=1,
        thickness=2e-9,
        cellSurface=1e-18,
        demagTensor=demag,
        damping=0.01,
    )

    # Create drivers
    primary_fl_driver = cmtj.ScalarDriver.getConstantDriver(0.1)
    primary_dl_driver = cmtj.ScalarDriver.getConstantDriver(0.2)
    secondary_fl_driver = cmtj.ScalarDriver.getConstantDriver(0.3)
    secondary_dl_driver = cmtj.ScalarDriver.getConstantDriver(0.4)

    # Set primary and secondary torque drivers
    layer.setPrimaryTorqueDrivers(primary_fl_driver, primary_dl_driver)
    layer.setSecondaryTorqueDrivers(secondary_fl_driver, secondary_dl_driver)

    # Set reference layers for torque effects
    reference = cmtj.CVector(1, 0, 0)
    secondary_reference = cmtj.CVector(0, 1, 0)
    layer.setReferenceLayer(reference)
    layer.setSecondaryReferenceLayer(secondary_reference)

    # Create junction for simulation
    junction = cmtj.Junction([layer])

    # Run a short simulation to confirm there's no error
    junction.runSimulation(1e-9, 1e-13, 1e-12)

    # Logs should have been collected
    log = junction.getLog()
    assert "time" in log
    assert "free_mx" in log
    assert "free_my" in log
    assert "free_mz" in log


def test_junction_torque_drivers():
    """Test setting torque drivers via junction interface"""

    # Create two basic SOT layers
    mag1 = cmtj.CVector(0, 0, 1)
    mag2 = cmtj.CVector(1, 0, 0)
    anis = cmtj.CVector(0, 0, 1)
    demag = [cmtj.CVector(0, 0, 0)] * 3

    layer1 = cmtj.Layer(
        id="free1",
        mag=mag1,
        anis=anis,
        Ms=1,
        thickness=2e-9,
        cellSurface=1e-18,
        demagTensor=demag,
        damping=0.01,
    )

    layer2 = cmtj.Layer(
        id="free2",
        mag=mag2,
        anis=anis,
        Ms=0.8,
        thickness=1e-9,
        cellSurface=1e-18,
        demagTensor=demag,
        damping=0.02,
    )

    # Create junction with both layers
    junction = cmtj.Junction([layer1, layer2])

    # Create drivers for both primary and secondary torques
    fl_driver = cmtj.ScalarDriver.getConstantDriver(0.1)
    dl_driver = cmtj.ScalarDriver.getConstantDriver(0.2)
    sec_fl_driver = cmtj.ScalarDriver.getConstantDriver(0.3)
    sec_dl_driver = cmtj.ScalarDriver.getConstantDriver(0.4)

    # Set primary torque drivers via junction interface
    junction.setLayerFieldLikeTorqueDriver("free1", fl_driver)
    junction.setLayerDampingLikeTorqueDriver("free1", dl_driver)

    # Set secondary torque drivers via junction interface
    junction.setLayerSecondaryFieldLikeTorqueDriver("free1", sec_fl_driver)
    junction.setLayerSecondaryDampingLikeTorqueDriver("free1", sec_dl_driver)

    # Set combined torque drivers for the second layer
    junction.setLayerPrimaryTorqueDrivers("free2", fl_driver, dl_driver)
    junction.setLayerSecondaryTorqueDrivers("free2", sec_fl_driver, sec_dl_driver)

    # Set reference layers
    reference = cmtj.CVector(1, 0, 0)
    junction.setLayerReferenceLayer("free1", reference)
    junction.setLayerReferenceLayer("free2", reference)

    # Set reference type for the second layer
    junction.setLayerReferenceType("free2", cmtj.Reference.fixed)

    # Run a short simulation to confirm there's no error
    junction.runSimulation(1e-9, 1e-13, 1e-12)

    # Logs should have been collected
    log = junction.getLog()
    assert "time" in log
    assert "free1_mx" in log
    assert "free1_my" in log
    assert "free1_mz" in log
    assert "free2_mx" in log
    assert "free2_my" in log
    assert "free2_mz" in log


def test_must_raise_if_sot_driver_set_for_sot_layer():
    """Test that an error is raised if a SOT driver is set for a SOT layer"""

    layer = cmtj.Layer.createSOTLayer(
        id="free",
        mag=cmtj.CVector(0, 0, 1),
        anis=cmtj.CVector(0, 0, 1),
        Ms=1,
        thickness=2e-9,
        cellSurface=1e-18,
        demagTensor=[cmtj.CVector(0, 0, 0)] * 3,
        damping=0.01,
    )

    # Set a SOT driver
    driver = cmtj.ScalarDriver.getConstantDriver(0.1)

    # This should raise an error
    with pytest.raises(RuntimeError):
        layer.setFieldLikeTorqueDriver(driver)

    # This should raise an error
    with pytest.raises(RuntimeError):
        layer.setDampingLikeTorqueDriver(driver)
