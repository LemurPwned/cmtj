import pytest
import cmtj


def test_gyromagnetic_ratio_runtime_change():
    """Test that gyromagnetic ratio can be changed during runtime and the change is properly reflected."""
    # Get the original value
    original_gyro = cmtj.constants.PhysicalConstants.gyromagnetic_ratio()
    
    # Should be the default value from constants.hpp
    expected_default = 220880.0
    assert original_gyro == expected_default, f"Expected default gyromagnetic ratio {expected_default}, got {original_gyro}"
    
    # Set a new value
    new_gyro = 300000.0
    cmtj.constants.PhysicalConstants.set_gyromagnetic_ratio(new_gyro)
    
    # Check that the value has changed
    current_gyro = cmtj.constants.PhysicalConstants.gyromagnetic_ratio()
    assert current_gyro == new_gyro, f"Expected gyromagnetic ratio {new_gyro}, got {current_gyro}"
    
    # Set another different value to ensure it's truly dynamic
    another_gyro = 150000.0
    cmtj.constants.PhysicalConstants.set_gyromagnetic_ratio(another_gyro)
    
    # Check that the value has changed again
    current_gyro = cmtj.constants.PhysicalConstants.gyromagnetic_ratio()
    assert current_gyro == another_gyro, f"Expected gyromagnetic ratio {another_gyro}, got {current_gyro}"
    
    # Reset to defaults and verify
    cmtj.constants.PhysicalConstants.resetToDefaults()
    reset_gyro = cmtj.constants.PhysicalConstants.gyromagnetic_ratio()
    assert reset_gyro == expected_default, f"Expected reset gyromagnetic ratio {expected_default}, got {reset_gyro}"


def test_all_physical_constants_runtime_changes():
    """Test that all physical constants can be changed during runtime."""
    # Store original values
    original_values = {
        'magnetic_permeability': cmtj.constants.PhysicalConstants.magnetic_permeability(),
        'gyromagnetic_ratio': cmtj.constants.PhysicalConstants.gyromagnetic_ratio(),
        'TtoAm': cmtj.constants.PhysicalConstants.TtoAm(),
        'hbar': cmtj.constants.PhysicalConstants.hbar(),
        'elementary_charge': cmtj.constants.PhysicalConstants.elementary_charge(),
        'boltzmann_constant': cmtj.constants.PhysicalConstants.boltzmann_constant()
    }
    
    # Define new test values (different from defaults)
    new_values = {
        'magnetic_permeability': 15e-7,
        'gyromagnetic_ratio': 250000.0,
        'TtoAm': 800000.0,
        'hbar': 1.1e-34,
        'elementary_charge': 1.7e-19,
        'boltzmann_constant': 1.5e-23
    }
    
    try:
        # Set new values
        cmtj.constants.PhysicalConstants.set_magnetic_permeability(new_values['magnetic_permeability'])
        cmtj.constants.PhysicalConstants.set_gyromagnetic_ratio(new_values['gyromagnetic_ratio'])
        cmtj.constants.PhysicalConstants.set_TtoAm(new_values['TtoAm'])
        cmtj.constants.PhysicalConstants.set_hbar(new_values['hbar'])
        cmtj.constants.PhysicalConstants.set_elementary_charge(new_values['elementary_charge'])
        cmtj.constants.PhysicalConstants.set_boltzmann_constant(new_values['boltzmann_constant'])
        
        # Verify all values have changed
        assert cmtj.constants.PhysicalConstants.magnetic_permeability() == new_values['magnetic_permeability']
        assert cmtj.constants.PhysicalConstants.gyromagnetic_ratio() == new_values['gyromagnetic_ratio']
        assert cmtj.constants.PhysicalConstants.TtoAm() == new_values['TtoAm']
        assert cmtj.constants.PhysicalConstants.hbar() == new_values['hbar']
        assert cmtj.constants.PhysicalConstants.elementary_charge() == new_values['elementary_charge']
        assert cmtj.constants.PhysicalConstants.boltzmann_constant() == new_values['boltzmann_constant']
        
    finally:
        # Always reset to defaults to avoid affecting other tests
        cmtj.constants.PhysicalConstants.resetToDefaults()
        
        # Verify reset worked
        assert cmtj.constants.PhysicalConstants.magnetic_permeability() == original_values['magnetic_permeability']
        assert cmtj.constants.PhysicalConstants.gyromagnetic_ratio() == original_values['gyromagnetic_ratio']
        assert cmtj.constants.PhysicalConstants.TtoAm() == original_values['TtoAm']
        assert cmtj.constants.PhysicalConstants.hbar() == original_values['hbar']
        assert cmtj.constants.PhysicalConstants.elementary_charge() == original_values['elementary_charge']
        assert cmtj.constants.PhysicalConstants.boltzmann_constant() == original_values['boltzmann_constant']


def test_gyromagnetic_ratio_persistence_across_operations():
    """Test that the changed gyromagnetic ratio persists across different operations."""
    # Set a custom value
    custom_gyro = 180000.0
    cmtj.constants.PhysicalConstants.set_gyromagnetic_ratio(custom_gyro)
    
    # Verify it's set
    assert cmtj.constants.PhysicalConstants.gyromagnetic_ratio() == custom_gyro
    
    # Perform some other constant operations
    original_mu = cmtj.constants.PhysicalConstants.magnetic_permeability()
    cmtj.constants.PhysicalConstants.set_magnetic_permeability(original_mu * 1.1)
    
    # Verify gyromagnetic ratio is still the custom value
    assert cmtj.constants.PhysicalConstants.gyromagnetic_ratio() == custom_gyro
    
    # Reset and verify
    cmtj.constants.PhysicalConstants.resetToDefaults()
    assert cmtj.constants.PhysicalConstants.gyromagnetic_ratio() == 220880.0


@pytest.fixture(autouse=True)
def reset_constants_after_each_test():
    """Ensure constants are reset to defaults after each test to avoid test interference."""
    yield  # This runs the test
    # Cleanup: reset to defaults after each test
    cmtj.constants.PhysicalConstants.resetToDefaults()
