#include <gtest/gtest.h>
#include "../core/2.0/fm.hpp"
#include <iostream>
#include <memory>

// Test fixture for Layer tests
class LayerTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Basic layer parameters
        mag = CVector<double>(0, 0, 1);
        anis = CVector<double>(0, 0, 1);
        demagTensor = {
            CVector<double>(0.1, 0, 0),
            CVector<double>(0, 0.1, 0),
            CVector<double>(0, 0, 0.8)
        };
    }

    CVector<double> mag;
    CVector<double> anis;
    std::vector<CVector<double>> demagTensor;
};

// Test basic construction
TEST_F(LayerTest, Construction) {
    // Test basic construction
    Layer<double> layer("test_layer", mag, anis, 1.3, 1e-9, 1e-12, demagTensor, 0.1);

    // Check basic properties
    EXPECT_EQ(layer.getMagnetisation(), mag);
    EXPECT_EQ(layer.Ms, 1.3);
    EXPECT_EQ(layer.thickness, 1e-9);
    EXPECT_EQ(layer.cellSurface, 1e-12);
    EXPECT_DOUBLE_EQ(layer.cellVolume, 1e-9*1e-12);
}

// Test setting and getting magnetization
TEST_F(LayerTest, MagnetizationSetGet) {
    Layer<double> layer("test_layer", mag, anis, 1e6, 1e-9, 1e-12, demagTensor, 0.1);

    // Test initial magnetization
    EXPECT_EQ(layer.getMagnetisation(), mag);

    // Set new magnetization
    CVector<double> newMag(1, 0, 0);
    layer.setMagnetisation(newMag);

    // Check if correctly set and normalized
    CVector<double> retrievedMag = layer.getMagnetisation();
    EXPECT_NEAR(retrievedMag[0], 1.0, 1e-10);
    EXPECT_NEAR(retrievedMag[1], 0.0, 1e-10);
    EXPECT_NEAR(retrievedMag[2], 0.0, 1e-10);
    EXPECT_NEAR(retrievedMag.length(), 1.0, 1e-10);

    // Check exception for zero vector
    CVector<double> zeroMag(0, 0, 0);
    EXPECT_THROW(layer.setMagnetisation(zeroMag), std::runtime_error);
}

// Test driver setting and evaluation
TEST_F(LayerTest, DriverSetGet) {
    Layer<double> layer("test_layer", mag, anis, 1e6, 1e-9, 1e-12, demagTensor, 0.1);

    // Test setting constant temperature driver
    auto tempDriver = ScalarDriver<double>::getConstantDriver(300.0);
    layer.setTemperatureDriver(tempDriver);

    // Test setting anisotropy driver
    auto anisDriver = ScalarDriver<double>::getConstantDriver(5e5);
    layer.setAnisotropyDriver(anisDriver);

    // Test setting external field driver
    AxialDriver<double> fieldDriver = AxialDriver<double>::getVectorAxialDriver(0.0, 0.0, 0.1);
    layer.setExternalFieldDriver(fieldDriver);

    // Calculate fields and verify they use the driver values
    double time = 0.0;
    CVector<double> stepMag = mag;
    CVector<double> bottomMag(0, 0, 0);
    CVector<double> topMag(0, 0, 0);

    // Calculate anisotropy field
    CVector<double> anis = layer.calculateAnisotropy(stepMag, time);

    // Calculate external field
    CVector<double> extField = layer.calculateExternalField(time);

    // Verify field values
    // Anisotropy field should point along z with magnitude related to anisotropy constant
    EXPECT_NEAR(anis[0], 0.0, 1e-10);
    EXPECT_NEAR(anis[1], 0.0, 1e-10);
    EXPECT_GT(anis[2], 0.0);  // Should be positive along z

    // External field should match our driver settings
    EXPECT_NEAR(extField[0], 0.0, 1e-10);
    EXPECT_NEAR(extField[1], 0.0, 1e-10);
    EXPECT_NEAR(extField[2], 0.1, 1e-10);
}

// Test LLG dynamics with different drivers
TEST_F(LayerTest, LLGDynamics) {
    const double Ms = 1;
    const double thickness = 1e-9;
    const double cellSurface = 1e-12;
    const double damping = 0.1;
    Layer<double> layer("test_layer", mag, anis, Ms, thickness, cellSurface, demagTensor, damping);

    // Set up drivers
    auto anisDriver = ScalarDriver<double>::getConstantDriver(5e5);
    layer.setAnisotropyDriver(anisDriver);

    // Set up external field along x to cause precession
    AxialDriver<double> fieldDriver = AxialDriver<double>::getVectorAxialDriver(0.1, 0.0, 0.0);
    layer.setExternalFieldDriver(fieldDriver);

    // Initial state
    CVector<double> initialMag = layer.getMagnetisation();
    EXPECT_NEAR(initialMag[0], 0.0, 1e-10);
    EXPECT_NEAR(initialMag[1], 0.0, 1e-10);
    EXPECT_NEAR(initialMag[2], 1.0, 1e-10);

    // Calculate LLG (dm/dt) at t=0
    double time = 0.0;
    double timeStep = 1e-13;
    CVector<double> bottomMag(0, 0, 0);
    CVector<double> topMag(0, 0, 0);

    // Calculate effective field
    CVector<double> Heff = layer.calculateHeff(time, timeStep, initialMag, bottomMag, topMag);

    // Calculate LLG derivative
    CVector<double> dmdt = layer.calculateLLG(time, timeStep, initialMag, bottomMag, topMag);

    // Since field is along x and mag along z, dmdt should be along -y (right-hand rule)
    EXPECT_LT(dmdt[1], 0.0);             // y component should be negative
    EXPECT_NEAR(dmdt[2], 0.0, 1e-8);     // z component should be near zero due to damping

    // Verify damping term
}

// Test SineDriver integration
TEST_F(LayerTest, SineDriverTest) {
    Layer<double> layer("test_layer", mag, anis, 1e6, 1e-9, 1e-12, demagTensor, 0.1);

    // Create sine driver: baseValue=0, amplitude=0.1, frequency=1GHz, phase=0
    auto sineDriver = ScalarDriver<double>::getSineDriver(0.0, 0.1, 1e9, 0.0);

    // Set external field with sine driver for z component
    AxialDriver<double> fieldDriver(
        ScalarDriver<double>::getConstantDriver(0.0),  // x component
        ScalarDriver<double>::getConstantDriver(0.0),  // y component
        sineDriver  // z component
    );
    layer.setExternalFieldDriver(fieldDriver);

    // Test field at different times
    double t1 = 0.0;  // sin(0) = 0
    double t2 = 2.5e-10;  // sin(pi/2) = 1.0
    double t3 = 5.0e-10;  // sin(pi) = 0
    double t4 = 7.5e-10;  // sin(3pi/2) = -1.0

    CVector<double> field1 = layer.calculateExternalField(t1);
    CVector<double> field2 = layer.calculateExternalField(t2);
    CVector<double> field3 = layer.calculateExternalField(t3);
    CVector<double> field4 = layer.calculateExternalField(t4);

    // Check z components follow sine wave
    EXPECT_NEAR(field1[2], 0.0, 1e-10);
    EXPECT_NEAR(field2[2], 0.1, 1e-10);
    EXPECT_NEAR(field3[2], 0.0, 1e-10);
    EXPECT_NEAR(field4[2], -0.1, 1e-10);
}

// Test StepDriver integration
TEST_F(LayerTest, StepDriverTest) {
    const double Ms = 1;
    const double thickness = 1e-9;
    const double cellSurface = 1e-12;
    const double damping = 0.1;
    Layer<double> layer("test_layer", mag, anis, Ms, thickness, cellSurface, demagTensor, damping);

    // Create step driver: baseValue=0, amplitude=0.1, timeStart=1ns, timeStop=2ns
    const double amplitude = 0.1;
    const double timeStart = 1e-9;
    const double timeStop = 2e-9;
    auto stepDriver = ScalarDriver<double>::getStepDriver(0.0, amplitude, timeStart, timeStop);
    layer.setAnisotropyDriver(stepDriver);

    // Test at different times
    double t1 = 0.5e-9;  // Before step
    double t2 = 1.5e-9;  // During step
    double t3 = 2.5e-9;  // After step

    // Get anisotropy with a unit vector along z
    CVector<double> testMag(0, 0, 1);

    CVector<double> anis1 = layer.calculateAnisotropy(testMag, t1);
    CVector<double> anis2 = layer.calculateAnisotropy(testMag, t2);
    CVector<double> anis3 = layer.calculateAnisotropy(testMag, t3);

    // Before and after the step, anisotropy should be zero
    // During the step, it should be 0.1 * magnetization direction
    EXPECT_NEAR(anis1[2], 0.0, 1e-10);
    EXPECT_NEAR(anis2[2], 2*amplitude/Ms, 1e-10);
    EXPECT_NEAR(anis3[2], 0.0, 1e-10);
}

// Test PulseDriver integration
TEST_F(LayerTest, PulseDriverTest) {
    Layer<double> layer("test_layer", mag, anis, 1, 1e-9, 1e-12, demagTensor, 0.1);

    // Create pulse driver: baseValue=0, amplitude=0.1, period=2ns, cycle=0.5 (50% duty cycle)
    auto pulseDriver = ScalarDriver<double>::getPulseDriver(0.0, 0.1, 2e-9, 0.5);

    // Set as z-component of external field
    AxialDriver<double> fieldDriver(
        ScalarDriver<double>::getConstantDriver(0.0),
        ScalarDriver<double>::getConstantDriver(0.0),
        pulseDriver
    );
    layer.setExternalFieldDriver(fieldDriver);

    // Test at different times
    double t1 = 0.5e-9;  // During first pulse
    double t2 = 1.5e-9;  // Between pulses
    double t3 = 2.5e-9;  // During second pulse
    double t4 = 3.5e-9;  // Between pulses

    CVector<double> field1 = layer.calculateExternalField(t1);
    CVector<double> field2 = layer.calculateExternalField(t2);
    CVector<double> field3 = layer.calculateExternalField(t3);
    CVector<double> field4 = layer.calculateExternalField(t4);

    // Check pulse pattern in z component
    EXPECT_NEAR(field1[2], 0.1, 1e-10);  // On
    EXPECT_NEAR(field2[2], 0.0, 1e-10);  // Off
    EXPECT_NEAR(field3[2], 0.1, 1e-10);  // On
    EXPECT_NEAR(field4[2], 0.0, 1e-10);  // Off
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
