#include <gtest/gtest.h>
#include "../core/2.0/fm.hpp"
#include <iostream>
#include <memory>

// Test fixture for STT and SOT layer tests
class STTSOTTest : public ::testing::Test {
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
        Ms = 1.0;
        thickness = 1e-9;
        cellSurface = 1e-12;
        damping = 0.1;
    }

    CVector<double> mag;
    CVector<double> anis;
    std::vector<CVector<double>> demagTensor;
    double Ms;
    double thickness;
    double cellSurface;
    double damping;
};

// Test STT Layer Construction
TEST_F(STTSOTTest, STTConstruction) {
    // Test with basic parameters
    double SlonczewskiParam = 0.2;
    double beta = 0.05;
    double spinPolarisation = 0.7;

    LayerSTT<double> sttLayer("test_stt", mag, anis, Ms, thickness, cellSurface,
                             demagTensor, damping, SlonczewskiParam, beta, spinPolarisation);

    // Check basic properties
    EXPECT_EQ(sttLayer.getMagnetisation(), mag);
    EXPECT_EQ(sttLayer.Ms, Ms);
    EXPECT_EQ(sttLayer.thickness, thickness);
    EXPECT_EQ(sttLayer.alternativeSTTSet, false);
}

// Test SOT Layer Construction
TEST_F(STTSOTTest, SOTConstruction) {
    // Test with basic parameters
    double fieldLikeTorque = 0.1;
    double dampingLikeTorque = 0.2;

    LayerSOT<double> sotLayer("test_sot", mag, anis, Ms, thickness, cellSurface,
                             demagTensor, damping, fieldLikeTorque, dampingLikeTorque);

    // Check basic properties
    EXPECT_EQ(sotLayer.getMagnetisation(), mag);
    EXPECT_EQ(sotLayer.Ms, Ms);
    EXPECT_EQ(sotLayer.thickness, thickness);
    EXPECT_EQ(sotLayer.getReferenceLayer(), CVector<double>(0, 0, 0));
    EXPECT_EQ(sotLayer.directTorques, false);
}

// Test SOT Layer with Drivers
TEST_F(STTSOTTest, SOTDrivers) {
    double fieldLikeTorque = 0.1;
    double dampingLikeTorque = 0.2;
    const std::vector<CVector<double>> demagTensor_zero = {
        CVector<double>(0.0, 0, 0),
        CVector<double>(0, 0.0, 0),
        CVector<double>(0, 0, 0.0)
    };
    LayerSOT<double> sotLayer("test_sot", mag, anis, Ms, thickness, cellSurface,
                             demagTensor_zero, damping, fieldLikeTorque, dampingLikeTorque, true);

    // Set up dynamic drivers for torques
    auto flDriver = ScalarDriver<double>::getSineDriver(0.0, 0.3, 1e9, 0.0);
    auto dlDriver = ScalarDriver<double>::getConstantDriver(0.4);

    // Set drivers
    sotLayer.setFieldLikeTorqueDriver(flDriver);
    sotLayer.setDampingLikeTorqueDriver(dlDriver);

    // Test LLG calculation with reference layer
    CVector<double> refLayer(1, 0, 0);  // Reference layer along x-axis
    sotLayer.setReferenceLayer(refLayer);
    // Calculate LLG at t=0
    double time = 0.0;
    double timeStep = 1e-13;
    CVector<double> bottomMag;
    CVector<double> topMag;

    // Calculate LLG derivative
    CVector<double> dmdt = sotLayer.calculateLLG(time, timeStep, mag, bottomMag, topMag);

    // With reference layer along x and magnetization along z,
    // SOT should produce torque components
    EXPECT_NE(dmdt[0], 0.0);
    EXPECT_NE(dmdt[1], 0.0);
}

// Test SOT Layer with Dynamic Torques
TEST_F(STTSOTTest, SOTDynamicTorques) {
    double fieldLikeTorque = 0.1;
    double dampingLikeTorque = 0.2;

    LayerSOT<double> sotLayer("test_sot", mag, anis, Ms, thickness, cellSurface,
                             demagTensor, damping, fieldLikeTorque, dampingLikeTorque, true);

    // Set reference layer
    CVector<double> refLayer(1, 0, 0);
    sotLayer.setReferenceLayer(refLayer);

    // Set sinusoidal driver for field-like torque
    auto flDriver = ScalarDriver<double>::getSineDriver(0.0, 0.3, 1e9, 0.0);
    sotLayer.setFieldLikeTorqueDriver(flDriver);

    // Test at different times
    double t1 = 0.0;            // sin(0) = 0
    double t2 = 2.5e-10;        // sin(π/2) = 1
    CVector<double> mz = CVector<double>(0, 0, 1);
    CVector<double> dmdt1 = sotLayer.calculateLLG(t1, 1e-13, mz, CVector<double>(), CVector<double>());
    CVector<double> dmdt2 = sotLayer.calculateLLG(t2, 1e-13, mz, CVector<double>(), CVector<double>());

    // The magnitude should be different at these times due to the sine driver
    double magnitude1 = dmdt1.length();
    double magnitude2 = dmdt2.length();

    // At t2, the sine driver produces maximum amplitude (0.3), so torque should be stronger
    EXPECT_NE(magnitude1, magnitude2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
