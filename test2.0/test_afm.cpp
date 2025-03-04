    #include <gtest/gtest.h>
    #include "../core/2.0/afm.hpp"
    #include <iostream>

    class AFMTest : public ::testing::Test {
    protected:
        void SetUp() override {
            // Common setup parameters for all tests
            Ms = 1400e3;
            thickness = 2e-9;
            cellSurface = 1e-14;
            damping = 0.01;

            // Initialize demagnetization tensor (simple diagonal)
            demagTensor = {
                DVector(1.0, 0.0, 0.0),
                DVector(0.0, 1.0, 0.0),
                DVector(0.0, 0.0, 1.0)
            };
        }

        double Ms, thickness, cellSurface, damping;
        std::vector<DVector> demagTensor;
    };

    // Test basic construction of the AFM layer
    TEST_F(AFMTest, Construction) {
        // Initial magnetizations - antiparallel
        DVector mag1(0.0, 0.0, 1.0);
        DVector mag2(0.0, 0.0, -1.0);
        DVector anis(0.0, 0.0, 1.0);

        // Negative J_AFM for antiferromagnetic coupling
        double J_AFM = -1e-20;  // Negative for AFM coupling

        LayerAFM<double> afmLayer("AFM1", mag1, mag2, anis, Ms, thickness,
                                cellSurface, demagTensor, damping, J_AFM);

        // Verify initialization
        EXPECT_EQ(afmLayer.getMagnetisation(0)[2], 1.0);
        EXPECT_EQ(afmLayer.getMagnetisation(1)[2], -1.0);
        EXPECT_DOUBLE_EQ(afmLayer.J_AFM->getCurrentScalarValue(0), J_AFM);
    }

    // Test the calculation of AFM exchange field
    TEST_F(AFMTest, ExchangeField) {
        DVector mag1(0.0, 0.0, 1.0);
        DVector mag2(0.0, 0.0, -1.0);
        DVector anis(0.0, 0.0, 1.0);

        double J_AFM = -1e-20;  // Strong AFM coupling

        LayerAFM<double> afmLayer("AFM1", mag1, mag2, anis, Ms, thickness,
                                cellSurface, demagTensor, damping, J_AFM);

        // Calculate exchange field from sublattice 2 to 1
        DVector Hafm = afmLayer.calculateAFMExchangeField(0, mag2);

        // Exchange field should point in the same direction as mag2, with negative J_AFM
        EXPECT_NEAR(Hafm[2], J_AFM / (Ms * thickness) * mag2[2], 1e-10);

        // The field should strengthen with stronger coupling
        double J_AFM_strong = -10e-20;
        LayerAFM<double> afmLayerStrong("AFM2", mag1, mag2, anis, Ms, thickness,
                                    cellSurface, demagTensor, damping, J_AFM_strong);

        DVector Hafm_strong = afmLayerStrong.calculateAFMExchangeField(0, mag2);
        EXPECT_NEAR(Hafm_strong[2], J_AFM_strong / (Ms * thickness) * mag2[2], 1e-10);
        EXPECT_TRUE(std::abs(Hafm_strong[2]) > std::abs(Hafm[2]));
    }

    // Test Néel vector and net magnetization calculation
    TEST_F(AFMTest, VectorCalculations) {
        DVector mag1(0.0, 0.0, 1.0);
        DVector mag2(0.0, 0.0, -1.0);
        DVector anis(0.0, 0.0, 1.0);

        double J_AFM = -1e-20;

        LayerAFM<double> afmLayer("AFM1", mag1, mag2, anis, Ms, thickness,
                                cellSurface, demagTensor, damping, J_AFM);

        // For perfect AFM, Néel vector should be 2.0 in the direction of mag1
        DVector neelVector = afmLayer.getNeelVector();
        EXPECT_NEAR(neelVector[2], 1.0, 1e-10);

        // Net magnetization should be near zero for perfect AFM
        DVector netMag = afmLayer.getNetMagnetisation();
        EXPECT_NEAR(netMag.length(), 0.0, 1e-10);

        // Now test with non-perfect AFM alignment
        DVector mag1Canted(0.866, 0.0, 0.5);    // 30 degrees from xy plane
        DVector mag2Canted(-0.866, 0.0, -0.5);  // 30 degrees from xy plane

        LayerAFM<double> afmLayerCanted("AFM2", mag1Canted, mag2Canted, anis, Ms, thickness,
                                    cellSurface, demagTensor, damping, J_AFM);

        // Check Néel vector is in the direction of mag1
        neelVector = afmLayerCanted.getNeelVector();
        EXPECT_NEAR(neelVector[0], 0.0, 1e-10);  // This should be close to zero
        EXPECT_NEAR(neelVector[1], 0.0, 1e-10);  // This should be close to zero
        EXPECT_NEAR(neelVector[2], 1.0, 1e-10);  // This should be close to zero

        // Net magnetization should still be near zero
        netMag = afmLayerCanted.getNetMagnetisation();
        EXPECT_NEAR(netMag.length(), 0.0, 1e-10);
    }

    // Test dynamics for a simple AFM layer under constant field
    TEST_F(AFMTest, BasicDynamics) {
        DVector mag1(0.0, 0.0, 1.0);
        DVector mag2(0.0, 0.0, -1.0);
        DVector anis(0.0, 0.0, 1.0);

        double J_AFM = -1e-21;  // Weak AFM coupling to see more movement

        LayerAFM<double> afmLayer("AFM1", mag1, mag2, anis, Ms, thickness,
                                cellSurface, demagTensor, damping, J_AFM);

        // Apply external field in x direction
        DVector externalField(1e5, 0, 0);  // 1000 Oe in x direction
        afmLayer.setExternalFieldDriver(AxialDriver<double>::getVectorAxialDriver(
            externalField[0], externalField[1], externalField[2]));

        // Step forward
        double time = 0;
        double timeStep = 1e-13;

        // Run a few steps to see dynamics
        for (int i = 0; i < 10; i++) {
            // The rk4_step method handles both sublattices
            afmLayer.rk4_step(time, timeStep, DVector(0,0,0), DVector(0,0,0));
            time += timeStep;
        }

        // Field in x should cause both sublattices to rotate in x-z plane
        // but maintain opposite z components
        EXPECT_GT(afmLayer.getMagnetisation(0)[0], 0);     // First sublattice should rotate +x
        EXPECT_GT(afmLayer.getMagnetisation(1)[0], 0);    // Second sublattice also +x
        EXPECT_GT(afmLayer.getMagnetisation(0)[2], 0);     // First keeps +z
        EXPECT_LT(afmLayer.getMagnetisation(1)[2], 0);    // Second keeps -z
    }

    // Test dynamics with changing external field
    TEST_F(AFMTest, ExternalFieldResponse) {
        DVector mag1(0.0, 0.0, 1.0);
        DVector mag2(0.0, 0.0, -1.0);
        DVector anis(0.0, 0.0, 0.0);  // No anisotropy for this test

        double J_AFM = -1e-21;  // Weak AFM coupling

        LayerAFM<double> afmLayer("AFM1", mag1, mag2, anis, Ms, thickness,
                                cellSurface, demagTensor, damping, J_AFM);

        // Run simulation with varying field
        double time = 0;
        double timeStep = 1e-13;
        double field_strength = 1e5;  // 1000 Oe

        // Store initial values
        DVector initial_m1 = afmLayer.getMagnetisation(0);
        DVector initial_m2 = afmLayer.getMagnetisation(1);

        // Apply field in x direction
        afmLayer.setExternalFieldDriver(AxialDriver<double>::getVectorAxialDriver(
            field_strength, 0, 0));

        // Run 100 steps to see significant motion
        for (int i = 0; i < 100; i++) {
            afmLayer.rk4_step(time, timeStep, DVector(0,0,0), DVector(0,0,0));
            time += timeStep;
        }

        // Store magnetizations after x-field
        DVector x_field_m1 = afmLayer.getMagnetisation(0);
        DVector x_field_m2 = afmLayer.getMagnetisation(1);

        // Verify x component increased for both sublattices
        EXPECT_GT(x_field_m1[0], initial_m1[0]);
        EXPECT_GT(x_field_m2[0], initial_m2[0]);

        // Now apply field in y direction
        afmLayer.setExternalFieldDriver(AxialDriver<double>::getVectorAxialDriver(
            0, field_strength, 0));

        // Run more steps
        for (int i = 0; i < 100; i++) {
            afmLayer.rk4_step(time, timeStep, DVector(0,0,0), DVector(0,0,0));
            time += timeStep;
        }

        // Verify y component increased for both sublattices
        EXPECT_GT(afmLayer.getMagnetisation(0)[1], x_field_m1[1]);
        EXPECT_GT(afmLayer.getMagnetisation(1)[1], x_field_m2[1]);
    }

    // Test the AFM relaxation process
    TEST_F(AFMTest, Relaxation) {
        // Start with slightly misaligned magnetizations
        DVector mag1(0.1, 0.0, 0.995);  // Slightly tilted from z
        DVector mag2(0.1, 0.0, -0.995); // Slightly tilted from -z
        mag1.normalize();
        mag2.normalize();

        DVector anis(0.0, 0.0, 1.0);  // z-axis anisotropy

        double J_AFM = -5e-21;  // Strong AFM coupling to see relaxation

        LayerAFM<double> afmLayer("AFM1", mag1, mag2, anis, Ms, thickness,
                                cellSurface, demagTensor, 0.5, J_AFM);  // Higher damping

        // No external fields, just let the system relax
        double time = 0;
        double timeStep = 1e-13;

        // Run for a few picoseconds to see relaxation
        for (int i = 0; i < 1000; i++) {
            afmLayer.rk4_step(time, timeStep, DVector(0,0,0), DVector(0,0,0));
            time += timeStep;
        }

        // Magnetizations should align closer to the z-axis due to anisotropy
        // and be more antiparallel due to AFM exchange
        EXPECT_NEAR(afmLayer.getMagnetisation(0)[0], 0.0, 0.05);  // x should go to zero
        EXPECT_NEAR(afmLayer.getMagnetisation(1)[0], 0.0, 0.05); // x should go to zero

        // z components should be closer to +/-1
        EXPECT_GT(afmLayer.getMagnetisation(0)[2], mag1[2]);
        EXPECT_LT(afmLayer.getMagnetisation(1)[2], mag2[2]);

        // Verify they remain antiparallel
        EXPECT_NEAR(afmLayer.getMagnetisation(0)[2], -afmLayer.getMagnetisation(1)[2], 0.01);
    }

    int main(int argc, char **argv) {
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }
