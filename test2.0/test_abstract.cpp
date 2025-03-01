#include <gtest/gtest.h>
#include "../core/2.0/llgb.hpp"


TEST(LLGBLayerTest, Construction) {
    // Test basic construction of LLGB layer
    CVector<double> mag(0, 0, 1);
    CVector<double> anis(0, 0, 1);
    std::vector<CVector<double>> demagTensor = {
        CVector<double>(0.1, 0, 0),
        CVector<double>(0, 0.1, 0),
        CVector<double>(0, 0, 0.8)
    };

    EXPECT_NO_THROW(LLGBLayer<double>("test", mag, anis, 1e6, 1e-9, 1e-12,
        demagTensor, 0.1, 800, 1.0, 0.8));

    // Test invalid construction
    EXPECT_THROW(LLGBLayer<double>("test", CVector<double>(0, 0, 0), anis, 1e6, 1e-9, 1e-12,
        demagTensor, 0.1, 800, 1.0, 0.8),
        std::runtime_error);

    EXPECT_THROW(LLGBLayer<double>("test", mag, CVector<double>(0, 0, 0), 1e6, 1e-9, 1e-12,
        demagTensor, 0.1, 800, 1.0, 0.8),
        std::runtime_error);

    EXPECT_THROW(LLGBLayer<double>("test", mag, anis, 1e6, 0, 1e-12,
        demagTensor, 0.1, 800, 1.0, 0.8),
        std::runtime_error);
}

TEST(LLGBLayerTest, DampingParameters) {
    CVector<double> mag(0, 0, 1);
    CVector<double> anis(0, 0, 1);
    std::vector<CVector<double>> demagTensor = {
        CVector<double>(0.1, 0, 0),
        CVector<double>(0, 0.1, 0),
        CVector<double>(0, 0, 0.8)
    };
    const double Tc = 800;
    double Tconst = 300;
    double damping = 0.1;
    LLGBLayer<double> layer("test", mag, anis, 1e6, 1e-9, 1e-12,
        demagTensor, damping, Tc, 1.0, 0.8);

    // Set temperature driver
    layer.setTemperatureDriver(ScalarDriver<double>::getConstantDriver(Tconst));

    double time = 0;

    // Test alpha parallel and perpendicular at T < Tc
    EXPECT_NEAR(layer.getAlphaParallel(time), damping * (Tconst / Tc) * 2 / 3, 1e-10);
    EXPECT_NEAR(layer.getAlphaPerpendicular(time), damping * (1 - (Tconst / Tc) * 1 / 3), 1e-10);

    // Test at T > Tc
    Tconst = 1000;
    layer.setTemperatureDriver(ScalarDriver<double>::getConstantDriver(Tconst));

    EXPECT_NEAR(layer.getAlphaParallel(time), damping * (Tconst / Tc) * 2 / 3, 1e-10);
    EXPECT_NEAR(layer.getAlphaPerpendicular(time), damping * (Tconst / Tc) * 2 / 3, 1e-10);
}
