#include <gtest/gtest.h>
#include "../core/llgb.hpp"

TEST(LLGBTest, Basic)
{
    const double Ms = 1.2;
    const double Tc = 400;
    const double susceptibility = 1.1;
    const double susceptibility_perp = 100;
    LLGBLayer<double> l1(
        "free",
        CVector<double>(1, 0, 0),
        1.2,
        1.1,
        100,
        0.1,
        0.1,
        0.9,
        1e-9,
        0.01
    );
    LLGBJunction<double> junction(std::vector < LLGBLayer<double>>{l1});
    junction.setLayerTemperatureDriver("all", ScalarDriver<double>::getConstantDriver(300));
    junction.runSimulation(10e-9, 1e-9, 1e-9);
    junction.saveLogs("test_llgb_basic.csv");
}