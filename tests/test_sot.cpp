#include <gtest/gtest.h>
#include "../core/junction.hpp"

TEST(SOT_TEST, BasicSOT_SOTLayer)
{
    const double time = 1e-9;
    const double sTime = 1e-11;
    // we want to create a new junction for this
    std::vector<CVector<double>> demagTensor = {
        {0.0405132594942993, -2.11556045217049e-13, -3.54967543157989e-14},
        {-2.11556045217049e-13, 0.0405132594942993, -3.55034146894455e-14},
        {-3.54967543157989e-14, -3.55034146894455e-14, 0.918973481008816}};

    Layer<double> l1 = Layer<double>::LayerSOT("free",                      // id
                                               CVector<double>(0., 0., 1.), // mag
                                               CVector<double>(0., 0., 1.), // anis
                                               0.5,                         // Ms
                                               1e-9,                        // thickness
                                               70e-9,                       // surface
                                               demagTensor,                 // demag
                                               0.03,                        // damping
                                               150e-8,
                                               150e-9);
    Junction<double> mtj({l1});

    ASSERT_THROW(mtj.setLayerFieldLikeTorqueDriver(
                     "free", ScalarDriver<double>::getSineDriver(0, 1, 1e9, 0)),
                 std::runtime_error);
    mtj.runSimulation(time, sTime, sTime);
}

TEST(SOT_TEST, BasicSOT_DynamicSOT)
{
    const double time = 1e-9;
    const double sTime = 1e-11;
    // we want to create a new junction for this
    std::vector<CVector<double>> demagTensor = {
        {0.0405132594942993, -2.11556045217049e-13, -3.54967543157989e-14},
        {-2.11556045217049e-13, 0.0405132594942993, -3.55034146894455e-14},
        {-3.54967543157989e-14, -3.55034146894455e-14, 0.918973481008816}};

    Layer<double> l1("free",                      // id
                     CVector<double>(0., 0., 1.), // mag
                     CVector<double>(0., 0., 1.), // anis
                     0.5,                         // Ms
                     1e-9,                        // thickness
                     70e-9,                       // surface
                     demagTensor,                 // demag
                     0.03                         // damping
    );
    Junction<double> mtj({l1});

    ASSERT_NO_THROW(mtj.setLayerFieldLikeTorqueDriver(
        "free", ScalarDriver<double>::getSineDriver(0, 1, 1e9, 0)));
    mtj.runSimulation(time, sTime, sTime);
}
