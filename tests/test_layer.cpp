#include "../core/junction.hpp"
#include <gtest/gtest.h>
#include <map>
class JunctionTest : public testing::Test
{
protected:
    void SetUp() override
    {

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
        mtj = new Junction<double>(
            {l1});
    }
    Junction<double> *mtj;
};

TEST_F(JunctionTest, IsLogEmpty)
{
    // check if the log is empty before the sim starts
    EXPECT_EQ(mtj->getLog().size(), 0);
}

TEST_F(JunctionTest, BasicSimulation)
{
    // check the steps in the simulation
    // check if the simulated mag is not zero.
    const double time = 1e-9;
    const double sTime = 1e-11;
    mtj->runSimulation(time, sTime, sTime);
    auto log = mtj->getLog();
    EXPECT_NO_THROW(log.at("time")) << "Time is not available in logs!";
    const unsigned int logSize = log["time"].size();
    const unsigned int expectedLogSize = time / sTime;
    ASSERT_EQ(expectedLogSize, logSize) << "Log size does not match the expected one!";
    ASSERT_NE(mtj->getLayerMagnetisation("free"), CVector<double>(0, 0, 0)) << "Magnetisation is zero!";
    // TODO:

    // try checking if the mag is set properly
    // check exclusive conditions on
    // -- temperature driver & solver selection EXPECT THROW
    // -- SOT driver vs layer definition EXPECT THROW
    // test ScalarDriver
};

TEST_F(JunctionTest, SetTest)
{
    // check if setting mags is correct
    ASSERT_EQ(mtj->getLayerMagnetisation("free"), CVector<double>(0, 0, 1));
    CVector<double> zeroVec(0, 0, 0);
    // this should throw an error, we cannot set a 0 vector!
    ASSERT_THROW(mtj->setLayerMagnetisation("free", zeroVec), std::runtime_error) << "Can set 0 mag!";

    CVector<double> magVec(1, 0, 0);
    mtj->setLayerMagnetisation("free", magVec);
    const auto vec = mtj->getLayerMagnetisation("free");
    ASSERT_EQ(vec, std::as_const(magVec)) << "Setting mag is ineffective!: " << vec.x << ":" << vec.y << ":" << vec.z;
}

TEST_F(JunctionTest, TemperatureSolver)
{
    mtj->setLayerTemperatureDriver(
        "free", ScalarDriver<double>::getConstantDriver(100));
    ASSERT_TRUE(mtj->layers[0].hasTemperature()); // this enforces the solver
}
