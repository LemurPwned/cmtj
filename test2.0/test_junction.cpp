#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "../core/2.0/junction.hpp"
#include "../core/2.0/fm.hpp"
#include "../core/2.0/driver.hpp"

// Helper function to check if two values are approximately equal
template <typename T>
bool approxEqual(T a, T b, T epsilon = 1e-5) {
    return std::abs(a - b) < epsilon;
}

class JunctionTest : public ::testing::Test {
protected:
    DVector mag1;
    DVector mag2;
    DVector anis;
    std::vector<DVector> demagTensor;
    std::shared_ptr<Layer<double>> layer1;
    std::shared_ptr<Layer<double>> layer2;

    void SetUp() override {
        mag1 = DVector(0, 0, 1);
        mag2 = DVector(1, 0, 0);
        anis = DVector(0, 0, 1);
        demagTensor = {
            DVector(1, 0, 0),
            DVector(0, 1, 0),
            DVector(0, 0, 1)
        };

        layer1 = std::make_shared<Layer<double>>("layer1", mag1, anis, 1400e3, 2e-9, 1e-14, demagTensor, 0.01);
        layer2 = std::make_shared<Layer<double>>("layer2", mag2, anis, 1400e3, 2e-9, 1e-14, demagTensor, 0.01);
    }
};

TEST_F(JunctionTest, JunctionCreation) {
    // Create junction
    std::vector<std::shared_ptr<AbstractLayer<double>>> layers = { layer1, layer2 };
    FMJunction<double> junction(layers);

    // Verify layer count
    EXPECT_EQ(junction.getLayerCount(), 2);

    // Verify layer IDs
    auto ids = junction.getLayerIds();
    EXPECT_EQ(ids.size(), 2);
    EXPECT_EQ(ids[0], "layer1");
    EXPECT_EQ(ids[1], "layer2");

    // Verify we can get layer by index
    auto retrievedLayer = junction.getLayer(0);
    EXPECT_EQ(retrievedLayer->getId(), "layer1");
}

TEST_F(JunctionTest, ClassicBilayerMagnetoresistance) {
    // Create junction
    std::vector<std::shared_ptr<AbstractLayer<double>>> layers = { layer1, layer2 };
    FMJunction<double> junction(layers);

    // Set MR mode and parameters
    junction.setMRMode(MRmode::CLASSIC);
    junction.setRp(100.0);
    junction.setRap(200.0);

    // Calculate MR - for perpendicular magnetizations, dot product is 0
    layer1->setMagnetisation(CVector<double>(1, 0, 0));
    layer2->setMagnetisation(CVector<double>(0, 0, 1));
    auto mr = junction.getMagnetoresistance();

    // Expected: Rp + ((Rap-Rp)/2) * (1-0) = 100 + 50 = 150
    EXPECT_NEAR(mr[0], 150.0, 1e-5);

    // Change to parallel alignment
    layer1->setMagnetisation(CVector<double>(1, 0, 0));
    layer2->setMagnetisation(CVector<double>(1, 0, 0));
    mr = junction.getMagnetoresistance();

    // Expected: Rp + ((Rap-Rp)/2) * (1-1) = 100 + 0 = 100
    EXPECT_NEAR(mr[0], 100.0, 1e-5);
}

TEST_F(JunctionTest, PinningLayerMagnetoresistance) {
    // Create test layer with required parameters
    DVector mag(0, 0, 1);
    auto layer = std::make_shared<Layer<double>>("layer", mag, anis, 1400e3, 2e-9, 1e-14, demagTensor, 0.01);

    // Set reference layer before creating junction
    layer->setReferenceLayer(CVector<double>(1, 0, 0));

    // Create junction
    std::vector<std::shared_ptr<AbstractLayer<double>>> layers = { layer };
    FMJunction<double> junction(layers);

    // Set MR mode and parameters
    junction.setMRMode(MRmode::CLASSIC);
    junction.setRp(100.0);
    junction.setRap(200.0);

    // Calculate MR - perpendicular, dot product is 0
    auto mr = junction.getMagnetoresistance();

    // Expected: Rp + ((Rap-Rp)/2) * (1-0) = 100 + 50 = 150
    EXPECT_TRUE(approxEqual(mr[0], 150.0));
}

TEST_F(JunctionTest, SettingDrivers) {
    // Create test layer with required parameters
    DVector mag(0, 0, 1);
    auto layer = std::make_shared<Layer<double>>("testLayer", mag, anis, 1400e3, 2e-9, 1e-14, demagTensor, 0.01);

    // Create junction
    std::vector<std::shared_ptr<AbstractLayer<double>>> layers = { layer };
    FMJunction<double> junction(layers);

    // Test setting temperature driver
    auto tempDriver = std::make_shared<NullDriver<double>>(300.0);
    EXPECT_NO_THROW(junction.setLayerTemperatureDriver("testLayer", tempDriver));

    // Test setting anisotropy driver
    auto anisDriver = std::make_shared<NullDriver<double>>(1e6);
    EXPECT_NO_THROW(junction.setLayerAnisotropyDriver("testLayer", anisDriver));

    // Test setting external field driver
    auto fieldDriver = std::make_shared<AxialDriver<double>>(CVector<double>(0, 0, 1e3));
    EXPECT_NO_THROW(junction.setLayerExternalFieldDriver("testLayer", fieldDriver));
}

TEST_F(JunctionTest, InvalidLayerId) {
    // Create test layer with required parameters
    DVector mag(0, 0, 1);
    auto layer = std::make_shared<Layer<double>>("validLayer", mag, anis, 1400e3, 2e-9, 1e-14, demagTensor, 0.01);

    // Create junction
    std::vector<std::shared_ptr<AbstractLayer<double>>> layers = { layer };
    FMJunction<double> junction(layers);

    // Test setting driver with invalid layer ID
    auto tempDriver = std::make_shared<NullDriver<double>>(300.0);

    EXPECT_THROW(junction.setLayerTemperatureDriver("invalidLayer", tempDriver), std::runtime_error);
}

TEST_F(JunctionTest, ShortSimulation) {
    // Create test layer with required parameters
    DVector mag(0, 0, 1);
    auto layer = std::make_shared<Layer<double>>("testLayer", mag, anis, 1400e3, 2e-9, 1e-14, demagTensor, 0.01);

    // Create junction
    std::vector<std::shared_ptr<AbstractLayer<double>>> layers = { layer };
    FMJunction<double> junction(layers);

    // Run simulation
    junction.runSimulation(1e-9, 1e-12, 1e-12);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
