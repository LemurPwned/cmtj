#include <gtest/gtest.h>
#include "../core/reservoir.hpp"
#include "../core/drivers.hpp"
#include <Eigen/Dense>

auto fitLeastSquares(Eigen::MatrixXd X, Eigen::MatrixXd Y)
{
    // using the ordinary least squares formulation
    // Eigen::Map<Eigen::VectorXd> X(xdata.data(), xdata.size());
    // Eigen::Map<Eigen::VectorXd> Y(ydata.data(), ydata.size());
    const auto result = (X.transpose() * X).inverse() * X.transpose() * Y;
    return result;
}

double sigmoidExpr(double x)
{
    return 1. / (1 + std::exp(-x));
}
double power2(double x)
{
    return pow(x, 2);
}

auto generateRandomBinaryInput(int size)
{
    // https://stackoverflow.com/questions/38333326/eigen-random-binary-vector-with-t-1s
    std::vector<int> esv(size, 0);
    std::fill_n(esv.begin(), (int)(size / 2), 1);
    Eigen::Map<Eigen::VectorXi> es(esv.data(), esv.size());
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(std::begin(esv), std::end(esv), g);
    return es;
}

auto generateTrainSignal(auto es, int delay, int size)
{
    // AND
    std::vector<int> trainY(es.size(), 0);
    for (int i = delay; i < es.size(); i++)
    {
        trainY[i] = es[i] * es[i - delay];
    }
    return trainY;
}

auto logisitcRegression(Eigen::MatrixXd X, Eigen::VectorXd Y, int k, int reservoirStates)
{
    const double eps = 1e-3;
    const double learningRate = 1e-3;
    // zero vector of parameters
    Eigen::VectorXd theta = Eigen::VectorXd::Constant(reservoirStates, 0);
    Eigen::VectorXd thetaNew = Eigen::VectorXd::Constant(reservoirStates, 0);
    const unsigned int maxSteps = 10000;
    unsigned int steps = 0;
    while (true)
    {
        theta = thetaNew;
        // theta.transpose() * X;
        thetaNew = theta + learningRate * X.transpose() * (Y - (X * theta).unaryExpr(&sigmoidExpr));
        // thetaNew = theta + learningRate * (Y - (theta.transpose() * X).unaryExpr(&sigmoidExpr)) * X.transpose();
        steps += 1;
        if ((steps > maxSteps) || (abs((theta - thetaNew).sum()) <= eps))
        {
            break;
        }
    }
    // try test error
    const auto pred = (X * theta).unaryExpr(&sigmoidExpr);
    const auto MSE = (Y - pred).unaryExpr(&power2).sum();
    std::cout << "MSE: " << MSE << std::endl;
    std::cout << "Optimised after: " << steps << " steps" << std::endl;
    return theta;
}

TEST(RESERVOIR_TEST, BasicReservoirSim)
{
    const double time = 4e-9;
    const double sTime = 4e-12;
    // we want to create a new junction for this
    std::vector<CVector<double>> demagTensor = {
        {0.0405132594942993, -2.11556045217049e-13, -3.54967543157989e-14},
        {-2.11556045217049e-13, 0.0405132594942993, -3.55034146894455e-14},
        {-3.54967543157989e-14, -3.55034146894455e-14, 0.998973481008816} };
    const double Ms = 1.63;
    const double gap = 60e-9;
    const double thickness = 1e-9;
    const double radius = 20e-9;
    const double surface = M_PI * pow(radius, 2);
    const double volume = surface * radius;
    const double Ku0 = 1e4;
    const double Hmax = 10e4;

    const int delay = 1;
    std::vector<double> trainingSignal = {
        0, 0, 1, 1, 0, 0, 1, 1, 1 };

    std::vector<double> ysignalAndN1 = {
        0, 0, 1, 0, 0, 0, 1, 1, 0 };
    ASSERT_EQ(trainingSignal.size(), ysignalAndN1.size());

    const std::vector<std::vector<DVector>>
        coordinateMatrix = {
                            {
                                DVector(0, 5 * gap, 0),
                                DVector(gap, 5 * gap, 0.),
                            },

                            {
                                DVector(0, 4 * gap, 0),
                                DVector(gap, 4 * gap, 0.),
                            },
                            {
                                DVector(0, 3 * gap, 0),
                                DVector(gap, 3 * gap, 0.),
                            },

                            {
                                DVector(0, 2 * gap, 0),
                                DVector(gap, 2 * gap, 0.),
                            },

                             {
                                DVector(0, gap, 0),
                                DVector(gap, gap, 0.),
                            },
                            {
                                DVector(0, 0, 0),
                                DVector(gap, 0, 0.),
                            } };
    // print coordinate matrix
    std::cout << "Coordinate matrix" << std::endl;
    for (const auto& row : coordinateMatrix)
    {
        for (auto& col : row)
        {
            std::cout << col << " ";
        }
        std::cout << std::endl;
    }
    Layer<double> l1("free",                      // id
        CVector<double>(0., 0., 1.), // mag
        CVector<double>(0., 0., 1.), // anis
        Ms,                          // Ms
        thickness,                   // thickness
        surface,                     // surface
        demagTensor,                 // demag
        0.5                          // damping
    );
    auto anisDriver = ScalarDriver<double>::getConstantDriver(
        Ku0);
    l1.setAnisotropyDriver(anisDriver);
    const std::vector<std::vector<Layer<double>>> layerMatrix = {
        // {l1, l1}, {l1, l1} };
     {l1, l1} , {l1, l1} , {l1, l1}, {l1, l1} };

    auto reservoir = Reservoir(coordinateMatrix, layerMatrix);
    const int k = trainingSignal.size();
    Eigen::MatrixXd mat(k, reservoir.noElements);
    std::cout << "Training over " << k << " steps" << std::endl;
    for (int step = 0; step < k; step++)
    {
        std::cout << "### STEP " << step << " ###" << std::endl;
        // STAGE I
        // pass the signal to the first group
        const double polariser = trainingSignal[step] == 1 ? 1 : -1;
        reservoir.setLayerExternalField(0,
            AxialDriver<double>(
                CVector<double>(0, 0, polariser * Hmax)));
        // STAGE 2
        // change the anisotropy to 0 for Group II and III
        std::cout << "STAGE 2" << std::endl;
        for (int i = 0; i < reservoir.noElements; i++)
        {
            if ((i % 6) < 2)
            {
                std::cout << "Setting anisotropy of layer: " << i << " to 0" << std::endl;
                // set it to 0
                reservoir.setLayerAnisotropy(i, ScalarDriver<double>::getConstantDriver(0));
            }
        }
        reservoir.runSimulation(time, sTime);
        // this should align the data from I to II and III

        // STAGE 3
        // The magnetisations are going to be fixed by changing their Ku back to their Ku0
        std::cout << "STAGE 3 " << std::endl;
        for (int i = 0; i < reservoir.noElements; i++)
        {
            if ((i % 6) < 2)
            {
                std::cout << "Setting anisotropy of layer: " << i << " to " << Ku0 << std::endl;
                // set it to 0
                reservoir.setLayerAnisotropy(i, ScalarDriver<double>::getConstantDriver(Ku0));
            }
        }
        reservoir.runSimulation(time, sTime);

        // Collection happens at stage 3
        std::cout << "Collecting states" << std::endl;
        std::vector<DVector> states = reservoir.collectFrozenMMatrix();
        std::vector<double> stateResult;
        std::cout << " Mag X:" << states[0].x << std::endl;
        // state transform
        std::transform(
            states.begin(), states.end(),
            std::back_inserter(stateResult), [](CVector<double> v) -> double
            { return 0.5 * v.x + 0.5; }); // adjust the mean of the output
        std::cout << states.size() << " " << mat.size() << ":" << mat.rows() << "," << mat.cols() << std::endl;
        for (int z = 0; z < states.size(); z++)
        {
            mat(step, z) = stateResult[z];
        }

        std::cout << "STAGE 4" << std::endl;
        // Group I has zeroed out anisotropy
        for (int i = 0; i < 2; i++)
        {
            std::cout << "Setting anisotropy of layer: " << i << " to 0" << std::endl;
            // set it to 0
            reservoir.setLayerAnisotropy(i, ScalarDriver<double>::getConstantDriver(0));
        }
        reservoir.runSimulation(time, sTime);
    }
    Eigen::Map<Eigen::VectorXd> Y(ysignalAndN1.data(), ysignalAndN1.size());
    std::cout << "X:\n"
        << mat << std::endl;
    std::cout << "Y: " << Y.transpose() << std::endl;
    std::cout << "Y: [" << Y.rows() << " " << Y.cols() << "]" << std::endl;
    Eigen::VectorXd weights = logisitcRegression(
        mat, Y,
        k, reservoir.noElements);
    std::cout << weights << std::endl;
    reservoir.saveLogs("reservoir.csv");
}
