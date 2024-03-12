#include "junction.hpp"
// import map
#include <map>

template <typename T = double>
class LLGBLayer
{
public:
    ScalarDriver<T> Ms_driver;
    ScalarDriver<T> temperatureDriver;
    AxialDriver<T> externalFieldDriver;
    T Tc;
    T susceptibility;
    T susceptibility_perp;
    T me;
    CVector<T> mag;
    T thickness;
    T damping;
    T alpha_perp_log, alpha_par_log;
    std::string id;
    LLGBLayer(
        const std::string& id,
        CVector<T> mag,
        T Ms,
        T Tc,
        T susceptibility,
        T susceptibility_perp,
        T me,
        T thickness,
        T damping) : id(id), mag(mag), Tc(Tc),
        susceptibility(susceptibility), susceptibility_perp(susceptibility_perp), me(me), thickness(thickness), damping(damping) {
        this->Ms_driver = ScalarDriver<T>::getConstantDriver(Ms);
    }

    T getAlphaParallel(T& time) {
        const T temp = this->temperatureDriver.getCurrentScalarValue(time);
        this->alpha_par_log = this->damping * (temp / this->Tc) * (2. / 3.);
        return this->alpha_par_log;
    }

    T getAlphaPerpendicular(T& time) {
        const T temp = this->temperatureDriver.getCurrentScalarValue(time);
        const T ratio = temp / this->Tc;
        if (temp >= this->Tc) {
            return this->damping * ratio * (2. / 3.);
        }
        this->alpha_perp_log = this->damping * (1 - ratio / 3.0);
        return this->alpha_perp_log;

    }

    CVector<T> getLongitudinal(T time, const CVector<T>& m) {
        const T temp = this->temperatureDriver.getCurrentScalarValue(time);
        const T ratio_susc = 1. / (2 * this->susceptibility);
        const T m2 = pow(m.length(), 2);
        if (temp >= this->Tc) {
            const T ratio_m = m2 / pow(this->me, 2);
            return ratio_susc * (1 - ratio_m) * m;
        }
        const T ratio_T = (this->Tc / (temp - this->Tc));
        const T ratio_T_adj = (3. / 5.) * ratio_T * m2 - 1.;
        return ratio_susc * ratio_T_adj * m;
    }

    CVector<T> getAnisotropyField(T time, const CVector<T>& m) {
        return (-1 / this->susceptibility_perp) * CVector<T>(m[0], m[1], 0);
    }


    const CVector<T> calculateHeff(T time, const CVector<T>& m) {
        const CVector<T> anis = this->getAnisotropyField(time, m);
        const CVector<T> hext = this->externalFieldDriver.getCurrentAxialDrivers(time);
        const CVector<T> temp = this->getLongitudinal(time, m);
        return anis + hext + temp;
    }

    CVector<T> calculateLLG(const T& time, const T& timeStep, const CVector<T>& m) {
        const CVector<T> heff = this->calculateHeff(time, m);
        return solveLLG(time, timeStep, m, heff);
    }

    const CVector<T> solveLLG(T time, T timeStep, const CVector<T>& m, const CVector<T>& heff) {
        const CVector<T> mxh = c_cross<T>(m, heff);
        const CVector<T> mxmxht = c_cross<T>(m, mxh);
        const CVector<T> llbTerm = c_dot(m, heff) * m;
        const T inv_mlen = pow(1. / m.length(), 2);
        const T gamma_p = GYRO / (1 + pow(this->damping, 2)); // LLGS -> LL form
        const CVector<T> dmdt = -1 * mxh - getAlphaPerpendicular(time) * mxmxht * inv_mlen + llbTerm * getAlphaParallel(time) * inv_mlen;
        return GYRO * dmdt;
    }


    // TODO: implement the stochastic part
    CVector<T> getOneFVector() {
        return CVector<T>();
    }
    CVector<T> getStochasticLangevinVector(T time, T timeStep) {
        return CVector<T>();
    }
    CVector<T> stochasticTorque(const CVector<T>& m, const CVector<T>& dW) {
        return CVector<T>();
    }
};

template <typename T = double>
class LLGBJunction
{
    friend class LLGBLayer<T>;
    const std::vector<std::string> vectorNames = { "x", "y", "z" };
    std::vector<LLGBLayer<T>> layers;
    std::unordered_map<std::string, std::vector<T>> log;
    unsigned int logLength = 0;
    unsigned int layerNo = 0;
public:
    explicit LLGBJunction(const std::vector<LLGBLayer<T>>& layers) {
        this->layers = layers;
        this->layerNo = layers.size();
    }

    void heunSolverStep(const T& t, const T& timeStep) {
        /*
            Heun method
            y'(t+1) = y(t) + dy(y, t)
            y(t+1) = y(t) + 0.5 * (dy(y, t) + dy(y'(t+1), t+1))
        */
        /*
            Stochastic Heun method
            y_np = y + g(y,t,dW)*dt
            g_sp = g(y_np,t+1,dW)
            y' = y_n + f_n * dt + g_n * dt
            f' = f(y, )
            y(t+1) = y + dt*f(y,t) + .5*(g(y,t,dW)+g_sp)*sqrt(dt)
        */
        std::vector<CVector<T>> fn(this->layerNo, CVector<T>());
        std::vector<CVector<T>> gn(this->layerNo, CVector<T>());
        std::vector<CVector<T>> dW(this->layerNo, CVector<T>());
        std::vector<CVector<T>> mNext(this->layerNo, CVector<T>());
        // first approximation

        // make sure that
        // 1. Thermal field is added if needed
        // 2. One/f noise is added if needed
        // 3. The timestep is correctly multiplied

        for (unsigned int i = 0; i < this->layerNo; i++)
        {
            fn[i] = this->layers[i].calculateLLG(
                t, timeStep, this->layers[i].mag);

            // draw the noise for each layer, dW
            dW[i] = this->layers[i].getStochasticLangevinVector(t, timeStep) + this->layers[i].getOneFVector();
            gn[i] = this->layers[i].stochasticTorque(this->layers[i].mag, dW[i]);

            mNext[i] = this->layers[i].mag + fn[i] * timeStep + gn[i] * sqrt(timeStep);
        }
        // second approximation
        for (unsigned int i = 0; i < this->layerNo; i++)
        {
            // first approximation is already multiplied by timeStep
            this->layers[i].mag = this->layers[i].mag + 0.5 * timeStep * (
                fn[i] + this->layers[i].calculateLLG(
                    t + timeStep, timeStep, mNext[i])
                ) + 0.5 * (gn[i] + this->layers[i].stochasticTorque(mNext[i], dW[i])) * sqrt(timeStep);
            // normalise only in classical
            // this->layers[i].mag.normalize();
        }
    }
    void eulerHeunSolverStep(const T& t, const T& timeStep) {
        /*
            Euler Heun method (stochastic heun)

            y_np = y + g(y,t,dW)*dt
            g_sp = g(y_np,t+1,dW)
            y(t+1) = y + dt*f(y,t) + .5*(g(y,t,dW)+g_sp)*sqrt(dt)

            with f being the non-stochastic part and g the stochastic part
        */
        // draw the noise for each layer, dW
        std::vector<CVector<T>> mPrime(this->layerNo, CVector<T>());
        for (unsigned int i = 0; i < this->layerNo; i++) {
            // todo: after you're done, double check the thermal magnitude and dt scaling there
            const CVector<T> dW = this->layers[i].getStochasticLangevinVector(t, timeStep) + this->layers[i].getOneFVector();

            const CVector<T> fnApprox = this->layers[i].calculateLLG(
                t, timeStep, this->layers[i].mag);
            const CVector<T> gnApprox = this->layers[i].stochasticTorque(this->layers[i].mag, dW);

            // theoretically we have 2 options
            // 1. calculate only the stochastic part with the second approximation
            // 2. calculate the second approximation of m with the stochastic and non-stochastic
            //    part and then use if for torque est.
            const CVector<T> mNext = this->layers[i].mag + gnApprox * sqrt(timeStep);
            const CVector<T> gnPrimeApprox = this->layers[i].stochasticTorque(mNext, dW);
            mPrime[i] = this->layers[i].mag + fnApprox * timeStep + 0.5 * (gnApprox + gnPrimeApprox) * sqrt(timeStep);
        }

        for (unsigned int i = 0; i < this->layerNo; i++) {
            this->layers[i].mag = mPrime[i];
            // this->layers[i].mag.normalize(); LLGB don't need to be normalised
        }
    }

    typedef void (LLGBJunction<T>::* runnerFn)(const T& t, const T& timeStep);
    std::tuple<runnerFn, SolverMode> getSolver(SolverMode mode) {
        auto runner = &LLGBJunction<T>::heunSolverStep;
        if (mode == HEUN)
            runner = &LLGBJunction<T>::heunSolverStep;
        else if (mode == EULER_HEUN)
            runner = &LLGBJunction<T>::eulerHeunSolverStep;

        return std::make_tuple(runner, mode);
    }

    /**
     * @brief Log computed layer parameters.
     * This function logs all the necessayr parameters of the layers.
     * @param t: current time
     * @param timeStep: timeStep of the simulation (unsued for now)
     * @param calculateEnergies: if true, also include fields for energy computation.
     */
    void logLayerParams(T& t, const T timeStep)
    {
        for (const auto& layer : this->layers)
        {
            const std::string lId = layer.id;
            // always save magnetisation
            for (int i = 0; i < 3; i++)
            {
                this->log[lId + "_m" + vectorNames[i]].emplace_back(layer.mag[i]);
            }
            this->log["alpha_parallel"].emplace_back(layer.alpha_par_log);
            this->log["alpha_perpendicular"].emplace_back(layer.alpha_perp_log);
        }
        this->log["time"].emplace_back(t);
        this->logLength++;
    }


    void
        saveLogs(std::string filename)
    {
        if (filename == "")
        {
            // if there's an empty fn, don't save
            throw std::runtime_error("The filename may not be empty!");
        }
        std::ofstream logFile;
        logFile.open(filename);
        for (const auto& keyPair : this->log)
        {
            logFile << keyPair.first << ";";
        }
        logFile << "\n";
        for (unsigned int i = 0; i < logLength; i++)
        {
            for (const auto& keyPair : this->log)
            {
                logFile << keyPair.second[i] << ";";
            }
            logFile << "\n";
        }
        logFile.close();
    }

    /**
     * Clears the simulation log.
     **/
    void clearLog()
    {
        this->log.clear();
        this->logLength = 0;
    }

    /**
     * Main run simulation function. Use it to run the simulation.
     * @param totalTime: total time of a simulation, give it in seconds. Typical length is in ~couple ns.
     * @param timeStep: the integration step of the RK45 method. Default is 1e-13
     * @param writeFrequency: how often is the log saved to? Must be no smaller than `timeStep`. Default is 1e-11.
     * @param persist: whether to save to the filename specified in the Junction constructor. Default is true
     * @param log: if you want some verbosity like timing the simulation. Default is false
     * @param mode: Solver mode EULER_HEUN, RK4 or DORMAND_PRICE
     */
    void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11,
        bool log = false,
        SolverMode mode = RK4)

    {
        if (timeStep > writeFrequency)
        {
            std::runtime_error("The time step cannot be larger than write frequency!");
        }
        const unsigned int totalIterations = (int)(totalTime / timeStep);
        const unsigned int writeEvery = (int)(writeFrequency / timeStep);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // pick a solver based on drivers
        auto [runner, _] = getSolver(mode);

        for (unsigned int i = 0; i < totalIterations; i++)
        {
            T t = i * timeStep;
            (*this.*runner)(t, timeStep);

            if (!(i % writeEvery))
            {
                logLayerParams(t, timeStep);
            }
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if (log)
        {
            std::cout << "Steps in simulation: " << totalIterations << std::endl;
            std::cout << "Write every: " << writeEvery << std::endl;
            std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
        }
    }
};
