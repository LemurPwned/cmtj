/**
 * @file junction.hpp
 * @author Jakub
 * @brief CMTJ Junction model
 * @version 1.0
 * @date 2022-03-22
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef CORE_JUNCTION_HPP_
#define CORE_JUNCTION_HPP_

#define _USE_MATH_DEFINES
#include <fstream>                  // for file save
#include <algorithm>                // for find_if
#include <functional>               // for bind, function
#include <array>                    // for array, array<>::value_type
#include <chrono>                   // for seconds, steady_clock, duration
#include <cmath>                    // for isnan, M_PI
#include <iostream>                 // for string, operator<<, basic_ostream
#include <random>                   // for mt19937, normal_distribution
#include <stdexcept>                // for runtime_error, invalid_argument
#include <string>                   // for operator+, operator==, basic_string
#include <type_traits>              // for enable_if<>::type
#include <unordered_map>            // for unordered_map
#include <vector>                   // for vector, __vector_base<>::value_type
#include "cvector.hpp"              // for CVector
#include "drivers.hpp"              // for ScalarDriver, AxialDriver
#include "noise.hpp"                // for OneFNoise

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 220880.0 // rad/Ts converted to m/As
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2. * M_PI)
#define ELECTRON_CHARGE 1.60217662e-19
#define BOLTZMANN_CONST 1.380649e-23

typedef CVector<double> DVector;
typedef CVector<float> FVector;

double operator"" _ns(unsigned long long timeUnit)
{
    return ((double)timeUnit) / 1e9;
}
double operator"" _ns(long double timeUnit)
{
    return ((double)timeUnit) / 1e9;
}

double operator"" _mT(unsigned long long tesla)
{
    return ((double)tesla) / 1000.0;
}

double operator"" _mT(long double tesla)
{
    return ((double)tesla) / 1000.0;
}

template <typename T>
inline CVector<T> calculate_tensor_interaction(const CVector<T>& m,
    const std::vector<CVector<T>>& tensor,
    const T& Ms)
{
    CVector<T> res(
        tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
        tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
        tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
    return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> calculate_tensor_interaction(const CVector<T>& m,
    const std::array<CVector<T>, 3>& tensor,
    const T& Ms)
{
    CVector<T> res(
        tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
        tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
        tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
    return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> c_cross(const CVector<T>& a, const CVector<T>& b)
{
    CVector<T> res(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);

    return res;
}

template <typename T>
inline T c_dot(const CVector<T>& a, const CVector<T>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
class EnergyDriver
{
public:
    static T calculateZeemanEnergy(CVector<T> mag, CVector<T> Hext, T cellVolume, T Ms)
    {
        return -MAGNETIC_PERMEABILITY * Ms * c_dot<T>(mag, Hext) * cellVolume;
    }

    static T calculateAnisotropyEnergy(CVector<T> mag, CVector<T> anis, T K, T cellVolume)
    {
        const T sinSq = 1.0 - pow(c_dot<T>(mag, anis) / (anis.length() * mag.length()), 2);
        return K * sinSq * cellVolume;
    }

    static T calculateIECEnergy(CVector<T> mag, CVector<T> other, T J, T cellSurface)
    {
        return -c_dot<T>(mag, other) * J * cellSurface;
    }

    static T calculateDemagEnergy(CVector<T> mag, CVector<T> Hdemag, T Ms, T cellVolume)
    {
        return -0.5 * MAGNETIC_PERMEABILITY * Ms * c_dot<T>(mag, Hdemag) * cellVolume;
    }
};

enum Reference
{
    NONE = 0,
    FIXED,
    TOP,
    BOTTOM
};

enum SolverMode
{
    EULER_HEUN = 0,
    RK4 = 1,
    DORMAND_PRICE = 2,
    HEUN = 3
};

// seems to be the faster so far.
// static std::mt19937 generator((std::random_device{}()));

template <typename T = double>
class Layer
{
private:

    OneFNoise<T>* ofn;

    ScalarDriver<T> temperatureDriver;
    ScalarDriver<T> IECDriverTop;
    ScalarDriver<T> IECDriverBottom;
    ScalarDriver<T> IECQuadDriverTop;
    ScalarDriver<T> IECQuadDriverBottom;

    ScalarDriver<T> currentDriver;
    ScalarDriver<T> anisotropyDriver;
    ScalarDriver<T> fieldLikeTorqueDriver;
    ScalarDriver<T> dampingLikeTorqueDriver;
    AxialDriver<T> externalFieldDriver;
    AxialDriver<T> HoeDriver;

    bool nonStochasticTempSet = false;
    bool nonStochasticOneFSet = true;
    bool temperatureSet = false;
    bool pinkNoiseSet = false;
    bool alternativeSTTSet = false;
    Reference referenceType = NONE;

    // the distribution is binded for faster generation
    // is also shared between 1f and Gaussian noise.
    std::function<T()> distribution = std::bind(std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));

    CVector<T> dWn, dWn2; // one for thermal, one for OneF
    Layer(
        const std::string& id,
        CVector<T> mag,
        CVector<T> anis,
        T Ms,
        T thickness,
        T cellSurface,
        const std::vector<CVector<T>>& demagTensor,
        T damping,
        T fieldLikeTorque,
        T dampingLikeTorque,
        T SlonczewskiSpacerLayerParameter,
        T beta,
        T spinPolarisation) : id(id),
        mag(mag),
        anis(anis),
        Ms(Ms),
        thickness(thickness),
        cellSurface(cellSurface),
        demagTensor(demagTensor),
        damping(damping),
        fieldLikeTorque(fieldLikeTorque),
        dampingLikeTorque(dampingLikeTorque),
        SlonczewskiSpacerLayerParameter(SlonczewskiSpacerLayerParameter),
        beta(beta),
        spinPolarisation(spinPolarisation)
    {
        if (mag.length() == 0)
        {
            throw std::runtime_error("Initial magnetisation was set to a zero vector!");
        }
        if (anis.length() == 0)
        {
            throw std::runtime_error("Anisotropy was set to a zero vector!");
        }
        // normalise magnetisation
        mag.normalize();
        dWn = CVector<T>(this->distribution);
        dWn.normalize();
        this->cellVolume = this->cellSurface * this->thickness;
        this->ofn = new OneFNoise<T>(0, 0., 0.);
    }

public:
    bool includeSTT = false;
    bool includeSOT = false;

    std::string id;
    T Ms = 0.0;

    // geometric parameters
    T thickness = 0.0;
    T cellVolume = 0.0, cellSurface = 0.0;

    CVector<T> H_log, Hoe_log, Hconst, mag, anis, referenceLayer;
    CVector<T> Hext, Hdipole, Hdemag, Hoe, HAnis, Hthermal, Hfluctuation;

    CVector<T> Hfl_v, Hdl_v;

    CVector<T> HIEC, HIECtop, HIECbottom;
    T Jbottom_log = 0.0, Jtop_log = 0.0;
    T J2bottom_log = 0.0, J2top_log = 0.0;
    T K_log = 0.0;
    T I_log = 0.0;

    // dipole and demag tensors
    std::vector<CVector<T>> demagTensor;
    std::vector<CVector<T>> dipoleBottom = std::vector<CVector<T>>{ CVector<T>(), CVector<T>(), CVector<T>() };
    std::vector<CVector<T>> dipoleTop = std::vector<CVector<T>>{ CVector<T>(), CVector<T>(), CVector<T>() };

    // LLG params
    T damping;

    // SOT params
    bool dynamicSOT = true;
    T fieldLikeTorque;
    T dampingLikeTorque;

    // STT params
    T SlonczewskiSpacerLayerParameter;
    T beta; // usually either set to 0 or to damping
    T kappa = 1; // for damping-like off -turning torque
    T spinPolarisation;

    T hopt = -1.0;

    Layer() {}
    explicit Layer(const std::string& id,
        const CVector<T>& mag,
        const CVector<T>& anis,
        T Ms,
        T thickness,
        T cellSurface,
        const std::vector<CVector<T>>& demagTensor,
        T damping) : Layer(id, mag, anis, Ms, thickness, cellSurface,
            demagTensor,
            damping, 0, 0, 0, 0, 0) {}

    /**
     * The basic structure is a magnetic layer.
     * Its parameters are defined by the constructor and may be altered
     * by the drivers during the simulation time.
     * If you want STT, remember to set the reference vector for the polarisation of the layer.
     * Use `setReferenceLayer` function to do that.
     * @param id: identifiable name for a layer -- e.g. "bottom" or "free".
     * @param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
     * @param anis: anisotropy of the layer. A normalised vector
     * @param Ms: magnetisation saturation. Unit: Tesla [T].
     * @param thickness: thickness of the layer. Unit: meter [m].
     * @param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
     * @param demagTensor: demagnetisation tensor of the layer.
     * @param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
     * @param fieldLikeTorque: [SOT] effective spin Hall angle (spin effectiveness) for Hfl.
     * @param dampingLikeTorque: [SOT] effective spin Hall angle (spin effectiveness) for Hdl.
     */
    explicit Layer(const std::string& id,
        const CVector<T>& mag,
        const CVector<T>& anis,
        T Ms,
        T thickness,
        T cellSurface,
        const std::vector<CVector<T>>& demagTensor,
        T damping,
        T fieldLikeTorque,
        T dampingLikeTorque) : Layer(id, mag, anis, Ms, thickness, cellSurface,
            demagTensor,
            damping,
            fieldLikeTorque,
            dampingLikeTorque, 0, 0, 0)
    {
        this->includeSTT = false;
        this->includeSOT = true;
        this->dynamicSOT = false;
    }

    /**
     * The basic structure is a magnetic layer.
     * Its parameters are defined by the constructor and may be altered
     * by the drivers during the simulation time.
     * If you want STT, remember to set the reference vector for the polarisation of the layer.
     * Use `setReferenceLayer` function to do that.
     * @param id: identifiable name for a layer -- e.g. "bottom" or "free".
     * @param mag: initial magnetisation. Must be normalised (norm of 1). Used for quicker convergence.
     * @param anis: anisotropy of the layer. A normalised vector
     * @param Ms: magnetisation saturation. Unit: Tesla [T].
     * @param thickness: thickness of the layer. Unit: meter [m].
     * @param cellSurface: surface of the layer, for volume calculation. Unit: meter^2 [m^2].
     * @param demagTensor: demagnetisation tensor of the layer.
     * @param damping: often marked as alpha in the LLG equation. Damping of the layer. Default 0.011. Dimensionless.
     * @param SlomczewskiSpacerLayerParameter: [STT] Slomczewski parameter. Default 1.0. Dimensionless.
     * @param beta: [STT] beta parameter for the STT. Default 0.0. Dimensionless.
     * @param spinPolarisation: [STT] polarisation ratio while passing through reference layer.
     */
    explicit Layer(const std::string& id,
        const CVector<T>& mag,
        const CVector<T>& anis,
        T Ms,
        T thickness,
        T cellSurface,
        const std::vector<CVector<T>>& demagTensor,
        T damping,
        T SlonczewskiSpacerLayerParameter,
        T beta,
        T spinPolarisation) : Layer(id, mag, anis, Ms, thickness, cellSurface,
            demagTensor,
            damping, 0, 0, SlonczewskiSpacerLayerParameter, beta, spinPolarisation)
    {
        this->includeSTT = true;
        this->includeSOT = false;
    }

    inline static Layer<T> LayerSTT(const std::string& id,
        const CVector<T>& mag,
        const CVector<T>& anis,
        T Ms,
        T thickness,
        T cellSurface,
        const std::vector<CVector<T>>& demagTensor,
        T damping,
        T SlonczewskiSpacerLayerParameter,
        T beta,
        T spinPolarisation)
    {
        return Layer<T>(
            id,
            mag,
            anis,
            Ms,
            thickness,
            cellSurface,
            demagTensor,
            damping,
            SlonczewskiSpacerLayerParameter,
            beta,
            spinPolarisation);
    }

    inline static Layer<T> LayerSOT(const std::string& id,
        const CVector<T>& mag,
        const CVector<T>& anis,
        T Ms,
        T thickness,
        T cellSurface,
        const std::vector<CVector<T>>& demagTensor,
        T damping,
        T fieldLikeTorque,
        T dampingLikeTorque)
    {
        return Layer<T>(id,
            mag,
            anis,
            Ms,
            thickness,
            cellSurface,
            demagTensor,
            damping,
            fieldLikeTorque,
            dampingLikeTorque);
    }

    /**
     * @brief Get the Id object
     *
     * @return const std::string
     */
    const std::string getId() const { return id; }

    /**
     * @brief Set the Alternative STT formulation
     *
     * @param alternativeSTT: True if you want to use the alternative STT formulation.
     */
    void setAlternativeSTT(bool alternativeSTT) { this->alternativeSTTSet = alternativeSTT; }
    void setKappa(T kappa) { this->kappa = kappa; }
    void setTopDipoleTensor(const std::vector<CVector<T>>& dipoleTensor)
    {
        this->dipoleTop = dipoleTensor;
    }

    void setBottomDipoleTensor(const std::vector<CVector<T>>& dipoleTensor)
    {
        this->dipoleBottom = dipoleTensor;
    }

    const bool hasTemperature()
    {
        return this->temperatureSet;
    }

    void setTemperatureDriver(const ScalarDriver<T>& driver)
    {
        this->temperatureDriver = driver;
        this->temperatureSet = true;
    }

    void setNonStochasticLangevinDriver(const ScalarDriver<T>& driver)
    {
        this->temperatureDriver = driver;
        // do not set the SDE flag here
        this->temperatureSet = false;
        this->nonStochasticTempSet = true;
    }

    void setOneFNoise(unsigned int sources, T bias, T scale) {
        this->ofn = new OneFNoise<T>(sources, bias, scale);
        this->pinkNoiseSet = true;
        this->nonStochasticOneFSet = true; // by default turn it on, but in the stochastic sims, we will have to turn it off
    }

    void setCurrentDriver(const ScalarDriver<T>& driver)
    {
        this->currentDriver = driver;
    }

    void setFieldLikeTorqueDriver(const ScalarDriver<T>& driver)
    {
        this->includeSOT = true;
        if (this->includeSTT)
            throw std::runtime_error("includeSTT was on and now setting SOT interaction!");
        if (!this->dynamicSOT)
            throw std::runtime_error("used a static SOT definition, now trying to set it dynamically!");
        this->fieldLikeTorqueDriver = driver;
    }

    void setDampingLikeTorqueDriver(const ScalarDriver<T>& driver)
    {
        this->includeSOT = true;
        if (this->includeSTT)
            throw std::runtime_error("includeSTT was on and now setting SOT interaction!");
        if (!this->dynamicSOT)
            throw std::runtime_error("used a static SOT definition, now trying to set it dynamically!");
        this->dampingLikeTorqueDriver = driver;
    }

    void setAnisotropyDriver(const ScalarDriver<T>& driver)
    {
        this->anisotropyDriver = driver;
    }

    void setExternalFieldDriver(const AxialDriver<T>& driver)
    {
        this->externalFieldDriver = driver;
    }
    void setOerstedFieldDriver(const AxialDriver<T>& driver)
    {
        this->HoeDriver = driver;
    }

    void setMagnetisation(CVector<T>& mag)
    {
        if (mag.length() == 0)
        {
            throw std::runtime_error("Initial magnetisation was set to a zero vector!");
        }
        this->mag = mag;
        this->mag.normalize();
    }

    void setIECDriverBottom(const ScalarDriver<T>& driver)
    {
        this->IECDriverBottom = driver;
    }

    void setIECDriverTop(const ScalarDriver<T>& driver)
    {
        this->IECDriverTop = driver;
    }

    void setQuadIECDriverTop(const ScalarDriver<T>& driver)
    {
        this->IECQuadDriverTop = driver;
    }

    void setQuadIECDriverBottom(const ScalarDriver<T>& driver)
    {
        this->IECQuadDriverBottom = driver;
    }

    /**
     * @brief Sets reference layer with a custom vector
     * Set reference layer parameter. This is for calculating the spin current
     * polarisation if `includeSTT` is true.
     * @param reference: CVector describing the reference layer.
     */
    void setReferenceLayer(const CVector<T>& reference)
    {
        this->referenceLayer = reference;
        this->referenceType = FIXED;
    }

    /**
     * @brief Set reference layer with enum
     * Can be used to refer to other layers in stack as reference
     * for this layer.
     * @param reference: an enum: FIXED, TOP, BOTTOM, or CUSTOM
     */
    void setReferenceLayer(Reference reference)
    {
        if ((reference == FIXED) && (!this->referenceLayer.length()))
        {
            throw std::runtime_error("Cannot set fixed polarisation layer to 0!"
                " Set reference to NONE to disable reference.");
        }
        this->referenceType = reference;
    }


    /**
     * @brief Get the Reference Layer object
     */
    CVector<T> getReferenceLayer()
    {
        // TODO: return other mags when the reference layer is not fixed.
        return this->referenceLayer;
    }

    /**
     * @brief Get the Reference Layer Type object (enum type is returned)
     */
    Reference getReferenceType()
    {
        return this->referenceType;
    }

    const CVector<T> calculateHeff(T time, T timeStep,
        const CVector<T>& stepMag, const CVector<T>& bottom, const CVector<T>& top,
        const CVector<T>& Hfluctuation = CVector<T>())
    {
        this->Hdipole = calculate_tensor_interaction(bottom, this->dipoleBottom, this->Ms) +
            calculate_tensor_interaction(top, this->dipoleTop, this->Ms);
        return calculateHeffDipoleInjection(time, timeStep, stepMag, bottom, top, this->Hdipole, Hfluctuation);
    }

    const CVector<T> calculateHeffDipoleInjection(T time, T timeStep,
        const CVector<T>& stepMag, const CVector<T>& bottom, const CVector<T>& top,
        const CVector<T>& dipole, const CVector<T>& Hfluctuation)
    {
        this->Hext = calculateExternalField(time);
        this->Hoe = calculateHOeField(time);

        this->Hdemag = calculate_tensor_interaction(stepMag, this->demagTensor, this->Ms);
        this->HIEC = calculateIEC(time, stepMag, bottom, top);
        this->HAnis = calculateAnisotropy(stepMag, time);
        const CVector<T> Heff = this->Hext    // external
            + this->HAnis // anistotropy
            + this->HIEC  // IEC
            + this->Hoe   // Oersted field
            + Hfluctuation
            // demag -- negative contribution
            - this->Hdemag
            // dipole -- negative contribution
            - dipole;
        return Heff;
    }

    CVector<T> calculateHOeField(const T& time)
    {
        this->Hoe_log = this->HoeDriver.getCurrentAxialDrivers(time);
        return this->Hoe_log;
    }

    CVector<T> calculateExternalField(const T& time)
    {
        this->H_log =
            this->externalFieldDriver.getCurrentAxialDrivers(time);
        return this->H_log;
    }

    CVector<T> calculateAnisotropy(const CVector<T>& stepMag, T& time)
    {
        this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
        const T nom = (2 * this->K_log) * c_dot<T>(this->anis, stepMag) / (this->Ms);
        return this->anis * nom;
    }

    CVector<T> calculateIEC_(const T J, const T J2, const CVector<T>& stepMag, const CVector<T>& coupledMag)
    {
        // below an alternative method for computing J -- it's here for reference only.
        // const T nom = J / (this->Ms * this->thickness);
        // return (coupledMag - stepMag) * nom; // alternative form
        // return (coupledMag + coupledMag * 2 * J2 * c_dot(coupledMag, stepMag)) * nom;
        return coupledMag * (J + 2 * J2 * c_dot(coupledMag, stepMag)) / (this->Ms * this->thickness);
    }

    CVector<T> calculateIEC(T time, const CVector<T>& stepMag, const CVector<T>& bottom, const CVector<T>& top)
    {
        this->Jbottom_log = this->IECDriverBottom.getCurrentScalarValue(time);
        this->Jtop_log = this->IECDriverTop.getCurrentScalarValue(time);

        this->J2bottom_log = this->IECQuadDriverBottom.getCurrentScalarValue(time);
        this->J2top_log = this->IECQuadDriverTop.getCurrentScalarValue(time);

        return calculateIEC_(this->Jbottom_log,
            this->J2bottom_log, stepMag, bottom) +
            calculateIEC_(this->Jtop_log, this->J2top_log, stepMag, top);
    }


    /**
     * @brief Main solver function. It is solver-independent (all solvers use this function).
     * This function is called by the solver to calculate the next step of the magnetisation.
     * It computes implicitly, all torques, given the current magnetisation and effective field.
     * @param time the time at which the solver is currently at.
     * @param m the current magnetisation (from the solver, may be a semi-step)
     * @param timeStep integration time
     * @param bottom magnetisation of the layer below
     * @param top magnetisation of the layer above
     * @param heff the effective field
     * @return const CVector<T> magnetisation after the step
     */
    const CVector<T> solveLLG(T time, const CVector<T>& m, T timeStep,
        const CVector<T>& bottom, const CVector<T>& top, const CVector<T>& heff)
    {
        const CVector<T> prod = c_cross<T>(m, heff);
        const CVector<T> prod2 = c_cross<T>(m, prod);
        const T convTerm = 1 / (1 + pow(this->damping, 2)); // LLGS -> LL form
        const CVector<T> dmdt = prod + prod2 * this->damping;
        CVector<T> reference;

        // decide what is to be the reference for (s)LLG-STT
        // dynamically substitute other active layers
        switch (this->referenceType)
        {
            // TODO: add the warning if reference layer is top/bottom and empty
        case FIXED:
            reference = this->referenceLayer;
            break;
        case TOP:
            reference = top;
            break;
        case BOTTOM:
            reference = bottom;
            break;
        default:
            break;
        }

        // extra terms
        if (this->includeSTT)
        {
            this->I_log = this->currentDriver.getCurrentScalarValue(time);
            // use standard STT formulation
            const T aJ = HBAR * this->I_log /
                (ELECTRON_CHARGE * this->Ms * this->thickness);
            // field like
            T eta = 0;
            if (this->alternativeSTTSet) {
                // this is simplified
                eta = (this->spinPolarisation) / (1 + this->SlonczewskiSpacerLayerParameter * c_dot<T>(m, reference));
            }
            else {
                // this is more complex model (classical STT)
                const T slonSq = pow(this->SlonczewskiSpacerLayerParameter, 2);
                eta = (this->spinPolarisation * slonSq) / (slonSq + 1 + (slonSq - 1) * c_dot<T>(m, reference));
            }
            const T sttTerm = GYRO * aJ * eta;
            const CVector<T> fieldLike = c_cross<T>(m, reference);
            // damping like
            const CVector<T> dampingLike = c_cross<T>(m, fieldLike);
            return (dmdt * -GYRO + dampingLike * -sttTerm * this->kappa + fieldLike * sttTerm * this->beta) * convTerm;
        }
        else if (this->includeSOT)
        {
            T Hdl, Hfl;
            // I log current density
            // use SOT formulation with effective DL and FL fields
            if (this->dynamicSOT)
            {
                // dynamic SOT is set when the driver is present
                Hdl = this->dampingLikeTorqueDriver.getCurrentScalarValue(time);
                Hfl = this->fieldLikeTorqueDriver.getCurrentScalarValue(time);
            }
            else
            {
                this->I_log = this->currentDriver.getCurrentScalarValue(time);
                Hdl = this->dampingLikeTorque * this->I_log;
                Hfl = this->fieldLikeTorque * this->I_log;
            }
            this->Hfl_v = reference * (Hfl - this->damping * Hdl);
            this->Hdl_v = reference * (Hdl + this->damping * Hfl);
            const CVector<T> cm = c_cross<T>(m, reference);
            const CVector<T> ccm = c_cross<T>(m, cm);
            const CVector<T> flTorque = cm * (Hfl - this->damping * Hdl);
            const CVector<T> dlTorque = ccm * (Hdl + this->damping * Hfl);
            return (dmdt + flTorque + dlTorque) * -GYRO * convTerm;
        }
        return dmdt * -GYRO * convTerm;
    }
    const CVector<T> calculateLLGWithFieldTorqueDipoleInjection(T time, const CVector<T>& m,
        const CVector<T>& bottom, const CVector<T>& top,
        const CVector<T>& dipole, T timeStep, const CVector<T>& Hfluctuation = CVector<T>())
    {
        // classic LLG first
        const CVector<T> heff = calculateHeffDipoleInjection(time, timeStep, m, bottom, top, dipole, Hfluctuation);
        return solveLLG(time, m, timeStep, bottom, top, heff);
    }

    /**
     * Compute the LLG time step. The efficient field vectors is calculated implicitly here.
     * Use the effective spin hall angles formulation for SOT interaction.
     * @param time: current simulation time.
     * @param m: current RK45 magnetisation.
     * @param bottom: layer below the current layer (current layer's magnetisation is m). For IEC interaction.
     * @param top: layer above the current layer (current layer's magnetisation is m). For IEC interaction.
     * @param timeStep: RK45 integration step.
     */
    const CVector<T> calculateLLGWithFieldTorque(T time, const CVector<T>& m, const CVector<T>& bottom,
        const CVector<T>& top, T timeStep, const CVector<T>& Hfluctuation = CVector<T>())
    {
        // classic LLG first
        const CVector<T> heff = calculateHeff(time, timeStep, m, bottom, top, Hfluctuation);
        return solveLLG(time, m, timeStep, bottom, top, heff);
    }


    /**
     * @brief RK4 step of the LLG equation.
     * Compute the LLG time step. The efficient field vectors is calculated implicitly here.
     * Use the effective spin hall angles formulation for SOT interaction.
     * @param time: current simulation time.
     * @param m: current RK45 magnetisation.
     * @param bottom: layer below the current layer (current layer's magnetisation is m). For IEC interaction.
     * @param top: layer above the current layer (current layer's magnetisation is m). For IEC interaction.
     * @param timeStep: RK45 integration step.
     */
    void rk4_step(T time, T timeStep, const CVector<T>& bottom, const CVector<T>& top)
    {
        CVector<T> m_t = this->mag;
        const CVector<T> k1 = calculateLLGWithFieldTorque(time, m_t, bottom, top, timeStep) * timeStep;
        const CVector<T> k2 = calculateLLGWithFieldTorque(time + 0.5 * timeStep, m_t + k1 * 0.5, bottom, top, timeStep) * timeStep;
        const CVector<T> k3 = calculateLLGWithFieldTorque(time + 0.5 * timeStep, m_t + k2 * 0.5, bottom, top, timeStep) * timeStep;
        const CVector<T> k4 = calculateLLGWithFieldTorque(time + timeStep, m_t + k3, bottom, top, timeStep) * timeStep;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        this->mag = m_t;
        if (isnan(this->mag.x))
        {
            throw std::runtime_error("NAN magnetisation");
        }
    }

    /**
     * @brief RK4 step of the LLG equation if dipole injection is present.
     * Compute the LLG time step. The efficient field vectors is calculated implicitly here.
     * Use the effective spin hall angles formulation for SOT interaction.
     * @param time: current simulation time.
     * @param m: current RK45 magnetisation.
     * @param bottom: layer below the current layer (current layer's magnetisation is m). For IEC interaction.
     * @param top: layer above the current layer (current layer's magnetisation is m). For IEC interaction.
     * @param timeStep: RK45 integration step.
     */
    void rk4_stepDipoleInjection(T time, T timeStep, const CVector<T>& bottom, const CVector<T>& top, const CVector<T>& dipole)
    {
        CVector<T> m_t = this->mag;
        const CVector<T> k1 = calculateLLGWithFieldTorqueDipoleInjection(time, m_t, bottom, top, dipole, timeStep) * timeStep;
        const CVector<T> k2 = calculateLLGWithFieldTorqueDipoleInjection(time + 0.5 * timeStep, m_t + k1 * 0.5, bottom, top, dipole, timeStep) * timeStep;
        const CVector<T> k3 = calculateLLGWithFieldTorqueDipoleInjection(time + 0.5 * timeStep, m_t + k2 * 0.5, bottom, top, dipole, timeStep) * timeStep;
        const CVector<T> k4 = calculateLLGWithFieldTorqueDipoleInjection(time + timeStep, m_t + k3, bottom, top, dipole, timeStep) * timeStep;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        this->mag = m_t;
    }

    void dormandPriceStep(T time, T timeStep, const CVector<T>& bottom, const CVector<T>& top)
    {
        CVector<T> m_t = this->mag;
        CVector<T> e_t;
        std::array<CVector<T>, 7> K;
        // if (this->hopt < 0)
        // {
        // }
        // this makes the step non-adaptive
        // we will deal with this problem later, this requires consistency in step
        // across all the layers
        this->hopt = timeStep;
        // define Buchter tableau below
        // there are redundant zeros, but for clarity, we'll leave them
        const std::array<double, 7> c = {
            0., 1. / 5., 3. / 10., 4. / 5., 8. / 9., 1., 1. };
        const std::array<double, 7> b = {
            35. / 384., 0, 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0. };
        const std::array<double, 7> b2 = {
            5179. / 57600., 0, 7571. / 16695., 393. / 640., -92097. / 339200., 187. / 2100., 1. / 40. };
        // extra braces are required even though the struct is 2D only
        const std::array<std::array<double, 7>, 7> aCoefs = { {{0, 0, 0, 0, 0, 0, 0},
                                                              {1. / 5., 0, 0, 0, 0, 0, 0},
                                                              {3. / 40., 9. / 40., 0, 0, 0, 0, 0},
                                                              {44. / 45., -56. / 15., 32. / 9., 0, 0, 0, 0},
                                                              {19372. / 6561., -25360. / 2187., 64448. / 6561., -212. / 729., 0, 0, 0},
                                                              {9017. / 3168., -355. / 33., 46732. / 5247., 49. / 176., -5103. / 18656., 0, 0},
                                                              {35. / 384., 0., 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0.}} };
        // compute the first
        K[0] = calculateLLGWithFieldTorque(time, this->mag, bottom, top, this->hopt);
        m_t = m_t + K[0] * b[0] * this->hopt;
        for (int i = 1; i < 7; i++)
        {
            CVector<T> kS;
            for (int j = 0; j < i; j++)
            {
                kS = kS + K[j] * aCoefs[i][j];
            }
            K[i] = calculateLLGWithFieldTorque(time + c[i] * this->hopt, this->mag + kS * this->hopt,
                bottom, top, this->hopt);
            m_t = m_t + K[i] * b[i] * this->hopt;           // this is function estimate
            e_t = e_t + K[i] * (b[i] - b2[i]) * this->hopt; // this is error estimate
        }
        // adapt the step size
        const T eps = 1e-3;
        const T s = pow((eps * timeStep) / (2 * e_t.length()), 1. / 5.);
        this->hopt = s * timeStep;
        m_t.normalize();
        this->mag = m_t;
        if (isnan(this->mag.x))
        {
            throw std::runtime_error("NAN magnetisation");
        }
    }

    CVector<T> non_stochastic_llg(const CVector<T>& cm, T time, T timeStep, const CVector<T>& bottom, const CVector<T>& top)
    {
        return calculateLLGWithFieldTorque(time, cm, bottom, top, timeStep);
    }

    CVector<T> stochastic_llg(const CVector<T>& cm, T time, T timeStep,
        const  CVector<T>& bottom, const CVector<T>& top, const CVector<T>& dW, const CVector<T>& dW2, const T& HoneF)
    {
        // compute the Langevin fluctuations -- this is the sigma
        const T convTerm = -GYRO / (1 + pow(this->damping, 2));
        const T Hthermal_temp = this->getLangevinStochasticStandardDeviation(time, timeStep);
        const CVector<T> thcross = c_cross(cm, dW);
        const CVector<T> thcross2 = c_cross(thcross, dW);
        const T scalingTh = Hthermal_temp * convTerm;

        // compute 1/f noise term
        const CVector<T> onefcross = c_cross(cm, dW2);
        const CVector<T> onefcross2 = c_cross(onefcross, dW2);
        const T scalingOneF = HoneF * convTerm;

        return (thcross + thcross2 * this->damping) * scalingTh + (onefcross + onefcross2 * this->damping) * scalingOneF;
    }

    const T getStochasticOneFNoise(T time) {
        if (!this->pinkNoiseSet)
            return 0;
        return this->ofn->tick();
    }

    T getLangevinStochasticStandardDeviation(T time, T timeStep)
    {
        if (this->cellVolume == 0.0)
            throw std::runtime_error("Cell surface cannot be 0 during temp. calculations!");
        const T currentTemp = this->temperatureDriver.getCurrentScalarValue(time);
        const T mainFactor = (2 * this->damping * MAGNETIC_PERMEABILITY * BOLTZMANN_CONST * currentTemp) / (this->Ms * this->cellVolume * timeStep);
        return sqrt(mainFactor);
    }

    CVector<T> nonStochasticLangevin(T time, T timeStep)
    {
        const T Hthermal_temp = this->getLangevinStochasticStandardDeviation(time, timeStep);
        const CVector<T> dW = CVector<T>(this->distribution);
        return dW * Hthermal_temp;
    }

    CVector<T> nonStochasticOneFNoise(T time, T timestep) {
        const T pinkNoise = this->ofn->tick();
        const CVector<T> dW2 = CVector<T>(this->distribution);
        return dW2 * pinkNoise;
    }

    /**
     * @brief Computes a single Euler-Heun step [DEPRECATED].
     * [DEPRECATED] This is the old Euler-Heun method, Heun is preferred.
     * Bottom and top are relative to the current layer.
     * They are used to compute interactions.
     * @param time: current time of the simulation
     * @param timeStep: integration time of the solver
     * @param bottom: bottom layer to the current layer
     * @param top: top layer to the current layer
     */
    void euler_heun_step(T time, T timeStep, const CVector<T>& bottom, const CVector<T>& top)
    {
        // we compute the two below in stochastic part, not non stochastic.
        this->nonStochasticTempSet = false;
        this->nonStochasticOneFSet = false;
        // this is Stratonovich integral
        if (isnan(this->mag.x))
        {
            throw std::runtime_error("NAN magnetisation");
        }
        // Brownian motion sample
        // Generate the noise from the Brownian motion
        // dW2 is used for 1/f noise generation
        CVector<T> dW = CVector<T>(this->distribution); // * sqrt(timeStep);
        CVector<T> dW2 = CVector<T>(this->distribution); //* sqrt(timeStep);
        // squared dW -- just utility
        dW.normalize();
        dW2.normalize();
        // f_n is the vector of non-stochastic part at step n
        // multiply by timeStep (h) here for convenience
        const T Honef = this->getStochasticOneFNoise(time);
        const CVector<T> f_n = non_stochastic_llg(this->mag, time, timeStep, bottom, top) * timeStep;
        // g_n is the stochastic part of the LLG at step n
        const CVector<T> g_n = stochastic_llg(this->mag, time, timeStep, bottom, top, dW, dW2, Honef) * timeStep;

        // actual solution
        // approximate next step ytilde
        const CVector<T> mapprox = this->mag + g_n;
        // calculate the approx g_n
        const CVector<T> g_n_approx = stochastic_llg(mapprox, time, timeStep, bottom, top, dW, dW2, Honef) * timeStep;
        // CVector<T> m_t = this->mag + f_n + g_n + (g_n_approx - g_n) * 0.5;
        CVector<T> m_t = this->mag + f_n + (g_n_approx + g_n) * 0.5;
        m_t.normalize();
        this->mag = m_t;
    }

    /**
     * @brief Computes a single Heun step.
     * This method is preferred over Euler-Heun method.
     * Bottom and top are relative to the current layer.
     * They are used to compute interactions.
     * @param time: current time of the simulation
     * @param timeStep: integration time of the solver
     * @param bottom: bottom layer to the current layer
     * @param top: top layer to the current layer
     */
    void heun_step(T time, T timeStep, const CVector<T>& bottom, const CVector<T>& top) {
        // we compute the two below in stochastic part, not non stochastic.
        this->nonStochasticTempSet = false;
        this->nonStochasticOneFSet = false;
        // this is Stratonovich integral
        if (isnan(this->mag.x))
        {
            throw std::runtime_error("NAN magnetisation");
        }
        // Brownian motion sample
        // Generate the noise from the Brownian motion
        const T Honef_scale = this->getStochasticOneFNoise(time);
        const T Hthermal_scale = this->getLangevinStochasticStandardDeviation(time, timeStep);
        const CVector<T> Hlangevin = dWn * Hthermal_scale;
        const CVector<T> Honef = dWn2 * Honef_scale;
        const CVector<T> m_t = this->mag;
        const CVector<T> f_n = this->calculateLLGWithFieldTorque(time, m_t, bottom, top, timeStep, Hlangevin + Honef);
        // immediate m approximation
        CVector<T> m_approx = m_t + f_n * timeStep;
        CVector<T> dW = CVector<T>(this->distribution);
        CVector<T> dW2 = CVector<T>(this->distribution);
        dW.normalize();
        dW2.normalize();
        m_approx.normalize();
        const CVector<T> Hlangevin_approx = dW * Hthermal_scale;
        const CVector<T> Honef_approx = dW2 * Honef_scale;
        const CVector<T> f_approx = this->calculateLLGWithFieldTorque(time + timeStep,
            m_approx, bottom, top, timeStep, Hlangevin_approx + Honef_approx);
        dWn = dW; // replace
        dWn2 = dW2;
        CVector<T> nm_t = this->mag + (f_n + f_approx) * 0.5 * timeStep;
        nm_t.normalize();
        this->mag = nm_t;
    }
};

template <typename T>
class Junction
{
    friend class Layer<T>;
    const std::vector<std::string> vectorNames = { "x", "y", "z" };

public:
    enum MRmode
    {
        NONE = 0,
        CLASSIC = 1,
        STRIP = 2
    };

    MRmode MR_mode;
    std::vector<Layer<T>> layers;
    T Rp, Rap = 0.0;

    std::vector<T> Rx0, Ry0, AMR_X, AMR_Y, SMR_X, SMR_Y, AHE;
    std::unordered_map<std::string, std::vector<T>> log;
    // std::string fileSave;
    unsigned int logLength = 0;
    unsigned int layerNo;
    std::string Rtag = "R";

    Junction() {}

    /**
     * @brief Create a plain junction.
     * No magnetoresistance is calculated.
     * @param layersToSet: layers that compose the junction
     */
    explicit Junction(const std::vector<Layer<T>>& layersToSet)
    {
        this->MR_mode = NONE;
        this->layers = layersToSet;
        this->layerNo = this->layers.size();
        if (this->layerNo == 0)
        {
            throw std::invalid_argument("Passed a zero length Layer vector!");
        }
    }
    explicit Junction(const std::vector<Layer<T>>& layersToSet, T Rp, T Rap) : Junction(
        layersToSet)
    {
        if (this->layerNo == 1)
        {
            // we need to check if this layer has a reference layer.
            if (!this->layers[0].referenceLayer.length())
            {
                throw std::invalid_argument("MTJ with a single layer must have"
                    " a pinning (referenceLayer) set!");
            }
        }
        if (this->layerNo > 2)
        {
            throw std::invalid_argument("This constructor supports only bilayers!"
                " Choose the other one with the strip resistance!");
        }
        this->Rp = Rp;
        this->Rap = Rap;
        this->MR_mode = CLASSIC;
        // A string representing the tag for the junction's resistance value.
        this->Rtag = "R_" + this->layers[0].id + "_" + this->layers[1].id;
    }

    /**
     * Creates a junction with a STRIP magnetoresistance.
     * Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
     * length of the layers passed (they directly correspond to each layer).
     * Calculates the magnetoresistance as per: __see reference__:
     * Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
     * @param Rx0
     * @param Ry0
     * @param AMR_X
     * @param AMR_Y
     * @param SMR_X
     * @param SMR_Y
     * @param AHE
     */
    explicit Junction(const std::vector<Layer<T>>& layersToSet,
        std::vector<T> Rx0,
        std::vector<T> Ry0,
        std::vector<T> AMR_X,
        std::vector<T> AMR_Y,
        std::vector<T> SMR_X,
        std::vector<T> SMR_Y,
        std::vector<T> AHE) : Rx0(std::move(Rx0)),
        Ry0(std::move(Ry0)),
        AMR_X(std::move(AMR_X)),
        AMR_Y(std::move(AMR_Y)),
        SMR_X(std::move(SMR_X)),
        SMR_Y(std::move(SMR_Y)),
        AHE(std::move(AHE))

    {
        this->layers = std::move(layersToSet);
        this->layerNo = this->layers.size();
        if (this->layerNo == 0)
        {
            throw std::invalid_argument("Passed a zero length Layer vector!");
        }
        if ((this->layerNo != (unsigned int)this->Rx0.size()) ||
            (this->layerNo != (unsigned int)this->Ry0.size()) ||
            (this->layerNo != (unsigned int)this->AMR_X.size()) ||
            (this->layerNo != (unsigned int)this->AMR_Y.size()) ||
            (this->layerNo != (unsigned int)this->AHE.size()) ||
            (this->layerNo != (unsigned int)this->SMR_X.size()) ||
            (this->layerNo != (unsigned int)this->SMR_Y.size()))
        {
            throw std::invalid_argument("Layers and Rx0, Ry, AMR, AMR and SMR must be of the same size!");
        }
        // this->fileSave = std::move(filename);
        this->MR_mode = STRIP;
    }

    /**
     * @brief Select a solver for the LLGS/sLLGS solutions
     * Available: DormandPrice, EulerHeun, RK4
     * @param solverMode one of EULER_HEUN, DORMAND_PRICE or RK4
     */
    const auto selectSolver(SolverMode solverMode)
    {
        switch (solverMode)
        {
        case EULER_HEUN:
            return &Layer<T>::euler_heun_step;

        case DORMAND_PRICE:
            return &Layer<T>::dormandPriceStep;

        case RK4:
            return &Layer<T>::rk4_step;

        case HEUN:
            return &Layer<T>::heun_step;

        default:
            return &Layer<T>::rk4_step;
        }
    }

    /**
     * @brief Get Ids of the layers in the junction.
     * @return vector of layer ids.
     */
    const std::vector<std::string> getLayerIds() const
    {
        std::vector<std::string> ids;
        std::transform(this->layers.begin(), this->layers.end(), std::back_inserter(ids),
            [](const Layer<T>& layer) { return layer.id; });
        return ids;
    }

    /**
     * Clears the simulation log.
     **/
    void clearLog()
    {
        this->log.clear();
        this->logLength = 0;
    }

    std::unordered_map<std::string, std::vector<T>>& getLog()
    {
        return this->log;
    }

    typedef void (Layer<T>::* scalarDriverSetter)(const ScalarDriver<T>& driver);
    typedef void (Layer<T>::* axialDriverSetter)(const AxialDriver<T>& driver);
    void scalarlayerSetter(const std::string& layerID, scalarDriverSetter functor, ScalarDriver<T> driver)
    {
        bool found = false;
        for (auto& l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                (l.*functor)(driver);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }
    void axiallayerSetter(const std::string& layerID, axialDriverSetter functor, AxialDriver<T> driver)
    {
        bool found = false;
        for (auto& l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                (l.*functor)(driver);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }
    void setLayerTemperatureDriver(const std::string& layerID, const ScalarDriver<T>& driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setTemperatureDriver, driver);
    }
    void setLayerNonStochasticLangevinDriver(const std::string& layerID, const ScalarDriver<T>& driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setNonStochasticLangevinDriver, driver);
    }
    void setLayerAnisotropyDriver(const std::string& layerID, const ScalarDriver<T>& driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setAnisotropyDriver, driver);
    }
    void setLayerExternalFieldDriver(const std::string& layerID, const AxialDriver<T>& driver)
    {
        axiallayerSetter(layerID, &Layer<T>::setExternalFieldDriver, driver);
    }
    void setLayerOerstedFieldDriver(const std::string& layerID, const AxialDriver<T>& driver)
    {
        axiallayerSetter(layerID, &Layer<T>::setOerstedFieldDriver, driver);
    }
    void setLayerCurrentDriver(const std::string& layerID, const ScalarDriver<T>& driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setCurrentDriver, driver);
    }
    void setLayerDampingLikeTorqueDriver(const std::string& layerID, const ScalarDriver<T>& driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setDampingLikeTorqueDriver, driver);
    }
    void setLayerFieldLikeTorqueDriver(const std::string& layerID, const ScalarDriver<T>& driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setFieldLikeTorqueDriver, driver);
    }

    void setLayerAlternativeSTT(const std::string& layerID, const bool alternative)
    {
        if (layerID == "all")
        {
            for (auto& l : this->layers)
            {
                l.setAlternativeSTT(alternative);
            }
        }
        else
            getLayer(layerID).setAlternativeSTT(alternative);
    }

    void setLayerOneFNoise(const std::string& layerID, unsigned int sources, T bias, T scale) {

        if (layerID == "all")
        {
            for (auto& l : this->layers)
            {
                l.setOneFNoise(sources, bias, scale);
            }
        }
        else
            getLayer(layerID).setOneFNoise(sources, bias, scale);
    }

    /**
     * Set IEC interaction between two layers.
     * The names of the params are only for convention. The IEC will be set
     * between bottomLayer or topLayer, order is irrelevant.
     * @param bottomLayer: the first layer id
     * @param topLayer: the second layer id
     */
    void setIECDriver(const std::string& bottomLayer, const std::string& topLayer, const ScalarDriver<T>& driver)
    {
        bool found = false;
        for (unsigned int i = 0; i < this->layerNo - 1; i++)
        {
            // check if the layer above is actually top layer the user specified
            if ((this->layers[i].id == bottomLayer) && (this->layers[i + 1].id == topLayer))
            {
                this->layers[i].setIECDriverTop(driver);
                this->layers[i + 1].setIECDriverBottom(driver);
                found = true;
                break;
            }
            else if ((this->layers[i].id == topLayer) && (this->layers[i + 1].id == bottomLayer))
            {
                this->layers[i].setIECDriverTop(driver);
                this->layers[i + 1].setIECDriverBottom(driver);
                found = true;
                break;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to match the layer order or find layer ids!");
        }
    }

    void setQuadIECDriver(const std::string& bottomLayer, const std::string& topLayer, const ScalarDriver<T>& driver)
    {
        bool found = false;
        for (unsigned int i = 0; i < this->layerNo - 1; i++)
        {
            // check if the layer above is actually top layer the user specified
            if ((this->layers[i].id == bottomLayer) && (this->layers[i + 1].id == topLayer))
            {
                this->layers[i].setQuadIECDriverTop(driver);
                this->layers[i + 1].setQuadIECDriverBottom(driver);
                found = true;
                break;
            }
            else if ((this->layers[i].id == topLayer) && (this->layers[i + 1].id == bottomLayer))
            {
                this->layers[i].setQuadIECDriverTop(driver);
                this->layers[i + 1].setQuadIECDriverBottom(driver);
                found = true;
                break;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to match the layer order or find layer ids!");
        }
    }

    void setLayerMagnetisation(const std::string& layerID, CVector<T>& mag)
    {
        bool found = false;
        for (auto& l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                l.setMagnetisation(mag);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }

    CVector<T> getLayerMagnetisation(const std::string& layerID)
    {
        return getLayer(layerID).mag;
    }

    Reference getLayerReferenceType(const std::string& layerID)
    {
        return getLayer(layerID).referenceType;
    }

    void setLayerReferenceLayer(const std::string& layerID, const CVector<T>& referenceLayer)
    {
        if (layerID == "all")
        {
            for (auto& l : this->layers)
            {
                l.setReferenceLayer(referenceLayer);
            }
        }
        else
            getLayer(layerID).setReferenceLayer(referenceLayer);
    }

    void setLayerReferenceType(const std::string& layerID, Reference referenceType)
    {
        if (layerID == "all")
        {
            for (auto& l : this->layers)
            {
                l.setReferenceLayer(referenceType);
            }
        }
        else
            getLayer(layerID).setReferenceLayer(referenceType);
    }

    Layer<T>& getLayer(const std::string& layerID)
    {
        const auto res = std::find_if(
            this->layers.begin(), this->layers.end(),
            [layerID](const auto& l) -> bool {return (l.id == layerID);}
        );
        if (res != this->layers.end()) {
            return *res;
        }
        throw std::runtime_error("Failed to find a layer with a given id!");
    }

    /**
     * @brief Log computed layer parameters.
     * This function logs all the necessayr parameters of the layers.
     * @param t: current time
     * @param timeStep: timeStep of the simulation (unsued for now)
     * @param calculateEnergies: if true, also include fields for energy computation.
     */
    void logLayerParams(T& t, T timeStep, bool calculateEnergies = false)
    {
        for (const auto& layer : this->layers)
        {
            const std::string lId = layer.id;

            if (calculateEnergies)
            {
                // TODO: avoid recomputation at a cost of a slight error
                // recompute the current Heff to avoid shadow persistence of the layer parameters
                // const CVector<T> heff = calculateHeff(t, timeStep, layer.m, layer.bottom, layer.top);
                this->log[lId + "_K"].emplace_back(layer.K_log);
                this->log[lId + "_Jbottom"].emplace_back(layer.Jbottom_log);
                this->log[lId + "_Jtop"].emplace_back(layer.Jtop_log);
                this->log[lId + "_I"].emplace_back(layer.I_log);
                for (int i = 0; i < 3; i++)
                {
                    this->log[lId + "_Hext" + vectorNames[i]].emplace_back(layer.Hext[i]);
                    this->log[lId + "_Hiec" + vectorNames[i]].emplace_back(layer.HIEC[i]);
                    this->log[lId + "_Hanis" + vectorNames[i]].emplace_back(layer.HAnis[i]);
                    this->log[lId + "_Hdemag" + vectorNames[i]].emplace_back(layer.Hdemag[i]);
                    this->log[lId + "_Hth" + vectorNames[i]].emplace_back(layer.Hfluctuation[i]);
                    if (layer.includeSOT)
                    {
                        this->log[lId + "_Hfl" + vectorNames[i]].emplace_back(layer.Hfl_v[i]);
                        this->log[lId + "_Hdl" + vectorNames[i]].emplace_back(layer.Hdl_v[i]);
                    }
                }
                if (layer.includeSTT | layer.includeSOT)
                    this->log[lId + "_I"].emplace_back(layer.I_log);
            }
            // always save magnetisation
            for (int i = 0; i < 3; i++)
            {
                this->log[lId + "_m" + vectorNames[i]].emplace_back(layer.mag[i]);
            }
        }
        if (this->MR_mode == CLASSIC && this->layerNo == 1)
        {
            this->log["R"].emplace_back(calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[0].referenceLayer)));
        }
        else if (MR_mode == CLASSIC && this->layerNo > 1)
        {
            const auto magnetoresistance = calculateMagnetoresistance(c_dot<T>(this->layers[0].mag,
                this->layers[1].mag));
            this->log[this->Rtag].emplace_back(magnetoresistance);
        }
        else if (MR_mode == STRIP)
        {
            const auto magnetoresistance = stripMagnetoResistance(this->Rx0,
                this->Ry0,
                this->AMR_X,
                this->SMR_X,
                this->AMR_Y,
                this->SMR_Y,
                this->AHE);
            this->log["Rx"].emplace_back(magnetoresistance[0]);
            this->log["Ry"].emplace_back(magnetoresistance[1]);
            this->log["Rz"].emplace_back(magnetoresistance[2]);
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

    typedef void (Layer<T>::* solverFn)(T t, T timeStep, const CVector<T>& bottom, const CVector<T>& top);

    /**
     * @brief Run Euler-Heun or RK4 method for a single layer.
     *
     * The Euler-Heun method should only be used
     * for stochastic simulations where the temperature
     * driver is set.
     * @param functor: solver function.
     * @param t: current time
     * @param timeStep: integration step
     */
    void runSingleLayerSolver(solverFn& functor, T& t, T& timeStep)
    {
        CVector<T> null;
        (this->layers[0].*functor)(
            t, timeStep, null, null);
    }

    /**
     * @brief Select a solver based on the setup.
     *
     * Multilayer layer solver iteration.
     * @param functor: solver function.
     * @param t: current time
     * @param timeStep: integration step
     * */
    void runMultiLayerSolver(solverFn& functor, T& t, T& timeStep)
    {
        // initialise with 0 CVectors
        std::vector<CVector<T>> magCopies(this->layerNo + 2, CVector<T>());
        // the first and the last layer get 0 vector coupled
        for (unsigned int i = 0; i < this->layerNo; i++)
        {
            magCopies[i + 1] = this->layers[i].mag;
        }

        for (unsigned int i = 0; i < this->layerNo; i++)
        {
            (this->layers[i].*functor)(
                t, timeStep, magCopies[i], magCopies[i + 2]);
        }
    }

    /**
     * @brief Calculate strip magnetoresistance for multilayer.
     *
     * Used when MR_MODE == STRIP
     * Magnetoresistance as per:
     * Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
     * Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
     * length of the layers passed (they directly correspond to each layer).
     * Calculates the magnetoresistance as per: __see reference__:
     * Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
     * @param Rx0
     * @param Ry0
     * @param AMR_X
     * @param AMR_Y
     * @param SMR_X
     * @param SMR_Y
     * @param AHE
     */
    std::vector<T> stripMagnetoResistance(const std::vector<T>& Rx0,
        const std::vector<T>& Ry0,
        const std::vector<T>& AMR_X,
        const std::vector<T>& SMR_X,
        const std::vector<T>& AMR_Y,
        const std::vector<T>& SMR_Y,
        const std::vector<T>& AHE)
    {
        T Rx_acc = 0.0;
        T Ry_acc = 0.0;

        for (unsigned int i = 0; i < this->layers.size(); i++)
        {
            const T Rx = Rx0[i] + AMR_X[i] * pow(this->layers[i].mag.x, 2) + SMR_X[i] * pow(this->layers[i].mag.y, 2);
            const T Ry = Ry0[i] + 0.5 * AHE[i] * this->layers[i].mag.z +
                (AMR_Y[i] + SMR_Y[i]) * this->layers[i].mag.x * this->layers[i].mag.y;
            Rx_acc += 1. / Rx;
            Ry_acc += 1. / Ry;
        }

        return { 1 / Rx_acc, 1 / Ry_acc, 0. };
    }

    /**
     * Calculate classic magnetoresistance.
     * Only for bilayer structures.
     * used when MR_MODE == CLASSIC
     * @param cosTheta: cosine between two layers.
     */
    T calculateMagnetoresistance(T cosTheta)
    {
        return this->Rp + (((this->Rap - this->Rp) / 2.0) * (1.0 - cosTheta));
    }

    std::vector<T> getMagnetoresistance()
    {
        // this is classical bilayer case
        if (this->MR_mode == CLASSIC && this->layerNo == 2)
        {
            return { calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[1].mag)) };
        }
        // this is the case when we use the pinning layer
        else if (this->MR_mode == CLASSIC && this->layerNo == 1)
        {
            return { calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[0].referenceLayer)) };
        }
        // this is strip magnetoresistance
        else if (this->MR_mode == STRIP)
        {
            return stripMagnetoResistance(this->Rx0,
                this->Ry0,
                this->AMR_X,
                this->SMR_X,
                this->AMR_Y,
                this->SMR_Y,
                this->AHE);
        }
        else
        {
            throw std::runtime_error("Magnetisation calculation is not supported for this structure!");
        }
    }

    /**
     * Main run simulation function. Use it to run the simulation.
     * @param totalTime: total time of a simulation, give it in seconds. Typical length is in ~couple ns.
     * @param timeStep: the integration step of the RK45 method. Default is 1e-13
     * @param writeFrequency: how often is the log saved to? Must be no smaller than `timeStep`. Default is 1e-11.
     * @param persist: whether to save to the filename specified in the Junction constructor. Default is true
     * @param log: if you want some verbosity like timing the simulation. Default is false
     * @param calculateEnergies: [WORK IN PROGRESS] log energy values to the log. Default is false.
     * @param mode: Solver mode EULER_HEUN, RK4 or DORMAND_PRICE
     */
    void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11,
        bool log = false, bool calculateEnergies = false,
        SolverMode mode = RK4)

    {
        if (timeStep > writeFrequency)
        {
            std::runtime_error("The time step cannot be larger than write frequency!");
        }
        const unsigned int totalIterations = (int)(totalTime / timeStep);
        const unsigned int writeEvery = (int)(writeFrequency / timeStep);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        auto solver = this->selectSolver(mode);
        // pick a solver based on drivers
        for (auto& l : this->layers)
        {
            if (l.hasTemperature())
            {
                if (mode != HEUN && mode != EULER_HEUN) {
                    std::cout << "[WARNING] Solver automatically changed to Heun for stochastic calculation." << std::endl;
                }
                // if at least one temp. driver is set
                // then use heun for consistency
                solver = this->selectSolver(HEUN);
                break;
            }
        }

        for (unsigned int i = 0; i < totalIterations; i++)
        {
            T t = i * timeStep;
            if (this->layerNo == 1)
            {
                runSingleLayerSolver(solver, t, timeStep);
            }
            else
            {
                runMultiLayerSolver(solver, t, timeStep);
            }

            if (!(i % writeEvery))
            {
                logLayerParams(t, timeStep, calculateEnergies);
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

#endif  // CORE_JUNCTION_HPP_
