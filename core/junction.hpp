#ifndef JUNCTION_H
#define JUNCTION_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <random>
#include <stdio.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <typeinfo>
#include <array>

#include "cvector.hpp"
#include "drivers.hpp"

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
inline CVector<T> calculate_tensor_interaction(CVector<T> &m,
                                               std::vector<CVector<T>> &tensor,
                                               T &Ms)
{
    CVector<T> res(
        tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
        tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
        tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
    return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> calculate_tensor_interaction(CVector<T> &m,
                                               std::array<CVector<T>, 3> &tensor,
                                               T &Ms)
{
    CVector<T> res(
        tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
        tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
        tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
    return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> c_cross(const CVector<T> &a, const CVector<T> &b)
{
    CVector<T> res(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);

    return res;
}

template <typename T>
inline T c_dot(CVector<T> &a, CVector<T> &b)
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

static std::default_random_engine generator;
template <typename T>
class Layer
{
private:
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

    bool temperatureSet = false;

    Reference referenceType = NONE;

    Layer(std::string id,
          CVector<T> mag,
          CVector<T> anis,
          T Ms,
          T thickness,
          T cellSurface,
          std::vector<CVector<T>> demagTensor,
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
        this->cellVolume = this->cellSurface * this->thickness;
        this->distribution = std::normal_distribution<T>(0, 1);
    }

public:
    bool includeSTT = false;
    bool includeSOT = false;

    std::string id;
    T Ms = 0.0;

    T thickness = 0.0;
    T cellVolume = 0.0, cellSurface = 0.0;

    CVector<T> H_log, Hoe_log, Hconst, mag, anisAxis, anis, referenceLayer;
    CVector<T> Hext, Hdipole, Hdemag, Hoe, HAnis, Hthermal, Hfl;

    CVector<T> HIEC, HIECtop, HIECbottom;
    T Jbottom_log = 0.0, Jtop_log = 0.0;
    T J2bottom_log = 0.0, J2top_log = 0.0;
    T K_log = 0.0;
    T I_log = 0.0;

    std::vector<CVector<T>> demagTensor;
    std::vector<CVector<T>> dipoleBottom = std::vector<CVector<T>>{CVector<T>(), CVector<T>(), CVector<T>()};
    std::vector<CVector<T>> dipoleTop = std::vector<CVector<T>>{CVector<T>(), CVector<T>(), CVector<T>()};

    T Hstart = 0.0, Hstop = 0.0, Hstep = 0.0;
    // LLG params
    T damping;

    // SOT params
    bool dynamicSOT = true;
    T fieldLikeTorque;
    T dampingLikeTorque;

    // STT params
    T SlonczewskiSpacerLayerParameter;
    T beta; // usually either set to 0 or to damping
    T spinPolarisation;

    std::normal_distribution<T> distribution;
    Layer() {}
    explicit Layer(std::string id,
                   CVector<T> mag,
                   CVector<T> anis,
                   T Ms,
                   T thickness,
                   T cellSurface,
                   std::vector<CVector<T>> demagTensor,
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
    explicit Layer(std::string id,
                   CVector<T> mag,
                   CVector<T> anis,
                   T Ms,
                   T thickness,
                   T cellSurface,
                   std::vector<CVector<T>> demagTensor,
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
    explicit Layer(std::string id,
                   CVector<T> mag,
                   CVector<T> anis,
                   T Ms,
                   T thickness,
                   T cellSurface,
                   std::vector<CVector<T>> demagTensor,
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

    inline static Layer<T> LayerSTT(std::string id,
                                    CVector<T> mag,
                                    CVector<T> anis,
                                    T Ms,
                                    T thickness,
                                    T cellSurface,
                                    std::vector<CVector<T>> demagTensor,
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

    inline static Layer<T> LayerSOT(std::string id,
                                    CVector<T> mag,
                                    CVector<T> anis,
                                    T Ms,
                                    T thickness,
                                    T cellSurface,
                                    std::vector<CVector<T>> demagTensor,
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

    void setTopDipoleTensor(std::vector<CVector<T>> dipoleTensor)
    {
        this->dipoleTop = dipoleTensor;
    }

    void setBottomDipoleTensor(std::vector<CVector<T>> dipoleTensor)
    {
        this->dipoleBottom = dipoleTensor;
    }

    const bool hasTemperature()
    {
        return this->temperatureSet;
    }

    void setTemperatureDriver(ScalarDriver<T> &driver)
    {
        this->temperatureDriver = driver;
        this->temperatureSet = true;
    }

    void setCurrentDriver(ScalarDriver<T> &driver)
    {
        this->currentDriver = driver;
    }

    void setFieldLikeTorqueDriver(ScalarDriver<T> &driver)
    {
        this->includeSOT = true;
        if (this->includeSTT)
            throw std::runtime_error("includeSTT was on and now setting SOT interaction!");
        if (!this->dynamicSOT)
            throw std::runtime_error("used a static SOT definition, now trying to set it dynamically!");
        this->fieldLikeTorqueDriver = driver;
    }

    void setDampingLikeTorqueDriver(ScalarDriver<T> &driver)
    {
        this->includeSOT = true;
        if (this->includeSTT)
            throw std::runtime_error("includeSTT was on and now setting SOT interaction!");
        if (!this->dynamicSOT)
            throw std::runtime_error("used a static SOT definition, now trying to set it dynamically!");
        this->dampingLikeTorqueDriver = driver;
    }

    void setAnisotropyDriver(ScalarDriver<T> &driver)
    {
        this->anisotropyDriver = driver;
    }

    void setExternalFieldDriver(AxialDriver<T> &driver)
    {
        this->externalFieldDriver = driver;
    }
    void setOerstedFieldDriver(AxialDriver<T> &driver)
    {
        this->HoeDriver = driver;
    }

    void setMagnetisation(CVector<T> &mag)
    {
        if (mag.length() == 0)
        {
            throw std::runtime_error("Initial magnetisation was set to a zero vector!");
        }
        this->mag = mag;
        this->mag.normalize();
    }

    void setIECDriverBottom(ScalarDriver<T> &driver)
    {
        this->IECDriverBottom = driver;
    }

    void setIECDriverTop(ScalarDriver<T> &driver)
    {
        this->IECDriverTop = driver;
    }

    void setQuadIECDriverTop(ScalarDriver<T> &driver)
    {
        this->IECQuadDriverTop = driver;
    }

    void setQuadIECDriverBottom(ScalarDriver<T> &driver)
    {
        this->IECQuadDriverBottom = driver;
    }

    /**
     * Set reference layer parameter. This is for calculating the spin current
     * polarisation if `includeSTT` is true.
     * @param reference: CVector describing the reference layer.
     */
    void setReferenceLayer(CVector<T> &reference)
    {
        this->referenceLayer = reference;
        this->referenceType = FIXED;
    }

    void setReferenceLayer(Reference reference)
    {
        if ((reference == FIXED) && (!this->referenceLayer.length()))
        {
            throw std::runtime_error("Cannot set fixed polarisation layer to 0! Set reference to NONE to disable reference.");
        }
        this->referenceType = reference;
    }

    CVector<T> getReferenceLayer()
    {
        // TODO: return other mags when the reference layer is not fixed.
        return this->referenceLayer;
    }

    Reference getReferenceType()
    {
        return this->referenceType;
    }

    const CVector<T> calculateHeff(T time, T timeStep, CVector<T> &stepMag, CVector<T> &bottom, CVector<T> &top)
    {
        this->Hdipole = calculate_tensor_interaction(bottom, this->dipoleBottom, this->Ms) +
                        calculate_tensor_interaction(top, this->dipoleTop, this->Ms);
        return calculateHeffDipoleInjection(time, timeStep, stepMag, bottom, top, this->Hdipole);
    }

    const CVector<T> calculateHeffDipoleInjection(T time, T timeStep, CVector<T> &stepMag, CVector<T> &bottom, CVector<T> &top, CVector<T> &dipole)
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
                                // demag -- negative contribution
                                - this->Hdemag
                                // dipole -- negative contribution
                                - dipole;

        return Heff;
    }

    CVector<T> calculateHOeField(T &time)
    {
        this->Hoe_log = this->HoeDriver.getCurrentAxialDrivers(time);
        return this->Hoe_log;
    }

    CVector<T> calculateExternalField(T &time)
    {
        this->H_log =
            this->externalFieldDriver.getCurrentAxialDrivers(time);
        return this->H_log;
    }

    CVector<T> calculateAnisotropy(CVector<T> &stepMag, T &time)
    {
        this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
        const T nom = (2 * this->K_log) * c_dot<T>(this->anis, stepMag) / (this->Ms);
        return this->anis * nom;
    }

    CVector<T> calculateIEC_(const T J, const T J2, CVector<T> stepMag, CVector<T> coupledMag)
    {
        const T nom = J / (this->Ms * this->thickness);
        // return (coupledMag - stepMag) * nom; // alternative form
        return (coupledMag + coupledMag * 2 * J2 * c_dot(coupledMag, stepMag)) * nom;
    }

    CVector<T> calculateIEC(T time, CVector<T> &stepMag, CVector<T> &bottom, CVector<T> &top)
    {
        this->Jbottom_log = this->IECDriverBottom.getCurrentScalarValue(time);
        this->Jtop_log = this->IECDriverTop.getCurrentScalarValue(time);

        this->J2bottom_log = this->IECQuadDriverBottom.getCurrentScalarValue(time);
        this->J2top_log = this->IECQuadDriverTop.getCurrentScalarValue(time);

        return calculateIEC_(this->Jbottom_log,
                             this->J2bottom_log, stepMag, bottom) +
               calculateIEC_(this->Jtop_log, this->J2top_log, stepMag, top);
    }

    const CVector<T> solveLLG(T time, CVector<T> m, T timeStep, CVector<T> bottom, CVector<T> top, CVector<T> heff)
    {
        const CVector<T> prod = c_cross<T>(m, heff);
        const CVector<T> prod2 = c_cross<T>(m, prod);
        const T convTerm = 1 / (1 + pow(this->damping, 2)); // LLGS -> LL form
        CVector<T> dmdt = prod + prod2 * this->damping;
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
            // this is more complex model
            // const T slonSq = pow(this->SlonczewskiSpacerLayerParameter, 2);
            // const T eta = (this->spinPolarisation * slonSq) / (slonSq + 1 + (slonSq - 1) * c_dot<T>(m, reference));
            // this is simplified
            const T eta = (this->spinPolarisation) / (1 + this->SlonczewskiSpacerLayerParameter * c_dot<T>(m, reference));
            const T sttTerm = GYRO * aJ * eta;
            const CVector<T> fieldLike = c_cross<T>(m, reference);
            // damping like
            const CVector<T> dampingLike = c_cross<T>(m, fieldLike);
            return (dmdt * -GYRO + dampingLike * -sttTerm + fieldLike * sttTerm * this->beta) * convTerm;
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
            const CVector<T> cm = c_cross<T>(m, reference);
            const CVector<T> ccm = c_cross<T>(m, cm);
            const CVector<T> flTorque = cm * (Hfl - this->damping * Hdl);
            const CVector<T> dlTorque = ccm * (Hdl + this->damping * Hfl);
            return (dmdt + flTorque + dlTorque) * -GYRO * convTerm;
        }
        return dmdt * -GYRO * convTerm;
    }
    const CVector<T> calculateLLGWithFieldTorqueDipoleInjection(T time, CVector<T> m, CVector<T> bottom, CVector<T> top, CVector<T> dipole, T timeStep)
    {
        // classic LLG first
        const CVector<T> heff = calculateHeffDipoleInjection(time, timeStep, m, bottom, top, dipole);
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
    const CVector<T> calculateLLGWithFieldTorque(T time, CVector<T> m, CVector<T> bottom, CVector<T> top, T timeStep)
    {
        // classic LLG first
        const CVector<T> heff = calculateHeff(time, timeStep, m, bottom, top);
        return solveLLG(time, m, timeStep, bottom, top, heff);
    }

    void rk4_step(T time, T timeStep, CVector<T> bottom, CVector<T> top)
    {
        CVector<T> m_t = this->mag;
        const CVector<T> k1 = calculateLLGWithFieldTorque(time, m_t, bottom, top, timeStep) * timeStep;
        const CVector<T> k2 = calculateLLGWithFieldTorque(time + 0.5 * timeStep, m_t + k1 * 0.5, bottom, top, timeStep) * timeStep;
        const CVector<T> k3 = calculateLLGWithFieldTorque(time + 0.5 * timeStep, m_t + k2 * 0.5, bottom, top, timeStep) * timeStep;
        const CVector<T> k4 = calculateLLGWithFieldTorque(time + timeStep, m_t + k3, bottom, top, timeStep) * timeStep;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        this->mag = m_t;
    }

    void rk4_stepDipoleInjection(T time, T timeStep, CVector<T> bottom, CVector<T> top, CVector<T> dipole)
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

    CVector<T> non_stochastic_llg(CVector<T> cm, T time, T timeStep, CVector<T> bottom, CVector<T> top)
    {
        return calculateLLGWithFieldTorque(time, cm, bottom, top, timeStep);
    }

    CVector<T> stochastic_llg(CVector<T> cm, T time, T timeStep, CVector<T> bottom, CVector<T> top, const CVector<T> &dW)
    {
        // compute the Langevin fluctuations -- this is the sigma
        const T Hthermal = this->getLangevinStochasticStandardDeviation(time);
        const CVector<T> thcross = c_cross(cm, dW);
        const CVector<T> thcross2 = c_cross(thcross, dW);
        const T scaling = -Hthermal * GYRO / (1 + pow(this->damping, 2));
        return (thcross + thcross2 * this->damping) * scaling;
    }

    T getLangevinStochasticStandardDeviation(T time)
    {
        if (this->cellVolume == 0.0)
            throw std::runtime_error("Cell surface cannot be 0 during temp. calculations!");
        const T currentTemp = this->temperatureDriver.getCurrentScalarValue(time);
        const T mainFactor = (2 * this->damping * BOLTZMANN_CONST * currentTemp) / (GYRO * this->Ms * this->cellVolume);
        return sqrt(mainFactor);
    }

    void euler_heun(T time, T timeStep, CVector<T> bottom, CVector<T> top)
    {
        // this is Stratonovich integral
        if (isnan(this->mag.x))
        {
            throw std::runtime_error("NAN magnetisation");
        }
        // Brownian motion sample
        // Generate the noise from the Brownian motion
        CVector<T> dW = CVector<T>(this->distribution, generator) * sqrt(timeStep);
        // squared dW -- just utility
        dW.normalize();
        // f_n is the vector of non-stochastic part at step n
        // multiply by timeStep (h) here for convenience
        const CVector<T> f_n = non_stochastic_llg(this->mag, time, timeStep, bottom, top) * timeStep;
        // g_n is the stochastic part of the LLG at step n
        const CVector<T> g_n = stochastic_llg(this->mag, time, timeStep, bottom, top, dW) * sqrt(timeStep);

        // actual solution
        // approximate next step ytilde
        const CVector<T> mapprox = this->mag + g_n;
        // calculate the approx g_n
        const CVector<T> g_n_approx = stochastic_llg(mapprox, time, timeStep, bottom, top, dW) * sqrt(timeStep);
        CVector<T> m_t = this->mag + f_n + g_n + (g_n_approx - g_n) * 0.5;
        m_t.normalize();
        this->mag = m_t;
    }
};

template <typename T>
class Junction
{
    friend class Layer<T>;
    const std::vector<std::string> vectorNames = {"x", "y", "z"};

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
    std::string fileSave;
    unsigned int logLength = 0;
    unsigned int layerNo;
    Junction() {}

    /**
     * Create a plain junction.
     * No magnetoresistance is calculated.
     */
    Junction(std::vector<Layer<T>> layersToSet, std::string filename = "")
    {
        this->MR_mode = NONE;
        this->layers = layersToSet;
        this->layerNo = this->layers.size();
        if (this->layerNo == 0)
        {
            throw std::invalid_argument("Passed a zero length Layer vector!");
        }
        this->fileSave = std::move(filename);
    }
    explicit Junction(std::vector<Layer<T>> layersToSet, std::string filename, T Rp, T Rap) : Junction(
                                                                                                  layersToSet, filename)
    {
        if (this->layerNo == 1)
        {
            // we need to check if this layer has a reference layer.
            if (!this->layers[0].referenceLayer.length())
            {
                throw std::invalid_argument("MTJ with a single layer must have a pinning (referenceLayer) set!");
            }
        }
        if (this->layerNo > 2)
        {
            throw std::invalid_argument("This constructor supports only bilayers! Choose the other one with the strip resistance!");
        }
        this->Rp = Rp;
        this->Rap = Rap;
        this->MR_mode = CLASSIC;
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
    explicit Junction(std::vector<Layer<T>> layersToSet,
                      std::string filename,
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
        this->fileSave = std::move(filename);
        this->MR_mode = STRIP;
    }

    /**
     * Clears the simulation log.
     **/
    void clearLog()
    {
        this->log.clear();
        this->logLength = 0;
    }

    std::unordered_map<std::string, std::vector<T>> &getLog()
    {
        return this->log;
    }

    typedef void (Layer<T>::*scalarDriverSetter)(ScalarDriver<T> &driver);
    typedef void (Layer<T>::*axialDriverSetter)(AxialDriver<T> &driver);
    void scalarlayerSetter(std::string &layerID, scalarDriverSetter functor, ScalarDriver<T> driver)
    {
        bool found = false;
        for (auto &l : this->layers)
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
    void axiallayerSetter(std::string &layerID, axialDriverSetter functor, AxialDriver<T> driver)
    {
        bool found = false;
        for (auto &l : this->layers)
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
    void setLayerTemperatureDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setTemperatureDriver, driver);
    }
    void setLayerAnisotropyDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setAnisotropyDriver, driver);
    }
    void setLayerExternalFieldDriver(std::string layerID, AxialDriver<T> driver)
    {
        axiallayerSetter(layerID, &Layer<T>::setExternalFieldDriver, driver);
    }
    void setLayerOerstedFieldDriver(std::string layerID, AxialDriver<T> driver)
    {
        axiallayerSetter(layerID, &Layer<T>::setOerstedFieldDriver, driver);
    }
    void setLayerCurrentDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setCurrentDriver, driver);
    }
    void setLayerDampingLikeTorqueDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setDampingLikeTorqueDriver, driver);
    }
    void setLayerFieldLikeTorqueDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setFieldLikeTorqueDriver, driver);
    }

    /**
     * Set IEC interaction between two layers.
     * The names of the params are only for convention. The IEC will be set
     * between bottomLayer or topLayer, order is irrelevant.
     * @param bottomLayer: the first layer id
     * @param topLayer: the second layer id
     */
    void setIECDriver(std::string bottomLayer, std::string topLayer, ScalarDriver<T> driver)
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

    void setQuadIECDriver(std::string bottomLayer, std::string topLayer, ScalarDriver<T> driver)
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

    void setLayerMagnetisation(std::string layerID, CVector<T> &mag)
    {
        bool found = false;
        for (auto &l : this->layers)
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

    CVector<T> getLayerMagnetisation(std::string layerID)
    {
        return getLayer(layerID).mag;
    }

    Reference getLayerReferenceType(std::string layerID)
    {
        return getLayer(layerID).referenceType;
    }

    void setLayerReferenceLayer(std::string layerID, CVector<T> referenceLayer)
    {
        if (layerID == "all")
        {
            for (auto &l : this->layers)
            {
                l.setReferenceLayer(referenceLayer);
            }
        }
        else
            getLayer(layerID).setReferenceLayer(referenceLayer);
    }

    void setLayerReferenceType(std::string layerID, Reference referenceType)
    {
        if (layerID == "all")
        {
            for (auto &l : this->layers)
            {
                l.setReferenceLayer(referenceType);
            }
        }
        else
            getLayer(layerID).setReferenceLayer(referenceType);
    }

    Layer<T> &getLayer(std::string layerID)
    {
        for (auto &l : this->layers)
        {
            if (l.id == layerID)
            {
                return l;
            }
        }
        throw std::runtime_error("Failed to find a layer with a given id!");
    }

    void logLayerParams(T &t, bool calculateEnergies = false)
    {
        for (const auto &layer : this->layers)
        {
            const std::string lId = layer.id;
            this->log[lId + "_K"].emplace_back(layer.K_log);
            for (int i = 0; i < 3; i++)
            {
                this->log[lId + "_m" + vectorNames[i]].emplace_back(layer.mag[i]);
                this->log[lId + "_Hext" + vectorNames[i]].emplace_back(layer.H_log[i]);
            }

            if (layer.includeSTT)
                this->log[lId + "_I"].emplace_back(layer.I_log);

            if (calculateEnergies)
            {
                this->log[lId + "_EZeeman"].emplace_back(EnergyDriver<T>::calculateZeemanEnergy(layer.mag,
                                                                                                layer.Hext,
                                                                                                layer.cellVolume,
                                                                                                layer.Ms));
                this->log[lId + "_EAnis"].push_back(EnergyDriver<T>::calculateAnisotropyEnergy(layer.mag,
                                                                                               layer.anis,
                                                                                               layer.K_log,
                                                                                               layer.cellVolume));
                this->log[lId + "_EDemag"].emplace_back(EnergyDriver<T>::calculateDemagEnergy(layer.mag,
                                                                                              layer.Hdemag,
                                                                                              layer.Ms,
                                                                                              layer.cellVolume));
                this->log[lId + "_EDipole"].emplace_back(EnergyDriver<T>::calculateDemagEnergy(layer.mag,
                                                                                               layer.Hdipole,
                                                                                               layer.Ms,
                                                                                               layer.cellVolume));
            }
        }
        if (MR_mode == CLASSIC)
        {
            const auto magnetoresistance = calculateMagnetoresistance(c_dot<T>(this->layers[0].mag,
                                                                               this->layers[1].mag));
            this->log["R_free_bottom"].emplace_back(magnetoresistance);
        }
        else if (this->MR_mode == CLASSIC && this->layerNo == 1)
        {
            this->log["R"].emplace_back(calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[0].referenceLayer)));
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
    saveLogs()
    {
        if (this->fileSave == "")
        {
            // if there's an empty fn, don't save
            std::cout << "Ignoring file save to an empty filename" << std::endl;
            return;
        }
        std::ofstream logFile;
        logFile.open(this->fileSave);
        for (const auto &keyPair : this->log)
        {
            logFile << keyPair.first << ";";
        }
        logFile << "\n";
        for (unsigned int i = 0; i < logLength; i++)
        {
            for (const auto &keyPair : this->log)
            {
                logFile << keyPair.second[i] << ";";
            }
            logFile << "\n";
        }
        logFile.close();
    }

    typedef void (Layer<T>::*solver)(T t, T timeStep, CVector<T> bottom, CVector<T> top);

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
    void runSingleLayerSolver(solver &functor, T &t, T &timeStep)
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
    void runMultiLayerSolver(solver &functor, T &t, T &timeStep)
    {
        std::vector<CVector<T>> magCopies(this->layerNo + 2, CVector<T>());
        // magCopies.resize(this->layerNo + 2)
        // the first and the last layer get 0 vector coupled
        // magCopies[0] = CVector<T>();
        // magCopies.at(this->layerNo + 1) = CVector<T>();
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
    std::vector<T> stripMagnetoResistance(std::vector<T> &Rx0,
                                          std::vector<T> &Ry0,
                                          std::vector<T> &AMR_X,
                                          std::vector<T> &SMR_X,
                                          std::vector<T> &AMR_Y,
                                          std::vector<T> &SMR_Y,
                                          std::vector<T> &AHE)
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

        return {1 / Rx_acc, 1 / Ry_acc, 0.};
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
            return {calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[1].mag))};
        }
        // this is the case when we use the pinning layer
        else if (this->MR_mode == CLASSIC && this->layerNo == 1)
        {
            return {calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[0].referenceLayer))};
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
     */
    void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11,
                       bool persist = true, bool log = false, bool calculateEnergies = false)

    {
        if (timeStep > writeFrequency)
        {
            std::runtime_error("The time step cannot be larger than write frequency!");
        }
        const unsigned int totalIterations = (int)(totalTime / timeStep);
        T t;
        const unsigned int writeEvery = (int)(writeFrequency / timeStep);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        // pick a solver based on drivers
        auto solv = &Layer<T>::rk4_step;
        for (auto &l : this->layers)
        {
            if (l.hasTemperature())
            {
                // if at least one temp. driver is set
                // then use euler_heun for consistency
                solv = &Layer<T>::euler_heun;
                break;
            }
        }

        for (unsigned int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;
            if (this->layerNo == 1)
            {
                runSingleLayerSolver(solv, t, timeStep);
            }
            else
            {
                runMultiLayerSolver(solv, t, timeStep);
            }

            if (!(i % writeEvery))
            {
                logLayerParams(t, calculateEnergies);
            }
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if (persist)
            saveLogs();
        if (log)
        {
            std::cout << "Steps in simulation: " << totalIterations << std::endl;
            std::cout << "Write every: " << writeEvery << std::endl;
            std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
        }
    }
};

#endif
