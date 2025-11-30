#ifndef CORE_CONSTANTS_HPP_
#define CORE_CONSTANTS_HPP_

#include <cmath>

class PhysicalConstants {
public:
  // Default values (matching current macros)
  static double MAGNETIC_PERMEABILITY;
  static double GYRO;
  static double TtoAm;
  static double HBAR;
  static double ELECTRON_CHARGE;
  static double BOLTZMANN_CONST;

  // Setters for Python interface
  static void setMagneticPermeability(double value) {
    MAGNETIC_PERMEABILITY = value;
  }
  static void setGyro(double value) { GYRO = value; }
  static void setTtoAm(double value) { TtoAm = value; }
  static void setHbar(double value) { HBAR = value; }
  static void setElectronCharge(double value) { ELECTRON_CHARGE = value; }
  static void setBoltzmannConst(double value) { BOLTZMANN_CONST = value; }

  // Getters
  static double getMagneticPermeability() { return MAGNETIC_PERMEABILITY; }
  static double getGyro() { return GYRO; }
  static double getTtoAm() { return TtoAm; }
  static double getHbar() { return HBAR; }
  static double getElectronCharge() { return ELECTRON_CHARGE; }
  static double getBoltzmannConst() { return BOLTZMANN_CONST; }

  // Reset to defaults
  static void resetToDefaults();
};

// Initialize static members with default values
double PhysicalConstants::MAGNETIC_PERMEABILITY = 12.566e-7;
double PhysicalConstants::GYRO = 220880.0; // rad/Ts converted to m/As
double PhysicalConstants::TtoAm = 795774.715459;
double PhysicalConstants::HBAR = 6.62607015e-34 / (2. * M_PI);
double PhysicalConstants::ELECTRON_CHARGE = 1.60217662e-19;
double PhysicalConstants::BOLTZMANN_CONST = 1.380649e-23;

void PhysicalConstants::resetToDefaults() {
  MAGNETIC_PERMEABILITY = 12.566e-7;
  GYRO = 220880.0;
  TtoAm = 795774.715459;
  HBAR = 6.62607015e-34 / (2. * M_PI);
  ELECTRON_CHARGE = 1.60217662e-19;
  BOLTZMANN_CONST = 1.380649e-23;
}

#endif // CORE_CONSTANTS_HPP_
