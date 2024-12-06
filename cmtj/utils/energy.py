import numpy as np


class EnergyCompute:
    """Energy density in [J/m^3] computing functions"""

    def __init__(self, cell_surface: float, thickness: float, log: dict[str, list[float]]) -> None:
        """Initialise energy computation class
        :param cell_surface: surface of the cell in [m^2]
        :param thickness: thickness of the cell in [m]
        :param log: log of the simulation (directly from .getLog())"""
        self.cell_surface = cell_surface
        self.thickness = thickness
        self.cell_volumne = self.cell_surface * thickness
        self.log = log

    def compute_from_log(self) -> dict[str, list[float]]:
        """
        Computes a log of energies over time and returns it
        in the same form of the
        """
        field_keys = list({k[:-1] for k in self.log if "_H" in k})
        mag_k = (k.replace("_mx", "") for k in self.log if "_mx" in k)
        mag_vectors = {k: np.asarray([self.log[f"{k}_mx"], self.log[f"{k}_my"], self.log[f"{k}_mz"]]) for k in mag_k}
        energy_data = {}
        for field_key in field_keys:
            if "J_" in field_key:
                eng_fn = self.calculate_energy_from_field_interfacial
            else:
                eng_fn = self.calculate_energy_from_field

            m_key = field_key.split("_")[0]  # get m key
            m = mag_vectors[m_key]
            field_series = np.asarray(
                [
                    self.log[f"{field_key}x"],
                    self.log[f"{field_key}y"],
                    self.log[f"{field_key}z"],
                ]
            )
            energy_data[f"energy_{field_key}"] = eng_fn(m, field_series)

        return energy_data

    def calculate_energy_from_field(self, m: np.ndarray, field_vector: np.ndarray) -> np.ndarray:
        """
        :param m: magnetisation
        :param field_vector: magnetic field vector (can be external, Oersted etc.)
        Compute generic energy density
        E = H * (mi0 Ms/V)
        where mi0 Ms is in [T], Ms in [A/m], H in [A/m]
        """
        return -np.sum(m * field_vector, axis=0) / self.cell_volumne

    def calculate_energy_from_field_interfacial(self, m: np.ndarray, field_vector: np.ndarray) -> np.ndarray:
        """
        :param m: magnetisation
        :param field_vector: magnetic field vector (can be IEC etc.)
        Compute generic energy density
        E = H * (mi0 Ms/A)
        where mi0 Ms is in [T], Ms in [A/m], H in [A/m]
        """
        return -np.sum(m * field_vector, axis=0) / self.cell_surface
