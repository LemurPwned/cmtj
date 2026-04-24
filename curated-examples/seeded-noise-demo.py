"""
Seeded thermal and 1/f noise demo.

This example runs a tiny stochastic LLG simulation twice with the same layer seed
and once with a different seed. The resulting chart compares identical-seed
repetitions against a different-seed run so the reproducibility is easy to see.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj import CVector, Junction, Layer, SolverMode, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # type: ignore[import-not-found]  # noqa: F401


def build_junction(seed: int | None = None) -> Junction:
    demag = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1.0)]
    layer = Layer(
        "free",
        CVector(1.0, 0.0, 0.0),
        CVector(0.0, 0.0, 1.0),
        1.0,
        1.0e-9,
        60e-9 * 60e-9,
        demag,
        damping=0.02,
    )
    junction = Junction([layer])
    if seed is not None:
        junction.setLayerSeed("free", seed)
    junction.setLayerTemperatureDriver("free", constantDriver(100.0))
    return junction


def run_simulation(seed: int | None = None):
    junction = build_junction(seed)
    junction.runSimulation(5e-10, 1e-12, 1e-12, solverMode=SolverMode.EulerHeun)
    print(f'Seed: {seed}')
    mx = junction.getLayerMagnetisation("free").x
    my = junction.getLayerMagnetisation("free").y
    mz = junction.getLayerMagnetisation("free").z
    print(f'Final m: {mx:.4f}, {my:.4f}, {mz:.4f}')
    return junction.getLog()


def plot_seed_comparison():
    same_seed_1 = run_simulation(1224)
    same_seed_2 = run_simulation(1224)
    different_seed = run_simulation(9876)

    time_ns = np.asarray(same_seed_1["time"]) * 1e9
    component = "free_mz"
    style = ["science", "no-latex"] if "science" in plt.style.available else ["default"]

    with plt.style.context(style):
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)

        axes[0].plot(time_ns, same_seed_1[component], color="navy", linewidth=2.0, label="Run 1")
        axes[0].plot(time_ns, same_seed_2[component], color="crimson", linewidth=1.2, linestyle="--", label="Run 2")
        axes[0].set_title("Same seed")
        axes[0].set_xlabel("t (ns)")
        axes[0].set_ylabel(r"$m_z$")
        axes[0].legend(frameon=False)

        axes[1].plot(time_ns, same_seed_1[component], color="navy", linewidth=2.0, label="Seed 1224")
        axes[1].plot(time_ns, different_seed[component], color="darkgreen", linewidth=1.6, label="Seed 9876")
        axes[1].set_title("Different seeds")
        axes[1].set_xlabel("t (ns)")
        axes[1].legend(frameon=False)

        for axis in axes:
            axis.grid(alpha=0.2)

        fig.suptitle("Seeded thermal noise reproducibility")
        fig.tight_layout()
        fig.savefig("./curated-examples/figures/seeded-noise-demo.png", dpi=300, bbox_inches="tight")


plot_seed_comparison()