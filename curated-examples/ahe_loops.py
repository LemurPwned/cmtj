import contextlib
import itertools as it
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from tqdm import tqdm

import cmtj.utils.plotting as plotting
from cmtj import AxialDriver, CVector, Junction, Layer, constantDriver
from cmtj.utils import FieldScan, OetoAm, calculate_resistance_series

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401


def AHE_simulate(
    thickness: float,
    J: float,
    Hscan: np.ndarray,
    Ms: float,
    Ku1: float,
    Ku2: float,
    Ku1_theta: float,
    Ku1_phi: float,
    Ku2_theta: float,
    Ku2_phi: float,
    Hdl: float,
    Hfl: float,
    asymmetry_DL: float = 1,
    asymmetry_FL: float = 1,
    AHE_asymmetry: float = 1,
    tstep: float = 1e-13,
    AHE_scaling: float = 1,
    alpha=5e-2,
):
    rho_f = 27e-8
    rho_h = 21e-8
    t_fm = 1e-9
    t_hm = thickness  # we only change the thickness of the HM
    w = 10e-6
    l = 80e-6
    area = w * l
    FM_R = rho_f * t_fm / area
    HM_R = rho_h * l / (w * t_hm)
    # parallel, current in plane
    T_R = 1.0 / FM_R + 1.0 / HM_R
    T_R = 1.0 / T_R
    demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1)]
    """
    Compute the current scan.
    The facing area is thickness x width
    """
    jden = 4.24e10
    HDL = Hdl / jden
    HFL = Hfl / jden

    Kdir1 = FieldScan.angle2vector(Ku1_theta, Ku1_phi)
    Kdir2 = FieldScan.angle2vector(Ku2_theta, Ku2_phi)
    pdir1 = CVector(0, 1.0, 0.0)
    pdir2 = CVector(0, 1.0, 0.0)
    layer_free = Layer.createSOTLayer(
        id="free",
        mag=Kdir1,
        anis=Kdir1,
        Ms=Ms,
        thickness=t_fm,
        cellSurface=area,
        demagTensor=demagTensor,
        damping=alpha,
        dampingLikeTorque=HDL,
        fieldLikeTorque=HFL,
    )

    layer_bottom = Layer.createSOTLayer(
        id="bottom",
        mag=Kdir2,
        anis=Kdir2,
        Ms=Ms,
        thickness=t_fm,
        cellSurface=area,
        demagTensor=demagTensor,
        damping=alpha,
        dampingLikeTorque=HDL * asymmetry_DL,
        fieldLikeTorque=HFL * asymmetry_FL,
    )

    j = Junction([layer_free, layer_bottom])
    j.setLayerAnisotropyDriver("free", constantDriver(Ku1))
    j.setLayerAnisotropyDriver("bottom", constantDriver(Ku2))
    j.setLayerReferenceLayer("free", pdir1)
    j.setLayerReferenceLayer("bottom", pdir2)
    j.setIECDriver("free", "bottom", constantDriver(J))
    """
    Run the simulation for range of H
    """
    j.runSimulation(4e-9, tstep, tstep)
    output = defaultdict(list)
    n_lay = 2

    for Hscan_ in (Hscan, Hscan[::-1]):
        for Hval in Hscan_:
            j.clearLog()
            j.setLayerExternalFieldDriver(
                "all",
                AxialDriver(
                    0,
                    0,
                    Hval,
                ),
            )

            j.runSimulation(15e-9, tstep, tstep, calculateEnergies=False)

            log = j.getLog()
            m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ("free", "bottom")])
            Rx0 = [100] * n_lay
            Ry0 = [0] * n_lay
            SMR = [5.11] * n_lay
            AMR = [0.41] * n_lay
            AHE = [AHE_scaling, AHE_scaling * AHE_asymmetry]

            _, Rxy = calculate_resistance_series(Rx0, Ry0, AMR=AMR, AHE=AHE, SMR=SMR, m=m, l=[l] * n_lay, w=[w] * n_lay)
            Rstable = Rxy[-100:].mean()
            output["Hscan"].append(Hval)
            output["R"].append(Rstable)

            theta_top, _, _ = FieldScan.vector2angle(log["free_mx"][-1], log["free_my"][-1], log["free_mz"][-1])
            theta_bottom, _, _ = FieldScan.vector2angle(
                log["bottom_mx"][-1], log["bottom_my"][-1], log["bottom_mz"][-1]
            )
            output["theta_top"].append(theta_top)
            output["theta_bottom"].append(theta_bottom)

    return output


big_results = defaultdict(list)
thfm = 1e-9
thicknesses = np.asarray([0.7, 1.88, 2.59, 3.48]) * 1e-9
K1_thetas = [1, 0.01, 0.1, 0.1]
AHE_asyms = [1.0, 1.0, 0.8, 0.5]
Ku1s = np.asarray([0.55, 1.0, 0.625, 0.625]) * 1e6
Ku2s = np.asarray([0.3, 0.27, 0.29, 0.27]) * 1e6
Ku2_thetas = [0.01, 0.01, 0.5, 0.5]
etas = [0.45, 3.5, 0.45, 10]
Js = np.asarray([8, 6, 3.5, 0.3]) * 1e-3
scales = [thfm / (th + thfm) for th in thicknesses]
Ms = [1.1, 1.18, 1.05, 1.15]
J = [j * scale for j, scale in zip(Js, scales)]
scanned_str = [r"$t_\mathrm{HM} (nm)$: " + f"{th:.2f}" for th in thicknesses / 1e-9]

HFL = np.asarray(list(it.repeat(90, len(thicknesses)))) * OetoAm
HDL = HFL * np.asarray(etas)

Hscan, _ = FieldScan.amplitude_scan(-7000 * OetoAm, 7000 * OetoAm, 30, 1, 0)

total_samples = len(thicknesses)
AHE_ouputs = []
for sample in tqdm(range(total_samples)):
    # if sample in (0, 2):
    # continue
    output = AHE_simulate(
        Hscan=Hscan,
        J=J[sample],
        thickness=thicknesses[sample],
        Ms=Ms[sample],
        Ku1=Ku1s[sample],
        Ku2=Ku2s[sample],
        Ku1_phi=0,
        Ku2_phi=0,
        Ku1_theta=K1_thetas[sample],
        Ku2_theta=Ku2_thetas[sample],
        Hdl=HDL[sample],
        Hfl=HFL[sample],
        AHE_scaling=5,
        AHE_asymmetry=AHE_asyms[sample],
        alpha=5e-2,
        tstep=1e-13,
    )
    AHE_ouputs.append(output)
    # break


with plt.style.context(["nature"]):
    fig, axs = plt.subplots(4, 2, figsize=(4.5, 3.5), sharex=True, dpi=400)
    regions = ("I", "II", "III", "IV")
    pt_th = list(reversed([1, 0.7, 0.5, 0.3]))
    angles_top = [0, 90, 90, 90]
    angles_bottom = [0, 90, 80, 0]
    s = 2
    lw = 1.2
    Hlim = 550
    arpops = "->, head_width=.08, head_length=.2"
    for i, (region, output, label) in enumerate(zip(regions, AHE_ouputs, thicknesses)):
        x_axis = np.asarray(output["Hscan"]) / 1e3
        rlabel = np.around(label * 1e9, 2)
        thtrue = f"{label * 1e9:.2f} nm"
        axs[i, 0].text(0.73, 0.7, thtrue, transform=axs[i, 0].transAxes)
        axs[i, 0].text(
            0.73,
            0.85,
            f"region {region}",
            transform=axs[i, 0].transAxes,
            color="royalblue",
        )

        axins2 = axs[i, 0].inset_axes([0.65, 0.17, 0.33, 0.52])
        plotting.create_stack(
            axins2,
            labels=["MgO", "Co", "Pt", "Co", "Ti"],
            colors=[
                "mediumturquoise",
                "lightgray",
                "lightblue",
                "lightgray",
                "peachpuff",
            ],
            angles=[None, angles_top[i], None, angles_bottom[i], None],
            heights=[0.3, 1.3, pt_th[i], 1.3, 0.5],
            text_fontsize=3.5,
            width=2.3,
            lw_arrow=0.1,
            ms=4,
            r=2,
            labelpad_left=0.49,
        )
        axs[i, 0].plot(x_axis, -1 * np.asarray(output["R"]), "-", lw=lw, color="k")
        axs[i, 0].axis("off")
        for _ax in (axs[i, 0], axs[i, 0]):
            _ax.set_ylabel(r"$\mathrm{R}_\mathrm{XY}$ ($\Omega$)")
        axs[i, 0].set_xlim([-Hlim, Hlim])
        if i != len(thicknesses) - 1:
            axs[i, 0].set_xlim([-Hlim, Hlim])

        ttop = np.asarray(output["theta_top"])
        tbottom = np.asarray(output["theta_bottom"])

        axs[i, 1].scatter(
            x_axis,
            ttop,
            s=s,
            label="top",
            color="forestgreen",
        )
        axs[i, 1].scatter(
            x_axis,
            tbottom,
            s=s,
            label="bottom",
            color="crimson",
        )
        axs[i, 1].set_xlim([-Hlim, Hlim])
        axs[i, 1].set_ylabel(r"$\theta$ (deg.)", labelpad=10)
        axs[i, 1].yaxis.tick_right()
        axs[i, 1].yaxis.set_label_position("right")
        # rotate the label
        axs[i, 1].yaxis.label.set_rotation(270)

        for j in range(len(x_axis) - 1):
            axs[i, 1].annotate(
                "",
                xy=(x_axis[j + 1], ttop[j + 1]),
                xytext=(x_axis[j], ttop[j]),
                arrowprops=dict(arrowstyle=arpops, color="forestgreen", lw=0.5),
            )
            axs[i, 1].annotate(
                "",
                xy=(x_axis[j + 1], tbottom[j + 1]),
                xytext=(x_axis[j], tbottom[j]),
                arrowprops=dict(arrowstyle=arpops, color="crimson", lw=0.5),
            )
    axs[-1, 0].set_xlabel(r"$\mathrm{H}_z$ ($\mathrm{kA/m}$)")
    axs[-1, 1].set_xlabel(r"$\mathrm{H}_z$ ($\mathrm{kA/m}$)")
    # move ylabels to the right

    axs[0, 1].legend()
    fig.align_ylabels(axs.flatten())
    fig.subplots_adjust(hspace=0, wspace=0.05)
    import matplotlib.transforms as mtransforms

    for label, ax in zip(("abcd"), axs[:, 0].flatten()):
        # label physical distance in and down:
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(
            0.0,
            1.0,
            f"({label})",
            transform=ax.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            bbox=dict(facecolor="none", edgecolor="none", pad=3.0),
        )
    for label, ax in zip(("efgh"), axs[:, 1].flatten()):
        # label physical distance in and down:
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(
            0.0,
            1.0,
            f"({label})",
            transform=ax.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            bbox=dict(facecolor="none", edgecolor="none", pad=3.0),
        )
    fig.savefig(
        "./curated-examples/figures/ahe-loops.png",
        dpi=350,
        bbox_inches="tight",
    )
