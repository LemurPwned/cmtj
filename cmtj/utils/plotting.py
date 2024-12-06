from itertools import permutations

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def get_sphere():
    r = 1
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0.0:pi:100j, 0.0 : 2.0 * pi : 100j]
    xs = r * sin(phi) * cos(theta)
    ys = r * sin(phi) * sin(theta)
    zs = r * cos(phi)
    return xs, ys, zs


def plot_trajectory_sphere(x, y, z, color="blue", alpha=1, ax=None):
    """Plot a trajectory in 3D. Normalises to unit sphere
    :param ax: matplotlib axis
    :param x: x-coordinates
    :param y: y-coordinates
    :param z: z-coordinates
    :param color: color of the trajectory
    :param alpha: alpha value of the trajectory
    """
    # Compute a unit sphere first
    xs, ys, zs = get_sphere()
    m = np.asarray([x, y, z])
    # make sure we are unit norm for m
    m = m / np.linalg.norm(m)
    if ax is None:
        with plt.style.context(["science", "nature"]):
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(1, 1, 1, projection="3d")
            ax.plot3D(m[0], m[1], m[2], color=color, alpha=alpha)
            ax.set_axis_off()
            ax.plot_surface(
                xs,
                ys,
                zs,
                rstride=2,
                cstride=2,
                color="azure",
                alpha=0.1,
                linewidth=0.1,
            )
            ax.scatter([0], [0], [1], color="crimson", alpha=1.0)
    else:
        ax.plot3D(m[0], m[1], m[2], color=color, alpha=alpha)
        ax.set_axis_off()
        ax.plot_surface(xs, ys, zs, rstride=2, cstride=2, color="azure", alpha=0.1, linewidth=0.1)
        ax.scatter([0], [0], [1], color="crimson", alpha=1.0)


def plot_coloured_trajectory(x, y, z, colormap="plasma", ax=None):
    """Plot a coloured trajectory in 3D. Normalises to unit sphere.
    Colour of the trajectory now designates the flow of time.
    :param ax: matplotlib axis
    :param x: x-coordinates
    :param y: y-coordinates
    :param z: z-coordinates
    :param colormap: colormap to use
    """
    import seaborn as sns

    xs, ys, zs = get_sphere()
    m = np.asarray([x, y, z])
    points = m.T.reshape(-1, 1, 3)
    segs = np.concatenate([points[:-1], points[1:]], axis=1)

    colors = sns.color_palette(colormap, len(segs))
    if ax is None:
        with plt.style.context(["science", "nature"]):
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(1, 1, 1, projection="3d")
            # plot the sphere firext
            ax.set_axis_off()
            ax.plot_surface(
                xs,
                ys,
                zs,
                rstride=2,
                cstride=2,
                color="azure",
                alpha=0.1,
                linewidth=0.1,
            )
            ax.add_collection(Line3DCollection(segs, colors=colors, alpha=1))
    else:
        ax.set_axis_off()
        ax.plot_surface(xs, ys, zs, rstride=2, cstride=2, color="azure", alpha=0.1, linewidth=0.1)
        ax.add_collection(Line3DCollection(segs, colors=colors, alpha=1))


def unpack_ndim_map(map, axes):
    """
    Unpack N-dimensional map into a list of 1-dimensional arrays
    :param map: N-dimensional map, each axis is separate parameter space.
    :param axes: list of axes to unpack.
    """
    # how long each one is
    sample_length = len(axes[0])
    perm_indx = permutations(range(sample_length), len(axes))

    ax_lists = [[] for _ in axes]
    value_list = []
    for indx in perm_indx:
        value_list.append(map[indx])
        for i, ax in enumerate(axes):
            ax_lists[i].append(ax[indx[i]])

    return ax_lists, value_list


def create_coordinates_plot(axes, ax_names, result_map, sample=0, alpha_black=0.01):
    """Create parallel coordinates plot for multidimensional parameter space.
    Modified from:
    https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib
    :param axes: N list of parameters
    :param ax_names: N list of parameter names
    :param result_map: map of values (N-dim)
    :param sample: if != 0, subsample the parameter space
    :param alpha_black: alpha value zero value
    """
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.patches as patches
    from matplotlib.path import Path

    with plt.style.context(["science", "nature"]):
        fig, host = plt.subplots(dpi=400)
        ax_lists, value_list = unpack_ndim_map(result_map, axes)

        norm = matplotlib.colors.Normalize(vmin=min(value_list), vmax=max(value_list), clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.magma)

        # organize the data
        ys = np.dstack([*ax_lists, value_list])[0]
        indices = np.arange(len(ys))
        if sample:
            indices = np.random.choice(indices, sample).ravel()
        ys = ys[indices]

        ymins = ys.min(axis=0)
        ymaxs = ys.max(axis=0)
        dys = ymaxs - ymins
        ymins -= dys * 0.05  # add 5% padding below and above
        ymaxs += dys * 0.05
        dys = ymaxs - ymins

        # transform all data to be compatible with the main axis
        zs = np.zeros_like(ys)
        zs[:, 0] = ys[:, 0]
        zs[:, 1:] = (ys[:, 1:] - ymins[1:]) / dys[1:] * dys[0] + ymins[0]

        axes = [host] + [host.twinx() for _ in range(ys.shape[1] - 1)]
        for i, ax in enumerate(axes):
            ax.set_ylim(ymins[i], ymaxs[i])
            ax.spines["top"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            if ax != host:
                ax.spines["left"].set_visible(False)
                ax.yaxis.set_ticks_position("right")
                ax.spines["right"].set_position(("axes", i / (ys.shape[1] - 1)))

        host.set_xlim(0, ys.shape[1] - 1)
        host.set_xticks(range(ys.shape[1]))
        host.set_xticklabels(ax_names, fontsize=8)
        host.tick_params(axis="x", which="major", pad=7)
        host.spines["right"].set_visible(False)
        host.xaxis.tick_top()
        host.set_title("Parallel Coordinates Plot")

        for j in range(ys.shape[0]):
            # create bezier curves
            # for each axis, there will a control vertex at the point itself, one at 1/3rd towards the previous and one
            #   at one third towards the next axis; the first and last axis have one less control vertex
            # x-coordinate of the control vertices: at each integer (for the axes) and two inbetween
            # y-coordinate: repeat every point three times, except the first and last only twice
            verts = list(
                zip(
                    list(np.linspace(0, len(ys) - 1, len(ys) * 3 - 2, endpoint=True)),
                    np.repeat(zs[j, :], 3)[1:-1],
                )
            )
            # for x,y in verts: host.plot(x, y, 'go') # to show the control points of the beziers
            codes = [Path.MOVETO] + [Path.CURVE4 for _ in range(len(verts) - 1)]
            path = Path(verts, codes)
            alpha = alpha_black if ys[j, -1] == 0 else 0.8
            patch = patches.PathPatch(
                path,
                facecolor="none",
                lw=0.5,
                edgecolor=mapper.to_rgba(ys[j, -1], alpha),
            )
            host.add_patch(patch)
        fig.tight_layout()


def rotation_matrix(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])


def create_stack(
    ax,
    colors,
    heights,
    angles,
    labels,
    width=2,
    labelpad_left=0.2,
    offset_x=0,
    offset_y=0,
    lw_arrow=1.5,
    ms=10,
    r=0.6,
    text_fontsize=4,
    reversed=True,
):
    """
    Create a material stack plot.
    If a given layer is to have no arrow, pass None.
    :param ax: matplotlib axis
    :param colors: list of colors
    :param heights: list of heights
    :param angles: list of angles
    :param labels: list of labels
    :param width: width of the bars
    :param labelpad_left: padding of the labels
    :param offset_x: offset of the patches in x direction
    :param offset_y: offset of the patches in y direction
    :param lw_arrow: linewidth of the arrows
    :param ms: mutation size of the arrows
    :param r: length of the arrows
    :param reversed: if True, the stack is reversed
    """
    [x, y] = [r, 0]
    first_offset = offset_y
    if reversed:
        heights = heights[::-1]
        colors = colors[::-1]
        angles = angles[::-1]
        labels = labels[::-1]
    for _i, (height, angle, color, label) in enumerate(zip(heights, angles, colors, labels)):
        ax.add_patch(patches.Rectangle((offset_x, offset_y), width, height, fill=True, color=color, zorder=10))
        ax.text(
            offset_x - labelpad_left,
            offset_y + height / 2,
            label,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=text_fontsize,
            zorder=11,
        )
        if angle is not None:
            [dx, dy] = np.dot(rotation_matrix(np.deg2rad(angle)), [x, y])
            x_mid = dx / 2
            y_mid = dy / 2
            centre_x = (offset_x + width) / 2 - x_mid
            centre_y = offset_y + height / 2 - y_mid
            ax.add_patch(
                patches.FancyArrowPatch(
                    (centre_x, centre_y),
                    (centre_x + dx, centre_y + dy),
                    mutation_scale=ms,
                    lw=lw_arrow,
                    color="black",
                    zorder=10,
                )
            )
        offset_y += height
    ax.set_ylim([first_offset - max(heights) / 2, offset_y + max(heights) / 2])
    ax.set_xlim([offset_x - width / 2, offset_x + width + width / 2])
    ax.axis("off")
    return ax
