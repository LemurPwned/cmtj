from itertools import permutations

import matplotlib.pyplot as plt
import numpy as np


def plot_trajectory_sphere(ax, x, y, z, color='blue', alpha=1):
    """Plot a trajectory in 3D. Normalises to unit sphere
    :param ax: matplotlib axis
    :param x: x-coordinates
    :param y: y-coordinates
    :param z: z-coordinates
    :param color: color of the trajectory
    :param alpha: alpha value of the trajectory
    """
    # Compute a unit sphere first
    r = 1
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0 * pi:100j]
    xs = r * sin(phi) * cos(theta)
    ys = r * sin(phi) * sin(theta)
    zs = r * cos(phi)

    with plt.style.context(['science', 'no-latex']):
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        m = np.asarray([x, y, z])
        # make sure we are unit norm for m
        m = m / np.linalg.norm(m)
        ax.plot3D(m[0], m[1], m[2], color=color, alpha=alpha)
        ax.set_axis_off()
        ax.plot_surface(xs,
                        ys,
                        zs,
                        rstride=2,
                        cstride=2,
                        color='c',
                        alpha=0.3,
                        linewidth=0.1)
        ax.scatter([0], [0], [1], color='crimson', alpha=1.0)


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


def create_coordinates_plot(axes,
                            ax_names,
                            result_map,
                            sample=0,
                            alpha_black=0.01):
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
    with plt.style.context(['science', 'nature']):
        fig, host = plt.subplots(dpi=400)
        ax_lists, value_list = unpack_ndim_map(result_map, axes)

        norm = matplotlib.colors.Normalize(vmin=min(value_list),
                                           vmax=max(value_list),
                                           clip=True)
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

        axes = [host] + [host.twinx() for i in range(ys.shape[1] - 1)]
        for i, ax in enumerate(axes):
            ax.set_ylim(ymins[i], ymaxs[i])
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            if ax != host:
                ax.spines['left'].set_visible(False)
                ax.yaxis.set_ticks_position('right')
                ax.spines["right"].set_position(
                    ("axes", i / (ys.shape[1] - 1)))

        host.set_xlim(0, ys.shape[1] - 1)
        host.set_xticks(range(ys.shape[1]))
        host.set_xticklabels(ax_names, fontsize=8)
        host.tick_params(axis='x', which='major', pad=7)
        host.spines['right'].set_visible(False)
        host.xaxis.tick_top()
        host.set_title('Parallel Coordinates Plot')

        for j in range(ys.shape[0]):
            # create bezier curves
            # for each axis, there will a control vertex at the point itself, one at 1/3rd towards the previous and one
            #   at one third towards the next axis; the first and last axis have one less control vertex
            # x-coordinate of the control vertices: at each integer (for the axes) and two inbetween
            # y-coordinate: repeat every point three times, except the first and last only twice
            verts = list(
                zip([
                    x for x in np.linspace(
                        0, len(ys) - 1, len(ys) * 3 - 2, endpoint=True)
                ],
                    np.repeat(zs[j, :], 3)[1:-1]))
            # for x,y in verts: host.plot(x, y, 'go') # to show the control points of the beziers
            codes = [Path.MOVETO
                     ] + [Path.CURVE4 for _ in range(len(verts) - 1)]
            path = Path(verts, codes)
            if ys[j, -1] == 0:
                alpha = alpha_black
            else:
                alpha = 0.8
            patch = patches.PathPatch(path,
                                      facecolor='none',
                                      lw=.5,
                                      edgecolor=mapper.to_rgba(
                                          ys[j, -1], alpha))
            host.add_patch(patch)
        fig.tight_layout()
