from itertools import permutations, product

import numpy as np
from multiprocess import Pool
from tqdm import tqdm

__all__ = ["distribute", "create_coordinates_plot"]


def distribute(simulation_fn, spaces, n_cores=8):
    """
    Distribute a function over a list of parameters in parallel.
    :param simulation_fn: function to be distributed
    :param spaces: list of lists of parameters
    :param n_cores: number of cores to use.
    """
    spaces = [np.asarray(space) for space in spaces]

    def _get_index(values):
        return [
            np.argwhere(space == values[i]).ravel()[0]
            for i, space in enumerate(spaces)
        ]

    iterables = list(product(*spaces))
    indexes = [_get_index(val) for val in iterables]

    def func_wrapper(iterable):
        return iterable, simulation_fn(*iterable)

    with Pool(processes=n_cores) as pool:
        for result in tqdm(pool.imap_unordered(func_wrapper, iterables),
                           total=len(iterables)):
            iterable, output = result
            indx = indexes[iterables.index(iterable)]
            yield indx, output


def unpack_ndim_map(map, axes):
    """
    Unpack N-dimensional map into a list of 1-dimensional arrays
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
    """
    Modified from
    https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib
    :param: axes: N llist of parameters
    :param result_map: map of values (N-dim)
    """
    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.path import Path
    with plt.style.context(['science', 'nature']):
        fig, host = plt.subplots(dpi=400)

        # create some dummy data
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
