from itertools import product
from typing import Callable, List

import numpy as np
from multiprocess import Pool
from tqdm import tqdm

__all__ = ["distribute", "create_coordinates_plot"]


def distribute(simulation_fn: Callable,
               spaces: List[List[float]],
               n_cores: int = None):
    """
    Distribute a function over a list of parameters in parallel.
    :param simulation_fn: function to be distributed
    :param spaces: list of lists of parameters
    :param n_cores: number of cores to use.
    :returns: index, simulation_fn output
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
