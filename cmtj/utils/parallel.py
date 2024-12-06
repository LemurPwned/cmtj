from itertools import product
from typing import Callable

import numpy as np
from multiprocess import Pool
from tqdm import tqdm

__all__ = ["distribute"]


def distribute(
    simulation_fn: Callable,
    spaces: list[list[float]],
    n_cores: int = None,
    shuffle: bool = False,
):
    """
    Distribute a function over a list of parameters in parallel.
    :param simulation_fn: function to be distributed
    :param spaces: list of lists of parameters
    :param n_cores: number of cores to use.
    :returns: tuple
        index (int): Index of the parameters in the spaces list, multiple dims.
        simulation_fn output (any): The output of the simulation function.
        index - index of the parameters in the spaces list, multiple dims.
    """
    spaces = [np.asarray(space) for space in spaces]

    def _get_index(values):
        return [np.argwhere(space == values[i]).ravel()[0] for i, space in enumerate(spaces)]

    iterables = list(product(*spaces))
    indexes = [_get_index(val) for val in iterables]
    # shuffle the indexes
    if shuffle:
        index_reshuffle = np.arange(len(indexes))
        np.random.shuffle(index_reshuffle)
        # reorder the indexes
        iterables = np.asarray(iterables)[index_reshuffle].tolist()
        indexes = np.asarray(indexes)[index_reshuffle].tolist()

    def func_wrapper(iterable):
        return iterable, simulation_fn(*iterable)

    with Pool(processes=n_cores) as pool:
        for result in tqdm(pool.imap_unordered(func_wrapper, iterables), total=len(iterables)):
            iterable, output = result
            indx = indexes[iterables.index(iterable)]
            yield indx, output
