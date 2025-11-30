from itertools import product
from typing import Callable

import numpy as np
from multiprocess import Pool
from tqdm import tqdm

from ..models.general_sb import LayerDynamic

__all__ = ["distribute", "parallel_vsd_sb_model"]


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


def parallel_vsd_sb_model(
    simulation_fn: Callable,
    frequencies: list[float],
    Hvecs: list[list[float]],
    layers: list[LayerDynamic],
    J1: list[float] = None,
    J2: list[float] = None,
    iDMI: list[float] = None,
    n_cores: int = None,
):
    """
    Parallelise the VSD SB model.
    :param simulation_fn: function to be distributed.
        This function must take a tuple of arguments, where the first argument is the
        frequency, then Hvectors, the list of layers and finally the list of J1 and J2 values.
    :param frequencies: list of frequencies
    :param Hvecs: list of Hvectors in cartesian coordinates
    :param layers: list of layers
    :param J1: list of J1 values
    :param J2: list of J2 values
    :param n_cores: number of cores to use.
    :returns: list of simulation_fn outputs for each frequency
    """
    if J1 is None:
        J1 = [0] * (len(layers) - 1)
    if J2 is None:
        J2 = [0] * (len(layers) - 1)
    if iDMI is None:
        iDMI = [0] * (len(layers) - 1)
    args = [(f, Hvecs, *layers, J1, J2, iDMI) for f in frequencies]
    with Pool(processes=n_cores) as pool:
        return list(tqdm(pool.imap(simulation_fn, args), total=len(frequencies)))
