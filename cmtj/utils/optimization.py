from typing import Callable, Dict

import numpy as np
from tqdm import tqdm


def coordinate_descent(operating_point: Dict[str, float],
                       fn: Callable,
                       best_mse: float = float("-inf"),
                       granularity: int = 10,
                       percentage: float = 0.05):
    """Performs coordinate descent on the operating point.
    :param operating_point: operating point to be optimised. Order of that dict matters.
    :param fn: function to be optimised
    :param best_mse: best mse so far
    :param granularity: granularity of the search
    :param percentage: percentage of the search
    :returns: best operating point, best mse
    """
    opt_params = operating_point
    for k, org_v in tqdm(operating_point.items(), desc="Coordinate descent"):
        new_params = operating_point.copy()
        for v in tqdm(np.linspace((1 - percentage) * org_v,
                                  (1 + percentage) * org_v, granularity),
                      desc=f"Optimising {k}",
                      leave=False):
            new_params[k] = v
            mse = fn(**new_params)
            if mse > best_mse:
                opt_params = new_params
                best_mse = mse
    return opt_params, best_mse
