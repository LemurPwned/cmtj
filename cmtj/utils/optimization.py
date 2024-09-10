from concurrent.futures import ProcessPoolExecutor
from typing import Callable

import numpy as np
from tqdm import tqdm


def coordinate_descent(
    operating_point: dict[str, float],
    fn: Callable,
    best_mse: float = float("-inf"),
    granularity: int = 10,
    percentage: float = 0.05,
):
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
        for v in tqdm(
            np.linspace((1 - percentage) * org_v, (1 + percentage) * org_v, granularity),
            desc=f"Optimising {k}",
            leave=False,
        ):
            new_params[k] = v
            mse = fn(**new_params)
            if mse > best_mse:
                opt_params = new_params
                best_mse = mse
    return opt_params, best_mse


def multiprocess_simulate(
    fn: Callable,
    error_fn: Callable,
    suggestions: list[dict],
    target: np.ndarray,
    fixed_parameters: dict,
):
    with ProcessPoolExecutor(max_workers=len(suggestions)) as executor:
        futures = [
            executor.submit(
                fn,
                **fixed_parameters,
                **suggestion,
            )
            for suggestion in suggestions
        ]
        errors = np.zeros(len(suggestions))
        for j, future in enumerate(futures):
            result = future.result()
            errors[j] = error_fn(target, result)
    return errors


def hebo_optimization_loop(
    cfg: dict,
    fn: Callable,
    error_fn: Callable,
    target: np.ndarray,
    fixed_parameters: dict,
    n_iters: int = 150,
    n_suggestions: int = 8,
):
    """Optimizes the parameters of a function using HEBO.
    See HEBO documentation for more details: https://github.com/huawei-noah/HEBO
    :param cfg: configuration of the design space
    :param fn: function to be optimised fn(**parameters, **fixed_parameters)
    :param error_fn: function to compute the error: error_fn(target, result)
    :param target: target data
    :param fixed_parameters: parameters that are fixed
    :param n_iters: number of iterations
    :param n_suggestions: number of suggestions per iteration
    """
    try:
        from hebo.design_space.design_space import DesignSpace
        from hebo.optimizers.hebo import HEBO
    except ImportError as e:
        raise ImportError("HEBO is not installed. Please install it with `pip install HEBO`") from e
    space = DesignSpace().parse(cfg)
    opt = HEBO(space)
    best_mse = float("inf")
    for i in tqdm(range(1, n_iters + 1), desc="HEBO optimization loop"):
        rec = opt.suggest(n_suggestions)
        errors = multiprocess_simulate(
            fn=fn,
            error_fn=error_fn,
            suggestions=rec.to_dict(orient="records"),
            target=target,
            fixed_parameters=fixed_parameters,
        )
        opt.observe(rec, errors)
        val = opt.y.min()
        if val < best_mse:
            best_mse = val
            best_params = opt.best_x.iloc[0].to_dict()
            print(f"iteration {i} best mse {best_mse}")
            print(best_params)
    return opt
