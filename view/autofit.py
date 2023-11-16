from itertools import groupby
from typing import List

import numpy as np
import streamlit as st
from bayes_opt import BayesianOptimization
from hebo.design_space.design_space import DesignSpace
from hebo.optimizers.hebo import HEBO
from helpers import read_data
from scipy.optimize import linear_sum_assignment
from simulation_fns import (compute_sb_mse, get_axis_angles,
                            prepare_simulation, simulate_sb)
from utils import extract_max_resonance_lines

from cmtj.utils import VectorObj, mu0
from cmtj.utils.procedures import PIMM_procedure


def hebo_fit(
    hvals: List[float],
    y: List[float],
    design_space: DesignSpace,
    N: int,
    n_suggestions: int = 1,
):
    def target_fn(yhat, y):
        return np.asarray([compute_sb_mse(target=y, data=yhat)])

    opt = HEBO(design_space)
    prog_bar = st.progress(0, text="Optimisation")
    try:
        for i in range(1, N + 1):
            rec = opt.suggest(n_suggestions=n_suggestions)
            param_values = rec.to_dict("records")
            for k, v in param_values[0].items():
                st.session_state[k] = v
            result_dictionary = simulate_sb(hvals=hvals)
            opt.observe(rec, target_fn(yhat=result_dictionary, y=y))
            try:
                progres_text = (
                    f"###\niteration {i}/{N}: {opt.y.min()}\n"
                    f"Params: {opt.best_x.iloc[0].to_dict()}"
                )
                prog_bar.progress(i, text=progres_text)
            except Exception as e:
                print(f"Error printing progress: {e}")
    except KeyboardInterrupt:
        print("Optimisation interrupted")
        return opt


def autofit(placeholder=None):
    N = st.session_state.N
    cfg = []
    bounds = {
        "Ms": (0.2, 2.0),
        "K": (10, 10e6),
    }
    for param_name in ("Ms", "K"):
        cfg.extend(
            {
                "name": f"{param_name}{i}",
                "type": "num",
                "lb": bounds[param_name][0],
                "ub": bounds[param_name][1],
            }
            for i in range(N)
        )
    cfg.extend(
        {
            "name": f"J{i}",
            "type": "num",
            "lb": -1e3,
            "ub": 1e3,
        }
        for i in range(N - 1)
    )
    try:
        target_data = read_data()
        h, f = target_data
        f = np.asarray(f).ravel().squeeze().tolist()
        target_data = (h, f)
        # TODO: sort
        # target_data = sorted(zip(*target_data))
    except AttributeError:
        st.write("Upload the data to start optimization!")
    else:
        hebo_fit(
            hvals=h,
            y=f,
            design_space=DesignSpace().parse(cfg),
            N=st.session_state.n_iters,
        )


def initialise_autofit(
    pbounds: dict,
    initial_values: dict,
    target_data: list,
    n_iters: int = 50,
):
    def optimise_sb(**kwargs):
        for k, v in kwargs.items():
            st.session_state[k] = v
        return simulate_sb(target_data[0], target_data[1])

    optimizer = BayesianOptimization(
        f=optimise_sb,
        pbounds=pbounds,
        # verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state=1,
    )
    optimizer.probe(initial_values)
    optimizer.maximize(
        init_points=5,
        n_iter=n_iters,
    )
    return optimizer


def initialise_autofit2(
    pbounds: dict,
    initial_values: dict,
    target_data: list,
    n_iters: int = 50,
):
    def mse_target(target_data, model_data):
        total_mse = 0
        for h, freqs in groupby(target_data, lambda x: x[0]):
            freqs = list(freqs)
            print(h, model_data)
            model_freqs = model_data[h]
            if isinstance(model_freqs, np.ndarray):
                model_freqs = model_freqs.tolist()
            cost_matrix = np.zeros((len(freqs), len(model_freqs)))
            for f in freqs:
                for mf in model_freqs:
                    cost = np.inf
                    if mf is not None:
                        cost = (f - mf) ** 2
                    print(cost)
                    cost_matrix[freqs.index(f), model_freqs.index(mf)] = cost
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            total_mse += cost_matrix[row_ind, col_ind].sum()
        return -total_mse

    def simulate(**kwargs):
        for k, v in kwargs.items():
            st.session_state[k] = v

        htheta, hphi = get_axis_angles(st.session_state.H_axis)

        Hvecs = [
            VectorObj.from_spherical(theta=htheta, phi=hphi, mag=mag)
            for mag in target_data[0]
        ]
        j, rparams = prepare_simulation()
        spec, freqs, _ = PIMM_procedure(
            j,
            Hvecs=Hvecs,
            int_step=st.session_state.int_step,
            resistance_params=rparams,
            max_frequency=45e9,
            simulation_duration=12e-9,
            disturbance=1e-6,
            wait_time=4e-9,
            disable_tqdm=True,
        )

        # model frequencies are the same for all H values
        model_data = extract_max_resonance_lines(
            spectrum=spec,
            h_vals=target_data[0],
            frequencies=freqs,
            N_layers=st.session_state.N,
        )
        return mse_target(target_data, model_data)

    optimizer = BayesianOptimization(
        f=simulate,
        pbounds=pbounds,
        # verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state=1,
    )
    optimizer.probe(initial_values)
    optimizer.maximize(
        init_points=5,
        n_iter=n_iters,
    )
    return optimizer
