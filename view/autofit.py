from typing import List

import numpy as np
import streamlit as st
from hebo.design_space.design_space import DesignSpace
from hebo.optimizers.hebo import HEBO
from helpers import read_data
from simulation_fns import compute_sb_mse, simulate_sb


def hebo_fit(
    hvals: List[float],
    y: List[float],
    design_space: DesignSpace,
    N: int,
    n_suggestions: int = 4,
):
    opt = HEBO(design_space)
    prog_bar = st.progress(0, text="Optimisation")
    try:
        for i in range(1, N + 1):
            rec = opt.suggest(n_suggestions=n_suggestions)
            param_values = rec.to_dict("records")
            errors = []
            # errors = multiprocess_simulate(
            #     fn=simulate_sb_wrapper,
            #     error_fn=compute_sb_mse,
            #     target=y,
            #     suggestions=param_values,
            #     fixed_parameters={
            #         "hvals": hvals,
            #         **get_fixed_arguments_from_state(),
            #     },
            # )

            """Multiprocess is not really stable, so we use a single process instead"""
            for param_set in param_values:
                for k, v in param_set.items():
                    st.session_state[k] = v
                result_dictionary = simulate_sb(hvals=hvals)
                errors.append(compute_sb_mse(target=y, data=result_dictionary))
            errors = np.asarray(errors)
            opt.observe(rec, errors)
            try:
                val = opt.y.min()
                progres_text = f"({i}/{N}) MSE: {val:.2f}"
                prog = int((i * 100) / N)
                prog_bar.progress(prog, text=progres_text)
            except ValueError as e:
                print(f"Error printing progress: {e}")
    except KeyboardInterrupt:
        print("Optimisation interrupted")
        return opt
    return opt


def autofit(placeholder=None):
    N = st.session_state.N
    cfg = []
    # keep the same units as in the GUI
    bounds = {
        "Ms": (0.2, 2.5),  # in T
        "K": (0.1, 500e3),  # in kJ/m^3
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
            "lb": -5e3,  # in uJ/m^2
            "ub": 5e3,
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
        result = hebo_fit(
            hvals=h,
            y=f,
            design_space=DesignSpace().parse(cfg),
            N=st.session_state.n_iters,
        )
        if len(result.y) == 0:
            placeholder.markdown("Optimisation failed!")
            return
        else:
            placeholder.markdown(
                f"""## OPTIMISATION COMPLETE\n
                Best MSE: {result.y.min():.2f}\n
                Best parameters: {result.best_x.iloc[0].to_dict()}
                """
            )
            for k, v in result.best_x.iloc[0].to_dict().items():
                st.session_state[k] = v
