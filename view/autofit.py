import numpy as np
import streamlit as st
from helpers import plot_optim, read_data
from simulation_fns import compute_sb_mse, simulate_sb_wrapper

is_optimiser_imported = False
try:
    from hebo.design_space.design_space import DesignSpace
    from hebo.optimizers.hebo import HEBO

    is_optimiser_imported = True
except (ImportError, AttributeError):
    is_optimiser_imported = False


def get_fixed_arguments_from_state():
    thickness = [st.session_state[f"thickness{i}"] for i in range(st.session_state.N)]
    anisotropy_axis = [st.session_state[f"anisotropy_axis{i}"] for i in range(st.session_state.N)]
    return {
        "thickness": thickness,
        "anisotropy_axis": anisotropy_axis,
        "H_axis": st.session_state.H_axis,
        "N": st.session_state.N,
    }


def hebo_fit(
    hvals: list[float],
    y: list[float],
):
    """Optimizes the parameters of a function using HEBO."""
    N = st.session_state.n_iters
    n_suggestions = st.session_state.n_suggestions
    cfg, fixed = get_config()
    all_fixed = fixed | get_fixed_arguments_from_state()
    design_space = DesignSpace().parse(cfg)
    opt = HEBO(design_space)
    prog_bar = st.progress(0, text="Optimisation in progress...")
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
                result_dictionary = simulate_sb_wrapper(hvals=hvals, **param_set, **all_fixed)
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

        try:
            opt_params = opt.best_x.iloc[0].to_dict()
            simulated_data = simulate_sb_wrapper(hvals=hvals, **opt_params, **all_fixed)
        except ValueError:
            simulated_data = None
    except KeyboardInterrupt:
        print("Optimisation interrupted")
    return opt, simulated_data


def get_config():
    N = st.session_state.N
    cfg = []
    # keep the same units as in the GUI
    fixed_parameters = {}
    for i in range(N):
        for param_name in ("Ms", "K", "J"):
            # check if fixed
            if param_name == "J" and i >= (N - 1):
                # we don't have a J
                continue
            if st.session_state[f"check_{param_name}{i}"]:
                # take lower bound
                fixed_parameters[f"{param_name}{i}"] = float(st.session_state[f"low_{param_name}{i}"])
            else:
                ub = float(st.session_state[f"up_{param_name}{i}"])
                lb = float(st.session_state[f"low_{param_name}{i}"])
                if ub == lb:
                    # fixed parameter
                    fixed_parameters[f"{param_name}{i}"] = ub
                else:
                    if ub < lb:
                        st.error(
                            f"Upper bound ({ub}) was smaller than lower bound ({lb}) for {param_name}{i}, swapping"
                        )
                        lb, ub = ub, lb
                    cfg.append(
                        {
                            "name": f"{param_name}{i}",
                            "type": "num",
                            "lb": float(st.session_state[f"low_{param_name}{i}"]),
                            "ub": float(st.session_state[f"up_{param_name}{i}"]),
                        }
                    )
    return cfg, fixed_parameters


def autofit(placeholder=None):
    if not is_optimiser_imported:
        st.toast("Hebo optimiser failed to import. Please report the issue to the developer!")
        return
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
        result, simulated_data = hebo_fit(hvals=h, y=f)
        if len(result.y) == 0:
            placeholder.markdown("Optimisation failed!")
            return
        else:
            # turn dict into list
            out_params = "\n".join(f"\t{k}: {v}" for k, v in result.best_x.iloc[0].to_dict().items())
            msg = f"""## OPTIMISATION COMPLETE\n
                Best MSE: {result.y.min():.2f}\n{out_params}
                """
            placeholder.markdown(msg)
            st.write(msg)
            plot_optim(h, f, simulated_data)

            for k, v in result.best_x.iloc[0].to_dict().items():
                st.session_state[k] = v
