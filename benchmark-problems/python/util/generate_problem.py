import importlib
import alpaqa.casadi_generator as cg
from sys import argv
import numpy as np

if len(argv) < 4 or argv[2] not in ["ocp", "ss", "ss2"]:
    print(f"Usage:    {argv[0]} <problem-name> ocp|ss <horizon>")
    exit(0)

name, formulation = argv[1:3]
N = int(argv[3])
module = importlib.import_module("nicholas_benchmark_problems.problems." + name)
name += "_" + formulation + "_" + str(N)
ocp = module.Problem(N=N)

txt_opts = {"delimiter": "\t", "newline": "\n"}

if formulation == "ocp":
    cg.generate_casadi_control_problem(
        f=ocp.f_dynamics,
        l=ocp.stage_cost,
        l_N=ocp.term_cost,
        c=ocp.stage_constr,
        c_N=ocp.term_constr,
        name=name,
    ).generate()

    with open(f"{name}.tsv", "w") as f:
        np.savetxt(f, ocp.input_constr_box, **txt_opts)
        np.savetxt(f, ocp.stage_constr_box, **txt_opts)
        np.savetxt(f, ocp.term_constr_box, **txt_opts)
        np.savetxt(f, [ocp.initial_state], **txt_opts)
        np.savetxt(f, [ocp.initial_guess], **txt_opts)
else:
    from nicholas_benchmark_problems.formulations.ss import ocp_to_ss

    ss = ocp_to_ss(ocp)
    cg.generate_casadi_problem(
        f=ss.cost,
        g=ss.constr,
        name=name,
        second_order=formulation == "ss2",
    )[0].generate()

    with open(f"{name}.tsv", "w") as f:
        np.savetxt(f, ss.C, **txt_opts)
        np.savetxt(f, ss.D, **txt_opts)
        np.savetxt(f, [ss.initial_state], **txt_opts)
        np.savetxt(f, [ss.initial_guess], **txt_opts)
