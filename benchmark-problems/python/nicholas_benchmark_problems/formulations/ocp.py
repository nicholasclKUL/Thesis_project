import numpy as np
import casadi as cs
from typing import List, Tuple
from dataclasses import dataclass

@dataclass
class OCProblem:
    N: int  # horizon
    Ts: float   # time step   
    state_var: cs.SX # SX: simbolyc expresion
    input_var: cs.SX
    state_names: List[str]
    input_names: List[str]
    initial_guess: np.ndarray 
    initial_state: np.ndarray
    f_dynamics: cs.Function
    stage_constr: cs.Function
    stage_constr_box: Tuple[np.ndarray, np.ndarray]
    term_constr: cs.Function
    term_constr_box: Tuple[np.ndarray, np.ndarray]
    input_constr_box: Tuple[np.ndarray, np.ndarray]
    stage_cost: cs.Function #?
    term_cost: cs.Function

    @property
    def nu(self) -> int:                                # nu: number of inputs -> integer (of course!) 
        assert self.input_var.shape[1] == 1
        return self.input_var.shape[0]

    @property
    def nx(self) -> int:                                # nx: number of states -> integer (of course!)
        assert self.state_var.shape[1] == 1
        return self.state_var.shape[0]

    @property
    def nc(self) -> int:                                # nc: number of constrains -> integer (of course!)
        return self.stage_constr_box[0].shape[0]

    @property
    def nc_N(self) -> int:                              # total number of constrains for the horizon N????
        return self.term_constr_box[0].shape[0]

    def simulate(self, uk):                             # uk: input at step k or all the inputs until step k???
        nx, nu = self.nx, self.nu
        N = len(uk) // nu
        result = np.empty((N * (nu + nx) + nx))         # creates a vector to store all the states and inputs for the full horizon
        result[:nx] = self.initial_state
        for i in range(N):
            result[i * (nx + nu) + nx:i * (nx + nu) + nx + nu] = uk[i * nu: i * nu + nu]
            result[(i + 1) * (nx + nu):(i + 1) * (nx + nu) + nx] = \
                self.f_dynamics(result[i * (nx + nu):i * (nx + nu) + nx], 
                                result[i * (nx + nu) + nx:i * (nx + nu) + nx + nu]).full().squeeze()
        return result

__all__ = [
    'OCProblem',
]
