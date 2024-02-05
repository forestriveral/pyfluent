import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from typing import Any, List, Dict, Union, Callable, Literal, Optional



def inflow_profile_plot():
    D = 0.15
    H = 0.125
    v_h = 4.88
    K = 0.4187

    params = {
        'u_star': 0.24,
        'z0': 0.00003,
    }

    bounds = {
        'u_star': (0.1, 0.5),
        'z0': (0.00001, 0.0001),
    }

    U_z = lambda z: params['u_star'] / K * np.log((z + params['z0']) / params['z0'])

    exp_vel_profile = np.loadtxt('inflow_vel_profile_exp.txt', skiprows=4)
    vel_ref_h = exp_vel_profile[:, 0]
    vel_ref = exp_vel_profile[:, 1]

    vel_fitting = np.vectorize(U_z)(vel_ref_h * D) / v_h

    print(vel_fitting)
    print(vel_ref)



if __name__ == '__main__':
    inflow_profile_plot()