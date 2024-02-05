import json
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as interp
from typing import Any, Dict, List, Literal, Optional, Union
from pathlib import Path



def field_data_load(path: str | Path) -> Dict[str, Any]:
    with open(path, "r+") as infile:
        wake_data = json.load(infile)
    return wake_data


def inflow_validation_plot(
    wake_data: Dict[str, Any],
    ) -> None:
    h_hub = 0.125
    D_rotor = 0.15
    v_hub = 4.88
    ti_hub = 0.07

    height, vel, turb, tke = [], [], [], []
    val_dist = [-3, 0, 3, 5, 7, 10, 15]
    val_name = [f'z_hub_{d}d' for d in val_dist]
    # val_data = ['coordinate', 'x-velocity', 'turb-intensity', 'turb-kinetic-energy']
    val_data = ['coordinate', 'x-velocity', 'turb-intensity']
    for name in val_name:
        for val_data_i, data_i in zip(val_data, [height, vel, turb, tke]):
            if val_data_i == 'coordinate':
                data_i.append(np.array(wake_data[name][val_data_i])[:, 2])
            else:
                data_i.append(np.array(wake_data[name][val_data_i]))

    inflow_height = np.array(wake_data['z_hub_inflow']['coordinate'])[:, 2]
    inflow_vel = np.array(wake_data['z_hub_inflow']['x-velocity'])
    inflow_turb = np.array(wake_data['z_hub_inflow']['turb-intensity'])
    # inflow_tke = np.array(wake_data['z_hub_inflow']['turb-kinetic-energy'])

    inflow_vel_interp = interp.interp1d(inflow_height, inflow_vel)
    inflow_turb_interp = interp.interp1d(inflow_height, inflow_turb)
    # inflow_tke_interp = interp.interp1d(inflow_height, inflow_tke)

    fig, ax = plt.subplots(1, 2, figsize=(20, 8), dpi=100)

    ax[0].plot(inflow_vel / v_hub, inflow_height / D_rotor, label=f'Inflow velocity',
               c="w", lw=0., markersize=7, marker="o", markeredgecolor='k', markeredgewidth=0.7,)
    for i, d in enumerate(val_dist):
        ax[0].plot(vel[i] / v_hub, height[i] / D_rotor, label=f'{d}D')
    ax[0].axvline(1.0, color='k', alpha=0.8, linestyle='--', linewidth=0.5)
    ax[0].axhline(h_hub / D_rotor + 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[0].axhline(h_hub / D_rotor, color='k', alpha=1., linestyle='-', linewidth=1.)
    ax[0].axhline(h_hub / D_rotor - 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[0].set_xlabel('velocity (m/s)')
    ax[0].set_xlim(0.5, 1.5)
    ax[0].set_ylabel('z/D (m)')
    ax[0].set_ylim(0, 3.)
    ax[0].legend(loc='best')

    ax[1].plot(inflow_turb, inflow_height / D_rotor, label=f'Inflow turbulence',
               c="w", lw=0., markersize=7, marker="o", markeredgecolor='k', markeredgewidth=0.7,)
    for i, d in enumerate(val_dist):
        ax[1].plot(turb[i], height[i] / D_rotor, label=f'{d}D')
    # ax[1].axvline(1.0, color='k', alpha=0.8, linestyle='--', linewidth=0.5)
    ax[1].axhline(h_hub / D_rotor + 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[1].axhline(h_hub / D_rotor, color='k', alpha=1., linestyle='-', linewidth=1.)
    ax[1].axhline(h_hub / D_rotor - 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[1].set_xlabel('Turbulence (%)')
    ax[1].set_xlim(0, 0.5)
    ax[1].set_ylabel('z/D (m)')
    ax[1].set_ylim(0, 3.)
    ax[1].legend(loc='best')

    plt.show()


if __name__ == "__main__":
    inflow_validation_plot(field_data_load('data/wake_data.txt'))