import os
import yaml
import copy
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import geatpy as ea
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from typing import Any, List, Dict, Union, Callable, Literal, Optional




class Loader(yaml.SafeLoader):
    def __init__(self, stream: str | Path) -> None:
        self._root = os.path.split(stream.name)[0]
        super().__init__(stream)

    def include(self, node: yaml.Node) -> Dict[Any, Any]:
        filename = os.path.join(self._root, self.construct_scalar(node))
        with open(filename, 'r') as f:
            return yaml.load(f, self.__class__)


Loader.add_constructor('!include', Loader.include)


class DotDict(dict):
    def __init__(self, *args, **kwargs) -> None:
        super(DotDict, self).__init__(*args, **kwargs)

    def __getattr__(self, key: str) -> Any:
        try:
            value = self[key]
            if isinstance(value, dict):
                value = DotDict(value)
            return value
        except:
            return super().__getattribute__(key)

    def __deepcopy__(self,
                     memo: Any,
                     _nil: List[Any] = []
                     ) -> Dict[Any, Any] | List[Any]:
        if memo is None:
            memo = {}
        d = id(self)
        y = memo.get(d, _nil)
        if y is not _nil:
            return y

        dict = DotDict()
        memo[d] = id(dict)
        for key in self.keys():
            dict.__setattr__(copy.deepcopy(key, memo),
                             copy.deepcopy(self.__getattr__(key), memo))
        return dict


class WakeConfig(DotDict):
    def __init__(self, input_config_path: str | Path) -> None:
        input_config = WakeConfig.load_yaml(Path(input_config_path).resolve())
        super().__init__(input_config)

    @classmethod
    def load_yaml(cls, filename: str, loader: yaml.SafeLoader = Loader) -> Dict[Any, Any]:
        with open(filename) as fid:
            return yaml.load(fid, loader)

    def to_yaml(self, output_file_path: str | Path) -> None:
        with open(output_file_path, "w+") as output_file:
            yaml.dump(
                self,
                output_file,
                sort_keys=False,
                default_flow_style=False
                )


def icem_mesh_generation(
    script_path: str | Path,
    log_name: Optional[str | Path] = 'meshing.log',
    write_mode: Literal['w+', 'a+'] = 'w+',
    ) -> None:
    split_line = '=' * 30
    log_path = Path(script_path).parent / log_name
    with open(log_path, write_mode) as output:
        output.write('\n' + split_line + ' MESH BEGIN ' + split_line + '\n') # split line for each mesh begin
        mesh_proc = subprocess.Popen(
            ['pwsh.exe', '-Command', 'icem-mesh', script_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        mesh_proc.wait()
        stdout, stderr = mesh_proc.communicate()
        output.write(stdout.decode('utf-8').replace('\r\n', '\n'))
        if stderr is not None:
            output.write(stderr.decode('utf-8').replace('\r\n', '\n'))
        output.write('\n' + split_line + ' MESH END ' + split_line + '\n')  # split line for each mesh end


def icem_script_modify(
    config: WakeConfig,
    new_script_name: Optional[str | Path] = None,
    ) -> str | Path:

    script_path = Path(config.general.mesh.path, config.general.mesh.script)
    domain = config.mesh.domain_scheme
    bl_max_len = round(domain.length.z * (1 - domain.grid.boundary_layer.len_ratio) / \
        (domain.grid.N_z - domain.grid.boundary_layer.number), 5)

    # Dictionary of modifications
    modifications = {
        "work_dir": '{' + f'{script_path.resolve().parent}' + '}',
        "proj_name": config.general.mesh.name,

        "l_x": str(domain.length.x),
        "l_y": str(domain.length.y),
        "l_z": str(domain.length.z),

        "lx_a": str(domain.axis.x.min * -1),
        "lx_b": str(domain.axis.x.max * 1),
        "ly_l": str(domain.axis.y.min * -1),
        "ly_r": str(domain.axis.y.max * 1),

        "N_x": str(domain.grid.N_x),
        "N_y": str(domain.grid.N_y),
        "N_z": str(domain.grid.N_z),
        "N_za": str(domain.grid.N_z - domain.grid.boundary_layer.number),
        "N_zb": str(domain.grid.boundary_layer.number),

        "max_len": str(domain.grid.max_len),
        "bl_ratio": str(domain.grid.boundary_layer.len_ratio),
        "bl_len": str(round(domain.length.z * domain.grid.boundary_layer.len_ratio, 5)),
        "bl_min_len": str(domain.grid.boundary_layer.min_len),
        "bl_max_len": str(bl_max_len),
        "bl_ratio_1": str(domain.grid.boundary_layer.ratio_1),
        "bl_ratio_2": str(domain.grid.boundary_layer.ratio_2),

        "tol": str(domain.grid.tolerance),
    }
    # print(modifications)

    # Python script to modify specific lines in a file based on a dictionary
    modified_lines = []
    with open(script_path, 'r+') as file:
        for line in file:
            if line.startswith('set'):
                parts = line.split()
                # Check if the line is a 'set' command and the parameter is in the dictionary
                if len(parts) >= 3 and parts[1] in modifications:
                    new_value = modifications[parts[1]]
                    modified_line = f"set {parts[1]} {new_value}\n"
                    modified_lines.append(modified_line)
                    continue
            modified_lines.append(line)

    if new_script_name is not None:
        new_script_path = Path(script_path).parent / 'modified_test_mesh.rpl'
    else:
        new_script_path = Path(script_path)

    # Write the modified content to a new file
    with open(new_script_path, 'w+') as new_file:
        new_file.writelines(modified_lines)

    return new_script_path


def turbine_udf_modify(
    config: WakeConfig,
    new_udf_name: Optional[str | Path | bool] = None,
    ) -> str | Path:

    udf_path = Path(config.udf_file.path, config.udf_file.turbine.name)

    modifications = {
        'inflow': '2.2',
        'theta': '0.0',
        'radius': '0.075',
        'pitch': '0.0',
        'tsr': '4.0',
        'rpm': '1120',
    }

    # Read the original file content
    with open(udf_path, 'r') as file:
        lines = file.readlines()

    # Modify the parameters as per the dictionary
    for i, line in enumerate(lines):
        if line.strip().startswith('#define'):
            parts = line.split()
            if len(parts) >= 3 and parts[1] in modifications:
                # Find the start and end index of the original value
                start_index = line.find(parts[2])
                end_index = start_index + len(parts[2])
                # Replace only the value part, keeping spacing and comments unchanged
                lines[i] = line[:start_index] + modifications[parts[1]] + line[end_index:]

    if new_udf_name is not None:
        new_udf_path = Path(udf_path).parent / f'modified_{config.udf_file.turbine.name}'
    else:
        new_udf_path = Path(udf_path)

    # Write the modified content back to the file
    with open(new_udf_path, 'w') as file:
        file.writelines(lines)

    return new_udf_path


def inflow_udf_modify(
    config: WakeConfig,
    new_udf_name: Optional[str | Path | bool] = None,
    ) -> str | Path:

    udf_path = Path(config.udf_file.path, config.udf_file.inflow.name)

    modifications = {
        'u0': '0.115',      # Friction velocity
        'k': '0.42',        # Karman constant
        'rho': '1.225',     # Air density
        'z0': '0.00003',    # Aerodynamic roughness length
        'h': '0.125',       # Hub height

        'Cmu': '0.09',
        'Cu1': '-0.412',    # First constant in velocity TKE and dissipation profile
        'Cu2': '3.77',      # Second constant in velocity TKE and dissipation profile

        'C1k': '0.005',     # First constant in k source equation
        'C2k': '0.10',      # Second constant in k source equation
        'sigma_k': '1.0',

        'C1w': '1.44',      # First constant in w source equation
        'C2w': '1.92',      # Second constant in w source equation
        'Cw': '0.48',
        'sigma_e': '1.50',

        'ht': '1.86'        # Possibly a height parameter
    }

    # Read the original file content
    with open(udf_path, 'r') as file:
        lines = file.readlines()

    # Modify the parameters as per the dictionary
    for i, line in enumerate(lines):
        if line.strip().startswith('#define'):
            parts = line.split()
            if len(parts) >= 3 and parts[1] in modifications:
                # Find the start and end index of the original value
                start_index = line.find(parts[2])
                end_index = start_index + len(parts[2])
                # Replace only the value part, keeping spacing and comments unchanged
                lines[i] = line[:start_index] + modifications[parts[1]] + line[end_index:]

    if new_udf_name is not None:
        new_udf_path = Path(udf_path).parent / f'modified_{config.udf_file.turbine.name}'
    else:
        new_udf_path = Path(udf_path)

    # Write the modified content back to the file
    with open(new_udf_path, 'w') as file:
        file.writelines(lines)

    return new_udf_path


def inflow_profile_plot(param=None, model='PengHu'):
    D = 0.15
    H = 0.125
    vel_hub = 4.88
    turb_hub = 0.07

    if model == 'PengHu':
        params = {
            'K': 0.4187,
            'rho': 1.225,
            'u_star': 0.3256,
            'z0': 0.000221,

            'Cu1': -0.5235,
            'Cu2': 4.6026,
            'C1k': 0.2720,
            'C2k': 0.7160,
            'C1w': 0.1560,
            'C4w': 0.0150,
        }
    elif model == 'YangYi':
        params = {
            'K': 0.4187,
            'rho': 1.225,
            'u_star': 0.3256,
            'z0': 0.000221,

            'Cmu': 0.4000,
            'C1': -0.9769,  # if < 0 abs should be small; if > 0 abs should be large; (better < 0)
            'C2': 7.2829,   # if C1 < 0 should be large; if C1 > 0 z always larger than negative (better > 0)
            'sig_e': 1.3,
            'C2e_1e': 0.88,
        }
        # repair_func = lambda x: -0.08985 * np.log(3.16246 * x) + 0.20843
        repair_func = lambda x:  0.3 * np.exp(-1.5 * x) + 0.0
    elif model == 'LiWang':
        params = {
            'K': 0.4187,
            'rho': 1.225,
            'u_star': 0.3256,
            'z0': 0.000221,

            'a': 0.8,
            'Cmu': 1.09,
            'C1': 0.288,
            'C2': 0.10,
        }
    elif model == 'ZiLong':
        params = {
            'K': 0.4187,
            'rho': 1.225,
            'u_star': 0.3256,
            'z0': 0.000221,
            'Cmu': 0.09,
        }
    else:
        raise ValueError(f'No such turbulence model: {model}')

    if param is not None:
        params.update(param)

    # ref_inflow_path = 'C:/Users/Li Hang/Documents/Projects/floris/floris/utils/dataset/baseline/2019_Lin_PorteAgel/Fig_3'
    ref_inflow_path = 'C:/Users/Li Hang/Documents/Projects/floris/floris/utils/dataset/baseline/2016_Bastankhah_PorteAgel/Fig_1'
    exp_vel_profile = np.loadtxt(f'{ref_inflow_path}/inflow_vel_profile_exp.txt', skiprows=4)
    exp_turb_profile = np.loadtxt(f'{ref_inflow_path}/inflow_turb_profile_exp.txt', skiprows=4)
    # print(exp_vel_profile.shape, exp_turb_profile.shape)

    turb_model = rans_inflow_turbulence_model(model, params)

    height = np.arange(0, 3.5, 0.1)
    velocity = np.vectorize(turb_model.U_z)(height * D)
    tke = np.vectorize(turb_model.k_z)(height * D)
    if hasattr(turb_model, 'C_z'):
        height_check = np.vectorize(turb_model.Z_z)(height * D)
        minus_check = np.vectorize(turb_model.C_z)(height * D)
        if np.any(minus_check < 0):
            tke = np.where(minus_check >= 0, tke, repair_func(height))
        # print('[z] <= :', - params['C2'] / params['C1'])
        # print('Height: ', height_check)
        print('Check: ', minus_check)
    turbulence = turb_model.tke_to_ti(tke, vel_hub)
    print('U at hub height: ', turb_model.U_z(H))
    print('TI at hub height: ', turb_model.tke_to_ti(turb_model.k_z(H), vel_hub) * 100)

    fig, ax = plt.subplots(1, 3, figsize=(18, 6), dpi=100)
    ax[0].plot(exp_vel_profile[:, 0],
               exp_vel_profile[:, 1],
               c="w", lw=0., label='Exp',
               markersize=7, marker="o",
               markeredgecolor='k',
               markeredgewidth=0.7,
               )
    ax[0].plot(velocity / vel_hub,
               height,
               c='r',
               lw=2.,
               ls='-',
               label='fitting',
               )
    ax[0].axvline(1.0, color='k', alpha=0.8, linestyle='--', linewidth=0.5)
    ax[0].axhline(H / D + 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[0].axhline(H / D, color='k', alpha=1., linestyle='-', linewidth=1.)
    ax[0].axhline(H / D - 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[0].set_xlabel('U [m/s]')
    ax[0].set_ylabel('z/d')
    ax[0].set_xlim([0.5, 1.5])
    ax[0].set_ylim([0, 3.5])
    ax[0].grid(alpha=0.5)
    ax[0].legend()

    ax[1].plot(exp_turb_profile[:, 0],
               exp_turb_profile[:, 1],
               c="w", lw=0., label='Exp',
               markersize=7, marker="o",
               markeredgecolor='k',
               markeredgewidth=0.7,
               )
    ax[1].plot(turbulence,
               height,
               c='r',
               lw=2.,
               ls='-',
               label='fitting',
               )
    ax[1].axvline(turb_hub, color='k', alpha=0.8, linestyle='--', linewidth=0.5)
    ax[1].axhline(H / D + 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[1].axhline(H / D, color='k', alpha=1., linestyle='-', linewidth=1.)
    ax[1].axhline(H / D - 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[1].set_xlabel('I [%]')
    ax[1].set_ylabel('z/d')
    ax[1].set_xlim([0, 0.15])
    ax[1].set_ylim([0, 3.5])
    ax[1].grid(alpha=0.5)
    ax[1].legend()

    diss_model = getattr(turb_model, 'e_z', None) or getattr(turb_model, 'w_z', None)
    diss = np.vectorize(diss_model)(height * D)
    ax[2].plot(diss,
               height,
               c='r',
               lw=2.,
               ls='-',
               label='fitting',
               )
    ax[2].axvline(turb_hub, color='k', alpha=0.8, linestyle='--', linewidth=0.5)
    ax[2].axhline(H / D + 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[2].axhline(H / D, color='k', alpha=1., linestyle='-', linewidth=1.)
    ax[2].axhline(H / D - 1 / 2, color='k', alpha=0.8, linestyle='--', linewidth=1.2)
    ax[2].set_xlabel('e [m2s-3]')
    ax[2].set_ylabel('z/d')
    ax[2].set_xlim([-10, 50.])
    ax[2].set_ylim([0, 3.5])
    ax[2].grid(alpha=0.5)
    ax[2].legend()

    plt.show()


def inflow_vel_profile_fitting(param=None):

    ref_inflow_path = 'C:/Users/Li Hang/Documents/Projects/floris/floris/utils/dataset/baseline/2016_Bastankhah_PorteAgel/Fig_1'
    exp_vel_data = np.loadtxt(f'{ref_inflow_path}/inflow_vel_profile_exp.txt', skiprows=4)
    ref_height = exp_vel_data[:, 1]; ref_vel = exp_vel_data[:, 0]

    fixed_params = {
        'D': 0.15,
        'H': 0.125,
        'v_h': 4.88,
        'K': 0.4187,
        }

    if param is not None:
        fixed_params.update(param)

    optimized_params = {
        'u_star': 0.3391,
        'z0': 0.000215,
        }

    bounds = {
        'u_star': (0.1, 0.5),
        'z0': (0.000001, 0.005),
    }

    initial = [optimized_params[key] for key in optimized_params.keys()]
    bound = [bounds[key] for key in optimized_params.keys()]

    # U_z function definition for fitting
    def U_z(z, u_star, z0):
        return u_star / fixed_params['K'] * np.log((z + z0) / z0)

    # Objective function to be minimized
    def obj_func(params, ref_z, ref_v):
        u_star, z0 = params
        vel = np.vectorize(U_z)(ref_z * fixed_params['D'], u_star, z0)
        vel_ref = ref_v * fixed_params['v_h']
        return np.mean((vel - vel_ref) ** 2)

    result = minimize(
        obj_func,
        initial,
        method='SLSQP',
        args=(ref_height, ref_vel),
        bounds=bound,
        )

    if result.success:
        print('[***] Optimization Results: ')
        print(f'    u_star: {result.x[0]:.4f}')
        print(f'    z0: {result.x[1]:.6f}')
        optimal_param = {'u_star': result.x[0], 'z0': result.x[1]}
        # inflow_profile_plot(param=optimal_param)

        return optimal_param  # Returns the fitted parameters
    else:
        raise RuntimeError(f"Optimization failed: {result.message}")


def inflow_turb_profile_fitting(param=None, model='LiWang'):

    ref_inflow_path = 'C:/Users/Li Hang/Documents/Projects/floris/floris/utils/dataset/baseline/2016_Bastankhah_PorteAgel/Fig_1'
    exp_turb_data = np.loadtxt(f'{ref_inflow_path}/inflow_turb_profile_exp.txt', skiprows=4)
    ref_height = exp_turb_data[:, 1]; ref_turb = exp_turb_data[:, 0]

    fixed_params = {
        'D': 0.15,
        'H': 0.125,
        'v_h': 4.88,

        'rho': 1.225,
        'K': 0.4187,
        'u_star': 0.3256,
        'z0': 0.000221,
        }

    if param is not None:
        fixed_params.update(param)

    if model == 'PengHu':
        optimized_params = {
            'Cu1': -0.33847,
            'Cu2': 2.61666,
            'C1k': 0.172,
            'C2k': 0.716,
            'C1w': 0.156,
            'C4w': 0.015,
            }
        bounds = {
            'Cu1': (-1.0, 0.0),
            'Cu2': (0.0, 10.0),
            'C1k': (0.0, 1.0),
            'C2k': (0.0, 1.0),
            'C1w': (0.0, 1.0),
            'C4w': (0.0, 1.0),
        }
    elif model == 'YangYi':
        optimized_params = {
            'Cmu': 0.9,
            'C1': -1.36,
            'C2': 14.3,
            'sig_e': 1.3,
            'C2e_1e': 0.88,
            }
        bounds = {
            'Cmu': (0.1, 1.0),
            'C1': (-10, 0.0),
            'C2': (0.0, 100.0),
            'sig_e': (0.1, 5.0),
            'C2e_1e': (0.1, 5.0),
        }
    elif model == 'LiWang':
        optimized_params = {
            'a': 0.8,
            'Cmu': 0.09,
            'C1': 3.088,
            'C2': 0.10,
            }
        bounds = {
            'a': (0.1, 1.0),
            'Cmu': (0.01, 0.1),
            'C1': (0.01, 10.0),
            'C2': (0.01, 10.0),
        }
    else:
        raise ValueError(f'No such turbulence model: {model}')

    initial = [optimized_params[key] for key in optimized_params.keys()]
    bound = [bounds[key] for key in optimized_params.keys()]

    # Objective function to be minimized
    def obj_func(params, ref_z, ref_v):
        opt_params = {key: params[i] for i, key in enumerate(optimized_params.keys())}
        opt_params.update(fixed_params)
        turb_model = rans_inflow_turbulence_model(model, opt_params)
        turb = turb_model.tke_to_ti(
            np.vectorize(turb_model.k_z)(ref_z * turb_model.D), turb_model.v_h)
        # loss = np.mean(np.sum((turb - ref_v) ** 2)) * 100.
        loss = np.mean(np.abs(turb - ref_v)) * 100.
        return loss

    result = minimize(
        obj_func,
        initial,
        method='SLSQP',
        args=(ref_height, ref_turb),
        bounds=bound,
        )

    if result.success:
        optimal_param = {p:result.x[i] for i, p in enumerate(optimized_params.keys())}
        print('[***] Optimization Results: ')
        for p, v in optimal_param.items():
            print(f'    {p}: {v:.4f}')
        optimal_param.update(fixed_params)
        inflow_profile_plot(param=optimal_param, model=model)

        return optimal_param  # Returns the fitted parameters
    else:
        raise RuntimeError(f"Optimization failed: {result.message}")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                MISCELLANEOUS                                 #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

class rans_inflow_turbulence_model():
    def __init__(self, model: str, param: Optional[Dict[str, Any]]) -> None:
        default_models = {
            'PengHu': rans_inflow_turbulence_model._PengHu_model,
            'YangYi': rans_inflow_turbulence_model._YangYi_model,
            'LiWang': rans_inflow_turbulence_model._LiWang_model,
            'ZiLong': rans_inflow_turbulence_model._ZiLong_model,
        }
        turb_model = default_models.get(model, None)
        if turb_model is None:
            raise ValueError(f'No such turbulence model: {model}')
        else:
            param.update(turb_model(param))
            for key, value in param.items():
                setattr(self, key, value)

    @staticmethod
    def tke_to_ti(k, v_h, factor=1.0):
        return np.sqrt(k * 2 / 3) / v_h
        # return np.sqrt(k * 1.098) / v_h

    @classmethod
    def _PengHu_model(cls, params: Dict[str, Any]) -> Dict[str, Any]:
        required_params = ['u_star', 'K', 'rho', 'z0', 'Cu1',
                           'Cu2', 'C1k', 'C2k', 'C1w', 'C4w']
        assert all([param in params.keys() for param in required_params]), \
            f'Not all required parameters are provided: {required_params}'

        U_z = lambda z: params['u_star'] / params['K'] * \
            np.log((z + params['z0']) / params['z0'])
        k_z = lambda z: params['u_star']**2 * (params['Cu1'] * \
            np.log((z + params['z0']) / params['z0']) + params['Cu2'])**2
        w_z = lambda z: params['u_star'] / (params['K'] * (z + params['z0'])) * \
            (params['Cu1'] * np.log((z + params['z0']) / params['z0']) + params['Cu2'])**2
        S_k = lambda z: params['rho'] * params['u_star']**3 / (z + params['z0']) * \
            (params['C1k'] * (params['Cu1'] * np.log((z + params['z0']) / params['z0']) + \
                params['Cu2'])**4 - params['C2k'])
        S_w = lambda z: params['rho'] * params['u_star']**2 / (z + params['z0'])**2 * \
            (params['C1w'] * (params['Cu1'] * np.log((z + params['z0']) / params['z0']) + \
                params['Cu2'])**4 - params['C4w'])

        return {'U_z': U_z, 'k_z': k_z, 'w_z': w_z, 'S_k': S_k, 'S_w': S_w}

    @classmethod
    def _YangYi_model(cls, params: Dict[str, Any]) -> Dict[str, Any]:
        required_params = ['u_star', 'K', 'rho', 'z0', 'Cmu',
                           'C1', 'C2', 'sig_e', 'C2e_1e']
        assert all([param in params.keys() for param in required_params]), \
            f'Not all required parameters are provided: {required_params}'

        U_z = lambda z: params['u_star'] / params['K'] * \
            np.log((z + params['z0']) / params['z0'])
        k_z = lambda z: params['u_star']**2 / np.sqrt(params['Cmu']) * np.sqrt(params['C1'] * \
            np.log((z + params['z0']) / params['z0']) + params['C2'])
        e_z = lambda z: params['u_star']**3 / (params['K'] * (z + params['z0'])) * \
            np.sqrt(params['C1'] * np.log((z + params['z0']) / params['z0']) + params['C2']) if \
                        (params['C1'] * np.log((z + params['z0']) / params['z0']) + params['C2']) >= 0 else \
                            params['u_star']**3 / (params['K'] * (z + params['z0']))
        S_k = lambda z: 0.0 * z
        S_e = lambda z: params['rho'] * params['u_star']**4 / (params['K'] * (z + params['z0']))**2 / \
            params['sig_e'] * (params['K']**2 * (1.5 * params['C1'] - params['C1'] * np.log((z + params['z0']) / \
                params['z0']) - params['C2']) + np.sqrt(params['Cmu']) * params['sig_e'] * params['C2e_1e'] * \
                    np.sqrt(params['C1'] * np.log((z + params['z0']) / params['z0']) + params['C2'])) if \
                        (params['C1'] * np.log((z + params['z0']) / params['z0']) + params['C2']) >= 0 else \
                            params['rho'] * params['u_star']**4 / (params['K'] * (z + params['z0']))**2 / \
                                params['sig_e'] * (params['K']**2 * (1.5 * params['C1'] - params['C1'] * np.log((z + params['z0']) / \
                                    params['z0']) - params['C2']))
        Z_z = lambda z: np.log((z + params['z0']) / params['z0'])
        C_z = lambda z: params['C1'] * np.log((z + params['z0']) / params['z0']) + params['C2']

        return {'U_z': U_z, 'k_z': k_z, 'e_z': e_z, 'S_k': S_k, 'S_e': S_e, 'Z_z':Z_z, 'C_z': C_z}

    @classmethod
    def _LiWang_model(cls, params: Dict[str, Any]) -> Dict[str, Any]:
        required_params = ['u_star', 'K', 'rho', 'z0', 'a', 'Cmu', 'C1', 'C2']
        assert all([param in params.keys() for param in required_params]), \
            f'Not all required parameters are provided: {required_params}'

        U_z = lambda z: params['u_star'] / params['K'] * \
            np.log((z + params['z0']) / params['z0'])
        k_z = lambda z: params['a'] * params['u_star']**2 * (params['C1'] * \
            np.log((z + params['z0']) / params['z0']) + params['C2'])**2
        e_z = lambda z: params['a']**2 * params['Cmu'] * params['u_star']**3 * \
            (params['C1'] * np.log((z + params['z0']) / params['z0']) + params['C2'])**4 \
                / (params['K'] * (z + params['z0']))
        S_k = lambda z: 0.0 * z
        S_e = lambda z: 0.0 * z

        return {'U_z': U_z, 'k_z': k_z, 'e_z': e_z, 'S_k': S_k, 'S_e': S_e}

    @classmethod
    def _ZiLong_model(cls, params: Dict[str, Any]) -> Dict[str, Any]:
        required_params = ['u_star', 'K', 'rho', 'z0', 'Cmu',]
        assert all([param in params.keys() for param in required_params]), \
            f'Not all required parameters are provided: {required_params}'

        U_z = lambda z: params['u_star'] / params['K'] * \
            np.log((z + params['z0']) / params['z0'])
        k_z = lambda z: params['u_star']**2 / np.sqrt(params['Cmu'])
        e_z = lambda z: params['u_star']**3 / (params['K'] * (z + params['z0']))
        S_k = lambda z: 0.0 * z
        S_e = lambda z: 0.0 * z

        return {'U_z': U_z, 'k_z': k_z, 'e_z': e_z, 'S_k': S_k, 'S_e': S_e}





if __name__ == "__main__":
    # config = WakeConfig('config_test.yaml')

    # icem_mesh_generation('mesh/test_mesh_eq.rpl', 'meshing_eq.log')

    # icem_script_modify(config)
    # turbine_udf_modify(config)
    # inflow_udf_modify(config)

    inflow_profile_plot()

    # inflow_turb_profile_fitting(inflow_vel_profile_fitting())
    # inflow_turb_profile_fitting()