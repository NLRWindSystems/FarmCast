from farmcast.generate_cases import generate_cases
from farmcast.generate_slurm_files import create_slurm_ff_files, create_slurm_ts_files
import os
import numpy as np

run_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.dirname(run_dir)

# Set the username for Kestrel
hpc_email = 'pbortolo@nrel.gov'
path2turbsim = '/projects/windse/cbay/solvers/turbsim'
path2fastfarm = '/projects/windse/cbay/solvers/FAST.Farm'
path2controller = '/home/pbortolo/ROSCO/ROSCO_v2p10p1/rosco/controller/build/libdiscon.so'

# Set the output directory for the generated files
set = 0 #0 debug, 1 full, 2 scan at 8m/s
# output_dir = os.path.join(os.path.dirname(base_dir), "FarmCast_runs_torque2026_set"+str(set))
output_dir = "/scratch/pbortolo/FarmCast_runs_torque2026_set"+str(set)

# Turbines in the farm
n_turbines = 3
model = "IEA-3.4-130-RWT"
rotor_diameter = 130.0
hub_height = 110.0


# Array of wind speeds in m/s
if set == 0:
    ws = [14.]
elif set == 1:
    ws = [6., 10., 14.]
elif set == 2:
    ws = [8.]
# Number of seeds
if set == 0:
    n_seeds = 1
else:
    n_seeds = 12
 # Array of turbulence intensities
if set == 0 or set == 2:
    TI = [0.1]
else:
    TI = [0.1, 0.2]
# Array of shear coefficients
if set == 0 or set == 2:
    shear = [0.1]
else:
    shear = [0.1, 0.2]
# Array of turbine spacing in rotor diameters
if set == 0 or set == 2:
    spacing = [6.]
else:
    spacing = [4., 6., 8.]
# Array of wind directions in degrees
if set == 0:
    wind_direction = [8.]
elif set == 1:
    wind_direction = [-8., 0., 8.]
elif set == 2:
    wind_direction = [-8., -4., 0., 4., 8.]
# Array of yaw misalignments for the upstream turbine (T1) in degrees
if set == 0:
    T1_yaw_misalignment = [30.]
elif set == 1:
    T1_yaw_misalignment = [-30., 0., 30.]
elif set == 2:
    T1_yaw_misalignment = [-30., -20., -10., 0., 10., 20., 30.]
# Array of yaw misalignments for the middle turbine (T2) in degrees
if set == 0:
    T2_yaw_misalignment = [-20.]
elif set == 1:
    T2_yaw_misalignment = [-20., 0., 20.]
elif set == 2:
    T2_yaw_misalignment = [-20., -10., 0., 10., 20.]
# Array of curtailment values for T1 and T2 in percentage
if set == 0:
    curtailment_T1T2 = [60.]
elif set == 1:
    curtailment_T1T2 = [60., 100.]
elif set == 2:
    curtailment_T1T2 = [60., 70., 80., 90., 100.]
# Wake model
Mod_Wake = 2 # Curled, 1 Polar, 3 Cartesian

# Estimate the total number of cases
n_cases = (
    len(ws)
    * n_seeds
    * len(TI)
    * len(shear)
    * len(spacing)
    * len(wind_direction)
    * len(T1_yaw_misalignment)
    * len(T2_yaw_misalignment)
    * len(curtailment_T1T2)
)
print(f"Total number of cases: {n_cases}")

# Generate the cases
turbsim_lr, turbsim_hr = generate_cases(
    n_turbines=n_turbines,
    model=model,
    rotor_diameter=rotor_diameter,
    hub_height=hub_height,
    ws=ws,
    n_seeds=n_seeds,
    TI=TI,
    shear=shear,
    spacing=spacing,
    wind_direction=wind_direction,
    T1_yaw_misalignment=T1_yaw_misalignment,
    T2_yaw_misalignment=T2_yaw_misalignment,
    curtailment_T1T2=curtailment_T1T2,
    domain_edge_LR=[2., 20., 3., 3., 3.],
    output_dir=output_dir,
    Mod_Wake = Mod_Wake,
    path2controller=path2controller,
)

# Create the slurm files for low res turbsim
slurm_dir = os.path.join(output_dir, "slurm_files", "turbsim_lr")
create_slurm_ts_files(turbsim_lr, slurm_dir, slurm_email = hpc_email, path2turbsim = path2turbsim)
slurm_dir = os.path.join(output_dir, "slurm_files", "turbsim_hr")
create_slurm_ts_files(turbsim_hr, slurm_dir, slurm_email = hpc_email, path2turbsim = path2turbsim)

# Create the slurm files for each case
create_slurm_ff_files(n_cases, n_turbines, output_dir, slurm_email = hpc_email, path2fastfarm = path2fastfarm)

print(f"All {n_cases} successfully generated in {output_dir}.")