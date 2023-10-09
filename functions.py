# declaration of functions for the LLM to use.
# populate functions with the functions you want to use in the LLM
# contains:
#   - build functions
#   - Analysis and Vis functions:
#       - UMAPS and SOAP
#       - Plotting convergence
#   - MLP functions:
#       - NEQUIP, *ACE, *CASTEP-GAP
#   - MD functions

import numpy as np
import os, sys, glob
import matplotlib.pyplot as plt
import seaborn as sns


# example function where search(query) is the name of the function
def search(query):
    return None
functions = [
    {
        "name": "search",
        "description": "google search for relevant information",
        "parameters": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Google search query",
                }
            },
            "required": ["query"],
        },
    },  
]



####################################################################################################
# Build functions and data
####################################################################################################
# build class: https://git.chem.ox.ac.uk/vld-group/theme-amorphous/c-mace-23/evolutionary_committee_models/-/blob/main/ase_build/build_structures.ipynb
# load-atoms (John G)
# 


####################################################################################################
# Analysis and Vis functions
####################################################################################################
'''
UMAPS and SOAP (Zak) : https://git.chem.ox.ac.uk/vld-group/theme-amorphous/c-mace-23/evolutionary_committee_models/-/blob/main/aG_data/structure_maps/make_structure_maps.ipynb
... and : https://git.chem.ox.ac.uk/vld-group/theme-amorphous/zif-umap-projection
Plotting convergence (Stephen Hoy) : https://git.chem.ox.ac.uk/vld-group/theme-amorphous/lipon-ace/-/blob/main/util/plot_convergence.py
losses plots NequIP / DFT plots (John) : https://git.chem.ox.ac.uk/vld-group/theme-ml/force-metrics/-/blob/main/visualise.ipynb
X-ray scattering and etc : https://git.chem.ox.ac.uk/vld-group/theme-amorphous/utilities/-/blob/main/examples.ipynb

TODO: add OVITO functionality and/or agent
'''

def rdf_analysis(atoms, start_idx=5, end_idx=-1):
    from ase.geometry.analysis import Analysis 
    geo_analysis_start = Analysis(atoms[start_idx])
    geo_analysis_end= Analysis(atoms[end_idx])

    rdf_start = geo_analysis_start.get_rdf(rmax=4, nbins=50, return_dists=True)[0]
    rdf_end = geo_analysis_end.get_rdf(rmax=4, nbins=50, return_dists=True)[0]
    
    # Plotting
    fig, axs = plt.subplots(1, 2, figsize=(12,4))

    axs[0].plot(rdf_start[1], rdf_start[0], label='liquid')
    axs[1].plot(rdf_end[1], rdf_end[0], label='amorphous')

    axs[0].set(xlabel='Distance (Å)', ylabel='g(r)', ylim=(-1, 5))
    axs[1].set(xlabel='Distance (Å)', ylim=(-1, 5))

    axs[0].legend(frameon=False, loc='upper right')
    axs[1].legend(frameon=False, loc='upper right')

functions.append({
    "name": "rdf_analysis",
    "description": "Analyse the radial distribution function of a trajectory",
    "parameters": {
        "atoms": {
            "type": "object",
            "description": "list of ASE atoms objects from a trajectory",
        }, "start_idx": {
            "type": "int",
            "description": "index of the starting structure",
        }, "end_idx": {
            "type": "int",
            "description": "index of the ending structure",
        },
    },
})

####################################################################################################
# MLP functions
####################################################################################################
'''
# nequip.ase.NequIPCalculator for DTF?
# NequIP central repo: https://github.com/mir-group/nequip
#
# TODO: $ nequip-train example.yaml --train-dir /path/to/training/session/ and document how to yaml
# TODO: $ nequip-evaluate --train-dir /path/to/training/session/
# TODO: nequip-deploy build --train-dir path/to/training/session/ where/to/put/deployed_model.pth
'''


functions.append({
    
})

####################################################################################################
# MD functions
####################################################################################################

def write_LAMMPS_in_file():
    rundir = 'lammps_runs/run_nvt'   # directory for the output files
    data_file = 'structures/random_1.data'  # random initial structure - density is determined by this
    random_seed = 2023        # random seed for lammps

    melt_timesteps = 1000     # Timesteps in fs for melting
    cool_timesteps = 1000     # Timesteps in fs for quenching (more = slower quench)
    melt_temperature = 9000   # in Kelvin
    final_temperature = 300   # in Kelvin

    # input file for LAMMPS
    lammps_nvt_input = f'''
    log {rundir}/log_nvt_C.dat append  # set output log file

    units metal  # set units for simulation (fs, K, eV, A etc.)
    atom_style atomic # set LAMMPS atom style (not molecule, or angle style etc.)

    read_data {data_file}  # read initial structure

    mass 1 12.011  # mass of carbon
    pair_style quip
    # carbon GAP potential from https://doi.org/10.1103/PhysRevB.95.094203
    pair_coeff * * potentials/carbon.xml \"\" 6   # arguments: apply to all types of atom with * *, potential file name, empty string, and atomic number

    neighbor 2.0 bin 
    neigh_modify every 1 delay 0 check yes   # build neighbour lists

    group carbon type 1  # just labelling
    timestep 0.001  # set timestep to 1 fs
    fix removeMomentum all momentum 1 linear 1 1 1 

    # setting these variables instructs lammps to calculate them so they can be used for fixes
    variable nAtoms equal atoms 
    compute T all temp 
    variable P equal press
    variable v equal vol
    compute MSD all msd
    compute pe_at all pe/atom
    compute PE all pe pair 
    variable PE_Atom equal c_PE/atoms

    # defined 'fixes' for calculating average properties over 100 fs here
    fix TempAve all ave/time 100 1 100 c_T 
    fix PressAve all ave/time 100 1 100 v_P 
    fix vAve all ave/time 100 1 100 v_v
    fix PEAve_Atom all ave/time 100 1 100 v_PE_Atom

    # Set up output files
    thermo_style custom step cpu temp f_TempAve press f_PressAve f_PEAve_Atom vol f_vAve c_MSD[4] 
    thermo_modify flush yes
    thermo 100
    dump traj all cfg 100 {rundir}/NVT/dump_nvt_C.*.cfg mass type xs ys zs id c_pe_at
    # ensure that carbon atoms are always ordered by id instead of random and add padding to filenames so they appear sorted
    dump_modify traj sort id element C pad 8
    # write binary restart files so you can restart the calculation if it crashes
    restart 100 {rundir}/restart_nvt_C.*

    # Set up NVT run

    # define variables for cumulative timesteps after which to switch temperatures
    variable Nrun1 equal {melt_timesteps}
    variable Nrun2 equal (${{Nrun1}}+{cool_timesteps})
    if "{melt_timesteps} > 0 && $(step) == 0" then "velocity all create {melt_temperature} {random_seed}"

    # (1) melting 
    run 0  # do pre-calcuations, but no actual timesteps
    fix integrate all nvt temp {melt_temperature} {melt_temperature} 0.1
    if "$(step) < ${{Nrun1}}" then "run ${{Nrun1}} upto"
    unfix integrate

    # (2) quenching
    fix integrate all nvt temp {melt_temperature} {final_temperature} 0.1
    if "$(step) < ${{Nrun2}}" then "run ${{Nrun2}} upto start ${{Nrun1}} stop ${{Nrun2}}"
    unfix integrate

    write_data {rundir}/out_data_quench_C.data # write last snapshot to lammps data file
    '''

    with open("lammps_runs/nvt.in", "w") as f:
        f.write(lammps_nvt_input)