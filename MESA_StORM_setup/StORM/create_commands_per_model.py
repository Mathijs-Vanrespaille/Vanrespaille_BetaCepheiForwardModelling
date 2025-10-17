import numpy as np
import h5py
import os


def split_line(line, sep) :
    head, sep_, tail = line.partition(sep)
    assert sep_ == sep
    return head, tail


def substring(line, sep_first, sep_second) :
    head, tail = split_line(line, sep = sep_first)
    if sep_second not in tail:
        return tail
    head, tail = split_line(tail, sep =sep_second)
    return head


def compute_f_K(histfile) :
    G = 6.6743e-11
    Msol = 1.988416e30
    Rsol = 6.96340e8
    header = np.loadtxt(histfile, skiprows=5, max_rows=1, dtype=str)
    data = np.loadtxt(histfile, skiprows=6, unpack=True)
    Xc = data[header=="center_h1"][0]
    mns = np.array(data[header=="model_number"][0], dtype=int)
    mass = data[header=="star_mass"][0]*Msol
    radius = 10**data[header=="log_R"][0]*Rsol
    all_f_K = np.sqrt(G*mass/radius**3)*86400/(2*np.pi)
    return mns, all_f_K


def convert_MESA_to_GSM(MESA_file, GSM_file) :
    header = np.loadtxt(MESA_file, max_rows=1)
    data = np.loadtxt(MESA_file, usecols=(1,2,4,6,8,9,18), skiprows=1, unpack=True)
    # CHECKING FOR DOUBLE POINTS
    r = data[0]
    I_toKeep = [0]
    for i in range(1, len(r)) :
        if not r[i-1] == r[i] :
            I_toKeep.append(i)
    M_r, rho, P, Gamma_1, N2, Omega_rot = data[1][I_toKeep], data[3][I_toKeep], data[2][I_toKeep], data[5][I_toKeep], data[4][I_toKeep], data[6][I_toKeep]
    r = r[I_toKeep]
    n = len(r)
    with h5py.File(GSM_file, "w") as f :
        f.create_dataset("r", (n,), dtype="f", data=r)
        f.create_dataset("M_r", (n,), dtype="f", data=M_r)
        f.create_dataset("rho", (n,), dtype="f", data=rho)
        f.create_dataset("P", (n,), dtype="f", data=P)
        f.create_dataset("Gamma_1", (n,), dtype="f", data=Gamma_1)
        f.create_dataset("N2", (n,), dtype="f", data=N2)
        f.create_dataset("Omega_rot", (n,), dtype="f", data=Omega_rot)
        f.attrs["n"] = n
        f.attrs["M_star"] = header[1]
        f.attrs["R_star"] = header[2]


def main():

    from pathlib import Path
    import os
    import sys
    frot_file = sys.argv[1] # As a fraction of the Keplerian critical rotation rate
    MESA_folder = sys.argv[2] # directory of a particular MESA run

    f_rots = np.loadtxt(frot_file)

    f_min = 2. # Minimum frequency in cycles per day
    f_max = 15.
    f_density = 40. # Density of modes per cycle per day
    ls = [0,1,1,1,2,2,2,2,2]
    ms = [0,-1,0,1,-2,-1,0,1,2]

    this_dir = f'{MESA_folder}/storm'
    MESA_dir = f'{MESA_folder}/gyre/input_models' # directory containing all the .GYRE file. Trying to change this to StORM after already having converted everything
    MESA_history = f'{MESA_folder}/{os.path.basename(MESA_folder)}.hist'
    mns, all_f_K = compute_f_K(MESA_history)

    if not os.path.isdir(f'{this_dir}/commands/') :
        os.mkdir(f'{this_dir}/commands/')
    if not os.path.isdir(f'{this_dir}/output/') :
        os.mkdir(f'{this_dir}/output/')
    if not os.path.isdir(f'{this_dir}/input_models/') :
        os.mkdir(f'{this_dir}/input_models/')

    partial_commands = f'{this_dir}/commands_template.txt'
    for MESA_file in os.listdir(MESA_dir) :
        all_lines = []

        mn = int(substring(MESA_file, "_mn", ".GYRE"))
        f_K = all_f_K[mns==mn][0]

        Xc = substring(MESA_file, "c", "_mn")
        completed_commands = f'{this_dir}/commands/commands_Xc{Xc}_mn{mn}.txt'

        # Check if the MESA file has already been converted into a GSM. Convert if not.
        GSM_file=f'{this_dir}/input_models/{MESA_file[:-5]}.HDF5'
        print(f"converting {MESA_dir}/{MESA_file}")
        convert_MESA_to_GSM(f"{MESA_dir}/{MESA_file}", GSM_file)

        # Open up the template to be completed. The first line of this template should be reading in the input
        with open(partial_commands, 'r') as f :
            lines = f.readlines()
        first_line = lines[0]
        all_lines.append(first_line.replace('INFILE', f"\'{GSM_file}\'"))

        for f_rot in f_rots :
            replacements = {'OUTFILE' : f"\'{this_dir}/output/storm_Xc{Xc}_mn{mn}_frot{f_rot:.2f}.HDF\'".format(f_rot=f_rot),
                            'FROT' : f"{f_rot}"}
            for line in lines[1:] :
                new_line = line
                for key in replacements:
                    if (replacements[key] != '') :
                        new_line = new_line.replace(key, replacements[key])
                all_lines.append(new_line)
            all_lines.append("\n")

        with open(completed_commands, 'w') as f :
            f.writelines(all_lines)

if __name__ == "__main__":
    main()
