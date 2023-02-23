# script for extracting parameters from multiwfn
# (of metal: mulliken spin, mulliken charge, total valence, free valence)
# currently specifically made for Fe, but that can be generalized
​
import shutil
import subprocess
import sys
import os
import numpy as np
​
owd = os.getcwd() # old working directory
​
list_of_file = ['YIDLOP'] # populate with names of folders of DFT calculations that were sent to job manager
dict_of_calcs = {'Hirshfeld': '1', 'Voronoi':'2', 'Mulliken': '5', 'ADCH': '11'}
​
temp_dict = {} # store extracted multiwfn parameters
​
for f in list_of_file:
    os.chdir(f + '/scr')
    subprocess.call("module load multiwfn/noGUI", shell=True)
    command_Z = "module load multiwfn/noGUI"
    command_A = '/opt/Multiwfn_3.7_bin_Linux_noGUI/Multiwfn '+f+ '.molden'
    final_dest = '/home/manets12/AuNanocage/solvated_guests/charge_files/'
    for key in dict_of_calcs:
        print('Current Key: '+key)
        proc = subprocess.Popen(command_A, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        calc_command = dict_of_calcs[key]
        commands = ['7', calc_command, '1', 'y', '0', 'q'] # for atomic charge type corresponding to dict key
        output = proc.communicate("\n".join(commands).encode())
        lines = str(output[0]).split('\\n')
        start = False
        counter = 0
        for num, line in enumerate(lines):
            print(line)
            if ('Total valences and free valences' in line):
                start = True
                continue
            if start and "Fe" in line:
                temp_dict[d] = {"Total valence": float(line.split()[-2]),
                            "Free valence": float(line.split()[-1])}
                counter += 1
        #rename the file as the new, desired name
        new_name = f+'_' +key+'.txt'
        os.rename(f+'.chg', new_name)
        shutil.copy(new_name, final_dest+new_name)   
    os.chdir(owd)
​
​
    #os.chdir('pdb'+d+'/scr')
    #subprocess.call("module load multiwfn/noGUI", shell=True)
    #proc = subprocess.Popen("Multiwfn *molden", stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    #commands = ['7', '5', '1', 'y', 'n']
    #output = proc.communicate("\n".join(commands).encode())
    #lines = str(output[0]).split('\\n')
    #start = False
    #counter = 0
    #for num, line in enumerate(lines):
    #    if ('Population of atoms' in line):
    #        start = True
    #        continue
    #    if start and "Fe" in line:
    #        temp_dict[d]["Mulliken spin"] = float(line.split()[-2])
    #        temp_dict[d]["Mulliken charge"] = float(line.split()[-1])
    #        counter += 1
    #os.chdir(owd)
#prin#t(temp_dict)