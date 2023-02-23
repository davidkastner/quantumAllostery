import pandas as pd

df = pd.read_csv("gas_phaseYIDLOP_Voronoi.txt", sep='\s+', names=["Atom",'x', 'y', 'z', "charge"])
k = 8.987551*(10**9)  #Coulombic constant in kg*m**3/(s**4*A**2)

#convert each column to list for quicker indexing
atoms = df['Atom']
charges = df['charge']
xs = df['x']
ys = df['y']
zs = df['z']

#pick the index of the atom at which the esp should be calculated
idx_atom = 0

#determine position and charge of the target atom
xo = xs[idx_atom]
yo = ys[idx_atom]
zo = zs[idx_atom]
chargeo = charges[0]
total_esp = 0

#unit conversion
A_to_m = 10**(-10)
KJ_J = 10**-3
faraday = 23.06   #kcal/(mol*V)
C_e = 1.6023*(10**-19)
one_mol = 6.02*(10**23)
cal_J = 4.184

for idx in range(0, len(atoms)):
    if idx == idx_atom:
        continue
    else:
        #Calculate esp and convert to units (A to m)
        r = (((xs[idx] - xo)*A_to_m)**2 + ((ys[idx] - yo)*A_to_m)**2 + ((zs[idx] - zo)*A_to_m)**2)**(0.5)
        total_esp = total_esp + (charges[idx]/r)

final_esp = k*total_esp*((C_e))*cal_J*faraday   #note that cal/kcal * kJ/J gives 1
print(str(final_esp) + ' kJ/(mol*e)')
    
