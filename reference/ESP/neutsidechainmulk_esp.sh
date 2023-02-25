import sys
import os
import glob
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from numpy import linalg as LA
from scipy.stats import norm


def open_file(filename_out):
   filein=open(filename_out,'r')
   lines=filein.readlines()
   charge=[]
   coordinate=[]
   for k in range(0,len(lines)):
    atomn=lines[k].split()[2]
    resi=lines[k].split()[4]
    resn=lines[k].split()[3]
    if atomn not in ['N', 'H', 'HN', 'CA', 'HA', 'C', 'O']:
		#if resn in ['ARG', 'LYS', 'ZN1']:
        if resn in ['THR', 'CYS', 'SER', 'ASN', 'GLN', 'GLY', 'TYR', 'MO1']:
		#if resn not in ['ARG', 'LYS', 'CY1', 'GLU']:
         charge.append(float(lines[k].split()[8]))
         coordinate.append(float(lines[k].split()[5]))
         coordinate.append(float(lines[k].split()[6]))
         coordinate.append(float(lines[k].split()[7]))                
   return charge,coordinate

def get_coordinate(indexs,coordinate,charge):
   sele_coord=[]
   sele_charge=[]
   for i in range(0,len(indexs)): 
     sele_coord.append(coordinate[indexs[i]-1])
     sele_charge.append(charge[indexs[i]-1])     
   return sele_coord,sele_charge

def electric_field(skiplines,indexs_coord,charge,coordinate):
    coord_trunc=np.concatenate((coordinate[:skiplines[0]],coordinate[skiplines[1]:skiplines[2]],coordinate[skiplines[3]:]),axis=0) 
    charge_trunc=np.concatenate((charge[:skiplines[0]],charge[skiplines[1]:skiplines[2]],charge[skiplines[3]:]),axis=0)
   # print indexs_coord[0]
    vector=coord_trunc-indexs_coord[0]
    vector_norm=LA.norm(vector, axis=1)
    electric_potential=sum(10**10/(vector_norm)*8.9875517873681*10**(9)*1.6*10**(-19)*charge_trunc)*23.06*4.184
    return electric_potential

def main():
  # we are interested in the electric field at C, and the projection along C-S and C-C direction
  indexs=[1,1,1]
  # The way of estimating lines to skip is: first write down the natural indexes that are removed, like the second aand the third gives [2,3], and then reduce 1 from the 
#first number gives [1,3]
  skiplines=[0,1,1,1]
  #print "E_Potential_Zn, E_Potential_X1, E_Potential_X2, E_Field, E_field_Zn-X1, E_field_Zn-X2, charge_Zn, charge_X1, charge_X2"
  charge,coordinate=open_file('mask_qm.xyz')
          #n_atoms,charge,coordinate=open_file(filepath+label+str(i)+'.out',filepath+'ptchrg'+label+str(i)+'.xyz')
  charge=np.array(charge)
  coordinate=np.array(coordinate)
  coordinate=np.reshape(coordinate, (-1, 3))
  indexs_coord,charge_sele=get_coordinate(indexs,coordinate,charge)
  EP=electric_field(skiplines,indexs_coord,charge,coordinate)
  print (EP)


	
if __name__ == "__main__":
  main()



