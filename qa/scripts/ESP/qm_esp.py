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
   filein=open(filename_out,'r')
   lines=filein.readlines()
   for k in range(0,len(lines)):
    if len(lines[k].split())==4:
      charge.append(float(lines[k].split()[0]))
      coordinate.append(float(lines[k].split()[1]))
      coordinate.append(float(lines[k].split()[2]))
      coordinate.append(float(lines[k].split()[3]))                
   #filein2=open(filename_xyz,'r')
   #lines2=filein2.readlines()
   #for j in range(0,len(lines2)):
    #if len(lines2[j].split())==4:
      #charge.append(float(lines2[j].split()[0]))
      #coordinate.append(float(lines2[j].split()[1]))
      #coordinate.append(float(lines2[j].split()[2]))
      #coordinate.append(float(lines2[j].split()[3]))
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
    vector=coord_trunc-indexs_coord[0]
    vector_norm=LA.norm(vector, axis=1)
    vector_unit=vector/np.reshape(vector_norm,(-1,1))
    electric_strength=10**20/(vector_norm*vector_norm)*89.875517873681*1.6*10**(-19)*charge_trunc
    electric_potential=sum(10**10/(vector_norm)*8.9875517873681*10**(9)*1.6*10**(-19)*charge_trunc)*23.06*4.184
    vector2=coord_trunc-indexs_coord[1]
    vector_norm2=LA.norm(vector2, axis=1)
    electric_potential2=sum(10**10/(vector_norm2)*8.9875517873681*10**(9)*1.6*10**(-19)*charge_trunc)*23.06*4.184
    vector3=coord_trunc-indexs_coord[2]
    vector_norm3=LA.norm(vector3, axis=1)
    electric_potential3=sum(10**10/(vector_norm3)*8.9875517873681*10**(9)*1.6*10**(-19)*charge_trunc)*23.06*4.184    
    total=(np.reshape(electric_strength,(-1,1))*vector_unit).sum(axis=0)
    strength=LA.norm(total)
    direction=total/strength   
    return electric_potential,electric_potential2,electric_potential3,strength,direction

def main():
  # we are interested in the electric field at C, and the projection along C-S and C-C direction
  indexs=[1,1,1]
  # The way of estimating lines to skip is: first write down the natural indexes that are removed, like the second aand the third gives [2,3], and then reduce 1 from the 
#first number gives [1,3]
  skiplines=[0,1,1,1]
  label='pea_smd'
  filepath='/home/zhongyue/MTCT/PEA_run-detail/'
  #print "E_Potential_Zn, E_Potential_X1, E_Potential_X2, E_Field, E_field_Zn-X1, E_field_Zn-X2, charge_Zn, charge_X1, charge_X2"
  charge,coordinate=open_file('qm.xyz')
          #n_atoms,charge,coordinate=open_file(filepath+label+str(i)+'.out',filepath+'ptchrg'+label+str(i)+'.xyz')
  charge=np.array(charge)
  coordinate=np.array(coordinate)
  coordinate=np.reshape(coordinate, (-1, 3))
  indexs_coord,charge_sele=get_coordinate(indexs,coordinate,charge)
  EP,EP2,EP3,strength,direction=electric_field(skiplines,indexs_coord,charge,coordinate)
  #vector 1 point from S to C, and 2 from CH3 to Cprod
  #vector1=indexs_coord[0]-indexs_coord[1]
  #vector1_unit=vector1/LA.norm(vector1)
  #vector2=indexs_coord[2]-indexs_coord[1]
  #vector2_unit=vector2/LA.norm(vector2)
  #strength1=sum(strength*direction*vector1_unit)
  #strength2=sum(strength*direction*vector2_unit)
  #print EP,EP2,EP3,strength, strength1, strength2, charge_sele[0], charge_sele[1], charge_sele[2]
  print(EP)


	
if __name__ == "__main__":
  main()


