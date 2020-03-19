import numpy as np
from shutil import copyfile
import os
import glob
import sys

"""
    This script is designed to automate the process of parsing transWORHP-results into the format that is needed by the a*-search algorithm.
    The process of parsing involves these steps, these are:
        - creating the correct file structure to perform autonomous test of motion primitive sets.
        - creating the necessary file name for a maneuver as described in primpathsearch.cpp -> NEW_PrimGridPathSearch3d::SetPrimitives
        - replacing the default header of trajectory files created by transWORHP by the header needed for the a* search 
"""

#output file from transWORHP
src = 'temp1.m'

folder_id = str(sys.argv[1])



#the first block used to determine the time scale factor and to determine whether it is a maneuver (n_discr>1) or a trim (n_discr ==1)
f = open(src, 'r')

line_counter = 0
for line in f:
    # the second line of a transWORHP file contains the time scale factor
    if line_counter == 1:
        time_scale = float(filter(None, line.split(" "))[1])
    line_counter = line_counter +1
n_discr = line_counter -2
f.close()

####################

#this block prepares the needed structure of a motion primitive, which is in case of a maneuver man_#id_XX_XX_describtion.txt where XX_XX contains which trims are connected, in our setting this is not used currently.
dst = ""
describtion = "describtion"
path2prims = "./primDirs/motion_prims_"+folder_id+"/"

if not os.path.isdir(path2prims):
    os.makedirs(path2prims)
    print ("The directory: " + path2prims + " has been created.")

dst = path2prims+"placeholder.txt"
copyfile(src, dst)
f1 = open(dst, 'r+')
line_counter = 0
cost = "0"
lines = f1.read().splitlines()
f1.close()

################

# this block goes through the directory and counts the number of maneuvers and trims, this is needed to create the correct primitive file

list_of_man = glob.glob(path2prims+"man_"+"*.txt")
list_of_man_id = []
for maneuver in list_of_man:
    list_of_man_id.append(int(maneuver[len(path2prims)+4:len(path2prims)+6]))

if len(list_of_man_id)>0:
    man_nmbr = max(list_of_man_id)+1
else:
    man_nmbr = 1

list_of_tp = glob.glob(path2prims+"tp_"+"*.txt")
list_of_tp_id = []
for trim in list_of_tp:
    list_of_tp_id.append(int(trim[len(path2prims)+3:len(path2prims)+5]))

if len(list_of_tp_id)>0:
    tp_nmbr = max(list_of_tp_id)+1
else:
    tp_nmbr = 1

################

#this block manipulates the lines of src. The first line is set to the prim id and a following zero. The second line contains the number od discretisation points and for all following lines the time collum is multiplied with the time scale factor. finally these lines are written in the placeholer file.

f1 = open(dst, 'w')
for i in range(len(lines)):
    if i == 0:
        if n_discr>0:
            lines[i] = str(man_nmbr).zfill(2) + " " + str(0)   
        else:
             lines[i] = str(tp_nmbr).zfill(2) + " " + str(0)   
    if i == 1:
        lines[i] = str(n_discr) + " " + cost
    if i > 1:
        line_list = filter(None, lines[i].split(" "))
        time = "{:0.16f}".format(float(line_list[0])*time_scale)
        line_list[0]  = "{:0.16f}".format(float(line_list[0])*time_scale)
        lines[i] = "     ".join(line_list)

for item in lines:
  f1.write("%s\n" % item)

f1.close()

################################

# in this final block  we modify the name of the file and store all the needed information in the name. The distiguish between tp and man follows from n_discr. assuming that a maneuver realzies a turning amneuver, the final heading is stored in the name of a maneuver as well. Finally placeholder is replaced by the correct name and the path to the file is printed to the console.


if n_discr == 1:
    dst2 = path2prims + "tp_"+str(nmbr).zfill(2)+"_"+describtion+".txt"
else:
    if np.sign(np.round(float(line_list[3])*180/np.pi))==-1:
        describtion = str(int(np.abs(np.round(float(line_list[3])*180/np.pi)))) + "_deg_" + "m"
    else:
        describtion = str(int(np.abs(np.round(float(line_list[3])*180/np.pi)))) + "_deg_" + "p"
    dst2 = path2prims + "man_"+str(man_nmbr).zfill(2)+"_XX_XX_"+describtion+".txt"
os.rename(dst, dst2)

print("The file: " + dst2 + " has been created.")

