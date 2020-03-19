import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from termcolor import colored
import sys
import os

class Colored_string:
    def __init__(self, name, color):
        self.length  = len(name)
        self.name = colored(name,color)
        self.clr  = color
    
dirName_lst = []

for i in range(len(sys.argv)):
    if sys.argv[i] is not None:
        primdir= sys.argv[i]
        prim_lst= primdir.split('/')
        for elem in prim_lst:
            if elem.find("motion_prims")>-1:
                dirName_lst.append(elem)
    else:
        dirName_lst.append("unspecified")

nmbrOfPrimSets = len(sys.argv)-1
lenOfEntry = 20
df_list = []
succrate_lst= []
time_div_edges_lst = []
time_lst = []
avrg_init_CV =[]
avrg_init_obj =[]
avrg_TW_iter = []
median_init_CV =[]
median_init_obj =[]
median_TW_iter = []
nmbr_of_trajectories = []

for i in range(nmbrOfPrimSets):
    resultfile = "Evaluation/result" + str(i+1)+ ".txt"
    df_list.append(pd.read_csv(resultfile,delim_whitespace=True,names=['dist', 'NoP', 'time', 'edges', 'vrtx', 'init_CV', 'init_obj', 'TWiteration', 'end_CV', 'end_obj', 'succ']))
    succrate = df_list[i]['succ'].sum()/float(df_list[i].shape[0])
    succrate_lst.append(str(succrate))
    time_lst.append(df_list[i]['time'].sum())
	
    avrg_init_CV.append(df_list[i]['init_CV'].sum()/float(df_list[i].shape[0]))
    avrg_init_obj.append(df_list[i]['init_obj'].sum()/float(df_list[i].shape[0]))
    avrg_TW_iter.append(df_list[i]['TWiteration'].sum()/float(df_list[i].shape[0]))

    median_init_CV.append(df_list[i]['init_CV'].median())
    median_init_obj.append(df_list[i]['init_obj'].median())
    median_TW_iter.append(df_list[i]['TWiteration'].median())

    nmbr_of_trajectories.append(df_list[i].shape[0])

displayDict={"Summary": dirName_lst, "Successrate": succrate_lst, "Total time": time_lst, "Average Initial CV": avrg_init_CV, "Average Initial Obj": avrg_init_obj, "Average TW Iterations": avrg_TW_iter, "Median Initial CV": median_init_CV, "Median Initial Obj": median_init_obj, "Median TW Iterations": median_TW_iter, "Number of Planning Problems": nmbr_of_trajectories}

def set_output_witdh(lenOfEntry):
    opener = "# "
    closer = " #"
    val_str = " | ".join(map(lambda x: str(x).rjust(lenOfEntry),displayDict['Summary']))
    max_width = len(val_str)    
    for i in range(len(displayDict)):
        if displayDict.keys()[i] is not "Summary":
            val_str = " | ".join(map(lambda x: str(x).rjust(lenOfEntry),displayDict.values()[i]))            
            if max_width<len(val_str):
                max_width = len(val_str)
    return max_width
            
def output_format(lenOfoutput):
    print "#"*lenOfoutput
    opener = "# "
    closer = " #"
    val_str = " | ".join(map(lambda x: str(x).rjust(lenOfEntry),displayDict['Summary']))
    nmbrSpaces = (lenOfoutput-len(opener)-len(closer)-len('Summary')-len(val_str))
    print opener + 'Summary' +  " "*nmbrSpaces+ val_str + closer
    for i in range(len(displayDict)):
        if displayDict.keys()[i] is not "Summary":
            print opener + "-"*(lenOfoutput-4) +closer
            nmbrSpaces = (lenOfoutput-len(opener)-len(closer)-len(displayDict.keys()[i])-len(val_str))
            val_str = " | ".join(map(lambda x: str(x).rjust(lenOfEntry),displayDict.values()[i]))            
            
            print opener + displayDict.keys()[i]  +  " "*nmbrSpaces+ val_str + closer
    print "#"*lenOfoutput


output_format(set_output_witdh(14))


