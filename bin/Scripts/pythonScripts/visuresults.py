import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from termcolor import colored
import sys

"""
    This python script is designed to generate a summary of the motion primitive evaluation performed by the bash script "test_main". It's
    essential purpose is to produce key performace indicators from the generated data during the evaluation process. Therefore, the data is
    stored in textfiles is converted into pandas own data structure DataFrame, which allows us to perform common data science operations on
    the collected data.

    To use this script as designed it has to be called with the following inputs:

        - the names of all tested motion primitive sets seperated by a space

"""

"""
    This first for-loop is used to cut the path to the motion primitive directories 
"""
    
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
a_star_succ = []

for i in range(nmbrOfPrimSets):
    resultfile = "Evaluation/result" + str(i+1)+ ".txt"
    df_list.append(pd.read_csv(resultfile,delim_whitespace=True,names=['dist', 'NoP', 'time', 'edges', 'vrtx', 'init_CV', 'init_obj', 'TWiteration', 'end_CV', 'end_obj', 'succ']))
    time_lst.append(df_list[i]['time'].sum()/float(df_list[i].shape[0]))
	
    avrg_init_CV.append(df_list[i]['init_CV'].sum()/float(df_list[i].shape[0]))
    avrg_init_obj.append(df_list[i]['init_obj'].sum()/float(df_list[i].shape[0]))
    avrg_TW_iter.append(df_list[i]['TWiteration'].sum()/float(df_list[i].shape[0]))

    median_init_CV.append(df_list[i]['init_CV'].median())
    median_init_obj.append(df_list[i]['init_obj'].median())
    median_TW_iter.append(df_list[i]['TWiteration'].median())

    nmbr_of_trajectories.append(df_list[i].shape[0])
    
    # computation of the successrates of A* and TW-reoptimization
    a_star_succ_df=df_list[i][df_list[i]['dist']>0.0]
    TW_succ    =a_star_succ_df[a_star_succ_df['succ']>0]
    a_star_succ.append(a_star_succ_df.shape[0]/float(df_list[i].shape[0]))
    succrate_lst.append(TW_succ.shape[0]/float(a_star_succ_df.shape[0]))
      
        

displayDict={"Summary": dirName_lst, "TW Successrate": succrate_lst, "Average time": time_lst, "Average Initial CV": avrg_init_CV, "Average Initial Obj": avrg_init_obj, "Average TW Iterations": avrg_TW_iter, "Median Initial CV": median_init_CV, "Median Initial Obj": median_init_obj, "Median TW Iterations": median_TW_iter, "Number of Trajectories": nmbr_of_trajectories, "A* Successrate": a_star_succ}

displayDict_header ={"Summary": dirName_lst}
displayDict_general={"Number of Trajectories": nmbr_of_trajectories}

displayDict_A={"Total time": time_lst, "A* Successrate": a_star_succ}

displayDict_TW={"TW Successrate": succrate_lst,"Average Initial CV": avrg_init_CV, "Average Initial Obj": avrg_init_obj, "Average TW Iterations": avrg_TW_iter, "Median Initial CV": median_init_CV, "Median Initial Obj": median_init_obj, "Median TW Iterations": median_TW_iter}

def set_output_witdh(lenOfEntry,dpDict):
    opener = "# "
    closer = " #"
    max_width = 0 
    for i in range(len(dpDict)):
        val_str = " | ".join(map(lambda x: str(x).rjust(lenOfEntry),dpDict.values()[i]))
        line_len = len(opener + dpDict.keys()[i]  + val_str + closer)
        if max_width< line_len:
            max_width = line_len
    return max_width
            
def output_format(lenOfoutput,dpDict, is_not_header):
    lenOfValues = 14
    opener = "# "
    closer = " #"
    for i in range(len(dpDict)):
        if is_not_header:
            print opener + "-"*(lenOfoutput-len(opener+closer)) +closer
        val_str = " | ".join(map(lambda x: str(x).rjust(lenOfValues),dpDict.values()[i]))            
        nmbrSpaces = (lenOfoutput-len(opener)-len(closer)-len(dpDict.keys()[i])-len(val_str))
            
        print opener + dpDict.keys()[i]  +  " "*nmbrSpaces+ val_str + closer
    
output_width = max(set_output_witdh(14,displayDict_header),
                   set_output_witdh(14,displayDict_general),
                   set_output_witdh(14,displayDict_A),
                   set_output_witdh(14,displayDict_TW))+10

print "#"*output_width
output_format(output_width,displayDict_header,False)
output_format(output_width,displayDict_general,True)
output_format(output_width,displayDict_A,True)
output_format(output_width,displayDict_TW,True)
print "#"*output_width

#print displayDict.values()
#print type(displayDict.values())
#with open("Evaluation/ai_file1.txt", 'w') as f:
#    f.write(" ".join(map(lambda x:str(x).replace('[','').replace(']','').replace('\'',''),displayDict.values())))
#    f.close()

    

