import tempfile
import os
from subprocess import call
import matplotlib.pyplot as plt
from f90nml import *
from array import *
import numpy as np
from math import *
import pandas as pd 
import seaborn as sns
from statistics import *
import shutil
import sys
from matplotlib.ticker import FormatStrFormatter


def readOutput():
    #instantiate output object
    data = {} #dictionary containing all critical output
    data['time'] = []
    data['temp'] = []
    data['pco2'] = []
    data['pop']  = []
    finalavgtemp=0;
    
    #read in population/avgtemp data
    output = open("output.dat","r")
    next(output) #skip the first line (of headers)
    
    for line in output: #iterates as many years as the program runs
        values = line.split()
        data['time'].append(float(values[0]))
        data['temp'].append(float(values[1]))
        data['pco2'].append(float(values[2]))
        data['pop'].append(float(values[3]))
    
    finalavgtemp=data['temp'][len(data['temp'])-1] # determine the final average temp
    
    output.close() # close output file
    
    df = pd.DataFrame(data)
    
    return df, finalavgtemp