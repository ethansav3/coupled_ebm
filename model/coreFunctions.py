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

notePath = os.getcwd()

def analyzeRun(dfModel,nameList,verbose):
    counter = 0
    maxima = dfModel.max();#find maxima from all columns in df
    maxPop = maxima[3];#find maxima in population column, peak popultion
    #print(maxPop)#peak population
    maxPopIndex = dfModel.loc[dfModel['pop']==maxPop].index#search rows for index of max pop
    maxPopTime=dfModel.iloc[maxPopIndex]['time_yrs'];#search column for time until peak pop is reached

    halfPop=maxPop/2;
    halfPopIndex = dfModel.loc[dfModel['pop']==halfPop].index#search rows for index of max pop
    halfPopTime=dfModel.iloc[maxPopIndex]['time_yrs'];#search column for time until peak pop is reached
    #while newDF.shape[0]>=2:
    #     newDF= dfModel.loc[(dfModel['pop'] > (halfPop-dP)) & (dfModel['pop'] < (halfPop+dP))]
    #     dT+=1
    #     print("dT: "+ str(dT))
    LhalfPop=0;
    UhalfPop=0;
    dP=50
    while True:
        counter += 1
        newDF= dfModel.loc[(dfModel['pop'] > (halfPop-dP)) & (dfModel['pop'] < (halfPop+dP))]
        if(newDF.shape[0]>=2):
            break
        if counter > 50:
            break
        dP+=50
#        print(str(dP)+": "+str(newDF.shape[0]))
#        print(newDF)
    try:
        if(nameList['ebm']['coupled']):
            LhalfPop = newDF.time_yrs.iloc[0]  #new dataframe, time column, first index
            UhalfPop = newDF.time_yrs.iloc[-1]      #new dataframe, time column, last row
    except TypeError:
        print('')

    #maxPopPlot=40

    #dictionary of population statistics
    popStats={'maxPop' : (maxPop/1000), 'maxTime': maxPopTime.mean(),
              'halfPop': (halfPop/1000),'LhalfTime': LhalfPop, 'UhalfTime': UhalfPop,
              'maxPopPlot': maxPop/1000}#maximum plotting range for population
    if(verbose):
        for k,v in popStats.items():
            print(k + " = " +str(v))
    return popStats

def habitableZone(nameList,newPco2,newA,runTime,dA):
    fullMaxPop=0;
    fullMaxPopA=0;
    coupled=True
    count=0
    minA=0;
    maxA=0;
    maxPop=0;
    while(True):
        relsolcon=newA**-2
        nameList['ebm']['coupled']=coupled
        nameList['ebm']['pco20']=newPco2/10**6#convert pco2 to bars
        nameList['ebm']['relsolcon']=newA**-2 #inverse square law for solar flux
        nameList['ebm']['runTime'] = runTime#change runtime
        dfModel, finalavgtemp, eqTime, eqTemp, equilibrium = runProgram("driver.exe",nameList,False)#False=no output
        life = (equilibrium) and (eqTemp<=373.15) and (eqTemp>=273.15)#determine habitability
        if(life and (count==0)):#if first habitable distance, make it minA
            minA = newA
            count+=1
        if(life and (count > 0)):#in the habitable zone, make it maxA
            if coupled:
                popStats = analyzeRun(dfModel,nameList,False)
                maxPop = popStats["maxPop"]
                if(maxPop >= fullMaxPop):
                    fullMaxPop = maxPop
                    fullMaxPopA = newA
            maxA = newA
        if((not life) and (count>0)):#if out of habitable zone, break out of loop
            break
        newA += dA
    return round(minA,2), round(maxA,2), fullMaxPop, fullMaxPopA

def printFolder():
    for root, dirs, files in os.walk('.'):
        for filename in files:
            print(filename)
        for dirname in dirs:
            print(dirname)
            
def runProgram(driver,nameList,output): #run the program with the given name  
    #make temporary directory to run in
    with tempfile.TemporaryDirectory() as dirpath:
        runFolder = newFolder(nameList,dirpath) #make the temporary folder
        
        call("./"+driver,shell=True) #run the driver
        
        dfModel, finalavgtemp, eqTime, eqTemp, equilibrium = readOutput() #read output into datframe
        
        os.chdir(notePath)
        call("rm -rf tmp*", shell=True)#delete the temporary folder and unlink it's contents
        if(output):
            if(equilibrium):
                print('Equilibrium Reached at Temp=' + str(eqTemp)+". At time="+str(eqTime)) 
                print('Final Temp(K): ' + str(finalavgtemp));
                print('Final Temp(F): ' + str(round((finalavgtemp-273.15)*(9/5)+32, 2)));
                print('')
            else:
                print("Equilibrium was not reached")
                print('')
     #   print('Final Temp(C): ' + str(round(finalavgtemp-273.15)));
        call("echo   ", shell=True)
    return dfModel, finalavgtemp, eqTime, eqTemp, equilibrium
            
def makeDefNamelist():
    nml = {
        'ebm': {
            'seasons' : True,
            'snowball' : False,
            'tend' : 3.154e7, #calculation length (s)
            'dt' : 8.64e4, #time step (s) default = 1 day
            'rot' : 7.27e-5, #earth's present rotation rate
            'a' : 1.49597892E13, #earths orbital semi-major axis (cm)
            'ecc' : 0.0167, #earths orbital eccentricity 
            'peri' : 76.25, #longitude of perihelion wrt vernal equinox (degrees)
            'obl' : 23.5, #earth's obliquity
            'cloudir' : -9.5, #reduction of outgoing infrared by clouds
            'pco20' : 2.84e-4,
            'ocean' : 0.7,
            'igeog' : 1,
            'groundalb' : 0.291,
            'relsolcon' : 1.0,
            'landsnowfrac' : 1.0,
            'fcloud' : 0.50,
            'd0' : 0.58, #initial diffusion coefficient
            'N0' : 1129,
            'Nmax' : 10200,
            'rBirth0' : 0.02,
            'rBirthMax' : 0.83,
            'rDeath0' : 0.015,
            'pBirth' : 1.00,
            'pDeath' : 1.00,
            'rTech' : 1.00,
            'dTpop' : 10,
            'opT' : 290.5,
            'rco2' : 1.9e-4,#3.459e-4,
            'En0' : 1.00,
            'fragility' : 1.00,
            'coupled' : True,
            'lverbose' : True,
            'runTime' : 100
        }
    }
    return nml

def newFolder(nml,dirpath): #delete old runfile and make a new one    
    write(nml,'../input.nml',force=True, sort=True) #write the updated namelist
    
    #symlink all relevant files into the temporary directory
    call("ln -s "+notePath+"/driver"+dirpath,shell=True) 
    call("ln -s "+notePath+"/output.dat "+dirpath,shell=True)
    call("ln -s "+notePath+"../input.nml "+dirpath,shell=True)
    return dirpath

def deleteFolder():
    #unlink all symbolic links
    call("unlink driver",shell=True)
    call("unlink output.dat",shell=True)
    call("unlink ../input.nml",shell=True)

    os.chdir(os.pardir)    #move one directory up, out of the temporary folder
    
    call("rm -rf tmp*",shell=True) #then remove the temporary folder


def readOutput():
    #instantiate output object
    data = {} #dictionary containing all critical output
    data['time'] = []
    data['temp'] = []
    data['pco2'] = []
    data['pop']  = []
    data['rBirth'] = []
    data['rDeath'] = []
    data['rGrowth'] = []
    finalavgtemp=0;
    equilibrium = None 

    if(os.path.getsize('output.dat')==3):
        equilibrium = False
    else:
        equilibrium = True
        
    #read in population/avgtemp data
    output = open("output.dat","r")
    next(output) #skip the first line (of headers)
    for line in output: #iterates as many years as the program runs
        values = line.split()
        data['time'].append(float(values[0]))
        data['temp'].append(float(values[1]))
        data['pco2'].append(float(values[2]))
        data['pop'].append(float(values[3]))
        data['rBirth'].append(float(values[4]))
        data['rDeath'].append(float(values[5]))
        data['rGrowth'].append(float(values[6]))
   
    try:
        finalavgtemp=data['temp'][len(data['temp'])-1] # determine the final average tem
        equilibrium=True
    except IndexError:
        finalavgtemp=np.NaN
#        print("Temperatures Exceeded 450 before equilibrium was reached")
        equilibrium = False
    
    if(equilibrium):
        #read in equilibrium conditions
        eqTemp=0;
        eqTime=0;
        equilibrium=open("equilibrium.dat","r")
        for line in equilibrium:
            values = line.split()
            eqTime = float(values[0])
            eqTemp = float(values[1])
    else:
        eqTemp = np.NaN
        eqTime = np.NaN
           
    output.close() # close output file
       
    df = pd.DataFrame(data)
 
    df['time_yrs'] = df['time']/60/60/24/365.25
    df['pco2_ppm'] = df['pco2']*10**6
 
    return df, finalavgtemp, eqTime,eqTemp,equilibrium
