import tempfile
import os
import math
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
from plottingFunctions import *

notePath = os.getcwd()
fullArr=[]
fullArrTr = []
fullMaxPop=0
def runModel(nameList, coupled, runTime, plot, save, analyze, driverName,maxPopList,distList,showInputs):
    saveName= 0
    count   = 0
    numCols = 4
    numRows = 4
    tempData = []
    pco2Data = []
#    allData = {}
#    allData['pco2'] = []
#    allData['temp'] = []
#    for j in np.linspace(10,100,numCols):
    for j in [nameList['ebm']['Nmax']]:#pop loops
        newMaxPop=int(j)#change max pop
#        newA = .9#change distance (AU)
        dA=.01
        popStats = {}
#        minA, maxA, fullMaxPop, fullMaxPopA = habitableZone(nameList,newMaxPop,newA,runTime,dA)
#        print("The Habitable Zone for maxPop="+str(round(newMaxPop/1000)) + " billion is:")
#        print("Minimum: "+str(minA) + " AU")
#        print("Maximum: "+str(maxA)+" AU")
#        print("Max Pop Reached: " + str(fullMaxPop)+ " billion" + ", at Distance: " + str(round(fullMaxPopA,3))+"\n")        
#        for i in np.linspace(minA,maxA,numRows):#distance loops
#            count += 1
#            print(str( (count/(numCols*numRows))*100 ) + "% Until Completion" )
        dictdTdP= {
        0.94 : 7.967*10**-2,
        0.97 : 4.94*10**-3,
        1.00 : 3.378*10**-3,
        1.03 : 1.959*10**-2,
        1.06 : 6.554*10**-4
        }
        for i in [nameList['ebm']['relsolcon']**-(1/2)]:#distance 
#            dimVar = (nameList['ebm']['rBirth0']*nameList['ebm']['dtemp'])/(newMaxPop*nameList['ebm']['rco2']*dictdTdP[i])
            #dimVar=(newMaxPop*nameList['ebm']['rco2']*dictdTdP[i])/(nameList['ebm']['rBirth0']*nameList['ebm']['dtemp'])
            dimVar=0
            if(nameList['ebm']['lverbose']==True):
                print("rBirth0: " + str(nameList['ebm']['rBirth0']))
                print("dT: " + str(nameList['ebm']['dtemp']))
                print("rco2: " + str(nameList['ebm']['rco2']))
                print("maxPop: " + str(newMaxPop))
                #print("dTdP: " +str(dictdTdP[i]))
                print("Dimensionless Variable: " + str(dimVar))

            #calculate inputs
            newA = round(i,3) #change distance (AU)
            nameList['ebm']['coupled']=coupled
            nameList['ebm']['relsolcon']=newA**-2 #inverse square law for solar flux
            nameList['ebm']['runTime'] = runTime#change runtime
            nameList['ebm']['Nmax'] = newMaxPop #change carrying capacity
            dtemp = nameList['ebm']['dtemp'] 
            #run program
            dfModel, finalavgtemp, eqTime, eqTemp, equilibrium = runProgram(driverName,nameList,True,showInputs)#False=no output
            pco2Data += dfModel['pco2'].tolist()
            tempData += dfModel['temp'].tolist()
            #plot the results
            if coupled and equilibrium:
                inputRunTime = runTime #how long I originally specified
                outputRunTime = math.ceil( dfModel['time_yrs'].iloc[-1] )#how long it actually ran
                if ((inputRunTime - outputRunTime) <= 3) and analyze:
                    popStats = analyzeRun(dfModel,nameList, False)#if True, then Print Dictionary Values
                    fullMaxPop=popStats["maxPop"]
                    popStats['maxPopPlot']=fullMaxPop+(3/100)*fullMaxPop#maximum population range
                inputs=[newA,newMaxPop,runTime,dtemp]
                if equilibrium and plot:
                    plotModelOutput(dfModel,inputs,eqTime,eqTemp,popStats,save,saveName,dimVar)#plot the output of our model, colored by pco2 
                if save: saveName+=1
                #partitition dataframe into numpy arrays
                popArr=np.asarray(dfModel["pop"])
                pop=np.ndarray.tolist(popArr/1000)
                pco2=dfModel["pco2_ppm"].tolist()
                time=dfModel["time_yrs"].tolist()
                temp=dfModel["temp"].tolist()
                runArr= [pop,temp,pco2,time]#make list of lists
                fullArrTr = zip(pop,temp,pco2,time)#make list of lists
                fullArr.append(runArr)#add the data from the run into the 2dList
    dfData = pd.DataFrame({'temp':tempData, 'pco2':pco2Data})
    return dfModel, dfData, equilibrium, eqTemp, eqTime
        #end = time.time()
    #print( "Elapsed Time: " + str(end-start))
    
def analyzeRun(dfModel,nameList,verbose):
    counter = 0
    maxima = dfModel.max();#find maxima from all columns in df
    maxPop = maxima[3];#find maxima in population column, peak popultion
    finalPop = dfModel['pop'][dfModel.index[-1]]
 #   print(maxPop)#peak population
 #   print(finalPop)
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
        dP+=50
       # print(str(dP)+": "+str(newDF.shape[0]))
       # print(newDF)
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
              'maxPopPlot': maxPop/1000, 'finalPop': (finalPop/1000)}#maximum plotting range for population
    if(verbose):
        for k,v in popStats.items():
            print(k + " = " +str(v))
    return popStats

def habitableZone(nameList,newMaxPop,newA,runTime,dA):
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
        nameList['ebm']['Nmax']=newMaxPop#convert pco2 to bars
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
            
def runProgram(driver,nameList,output,showInputs): #run the program with the given name  
    #make temporary directory to run in
    with tempfile.TemporaryDirectory() as dirpath:
        runFolder = newFolder(nameList,dirpath) #make the temporary folder
        
        call("./"+driver,shell=True) #run the driver

        if(showInputs):
            for i in nameList['ebm']:  print(i,": ",nameList['ebm'][i])
        dfModel, finalavgtemp, eqTime, eqTemp, equilibrium = readOutput() #read output into datframe
        
        os.chdir(notePath)
        call("rm -rf tmp*", shell=True)#delete the temporary folder and unlink it's contents
        if(output):
            if(equilibrium):
                if(nameList['ebm']['lverbose']):
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
            'rTech' : 1.00,
            'opT' : 290.5,
            'rco2' : 1.9e-4,#3.459e-4,
            'En0' : 1.00,
            'fragility' : 1.00,
            'coupled' : True,
            'lverbose' : True,
            'runTime' : 100,
            'dtemp' : 1.73,
            'dpco2' : 6.3e-5
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
