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
import pdb
from scipy.stats import linregress
from statistics import *
import shutil
import sys
from matplotlib.ticker import FormatStrFormatter
from plottingFunctions import *

notePath = os.getcwd()
fullArr=[]
fullArrTr = []
fullMaxPop=0
def runModel(nameList, coupled, runTime, plot, save, analyze, driverName,maxPopList,distList,showInputs,experiment=1,analyzeVerbose=False, printOutput=False, scaleInitPop=False):
    saveName= 0
    count   = 0
    numCols = 4
    numRows = 4
    tempData = []
    pco2Data = []
    lverbose=nameList['ebm']['lverbose']
#    allData = {}
#    allData['pco2'] = []
#    allData['temp'] = []
#    for j in np.linspace(10,100,numCols):
#    for j in [nameList['ebm']['Nmax']]:#pop loops
    eqList = []
    eqTempList = []
    eqTimeList = []
    popDeathList = []       
    for j in maxPopList:#pop loops
        newMaxPop=j#change max pop
#        newA = .9#change distance (AU)
        dA=.01
        popStats = {}
        if(experiment == 1):
            dictdTdP= {
            0.94 : 6.178*10**-2,
            0.96 : 1.655*10**-2,
            0.98 : 1.044*10**-2,
            1.00 : 8.838*10**-3,
            1.02 : 9.067*10**-3
            }
        if(experiment == 2):
            dictdTdP= {
            0.975 : 2.649*10**-1,
            1.0075 : 3.689*10**-3,
            1.040 : 6.846*10**-4,
            1.0725 : 3.130*10**-4,
            1.105 : 1.678*10**-4
            }
            dictPco20= {
            0.975 : 20.806,
            1.0075 : 1054.473,
            1.040 : 13150.274,
            1.0725 : 42286.411,
            1.105 : 93567.242
            }
#        for i in [nameList['ebm']['relsolcon']**-(1/2)]:#distance 
        for i in distList:#distance 
            #            dimVar = (nameList['ebm']['rBirth0']*nameList['ebm']['dtemp'])/(newMaxPop*nameList['ebm']['rco2']*dictdTdP[i]):
            initPop = nameList['ebm']['N0']
            anthroPop = 0 # number of people required to force the climate out of equilibrium in a generation
            if coupled: 
                if experiment==1 or experiment==2: 
                    dimVar=(j*nameList['ebm']['rco2']*dictdTdP[i])/(nameList['ebm']['rBirth0']*nameList['ebm']['dtemp'])
                    anthroPop = (nameList['ebm']['rBirth0']*nameList['ebm']['dtemp'])/(nameList['ebm']['rco2']*dictdTdP[i])
                else:
                    dimVar = 0
                    anthroPop = 0
                if scaleInitPop: 
                    nameList['ebm']['N0']=anthroPop/1000
                    if lverbose: print(f"Initial Pop: {initPop} million,  AnthroPop: {anthroPop:.3f} million,  AnthroPop/1000: {anthroPop/10**3:.3f} million")
            if not coupled:
                dimVar=0
            if experiment ==2: nameList['ebm']['pco20'] = dictPco20[i]*10**-6
            if(lverbose and showInputs):
                print("rBirth0: " + str(nameList['ebm']['rBirth0']))
                print("dT: " + str(nameList['ebm']['dtemp']))
                print("rco2: " + str(nameList['ebm']['rco2']))
                print("maxPop: " + str(j/1000),"billion")
                if coupled: print("dTdP: ",dictdTdP[i])
                print("Distance: " +str(round(i,4)))
                print("Dimensionless Variable: " + str(round(dimVar,5)))
            #calculate inputs
            newA = float(round(i,4)) #change distance (AU)
            nameList['ebm']['coupled']=coupled
            nameList['ebm']['relsolcon']=newA**-2 #inverse square law for solar flux
            nameList['ebm']['runTime'] = runTime#change runtime
            nameList['ebm']['Nmax'] = newMaxPop #change carrying capacity
            dtemp = nameList['ebm']['dtemp'] 
            #run program
            dfModel, finalavgtemp, eqTime, eqTemp, equilibrium, popDeath = runProgram(driverName,nameList,True,showInputs,dimVar, printOutput)#False=no output
            eqList.append(equilibrium)
            eqTempList.append(eqTemp)
            eqTimeList.append(eqTime)
            popDeathList.append(popDeath)
            pco2Data += dfModel['pco2'].tolist()
            tempData += dfModel['temp'].tolist()
            #plot the results
            if coupled and equilibrium:
                inputRunTime = runTime #how long I originally specified
                outputRunTime = math.ceil( dfModel['time_yrs'].iloc[-1] )#how long it actually 
                if analyze:
                    popStats = analyzeRun(dfModel,nameList, anthroPop, analyzeVerbose)#if True, then Print Dictionary Values 
                    maxTime = popStats['maxTime']                  
#                     if outputRunTime - maxTime < 50:
#                         print("problem occured")
                    fullMaxPop=popStats["maxPop"]
                    popStats['maxPopPlot']=fullMaxPop+(3/100)*fullMaxPop#maximum population range
                inputs=[newA,newMaxPop,runTime,dtemp]
                if equilibrium and plot:
                    plotModelOutput(dfModel,inputs,eqTime,eqTemp,popStats,save,saveName,dimVar, exp=experiment)#plot the output of our model, colored by pco2 
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
            if coupled and (not equilibrium):
                if lverbose and printOutput: print("Coupled, but no Equilibrium was reached")
            if (not coupled) and equilibrium:
                if lverbose and printOutput: print("Not coupled, but Equilibrium was reached")
    dfData = pd.DataFrame({'temp':tempData, 'pco2':pco2Data})
    return dfModel, dfData, eqList, eqTempList, eqTimeList, popDeathList

    
def habitableZoneFinder_exp1(nameList, lverbose, showOutput=False):
    '''Given a nameList of values, this function outputs the
    temperature defined habitable zone'''
    minD = 0
    maxD = 0
    pco20 = 284                           
    nameList['ebm']['pco20']=pco20*10**-6  #Convert initial pCO2 to Bars
    newA = 1       
    dA = .001
    while(True):    
        distList = [newA]
        nameList['ebm']['relsolcon']=newA**-2 #inverse square law for solar flux
        dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, False, False, False,"driver.exe",[10000],distList,False, printOutput=False)
        if(not eq):
            minD = newA
            print("Minimum Distance: ", minD)
            break
        else:
            newA -= dA
            if lverbose: print("Distance: ",newA,",    Equilibrium Temp: ", eqTemp )
    pco20 = 284                        
    nameList['ebm']['pco20']=pco20*10**-6  #Convert initial pCO2 to Bars
    newA = 1            
    while(True):    
        distList = [newA]
        nameList['ebm']['relsolcon']=newA**-2 #inverse square law for solar flux
        dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, False, False, False,"driver.exe",[10000],distList,False, printOutput=False)
        if(eqTemp < 273.15):
            maxD = newA
            print("Maximum Distance", maxD)
            break
        else:
            if lverbose: print("Distance: ",newA,",    Equilibrium Temp: ",eqTemp )
            newA += dA
    return minD, maxD
    
def habitableZoneFinder_exp2(namelist, lverbose=False, showOutput=False):
    '''Given a nameList of values, this function outputs the
    pco2 defined habitable zone'''
    minA, maxA = 0,0
    newA = 1
    dA = 0.001
    while(True):#find minimum distance
        pco2Earth = pco2Finder(287.09, nameList, newA, lverbose=False)
        if lverbose: print(f"Distance: {newA},  pCO2: {pco2Earth}")
        if(np.isnan(pco2Earth)): 
            minA = newA
            print(f"Min: {minA}")
            break
        else:
            newA -= dA
    newA=1
    while(True):#find maximum distance
        pco2Earth = pco2Finder(287.09, nameList, newA, lverbose=False)
        if lverbose: print(f"Distance: {newA},  pCO2: {pco2Earth}")
        if(pco2Earth >= 10**5): 
            maxA = newA
            print(f"Max: {maxA}")
            break
        else:
            newA += dA    
    return minA, maxA
    
def roomPco2(nameList,distList,name="driver.exe"):
    '''
    Returns value of pCO2 that will make the equilibrium
    temperature approximately room temperature.
    '''
    pco20 = 10
    pco2Room=0
    count=0
    while(True): #find initial pco2 such that temp is > 295 (about 70 fahrenheit)
        nameList['ebm']['pco20']=pco20/10**6
        dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, coupled, runTime, plot, save, analyze,name,maxPopList,distList,showInputs,experiment=2)
        if(eqTemp >= 295):
            pco2Room = pco20 #the value of pco2 which makes temp be habitable for humans
            print("\n","Distance: ",round(distList[0],4),"\n Room pCO2: ", round(pco2Room,4))
            print(" Number of Iterations: ", count)
            print(" Room Temp: ", eqTemp," Kelvin")
            print(" Room Temp: ", round((eqTemp-273.15)*(9/5)+32,2)," Fahrenheight","\n")
            break
        else:
            pco20 *= 1.5
            count += 1
    return roomPco2


def outlierFinder(tempList,eqTemp,inCount):    
    '''Given current temperature array, size of
    current temperature array, and a temperature,
    returns true if the current temperature is
    NOT an outlier.  If it IS an outlier, return False.'''
    if(inCount > 5):
        data_mean, data_width = np.mean(tempList), np.std(tempList)
        cut_off = data_width*4 #find cutoff for outliers
        lower, upper = data_mean - cut_off, data_mean + cut_off
        if((eqTemp < lower) or (eqTemp > upper)):
   #         print("outlier detected!!! \n Mean Temp: ", data_mean, "\n Equilibrium Temp: ", eqTemp )
            return False
        else:
            return True
    else:
        return True
def deathTimeFinder(dfModel):
    '''Returns a list containing some death statistics.  
    Index 0 gives the time it took for the population to decline by 10% after it peaked
    Index 1 gives the ratio of peak population to carrying capacity'''
    deathStatsList = []
    maxima = dfModel.max()
    maxPop = maxima[3]
    maxPopIndex = dfModel.loc[dfModel['pop']==maxPop].index#search rows for index of max pop
    maxPopTime=dfModel.iloc[maxPopIndex]['time_yrs'];#search column for time until peak pop is reached
    tenPerMaxPop = maxPop - 0.1*maxPop
    #tenPerMaxPopIndex = dfModel.loc[dfModel['pop']==tenPerMaxPop]
    finalPop = dfModel['pop'][dfModel.index[-1]]
    finalTime = dfModel['time_yrs'][dfModel.index[-1]]
    if(finalPop <= tenPerMaxPop):
        finalIndex = dfModel.index[-1]
        tenPerMaxPopIndex = 0
        maxIndex = int(maxPopIndex[0])
        for i in np.arange(maxIndex, finalIndex, 1):
            pop = dfModel["pop"][i]
            if(pop <= tenPerMaxPop):
                tenPerMaxPopIndex = i
                break
        tenPerMaxPopTime = dfModel.iloc[tenPerMaxPopIndex]['time_yrs'];
        dt = float(tenPerMaxPopTime) - float(maxPopTime)
        return dt
    if(finalPop > tenPerMaxPop):
        gen20 = float(finalTime-maxPopTime[maxPopIndex[0]])#500 years
        percentDeadByEnd= (1-finalPop/maxPop)*100
        return gen20*(10/percentDeadByEnd)


def analyzeRun(dfModel, nameList, anthroPop, verbose):
    counter = 0
    maxima = dfModel.max();#find maxima from all columns in df
    maxPop = maxima[3];#find maxima in population column, peak popultion
    finalPop = dfModel['pop'][dfModel.index[-1]]
   # print(maxPop)#peak population
   # print(finalPop)
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
    tDecline = deathTimeFinder(dfModel)
    popStats={'maxPop' : (maxPop/1000), 'maxTime': maxPopTime.mean(),
              'halfPop': (halfPop/1000),'LhalfTime': LhalfPop, 'UhalfTime': UhalfPop,
              'maxPopPlot': maxPop/1000, 'finalPop': (finalPop/1000), 'initPop' : nameList['ebm']['N0'],
              'anthroPop' :anthroPop, 'tDecline': tDecline }#maximum plotting range for population
    if(verbose):
        for k,v in popStats.items():
            print(k + " = " +str(v))
    return popStats


def regionSweep(distances, nameList, verbose=False, exp=1):
    '''Does a parameter sweep of the region of parameter
    space occupied by some distance (or list of distances),
    returns a dictionary containing all relevant data'''
    # for i in :
    iniTemp = 0
    prevTemp = 0
    eqTemp = 0
    outCount=0
    dictPco20= {
    0.975 : 20.806,
    1.0075 : 1054.473,
    1.040 : 13150.274,
    1.0725 : 42286.411,
    1.105 : 93567.242
    }
#     dictPco20= {
#     0.97 : 20,
#     1.0025 : 2750,
#     1.035 : 21500,
#     1.0675 : 62500,
#     1.1 : 130500
#     }
    dictData = {}
    for i in distances:
        dP = 5
        if exp==1: 
            pco20 = 284
            dP = 5
        if exp==2: 
            pco20 = dictPco20[i]
            dP = (1/100)*pco20
        inCount = 0
        pco2List = []
        tempList = []
        distList = [i]
        maxPopList = [10000]
        nameList['ebm']['pco20']=pco20/10**6
        if(inCount == 0):# if first loop, set final temp to initial temp
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, False, False, False,"driver.exe",[10000],distList,False, printOutput=False)
            finalTemp = dfModel['temp'][dfModel.index[-1]] #find final temperature
            iniTemp = eqTemp[0]
            if verbose: print(inCount,")","  dT: ", round(abs(eqTemp[0] - iniTemp),4) ,", pco20: ",pco20,", distance: ",i)
            if verbose: print("Initial Temp: ", iniTemp)
            if verbose: print("\n")
            pco20 -= dP
        while(abs(eqTemp[0] - iniTemp) < 2): #increment pco2s until dT >= 5
            nameList['ebm']['pco20']=pco20/10**6
            prevTemp = eqTemp[0]
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, False, False, False,"driver.exe",[10000],distList,False, printOutput=False)
            if eq:
                finalTemp = dfModel['temp'][dfModel.index[-1]] #find final temperature
                if(outlierFinder(tempList,eqTemp[0],inCount)):
                    pco2List.append(pco20)
                    tempList.append(eqTemp[0])
                    outCount += 2 #counter used to index arrays
                    inCount += 1 #counter used to specify initial temp
                    if verbose: print(inCount,")","  dT: ", round(abs(eqTemp[0] - iniTemp),4) ,", pco20: ",pco20,", distance: ",i)
                    if verbose: print("Equilibrium Temp: ",round(eqTemp[0],2))
                    if verbose: print("\n")
            pco20 -= dP
        dP = 5
        if exp==1: 
            pco20 = 284
            dP = 5
        if exp==2: 
            pco20 = dictPco20[i]
            dP = (2/100)*pco20             
        inCount = 0
        nameList['ebm']['pco20']=pco20/10**6
        if(inCount == 0):# if first loop, set final temp to initial temp
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, False, False, False,"driver.exe",[10000],distList,False, printOutput=False)
            iniTemp = eqTemp[0]
            if verbose: print(inCount,")","  dT: ", round(abs(eqTemp[0] - iniTemp),4) ,", pco20: ",pco20,", distance: ",i)
            if verbose: print("Initial Temp: ", iniTemp)
            if verbose: print("\n")
            pco20 += dP
        while(abs(eqTemp[0] - iniTemp) < 2): #increment pco2s until dT >= 5
            nameList['ebm']['pco20']=pco20/10**6
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, False, False, False,"driver.exe",[10000],distList,False, printOutput=False)
            if eq:
                finalTemp = dfModel['temp'][dfModel.index[-1]] #find final temperature
                if(outlierFinder(tempList,eqTemp[0],inCount)): #if not an outlier
                    pco2List.append(pco20)
                    tempList.append(eqTemp[0])
                    outCount += 2 #counter used to index arrays
                    inCount += 1 #counter used to specify initial temp
                    if verbose: print(inCount,")","  dT: ", round(abs(eqTemp[0] - iniTemp),4) ,", pco20: ",pco20,", distance: ",i)
                    if verbose: print("Equilibrium Temp: ",round(eqTemp[0],2))
                    if verbose: print("\n")
            pco20 += dP
        dictData[outCount-1] = pco2List
        dictData[outCount] = tempList
    return dictData

def pco2Finder(goalEqTemp, nameList, distList, maxPopList,lverbose=False, driver='driver.exe'):
    currEqTemp = 0
    pco20    = 10.3
    exp=0
    nameList['ebm']['pco20'] = pco20*10**-6
    goalPco2 = 0
    plot = False
    save=False
    showInputs=False
    analyze=True
    popDeath = []
    dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, plot, save, analyze,driver,maxPopList,distList,showInputs,analyzeVerbose=False, experiment=0, scaleInitPop=True)
    if(eqTemp[0]>= goalEqTemp):
        return np.nan, np.nan
 #   print(f"Equilbrium Reached at temp: {eqTemp:.2f}, and time: {eqTime:.0f}","\n")
    while(True):
        currEqTemp = eqTemp[0]
        call("rm -rf tmp*",shell=True)
        if(currEqTemp >= goalEqTemp):
            goalPco2 = pco20
            break
        if(goalEqTemp-currEqTemp > 20):
            if lverbose: print("\n >10 \n")
            pco20 *= 2 
            nameList['ebm']['pco20'] = pco20*10**-6
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, plot, save, analyze,driver,maxPopList,distList,showInputs, analyzeVerbose=False,experiment=exp, scaleInitPop=True)
            if np.isnan(eqTemp[0]): return np.nan, np.nan
            currEqTemp = eqTemp[0]
            if lverbose: print(f"pCO20:{pco20:.3f},  Equilbrium Reached at temp: {eqTemp[0]:.2f}, and time: {eqTime[0]:.0f},  goalTemp: {goalEqTemp:.0f}","\n")
        
        if((goalEqTemp-currEqTemp <= 20) and (goalEqTemp-currEqTemp > 10)):
            if lverbose: print("\n <=10,  >=5 \n")
            pco20 += 0.5*pco20
            nameList['ebm']['pco20'] = pco20*10**-6
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 3000, plot, save, analyze,driver,maxPopList,distList,showInputs,analyzeVerbose=False, experiment=exp, scaleInitPop=True)
            if np.isnan(eqTemp[0]): return np.nan, np.nan
            currEqTemp = eqTemp[0]
            if lverbose: print(f"pCO20:{pco20:.3f},  Equilbrium Reached at temp: {eqTemp[0]:.2f}, and time: {eqTime[0]:.0f},  goalTemp: {goalEqTemp:.0f}","\n")
        
        if((goalEqTemp-currEqTemp <= 10) and (goalEqTemp-currEqTemp > 2)):
            if lverbose: print("\n <=5,  >=1 \n")
            pco20 += 0.1*pco20
            nameList['ebm']['pco20'] = pco20*10**-6
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 3000, plot, save, analyze,driver,maxPopList,distList,showInputs, analyzeVerbose=False,experiment=exp, scaleInitPop=True)
            if np.isnan(eqTemp[0]): return np.nan, np.nan
            currEqTemp = eqTemp[0]
            if lverbose: print(f"pCO20:{pco20:.3f},  Equilbrium Reached at temp: {eqTemp[0]:.2f}, and time: {eqTime[0]:.0f},  goalTemp: {goalEqTemp:.0f}","\n")
        
        if(goalEqTemp-currEqTemp <= 2):
            if lverbose: print("\n <= 1 \n")
            pco20 += 0.05*pco20
            nameList['ebm']['pco20'] = pco20*10**-6
            dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 3000, plot, save, analyze,driver,maxPopList,distList,showInputs, analyzeVerbose=False,experiment=exp, scaleInitPop=True)
            if np.isnan(eqTemp[0]): return np.nan, np.nan
            currEqTemp = eqTemp[0]
            if lverbose: print(f"pCO20:{pco20:.3f},  Equilbrium Reached at temp: {eqTemp[0]:.2f}, and time: {eqTime[0]:.0f},  goalTemp: {goalEqTemp:.0f}","\n")
    dT = float(abs(goalEqTemp - currEqTemp))
    return round(goalPco2,3), dT

def timeRemaining(percentDone, timeMinutes):
    timePerPercent = timeMinutes/percentDone
    percentLeft=100-percentDone
    return percentLeft * timePerPercent

def dTdPFinder(goalPco2, distance, nameList, maxPopList, driver="driver.exe"):
    plot=False
    save=False
    show=False
    showInputs=False
    exp=0
    lverbose = False
    analyze=False
    if not np.isnan(goalPco2):
        distList  = [distance]
        #find mean
        nameList['ebm']['pco20'] = goalPco2*10**-6
        dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, plot, save, analyze,driver,maxPopList,distList,showInputs, experiment=exp, printOutput=lverbose, scaleInitPop=True)
        eqTempMean = eqTemp[0]
        #find upper bound
        pco20Plus = goalPco2 + 0.01*goalPco2 #increment pco2 by 1% of its value
        nameList['ebm']['pco20'] = pco20Plus*10**-6
        dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, plot, save, analyze,driver,maxPopList,distList,showInputs, experiment=exp, printOutput=lverbose, scaleInitPop=True)
        eqTempPlus = eqTemp[0]
        #find lower bound
        pco20Minus = goalPco2 - 0.01*goalPco2 #increment pco2 by 1% of its value
        nameList['ebm']['pco20'] = pco20Minus*10**-6
        dfModel,dfData,eq, eqTemp, eqTime, popDeath = runModel(nameList, False, 500, plot, save, analyze,driver,maxPopList,distList,showInputs, experiment=exp,printOutput=lverbose, scaleInitPop=True)
        eqTempMinus = eqTemp[0]
        m1 = (eqTempMean - eqTempMinus)/(goalPco2 - pco20Minus)
        m2 = (eqTempPlus - eqTempMean)/(pco20Plus - goalPco2)
        m3 = (eqTempPlus - eqTempMinus)/(pco20Plus - pco20Minus)
#        print(f"m1: {m1:.3e}, m2: {m2:.3e}, m3: {m3:.3e}")
        return (m1+m2)/2
    else:
        return np.nan

def linearRegressions(dictData, distLin, plotSlopes=False, plotData=False, saveName="regress"):
    '''Finds slopes, given data and distances'''
    slopes = []
    distances = []
    fitErrors = [] 
    dictDistances = {}
    #----------------------------------------------dTdP--------------------------------------------
    for i,v in enumerate(distLin):
        dictDistances[i] = v
    distNum = 0
    count = 0
    for k,v in dictData.items():
        count +=1
        if(count==1): 
            pco2Listing = v; #list of pco2 values
            print("Distance: ",dictDistances[distNum])
            distances.append(dictDistances[distNum])
        if(count==2): 
            tempListing = v; #list of temp values
            print("pco2 length: ",len(pco2Listing))
            print("temp length: ",len(tempListing))
            m = linregress(pco2Listing,tempListing).slope
            slopes.append(m)#-------------------------------------------add to slope array
            b = linregress(pco2Listing,tempListing).intercept
            error = linregress(pco2Listing,tempListing).stderr
            fitErrors.append(error)#------------------------------------add to error array
            print("Slope: " + str("{:.3e}".format(m)))
            print("Intercept: ",round(b))
            if plotData:
                plt.scatter(pco2Listing,tempListing,label='Model Runs',alpha=.5,s=10)
                plt.plot(pco2Listing,np.asarray(pco2Listing)*m+b,label='linear fit',color="orange")
                plt.legend(loc='best')
                plt.show()
            print("\n")
            count = 0
            tempListing = []
            pco2Listing = []
            distNum += 1
    if plotSlopes:
        sns.set_context("paper")
        sns.set_style("darkgrid")
#         plt.figure(1)
        plt.plot(distances,slopes)# from scipy.stats import linregress
       
        plt.scatter(distances,slopes,c="orange")# from scipy.stats import linregress
       #plt.scatter(distances[1],slopes[1])# from scipy.stats import linregress
        # plt.errorbar(distances, slopes, yerr = fitErrors,fmt='o',capthick=3,ms=5,uplims=True,lolims=True)
        plt.errorbar(distances, slopes, color="orange",yerr = fitErrors,fmt='o',capthick=3,ms=5)
        ax = plt.gca()
        ax.set_xticks(distLin)
        plt.title("dT/dP vs Orbital Distance");
        plt.xlabel("Distance (AU)");
        plt.ylabel("dTdP (Kelvin/ppm)");
    #    plt.legend()
        print("save dTdP")
        plt.savefig("dTdP_"+saveName+".png")
        plt.show()
        plt.close('all')
    #----------------------------------------------dPdT--------------------------------------------
    slopes = []
    distances = []
    fitErrors = [] 
    dictDistances = {}
    for i,v in enumerate(distLin):
        dictDistances[i] = v
    distNum = 0
    count = 0
    for k,v in dictData.items():
        count +=1
        if(count==1): 
            pco2Listing = v; #list of pco2 values
     #       print("Distance: ",dictDistances[distNum])
            distances.append(dictDistances[distNum])
        if(count==2): 
            tempListing = v; #list of temp values
            m = linregress(tempListing,pco2Listing).slope
            slopes.append(m)#-------------------------------------------add to slope array
            b = linregress(tempListing,pco2Listing).intercept
            error = linregress(tempListing,pco2Listing).stderr
            fitErrors.append(error)#------------------------------------add to error array
            count = 0
            tempListing = []
            pco2Listing = []
            distNum += 1
    if plotSlopes:
        sns.set_context("paper")
        sns.set_style("darkgrid")
  #      plt.figure(2)
        plt.plot(distances,slopes)# from scipy.stats import linregress
        #plt.scatter(distances[1],slopes[1])# from scipy.stats import linregress
        # plt.errorbar(distances, slopes, yerr = fitErrors,fmt='o',capthick=3,ms=5,uplims=True,lolims=True)
        plt.errorbar(distances, slopes, yerr = fitErrors,color="orange",fmt='o',capthick=3,ms=5)
        ax = plt.gca()
        ax.set_xticks(distLin)
        plt.title("dP/dT vs Orbital Distance");
        plt.xlabel("Distance (AU)");
        plt.ylabel("dPdT (ppm/Kelvin)");
        print("savedPdT")
        plt.savefig("dPdT_"+saveName+".png")
        plt.show()



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
        dfModel, finalavgtemp, eqTime, eqTemp, equilibrium, popDeath = runProgram("driver.exe",nameList,False)#False=no output
        life = (equilibrium) and (eqTemp<=373.15) and (eqTemp>=273.15)#determine habitability
        if(life and (count==0)):#if first habitable distance, make it minA
            minA = newA
            count+=1
        if(life and (count > 0)):#in the habitable zone, make it maxA
            if coupled:
                popStats = analyzeRun(dfModel,nameList,True)
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
            
def runProgram(driver,nameList,output,showInputs,dimVar, printOutput=False): #run the program with the given name  
    #make temporary directory to run in
    with tempfile.TemporaryDirectory() as dirpath:
        runFolder = newFolder(nameList,dirpath) #make the temporary folder
        
        call("./"+driver,shell=True) #run the driver

        if(showInputs):
            for i in nameList['ebm']:  print(i,": ",nameList['ebm'][i])
        dfModel, finalavgtemp, eqTime, eqTemp, equilibrium, popDeath = readOutput() #read output into datframe
        
        os.chdir(notePath)
        call("rm -rf tmp*", shell=True)#delete the temporary folder and unlink it's contents
        if(output):
            if(equilibrium):
                if(nameList['ebm']['lverbose'] and printOutput):
#                    print("Dimensionless Variable: ",str(round(dimVar,5)))
                    print('Final Temp(K): ' + str(finalavgtemp));
                    print('Final Temp(F): ' + str(round((finalavgtemp-273.15)*(9/5)+32, 2)));
                    print('Equilibrium Reached at Temp=' + str(eqTemp)+". At time="+str(eqTime)) 
                    print('Equilibrium Temp(F): ' + str(round((eqTemp-273.15)*(9/5)+32, 2)));
                    print('')
            else:
                if printOutput: print("Equilibrium was not reached,\n")
     #   print('Final Temp(C): ' + str(round(finalavgtemp-273.15)));
        call("echo   ", shell=True)
    return dfModel, finalavgtemp, eqTime, eqTemp, equilibrium, popDeath
#---------------------------Make NameList------------------------
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
            'Nmax' : 10000,
            'rBirth0' : 0.04,
            'rBirthMax' : 0.83,
            'rDeath0' : 0.036,
            'rTech' : 1.00,
            'opT' : 290.5,
            'rco2' : 0.000275,#3.459e-4,
            'En0' : 1.00,
            'fragility' : 1.00,
            'coupled' : True,
            'lverbose' : True,
            'runTime' : 100,
            'dtemp' : 5,
            'dpco2' : 200e-6
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
    global popDeath #how many people died after 2 generations 
    if(equilibrium):
        #read in equilibrium conditions
        eqTemp=0;
        popDeath = 0;
        eqTime=0;
        equilibrium=open("equilibrium.dat","r")
        counter = 0
        for line in equilibrium:
            values = line.split()
            if counter==0:
                eqTime = float(values[0])
                eqTemp = float(values[1])
            else: 
                popDeath = float(values[0])
            counter += 1
    else:
        eqTemp = np.NaN
        eqTime = np.NaN
           
    output.close() # close output file
       
    df = pd.DataFrame(data)
 
    df['time_yrs'] = df['time']/60/60/24/365.25
    df['pco2_ppm'] = df['pco2']*10**6
 
    return df, finalavgtemp, eqTime,eqTemp,equilibrium, popDeath
