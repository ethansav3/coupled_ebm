import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import seaborn as sns
from statistics import *
from math import *
import shutil
import sys
from matplotlib.ticker import FormatStrFormatter

def plotModelOutput(df,inputs,eqTime,eqTemp,popStats):
    sns.set_style('darkgrid')
    sns.set_context('poster',rc={'font.size': 30.0,
     'axes.labelsize': 26.0,
     'axes.titlesize': 24.0,
     'xtick.labelsize': 26.0,
     'ytick.labelsize': 26.0,
     'legend.fontsize': 22.0})
    timer = np.asarray(df['time'])
    temp = np.asarray(df['temp'])
    pop = np.asarray(df['pop'])
    pco2=np.asarray(df['pco2'])
    
    timer = timer/60/60/24/365.25; #convert seconds to years
    timer = timer+1820;
    pop = pop/1000
    pco2 = pco2*10**6
    
    fig, (ax2, ax1) = plt.subplots(2,sharex=True,figsize=(24.5,11.8),dpi=200) #set up figure, share the x axis
    fig.suptitle("Distance: " + str(inputs[0]) +" AU,  $pCO_{2}$: " + str(round(inputs[1])) +"ppm,  $T_{eq}$: "+str(eqTemp)+"K",x=.41)
     #plot time vs temp (K)
    line1 = ax1.scatter(timer,temp,c=pco2,cmap='jet')
    ax1.set_title('Temperature vs Time')
    ax1.set(ylabel='Temperature (K)',xlabel='Time (years)')
    ax1.set_yticks(np.linspace(min(temp),max(temp),4))
    color='black'
    linestyle='--'
    alpha=.5
    
    ax1.set_xlim(min(timer),max(timer))

    
    sns.set_style('darkgrid')

    #plot time vs pop    
    line2 = ax2.scatter(timer,pop,c=pco2,cmap='jet')
    
    ax2.set(ylabel='Population (billions)')
    ax2.set_title("Population vs Time")
    ax2.set_yticks(np.linspace(min(pop),popStats['maxPopPlot'],4))
    sns.set_palette('colorblind') 
    #horizontal lines
    ax2.axhline(y=popStats['maxPop'],c='black',linestyle='--',label='$N_{max}=$'+str(round(popStats['maxPop'],1))+" billion")
    ax2.axhline(y=popStats['halfPop'],c='black',linestyle=':',label='$N_{max}=$'+str(round(popStats['halfPop'],1))+" billion")       
    #vertical lines
    ax2.axvline(x=(popStats['LhalfTime']+1820),linestyle='--',c=(0,0,.7))
    ax1.axvline(x=(popStats['LhalfTime']+1820),linestyle='--',c=(0,0,.7),label="$t_{1/2}^{-}=$"+str(round((popStats['LhalfTime']+1820),1)))
    ax2.axvline(x=(popStats['maxTime']+1820),c=(0,.7,0),linestyle='--') 
    ax1.axvline(x=(popStats['maxTime']+1820),c=(0,.7,0),linestyle='--',label="$t_{max}=$"+str(round((popStats['maxTime']+1820),1)))
    ax2.axvline(x=(popStats['UhalfTime']+1820),linestyle='--',c=(.7,0,0))
    ax1.axvline(x=(popStats['UhalfTime']+1820),linestyle='--',c=(.7,0,0),label="$t_{1/2}^{+}=$"+str(round((popStats['UhalfTime']+1820),1))) 
    ax2.set_xlim(min(timer),max(timer))
    ax2.set_ylim(min(pop),popStats['maxPopPlot'])
    fig.colorbar(line2,label='pCO2 (ppm)',ax=[ax1,ax2])
    ax2.legend(loc='best')
    ax1.legend(loc='lower right')
    plt.show()
    
def plotTruePopCo2(dfPopCo2):
    #true population and co2 data
    if dfPopCo2.population.mean()>15000: dfPopCo2.population = dfPopCo2.population/1000
    co2 = np.asarray(dfPopCo2.co2_ppm)
    population = np.asarray(dfPopCo2.population)
    timeP = np.asarray(dfPopCo2.time)
    
    sns.set(context='notebook', style='darkgrid')
    line=plt.scatter(np.log(population),co2,c=timeP,cmap='jet');
    plt.ylabel('pCO2')
    plt.xlabel('Log(Population)')
    #plt.scatter(np.log(modelPop),modelPco2);
    plt.axvline(x=np.log(population[12]),c='black',alpha=.5, label="Start of Anthropocene")
    plt.axhline(y=co2[12],c='black',alpha=.5)
    plt.title('True Population ($x10^{9}$)  vs pCO2 (ppm)')

    plt.legend(loc='best')
    plt.colorbar(line);
    
def compareModelOutput(modelData,dfTemp,dfPopCo2):
    #true population and co2 data
    if dfPopCo2.population.mean()>15000: 
        dfPopCo2.population = dfPopCo2.population/1000
    co2 = np.asarray(dfPopCo2.co2_ppm)
    population = np.asarray(dfPopCo2.population)
    timeP = np.asarray(dfPopCo2.time)
    
    dfTemp.columns = ['year','anomaly','smoothed']; #rename column headers

    #convert from anomaly series to temperature series with BASELINE = 14.5 C
    dfTemp['tempC']=dfTemp['anomaly']+14.5;
    dfTemp['tempC_smooth']=dfTemp['smoothed']+14.5;
    #then convert from Celsius to Kelvin by adding 273.15
    dfTemp['tempK']=dfTemp['tempC']+273.15;
    dfTemp['tempK_smooth']=dfTemp['tempC_smooth']+273.15;
    temp_smooth = dfTemp['tempK_smooth']
    timeT = np.asarray(dfTemp.year)
    tempK = np.asarray(dfTemp.tempK)
    tempK_smooth = np.asarray(dfTemp.tempK_smooth)

    modelTime = np.asarray(modelData['time'])
    modelTemp = np.asarray(modelData['temp'])
    modelPop = np.asarray(modelData['pop'])
    modelPco2 = np.asarray(modelData['pco2'])
    modelTime = modelTime/60/60/24/365.25; #convert secons to years
    modelTime = modelTime + 1820;
    modelPco2 = modelPco2*10**6

    sns.set(context='talk', style='darkgrid')
    size = 30

    fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True,figsize=(15,10));
    ax1.scatter(timeP,population,s=size, label='True Population');
    ax1.plot(modelTime,modelPop,c='black', label='Model Population');
    ax1.set_title("Population vs Time")
    ax1.legend(loc='best');

    ax2.scatter(timeP,co2,s=size, label='True pC02');
    ax2.plot(modelTime,modelPco2,c='black', label='Model pC02');
    ax2.set_title("pCO2 vs Time");
    ax2.legend(loc='best');

    ax3.scatter(timeT,tempK_smooth,label='True Temp', alpha=.6,s=15);
    ax3.plot(modelTime,modelTemp,c='black', label='Model Temp');
    ax3.set_title("Temp vs Time");
    ax3.legend(loc='best');
    
def plotModelInput(nameList):
    sns.set_style('darkgrid')
    sns.set_context('notebook',rc={'font.size': 30.0,
     'axes.labelsize': 26.0,
     'axes.titlesize': 24.0,
     'xtick.labelsize': 26.0,
     'ytick.labelsize': 26.0,
     'legend.fontsize': 22.0})
    tempList=[]
    gRateList=[]
    bRateList=[]
    dRateList=[]
    
    width=5#line thickness

    fig,ax=plt.subplots(figsize=(8,5))
    opT = nameList['ebm']['opT']
    dtPop = nameList['ebm']['dTpop']

    for ann_tempave in range(int(opT-25),int(opT+25)):
        tempList.append(ann_tempave)
        gRateList.append((1.00/35.90)*exp(-( (ann_tempave-opT)/(dtPop) )**2 )-(1.00/70.00)*exp(( (ann_tempave-opT)/(dtPop))**2 ))
        bRateList.append((1.00/35.90)*exp(-( (ann_tempave-opT)/(dtPop) )**2 ))
        dRateList.append(-(1.00/70.00)*exp(( (ann_tempave-opT)/(dtPop))**2 ))

    plt.plot(tempList,bRateList,label='$r_{birth}$',linewidth=width)
    plt.plot(tempList,gRateList,label='$r$',linewidth=width)
    plt.plot(tempList,dRateList,label='$r_{death}$',linewidth=width)

    plt.xlabel('Temperature (K)')


    plt.axhline(y=0,color='black')

    plt.ylabel('Population Rates ($yr^{-1}$)')
    plt.ylim(-.04,.04)
    plt.yticks(np.linspace(-.04,.04,4))
    
    plt.xlim(opT-25,opT+25)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0d'))
    plt.legend()
    plt.show()
    
