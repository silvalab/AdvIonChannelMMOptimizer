import subprocess, os, shlex, csv
import pandas as pd
from numpy import inf
import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import tkinter 
import math  
from pylab import *



def yline(start,end,Voltage):
    point1 = [start, Voltage]
    point2 = [end, Voltage]

    x_values = [point1[0], point2[0]]
    
    y_values = [point1[1], point2[1]]
    


    plt.plot(x_values, y_values,c = 'black',linewidth=5)

def xline(start,Voltage_start,Voltage_end):
    point1 = [start, Voltage_start]
    point2 = [start, Voltage_end]

    x_values = [point1[0], point2[0]]
    
    y_values = [point1[1], point2[1]]
    


    plt.plot(x_values, y_values,c = 'black',linewidth=5)
class Step:
    def __init__(self,T_tot,vm_prev,dt = None,vm=None,end=None):
        self.dt = dt
        self.T_tot = T_tot
        self.vm = vm
        self.vm_prev = vm_prev
        self.end = end
    def draw(self):
        
        if isinstance(self.vm, float) and isinstance(self.dt, float):
            yline(self.T_tot,self.T_tot+self.dt,self.vm)
            if(self.end):
                if self.vm > self.vm_prev:
                    xline(self.T_tot,self.vm_prev,self.vm)
                    xline(self.T_tot+self.dt,self.vm_prev,self.vm)
                else:
                    xline(self.T_tot,self.vm,self.vm_prev)
                    xline(self.T_tot+self.dt,self.vm,self.vm_prev)
            else:
                if isinstance(self.vm_prev, float):
                    if self.vm > self.vm_prev:
                        xline(self.T_tot,self.vm_prev,self.vm)
                        
                    else:
                        xline(self.T_tot,self.vm,self.vm_prev)
                else: ## self.vm_prev is a list
                        xline(self.T_tot,self.vm_prev[0],self.vm)
                        xline(self.T_tot,self.vm_prev[-1],self.vm)
        elif isinstance(self.dt, float): ### vm is a vector
            if(self.end):
                for V in self.vm:
                    yline(self.T_tot,self.dt,V)
                xline(self.T_tot,self.vm_prev,self.vm[0])
                xline(self.T_tot+self.dt,self.vm_prev,self.vm[0])
            else:
                for V in self.vm:
                    yline(self.T_tot,self.dt,V)
                if isinstance(self.vm_prev, float):
                    xline(self.T_tot,self.vm_prev,self.vm[0])
                    xline(self.T_tot,self.vm_prev,self.vm[-1])
                else: ## self.vm_prev is a list
                    xline(self.T_tot,self.vm_prev[0],self.vm[0])
                    xline(self.T_tot,self.vm_prev[1],self.vm[0])
                    xline(self.T_tot,self.vm_prev[1],self.vm[1])
                    xline(self.T_tot,self.vm_prev[0],self.vm[1])
        elif isinstance(self.vm, float): ##need to extract the vector of times
            if(self.end):
                for t in self.dt:
                    yline(self.T_tot,self.T_tot+t,self.vm)
                xline(self.T_tot,self.vm_prev,self.vm)
                xline(self.T_tot+self.dt[-1],self.vm_prev,self.vm)
            else:
                for t in self.dt:
                       
                        yline(self.T_tot,self.T_tot+t,self.vm)
                        xline(self.T_tot+t,self.vm,self.vm+50)
                if isinstance(self.vm_prev, float):
                    xline(self.T_tot,self.vm_prev,self.vm)
                else:
                    xline(self.T_tot,self.vm_prev[0],self.vm)
                    xline(self.T_tot,self.vm_prev[1],self.vm)
                
class HP:
    def __init__(self,vm):
        self.dt = -20
        self.vm = vm
    def draw(self):
        yline(self.dt,0,self.vm)

class Step_Info:
    def __init__(self,dt=None,vm=None):
        self.dt = dt
        self.vm = vm
        self.y = None
        self.error = None
    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)
    def find_missing(self):
        if self.dt is None:
            return 0
        elif self.vm is None:
            return 1
        else:
            return -1
def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list

def find_closest_biggest(list,index):
    #print(list)
    #print(index)
    if len(list) == 1:
        return list[0]
    mins = []
    for item in list:
        if(index < item):
            mins.append(item-index)
        else:
            mins.append(inf)
    #print(mins)
    min_value = min(mins)
    min_index = mins.index(min_value)
    return list[min_index]
    
    
    
    
    
def has_validation(nameProto):
    for idx,e in enumerate(nameProto):
        if 'has_validation' in e:
            x = e.find(":")
            has_v = e[x+2:]
            return has_v
            
def plot_data(protocol,index,nameProto,nameDat,steps_protocol):    
    for idx,e in enumerate(nameProto):
        #print("e",e)
        if 'name' in e:
            x = e.find(":")
            name = e[x+2:]
        if 'v0' in e:
            x = e.find(":")
            v0 = e[x+2:]
        if 'step {' in e:
            upper_lim = find_closest_biggest(index,idx)
            steps_protocol.append(Step_Info())
            #print("upper_lim",upper_lim)
            for line in range(idx,upper_lim,1):
                if 'vm' in nameProto[line]:
                    #print(nameProto[line])
                    x = nameProto[line].find(":")
                    vm = nameProto[line][x+2:]
                    steps_protocol[-1].vm = float(vm)
                    #print(steps_protocol[-1].vm)
                if 'dt' in nameProto[line]:
                    #print(nameProto[line])
                    x = nameProto[line].find(":")
                    dt = nameProto[line][x+2:]
                    steps_protocol[-1].dt = float(dt)
                    #print(steps_protocol[-1].dt)
    
    
    
    
    
    # name.dat
    
    
    x = []
    y = []
    error = []
    
    for e in nameDat:
       #print("e",e)
       txt = e.split()
       x.append(txt[0])
       y.append(txt[1])
       error.append(txt[2])
    x = [float(z) for z in x]
    y = [float(z) for z in y]
    error = [float(z) for z in error]
    
    
    #print(x)
    #print(y)
    #print(error)
    
    
    
    for s in steps_protocol:
        missing = s.find_missing()
        if (missing ==0): ##dt is missing
            s.dt = [float(z) for z in x]
            s.y = [float(z) for z in y]
            s.error = [float(z) for z in error]
        elif (missing == 1):
            s.vm = [float(z) for z in x]
            s.y = [float(z) for z in y]
            s.error = [float(z) for z in error]
    #call step draw function here with dt,vm, v0....
    #print(type(s.dt))
    #print(type(s.vm))
    vm_prev = float(v0)
    hp = HP(float(v0))
    hp.draw()
    T_tot = 0
    
    for s in steps_protocol:
        
        step = Step(T_tot,vm_prev,vm = s.vm,dt = s.dt,end=0)
        step.draw()
        
        if isinstance(s.dt, float):
            T_tot = T_tot + s.dt
        else:
            T_tot = T_tot + s.dt[-1]
        if isinstance(s.vm, float):
            vm_prev = None
            vm_prev = s.vm
        else:
            vm_prev = []
            vm_prev.append(s.vm[-1])
            vm_prev.append(s.vm[0])
    plt.title(name)
    plt.xlabel('Time (ms)',fontsize=20)
    plt.ylabel('Voltage (mV)',fontsize=20)
    ax = gca()
    for axis in ['left','bottom']:
        ax.spines[axis].set_linewidth(2)

    for axis in ['top','right']:
        ax.spines[axis].set_linewidth(0)
    fontsize = 20  
    for tick in ax.xaxis.get_major_ticks():
        #print(tick.properties())
        tick.label.set_fontsize(fontsize)
        tick.label.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    #print(plt.ylim())
    #plt.ylim((-100, 60))
    plt.show()
    #print(plt)
    
    
    plt.figure()
    plt.scatter(x,y,c = 'black',linewidth=5)
    plt.errorbar(x,y,error,ls='none',fmt='none',c='black',capsize=5)
    plt.title(name)
    plt.xlabel('Independent Experimental Data',fontsize=20)
    plt.ylabel('Dependent Experimental Data',fontsize=20)
    ax = gca()
    for axis in ['left','bottom']:
        ax.spines[axis].set_linewidth(2)

    for axis in ['top','right']:
        ax.spines[axis].set_linewidth(0)
    fontsize = 20  
    for tick in ax.xaxis.get_major_ticks():
        #print(tick.properties())
        tick.label.set_fontsize(fontsize)
        tick.label.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    plt.show()
        











def extract_data(protocol):
    nameProto = []
    nameDat = []
    # name_valProto = []
    # name_valDat = []

    Dat = []
    #DatVal = []

    


# find .dat column index
    userInput = open('sampleData.csv','r')

    reader = csv.reader(userInput)
    counter = -1
    for row in reader:
        for e in row:
            counter=counter+1
            if e == protocol+'.dat':
                break
        break
    userInput.close()

# find _val.dat column index
    userInput = open('sampleData.csv','r')
    reader = csv.reader(userInput)
    counterVal = -1
    for row in reader:
        for e in row:
            counterVal=counterVal+1
            if e == protocol+'_val.dat':
                break
        break
    userInput.close()




    userInput = open('sampleData.csv','r')
    data = csv.DictReader(userInput)
    for col in data:         
        nameProto.append(col[protocol+'.prototxt'])
           
        #name_valProto.append(col[protocol+'_val.prototxt'])


    userInput.close()


    userInput = open('sampleData.csv','r')
    reader = pd.read_csv(userInput,skipinitialspace=True)
    dat = reader.iloc[:,[counter,counter+1,counter+2]]
    dat = dat.dropna()

    numrow = len(dat)
    col = dat.head(numrow)
    clmn = list(col)

    for e in range(1,numrow):
        for i in clmn:
            Dat.append(col[i][e])

    for e in range(0,numrow-1):
        nameDat.append(Dat[e*3]+" "+Dat[e*3+1]+" "+Dat[e*3+2])
    userInput.close()

    # userInput = open('sampleData.csv','r')
    # reader = pd.read_csv(userInput,skipinitialspace=True)
    # dat = reader.iloc[:,[counterVal,counterVal+1,counterVal+2]]
    # dat = dat.dropna()

    # numrow = len(dat)
    # col = dat.head(numrow)
    # clmn = list(col)

    #for e in range(1,numrow):
        #for i in clmn:
            #DatVal.append(col[i][e])

    #for e in range(0,numrow-1):
        #name_valDat.append(DatVal[e*3]+" "+DatVal[e*3+1]+" "+DatVal[e*3+2])
    userInput.close()




            
    nameProto = list(filter(None, nameProto))
    #name_valProto = list(filter(None, name_valProto))
    #name_valDat = list(filter(None, name_valDat))
    nameDat = list(filter(None, nameDat))
    return nameProto,nameDat
    
def extract_data_val(protocol):
    
    name_valProto = []
    name_valDat = []

    
    DatVal = []

   
    userInput = open('sampleData.csv','r')
    reader = csv.reader(userInput)
    counterVal = -1
    for row in reader:
        for e in row:
            counterVal=counterVal+1
            if e == protocol+'_val.dat':
                break
        break
    userInput.close()




    userInput = open('sampleData.csv','r')
    data = csv.DictReader(userInput)
    for col in data:         
        name_valProto.append(col[protocol+'_val.prototxt'])


    userInput.close()




    userInput = open('sampleData.csv','r')
    reader = pd.read_csv(userInput,skipinitialspace=True)
    dat = reader.iloc[:,[counterVal,counterVal+1,counterVal+2]]
    dat = dat.dropna()

    numrow = len(dat)
    col = dat.head(numrow)
    clmn = list(col)

    for e in range(1,numrow):
        for i in clmn:
            DatVal.append(col[i][e])

    for e in range(0,numrow-1):
        name_valDat.append(DatVal[e*3]+" "+DatVal[e*3+1]+" "+DatVal[e*3+2])
    userInput.close()




            
    #nameProto = list(filter(None, nameProto))
    name_valProto = list(filter(None, name_valProto))
    name_valDat = list(filter(None, name_valDat))
    #nameDat = list(filter(None, nameDat))
    return name_valProto,name_valDat
    
    
    
    
###########MAIN ROUTINE HERE############################

userInput = open('sampleData.csv','r',encoding='utf-8-sig')
data = csv.DictReader(userInput)
#print(data.fieldnames)
protocols = []
validation = []
solver = []
inputs = []

# extract protocols
for col in data:
    protocols.append(col['protocols.lst'])
    validation.append(col['validation.lst'])
    solver.append(col['solver_settings'])
    inputs.append(col['Inputs'])
userInput.close()
# delete empty elements
protocols = list(filter(None, protocols))
validation = list(filter(None, validation))


for protocol in protocols:
    
    nameProto,nameDat = extract_data(protocol)
    index = get_index_positions(nameProto,'}')
    val = has_validation(nameProto)
    if(val):
        name_valProto,name_valDat = extract_data_val(protocol)
        steps_protocol_val = [] 
    steps_protocol = [] 
    plot_data(protocol,index,nameProto,nameDat,steps_protocol)
    if(val):
        index = get_index_positions(name_valProto,'}')
        plot_data(protocol,index,name_valProto,name_valDat,steps_protocol_val)

    
# validation only protocols
for protocol in validation:
    nameProto,nameDat = extract_data(protocol)
    steps_validation = []   
    index = get_index_positions(nameProto,'}')
    plot_data(protocol,index,nameProto,nameDat,steps_validation)

