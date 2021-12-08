import subprocess, os, shlex, csv
import pandas as pd



def has_validation(nameProto):
    for idx,e in enumerate(nameProto):
        if 'has_validation' in e:
            x = e.find(":")
            has_v = e[x+2:]
            return has_v


######Populate protocols, validation, solver, inputs)
cwd = os.getcwd()
isdir = os.path.isdir('Sample_Data')
if isdir == False: 
    os.mkdir('Sample_Data')

userInput = open('sampleData.csv','r')
data = csv.DictReader(userInput)

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

# write to solver
file = open("solver.txt","w")
counter = 0;
for e in solver:
    file.write(e+inputs[counter]+'\n')
    counter = counter+1
    
file.close()

# write protocol and validation to .lst files
file = open(cwd+"/Sample_Data"+"/protocols.lst","w")
for e in protocols:
    file.write(e+".prototxt\n")
file.close()

file = open(cwd+"/Sample_Data"+"/validation.lst","w")
for e in validation:
    file.write(e+".prototxt\n")
file.close()

# two .prototxt and two .dat files for each complete protocol
# create and write to each file
####

###
for protocol in protocols:
    nameProto = []
    nameDat = []
    Dat = []
    
    name_valProto = []
    name_valDat = []
    DatVal = []

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

#populate nameProto
    userInput = open('sampleData.csv','r')
    data = csv.DictReader(userInput)
    for col in data:         
        nameProto.append(col[protocol+'.prototxt'])
    userInput.close()
    
    val = has_validation(nameProto)
    if(val):
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
    
#populate nameDat.
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

           
    nameProto = list(filter(None, nameProto))
    name_valProto = list(filter(None, name_valProto))
    name_valDat = list(filter(None, name_valDat))
    nameDat = list(filter(None, nameDat))


    


    # name.prototxt
    file = open(cwd+"/Sample_Data/"+protocol+".prototxt","w")
    for e in nameProto:
        file.write(e+'\n')
    file.close()

    # name.dat
    file = open(cwd+"/Sample_Data/"+protocol+".dat","w")
    for e in nameDat:
        file.write(e+'\n')
    file.close()
    if(val):
        # name_val.prototxt
        file = open(cwd+"/Sample_Data/"+protocol+"_val.prototxt","w")
        for e in name_valProto:
            file.write(e+'\n')
        file.close()

        # name_val.dat
        file = open(cwd+"/Sample_Data/"+protocol+"_val.dat","w")
        for e in name_valDat:
            file.write(e+'\n')
        file.close()
    
# for protocol in protocols:
    # nameProto = []
    # nameDat = []
    # Dat = []
    
    # name_valProto = []
    # name_valDat = []
    # DatVal = []

# # find .dat column index
    # userInput = open('sampleData.csv','r')

    # reader = csv.reader(userInput)
    # counter = -1
    # for row in reader:
        # for e in row:
            # counter=counter+1
            # if e == protocol+'.dat':
                # break
        # break
    # userInput.close()

# # find _val.dat column index
    # userInput = open('sampleData.csv','r')
    # reader = csv.reader(userInput)
    # counterVal = -1
    # for row in reader:
        # for e in row:
            # counterVal=counterVal+1
            # if e == protocol+'_val.dat':
                # break
        # break
    # userInput.close()



# #populate nameProto and name_valProto
    # userInput = open('sampleData.csv','r')
    # data = csv.DictReader(userInput)
    # for col in data:         
        # nameProto.append(col[protocol+'.prototxt'])
           
        # name_valProto.append(col[protocol+'_val.prototxt'])
    # userInput.close()

# #populate nameDat.
    # userInput = open('sampleData.csv','r')
    # reader = pd.read_csv(userInput,skipinitialspace=True)
    # dat = reader.iloc[:,[counter,counter+1,counter+2]]
    # dat = dat.dropna()
    # numrow = len(dat)
    # col = dat.head(numrow)
    # clmn = list(col)
    # for e in range(1,numrow):
        # for i in clmn:
            # Dat.append(col[i][e])

    # for e in range(0,numrow-1):
        # nameDat.append(Dat[e*3]+" "+Dat[e*3+1]+" "+Dat[e*3+2])
    # userInput.close()

# #populate DatVal
    # userInput = open('sampleData.csv','r')
    # reader = pd.read_csv(userInput,skipinitialspace=True)
    # dat = reader.iloc[:,[counterVal,counterVal+1,counterVal+2]]
    # dat = dat.dropna()
    # numrow = len(dat)
    # col = dat.head(numrow)
    # clmn = list(col)
    # for e in range(1,numrow):
        # for i in clmn:
            # DatVal.append(col[i][e])

    # for e in range(0,numrow-1):
        # name_valDat.append(DatVal[e*3]+" "+DatVal[e*3+1]+" "+DatVal[e*3+2])
    # userInput.close()




            
    # nameProto = list(filter(None, nameProto))
    # name_valProto = list(filter(None, name_valProto))
    # name_valDat = list(filter(None, name_valDat))
    # nameDat = list(filter(None, nameDat))


    


    # # name.prototxt
    # file = open(cwd+"/Sample_Data/"+protocol+".prototxt","w")
    # for e in nameProto:
        # file.write(e+'\n')
    # file.close()

    # # name.dat
    # file = open(cwd+"/Sample_Data/"+protocol+".dat","w")
    # for e in nameDat:
        # file.write(e+'\n')
    # file.close()

    # # name_val.prototxt
    # file = open(cwd+"/Sample_Data/"+protocol+"_val.prototxt","w")
    # for e in name_valProto:
        # file.write(e+'\n')
    # file.close()

    # # name_val.dat
    # file = open(cwd+"/Sample_Data/"+protocol+"_val.dat","w")
    # for e in name_valDat:
        # file.write(e+'\n')
    # file.close()
    
# validation only protocols
for protocol in validation:
    nameProto = []
    nameDat = []

    Dat = []
    userInput = open('sampleData.csv','r')
    data = csv.DictReader(userInput)
    reader = csv.reader(userInput,skipinitialspace=True)
    counter = -1
    for row in reader:
        for e in row:
            counter=counter+1
            if e == protocol+'.dat':
                break
        break
    userInput.close()
    
    userInput = open('sampleData.csv','r')
    data = csv.DictReader(userInput)
    for col in data:         
        nameProto.append(col[protocol+'.prototxt'])
           
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
    
    nameProto = list(filter(None, nameProto))
    nameDat = list(filter(None, nameDat))
    file = open(cwd+"/Sample_Data/"+protocol+".prototxt","w")
    for e in nameProto:
        file.write(e+'\n')
    file.close()


    file = open(cwd+"/Sample_Data/"+protocol+".dat","w")
    for e in nameDat:
        file.write(e+'\n')
    file.close()

    userInput.close()



