# Advanced IonChannel Markov Model Optimizer


This is the main directory for the program with necessary files required run a sample K+ conductance voltage optimization. 

*NOTE: The software utilizes the Intel MKL libraries for matrix multiplication, addition, etc. Please be sure to download the libraries and amend the sample Makefile appropriately. See the Makefile section for more information*

Overview of the directories:

Sample_Data: Includes sample files: WTgv.dat, WTgv.prototxt, WTgv_val.dat, WTgv_val.prototxt, protocolsWT.lst, and validationWT.lst. 
  -WTgv.dat: experimental data to fit in the form of "Voltage" "Normalized Conductance" "SEM/SD of Measurement"
  -WTgv.prototxt:  file encoding the voltage protocol to be utilized in MarkovChannel.cc.  
          name: wtgv #name of the voltage protocol
          source: WTgv.dat #name of the experimental data file
          v0: -70.0 #holding potential (mV)
          normalize: 1 #yes, normalize the output
          FLAG: 0 #keep zero for now
          weight: 1 #protocol weight for in the overall cost function
          has_validation: 1 #indicates that WTgv_val.dat and WTgv_val.prototxt need to be included as well

          step { # GV voltage protocol has one step. For each voltage recorded in the .dat file after holding at v0, simulate the model                   for 50 ms, outputting the maximum open probabillity (PEAK)  every 1.0 ms.
           dt: 50.0 
           stype: PEAK
           stepsize: 1.0
           }
   -WTgv_val.dat: record the desired voltage/normalized conductances to use as additional validation points
   -WTgv_val.protoxt: corresponding .prototxt for the WTgv_val.dat
   -protocolsWT.lst: list of the protocols to train on, just WTgv.prototxt for now
   -validationWT.lst: empty for now, MUST INCLUDE, when fitting current traces this will eventually become populated
