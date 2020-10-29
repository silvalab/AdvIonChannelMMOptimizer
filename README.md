# Advanced Ion Channel Markov Model Optimizer


## Run a Containerized Sample K+ GV Optimization

  - Download and install Docker for your operating system. https://www.docker.com/get-started
  - To verify your installation was successful, execute: 
      - docker run hello-world
  
  -  If successfull, run:
      - docker pull mangoldk/silvalabwustl:k_sample_opt2
      - docker run mangoldk/silvalabwustl:k_sample_opt2


## Quickstart Guide for Running on Your Own Cluster
  - Download all files
  - Install Intel MKL libraries. Update paths in Makefile and job.sh.
  - run "make" 
  - run "./job.sh 3 3 10 test"
*NOTE: The software utilizes the Intel MKL libraries for matrix multiplication, addition, etc. Please be sure to download the libraries and amend the sample Makefile appropriately. See the Makefile section for more information*
   
## Folders/File Explanations:

- **BMKFiles:**
    
    Output .txt files for all rooted graphs by size. Please see https://users.cecs.anu.edu.au/~bdm/. These files serve as input
    to parse_convert.py to create a parsed list of rooted graphs with biophysical restrictions.
    
- **Sample_Data:** 

  Includes sample files: WTgv.dat, WTgv.prototxt, WTgv_val.dat, WTgv_val.prototxt, protocolsWT.lst, and validationWT.lst. 
  
  - WTgv.dat: experimental data to fit in the form of "Voltage" "Normalized Conductance" "SEM/SD of Measurement"
  
  - WTgv.prototxt:  file encoding the voltage protocol to be utilized in MarkovChannel.cc.         
      
    - name: wtgv #name of the voltage protocol
      
    - source: WTgv.dat #name of the experimental data file
      
    - v0: -70.0 #holding potential (mV)
          
    - normalize: 1 #yes, normalize the output
          
    - FLAG: 0 #keep zero for now
         
    - weight: 1 #protocol weight for in the overall cost function
          
    - has_validation: 1 #indicates that WTgv_val.dat and WTgv_val.prototxt need to be included as well

     step { #GV voltage protocol has one step. For each voltage recorded in the .dat file after holding at v0, simulate the model                   for 50 ms, outputting the maximum open probabillity (PEAK)  every 1.0 ms.
      
     dt: 50.0 
      
     stype: PEAK
      
     stepsize: 1.0
     
     }
   
   - WTgv_val.dat: record the desired voltage/normalized conductances to use as additional validation points
   
   - WTgv_val.protoxt: corresponding .prototxt for the WTgv_val.dat
   
   - protocolsWT.lst: list of the protocols to train on, just WTgv.prototxt for now
  
   - validationWT.lst: empty for now, MUST INCLUDE, when fitting current traces this will eventually become populated

- **State3:**

  - Sample directory for storing results of 3-state model optimizations. When the program runs and you indicate a "Version," the results    folder is stored here.
 
  - State3parsedDT4CLT4.txt: This is the output of parsing unique rooted models (see the Python supplelmental code one directory up for more information)
  - Lists the 3-state models by ID, state used as the "root" (open), the list of edges
  
- **includeH:**
 
     -  Headers to includes for the various classes in the program
  
- **src:**
  
  - .cc files of course!
  
- **Makefile**
  
  - Sample Makefile to get you started. The paths to the Intel MKL include/link libraries will need to updated for your installation.
  
- **Sobol License**
   - License to use the direction numbers for the Sobol sequence generation
  
- **job.sh**
  
   - Sample script to show how to call the program. The script takes four required arguments (Number of States) (ID of Model) (Sobol Starts Desired) (Version)
  
   - Thus to optimize the fully connected 3-state model on the sample data (with 10 Sobol starts) and save results to folder "test", enter: 
  
     - ./job.sh 3 3 10  test
  
- **joe-kuo-other-4.5600**
   - Direction numbers for the Sobol sequence
  
- **solver.txt**
  - File of the simulation input parameters. 

       - double init_mu: 0

       - double init_std: 3.0

       - double mut_prob: 0.5

       - double update_mu_rates: 0 #bounds for gaussian perturbation, rates

       - double update_std_rates: 1

       - double update_mu_args: 0 ##bounds for gaussian perturbation, voltage dependent arguments

       - double update_std_args: 3

       - double gamma: 0.001 #factor to increase temperature for each worse solution encountered in simulated annealing

      - double t0: 0.00001 #starting temperature

      - double rate_min: -6 #ln(rate) limits to initialize

      - double rate_max: 6

      - double arg_min: -60 #voltage dependent arguments min and max

      - double arg_max: 60

      - int k_max: 1000000 #how many iterations to run the simulation

      - int step: 1000 #can be used to manipulate the temperature

      - int display: 10000 #how often the current best cost function and value is computed along with overfitting metrics

      - int n_chains: 10 #number of noninteracting chains used in simulated annealing

      - int restart: 1

      - string AWS_S3_path: "s3://my-bucket/" #If you want to run the containerized optimization on AWS Batch, list your bucket name here for S3 saving

      - int snapshot: 50000 #how often .model and .txt files are written to show progress

  
