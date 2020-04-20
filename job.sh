#!/bin/bash


echo Your args are: "$@"
#update the contents of the LD_LIBRARY_PATH to make sure the Intel MKL libraries are available. Ammend the following statement for your specific installation
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export/intel/compilers_and_libraries_2017.6.256/linux/compiler/lib/intel64:/export/intel/compilers_and_libraries_2017.6.256/linux/mkl/lib/intel64
echo $LD_LIBRARY_PATH

#Standard Notation for calling the program is
# ./Evalulator N(Number of states) Model_Number(ID of the model of within N to optimize) Times (initial starts) Version

# Therefore to optimize the fully connected 3-state model (3rd listed model in State3/State3parsedDT4CLT4.txt) with 10 Sobol starts and save the results to folder "test"
./Evalulator 3 3 10 test

# To restart the optimization from the last saved progress point, append on restart
./Evalulator 3 3 10 test restart

if [ $# -eq 4 ] 
then
	echo "restart"
	./Evalulator solver.txt Sample_Data/protocolsWT.lst Sample_Data/validationWT.lst $1 State$1/State$1parsedDT4CLT4.txt $2 joe-kuo-other-4.5600 $3 $4 State$1/$4/Model$2/iter_state.txt
else
	
	./Evalulator solver.txt Sample_Data/protocolsWT.lst Sample_Data/validationWT.lst $1 State$1/State$1parsedDT4CLT4.txt $2 joe-kuo-other-4.5600 $3 $4
fi
