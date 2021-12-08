#!/bin/bash
echo $1
echo $2
echo $3 
echo $4

pip install pandas


python import_and_run.py
dos2unix solver.txt

	
if [ $# -eq 5 ] 
then
	echo "restart"
	./run solver.txt Sample_Data/protocols.lst Sample_Data/validation.lst $1 State$1/State$1parsedDT4CLT4.txt $2 joe-kuo-other-4.5600 $3 $4 State$1/$4/Model$2/iter_state.txt
else
	
	./run solver.txt Sample_Data/protocols.lst Sample_Data/validation.lst $1 State$1/State$1parsedDT4CLT4.txt $2 joe-kuo-other-4.5600 $3 $4
fi
