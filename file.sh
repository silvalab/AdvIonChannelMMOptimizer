#!/bin/bash
echo Your container args are: "$@"
echo $1
echo $2
echo $3 
echo $4

pip install pandas
cd /usr/src/dockerdeploy/


python import_and_run.py
dos2unix solver.txt

echo $AWS_ACCESS_KEY_ID
echo $AWS_SECRET_ACCESS_KEY
aws configure list
aws s3 ls --profile=default
	
if [ $# -eq 5 ] 
then
	echo "restart"
	./run solver.txt Sample_Data/protocols.lst Sample_Data/validation.lst $1 State$1/State$1parsedDT4CLT4.txt $2 joe-kuo-other-4.5600 $3 $4 State$1/$4/Model$2/iter_state.txt
else
	
	./run solver.txt Sample_Data/protocols.lst Sample_Data/validation.lst $1 State$1/State$1parsedDT4CLT4.txt $2 joe-kuo-other-4.5600 $3 $4
fi


