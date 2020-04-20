#!/bin/bash

# NOTE: #PBS comments request resources and have to come before any commands

# This option will send email when the job starts/finishes
# To receive this email, you'll need a .forward file in your home
# directory that contains the email address you wish to receive your
# email at.

# supposedly stop emailing me
#PBS -m n
# Here I request 1 node with 1 processors per node (ppn) for 1 hour of
# walltime.
#PBS -l nodes=1:ppn=1,walltime=24:00:00

echo Your args are: "$@"
echo $LD_LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/export/intel/compilers_and_libraries_2017.6.256/linux/compiler/lib/intel64:/export/intel/compilers_and_libraries_2017.6.256/linux/mkl/lib/intel64
echo $LD_LIBRARY_PATH
echo $N
echo $NUMBER
echo $TIMES
echo $VERSION
ulimit -c unlimited
ulimit -s unlimited
if [ $# -eq 1 ] 
then
	echo "restart"
	/home/mangoldk/GitHub/Evalulator /home/mangoldk/GitHub/solver.txt /home/mangoldk/GitHub/Master_Data2/protocolsWT.lst /home/mangoldk/GitHub/Master_Data2/validationWT.lst ${N} /home/mangoldk/GitHub/State${N}/State${N}parsedDT4CLT4.txt ${NUMBER} /home/mangoldk/GitHub/joe-kuo-other-4.5600 ${TIMES} ${VERSION} /home/mangoldk/GitHub/State${N}/${VERSION}/Model${NUMBER}/iter_state.txt
else
	
	/home/mangoldk/GitHub/Evalulator /home/mangoldk/GitHub/solver.txt /home/mangoldk/GitHub/Master_Data2/protocolsWT.lst /home/mangoldk/GitHub/Master_Data2/validationWT.lst ${N} /home/mangoldk/GitHub/State${N}/State${N}parsedDT4CLT4.txt ${NUMBER} /home/mangoldk/GitHub/joe-kuo-other-4.5600 ${TIMES} ${VERSION}
fi
