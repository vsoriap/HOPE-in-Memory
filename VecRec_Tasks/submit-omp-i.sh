#!/bin/bash

#SBATCH --job-name=submit-omp-i.sh
#SBATCH -D .
#SBATCH --output=submit-omp-i.sh.o%j
#SBATCH --error=submit-omp-i.sh.e%j

USAGE="\n USAGE: script.sh PROG seq1 seq2 n_th\n
        PROG   -> omp program name\n
        seq1   -> first sequence \n
	seq2   -> second sequence \n
	n_th     -> number of threads\n"

if (test $# -lt 4 || test $# -gt 4)
then
	echo -e $USAGE
        exit 0
fi

make $1


export OMP_NUM_THREADS=$4

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export LD_PRELOAD=${EXTRAE_HOME}/lib/libomptrace.so
./$1 $2 $3 
unset LD_PRELOAD

mpi2prv -f TRACE.mpits -o $1-$4-${HOST}.prv -e $1 -paraver
rm -rf  TRACE.mpits set-0 >& /dev/null
