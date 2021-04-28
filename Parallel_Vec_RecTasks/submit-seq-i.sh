#!/bin/bash

#SBATCH --job-name=submit-seq-i.sh
#SBATCH -D .
#SBATCH --output=submit-seq-i-o%j
#SBATCH --error=submit-seq-i-e%j

USAGE="\n USAGE: ./submit-seq-i.sh n p\n
	n           -> Size of the Matrix
	p           -> Preconditioner"


if (test $# -lt 2 || test $# -gt 2)
then
        echo -e $USAGE
        exit 0
fi

export OMP_NUM_THREADS=1
export PROG=hirschberg_i
export str1=$1
export str2=$2

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export LD_PRELOAD=${EXTRAE_HOME}/lib/libseqtrace.so
./$PROG $str1 $str2 
unset LD_PRELOAD

mpi2prv -f TRACE.mpits -o hirschberg_i-1-${HOST}.prv -e hirschberg_i -paraver
rm -rf  TRACE.mpits set-0 >& /dev/null
