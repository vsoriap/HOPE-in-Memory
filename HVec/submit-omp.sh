#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .
#SBATCH --output=submit-omp.sh.o%j
#SBATCH --error=submit-omp.sh.e%j


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

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export OMP_NUM_THREADS=$4

/usr/bin/time -o time-$1-$4-${HOST} ./$1 $2 $3
