#!/bin/bash

#SBATCH --job-name=submit-seq.sh
#SBATCH -D .
#SBATCH --output=submit-seq.sh.o%j
#SBATCH --error=submit-seq.sh.e%j


USAGE="\n USAGE: script.sh PROG size n_th\n
        PROG   -> sequential program name\n
        seq1   -> first sequence \n
	seq2   -> second sequence \n"

if (test $# -lt 3 || test $# -gt 3)
then
	echo -e $USAGE
        exit 0
fi

make $1

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

valgrind --tool=massif --stacks=yes ./$1 $2 $3
