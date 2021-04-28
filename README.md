# HOPE in Memory
This repo contains the source code of **H**irschberg algorithm implemented in **Op**enMP that is **E**fficient **in Memory**. The algorithm goal is to globally align sequences of characters of big sizes. The algorithm uses Levenshtein Distance and Hirschberg algorithm to find the best alignment, while using the smallest amount of memory. The algorithm has been succesfully tested with up to two input strings of 14MB each. Classical implementations of the Smith-Waterman-Gotoh algoritm cannot process these inputs, because it would need 196 Exabytes of RAM to allocate the entire score matrix.

This implementation does not only focus on memory efficiency, but also in performance. There are four versions of the code, each one introduce incremental steps to vectorize and parallelize the source code. The final version of the code is able to achieve 20 GCUPS (Million of Cell Updates per second) on a Intel Xeon CPU E5645 with 24 cores.

## Source Tree Organization
- HSeq: **H**irschberg **Sequential** is the baseline implementation that executes the algorithm without parallelization or vectorization.
- HVec: **H**irschberg **Vec**tor is a vectorized version that exploits anti-diagonal independency.
- VecRec_Tasks: **Vec**tor & **Rec**ursive **Tasks** is a vector and parallel version that incrementally adds OpenMP recursive tasks. The parallelism exploited by these tasks is the decomposition of the score matrix done by Hirchberg Algorithm.
- Parallel_Vec_RecTasks: This is the **final** version that adds fine grained parallelism using parallel blocks computing the anti-diagonals.
- data: This directory contains some examples to execute with the algorithm. It also contains some reference outputs (tiny.out)

## Compilation and run

To compile the different source codes use the make files with the two possible variants hirchberg and hirchberg_omp. For example:

    cd Parallel_Vec_RecTasks
    make hirschberg_omp
    ./hirschberg_omp ../data/seq1_tiny.txt ../data/seq2_tiny.txt

To launch the experiments use SLURM with the submit.sh files:
 
    sbatch submit-omp.sh hirschberg_omp ../data/seq1_tiny.txt ../data/seq2_tiny.txt $NUM_THREADS

To compare the output use diff against the reference output inside data. The different versions are named short.out, medium.out, large.out, double_large.out and quadra_large.out

    ./hirschberg_omp ../data/seq1_tiny.txt ../data/seq2_tiny.txt >tiny.out
    diff tiny.out ../data/tiny.out
