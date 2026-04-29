# ParMODec
from https://doi.org/10.1007/s10898-019-00778-x
A parallel objective space decomposition algorithm for  extraction of the exact Pareto front of Multiobjective Integer programming prob.

This implementation uses CPLEX UI together with "omp" ["MOBB_parallel.cpp" contains input preparation and calls parallel loops].

Uploaded versions are customized for 4 objective Assignment problem instances. Those who want to try more instances can use the "AssignmentTestInputs" folder.

Those who want to convert the problem structure , for instance into a Knapsack Problem, should modify  "root_solve.h". 
