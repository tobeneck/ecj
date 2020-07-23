# Implementation of the traceable evolutionary algorithm (T-EA)

This repository contains the implementation of the T-EA, used in the paper "T-EA: A Traceable Evolutionary Algorithm". It is based on the in java implemented an evolutionary computation framework ECJ. For more information on ECJ and how to use it refer to  [ECJ's official website](http://cs.gmu.edu/~eclab/projects/ecj/) or this [tutorials](https://parabon.com/dev-center/origin/ecj/index.html).

Currently, bit-vector and integer-vector representations are supported. Three traceable problems are already implemented, the Max Ones Problem, the 0/1 Knapsack Problem and the (Un)bound Knapsack Problem. They can be found in the ecj/src/main/java/ec/app/TraceableProblems folder, alongside the statistics classes. Besides these implementations, python scripts for generating many testruns at once are provided in the GenerateTests folder.
