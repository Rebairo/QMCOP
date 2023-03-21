# The Program Executes Structure of Mini-proteins from sequence:

## Script defs:

1. qmc.sh ---- Executes the following programs in order.
	1. qgen -- generates initial sampling by sobol quasi monte carlo.
	2. lmols -- Initial Model builder.
	3. pdbread -- file format converter.
	4. Clust -- clustering the files.
	5. diheds -- find the diheds of the model.
2. Baker-tranfomation:
    This program is an attempt to tranform the quasi-random deterministic variables generated from sobol quasi-random sequence into more uniform (i.i.d) variables by random shift operation in base 1 followed by a baker/tent transformation to improve scoring in the energy function.  
