This tutorial demonstrates how to link a simulation code with libTPS.  The tutorial uses a simple Ising model simulation code that I cooked up for the purpose of demonstration (i.e., it is not in any way optimized and may not be correct).  

The steps to using the library are as follows:

step1: Obtain a simulation code capable of generating trajectories. (demonstrated here in terms of the Ising model).

step2: Write a plugin that allows the code to talk with the TPS library.  This works best if the code is object oriented like the Ising code, but this is not a requirement.

step3: Do transition path sampling.  Sample and visualize the transition path ensemble.

step4: Sample and visualize the transition state ensemble

step5: Compute committer surface
