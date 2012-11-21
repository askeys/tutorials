This tutorial demonstrates how to use umbrella sampling with replica exchange.  The tutorial uses a simple Ising model simulation code that I cooked up for the purpose of demonstration (i.e., it is not in any way optimized and may not be correct).  

Instructions to accompany the the codes are given here:
https://sites.google.com/site/aaronskeys/resources/tutorials/umbrella-sampling

prerequisite: Obtain a simulation code capable of generating trajectories. (demonstrated here in terms of the Ising model).

The steps to using the library are as follows:

step0: Generate a nucleation trajectory by brute force
 
step1: Simple umbrella sampling with a single window and moving biasing potential.

step2: Naive replica exchange on a single processor.

step3: Replica exchange with MPI
