# Distributed Online Learning and Control (DOLC) Framework 
This code base helps learn a system dynamics from time-series data using different variants of Dynamic Mode Decomposition (DMD) that includes DMD, Distributed DMD (D-DMD), Sparse Distributed DMD (SD-DMD), Recursive DMD (R-DMD), Recursive Distributed DMD (RS-DMD), Recursive Sparse Distributed DMD (RSD-DMD). 


The codes are commented as per the need and easy to follow. We present three example codes to illustrate how to use the DOLC Framework. 
1. Oscillators_L.mlx - This Matlab code generates data for a network of oscillators and computes open loop Koopman matrix using different DMD variants such as DMD, Distributed DMD (D-DMD), Sparse Distributed DMD (SD-DMD), Recursive DMD (R-DMD), Recursive Distributed DMD (RS-DMD), Recursive Sparse Distributed DMD (RSD-DMD) and compares their learning rates, eigenvalues, and computational times. 
2. Oscillators_LC.mlx - This Matlab code considers a network of heterogenous oscillators and computes the Koopman matrix, input matrix using RD-DMD and then based on the stopping criterion designs an LQR control based on the learned system and controls it in real-time to achieve desired objectives. 
