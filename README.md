# Parallel-CDM

This repository contains a parallel computational damage mechanics (CDM) MATLAB code. We utilize parallelization techniques in order to significantly accelerate CDM simulations. This an updated parallel code, see the original https://github.com/roshanphilip/UAL-codes. The package can be used with UAL solver (https://doi.org/10.1007/s00466-024-02473-5) or Newton-Raphson, and for local or non-local gradient damage.

We invite you to try out different continuum damage mechanics problems using this parallel code.

Steps to run the code:
1. Download code and mesh files. To change the problem parameters, open the corresponding parameter file of the mesh and adjust. Create a folder to save the data files at each increment.
2. Open main script file - "FEM_Main_Script_2D.m"
3. Enter the file path and results storage folder path
4. Create parallel threads pool depending on your machine or leave default pool initialization.
5. Run main script file - "FEM_Main_Script_2D.m". The code outputs data files at each increment and a final plot of the damage contours.
