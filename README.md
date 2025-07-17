# Parallel-CDM

This repository contains a MATLAB-based code and the relevant documentation to solve Computational Damage Mechanics (CDM) problems in parallel. This is a significantly enhanced and updated version of a serial-based CDM code published [here](https://github.com/roshanphilip/UAL-codes). The repo contains the following folders:

**a)	Paper:** contains the article files  <br /> 
**b)	Code:** contains all the [code files](https://github.com/Habiba-Eldababy/Parallel-CDM/tree/main/Code)  <br /> 
**c)	Documentation:** contains a [document](https://github.com/Habiba-Eldababy/Parallel-CDM/blob/main/Documentation/ParallelCDM_Documentation.md) with a thorough explanation of the code functionalities, limitations, how-to-run instructions, and more. <br /> 
**d)	Test Problems:** contains the mesh details and hyperparameters for several [benchmark problems](https://github.com/Habiba-Eldababy/Parallel-CDM/tree/main/Test%20Problems)<br /> 

## Software Features

This code performs parallel‐based simulations of continuum damage mechanics problems. In its current formulation, the following features are available:

- The code models the performance of quasi‐brittle materials, and the user can select between two relevant damage models that capture such behavior: the Mazars [1] and the Geers [2] damage laws.
- The code performs quasi‐static simulations for 2D problems (plane‐strain or plane‐stress).
- The code supports the analysis of both structured and unstructured meshes, with quadrilateral elements that contain 4 Gauss points.
- Two numerical solvers are available, namely the Newton‐Raphson (NR) [3] and the Unified Arc‐Length (UAL) [4] methods.
- The code supports displacement‐driven analysis (Dirichlet BCs are imposed).
- The code can be run either in serial or parallel mode.

## Software Prerequisites

This code requires **MATLAB** and the **Parallel Computing Toolbox** to be installed. The code is compatible with MATLAB version **R2022b and later**. While we recommend using GMSH to generate a new mesh, the user can pick any relevant software of their choice. The generated mesh should be stored in a `.txt` file, modified by the user to match the format that is explained in the Documentation (see the “IV. Model files structure” section).

## Code structure

The code is a compilation of several `.m` files. The main script is `FEM_Main_Script_2D.m`, which is the only code‐related file the user needs to open and edit. In this file the user specifies the <ins>working/saving directories</ins> (line 12) and the <ins>model name</ins> (line 15). This file then calls the solver (Newton‐Raphson or UAL) and the analysis starts.

Each file named `func____.m` serves a specific purpose, and a brief description of its functionality is provided at the beginning of the file. The overall code follows a typical FEM workflow (element‐level procedures, assembly to global matrices, solution for unknown degrees of freedom).

The analysis generates output data at each load increment, which are stored in the `Saved Files` folder. It also prints contour plots for damage, local strain, and nonlocal strain at the last load increment, which are stored in the `Images` folder.

## Contributions and support

We invite you to explore our parallel-based code to substantially accelerate computationally intensive CDM simulations. We welcome contributions from the community – a list of potential aspects for improvement are listed in the Contributing.md file or associated Section X in the Documentation. For further instructions or inquiries, please contact:  <br /> 

* Habiba Eldababy (hed279@nyu.edu)  <br /> 
* Mostafa Mobasher (mm11504@nyu.edu) <br /> 
