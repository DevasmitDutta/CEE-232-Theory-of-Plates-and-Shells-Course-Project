### This repository contains the finite element simulation routines for generating the deflected shapes of plate using planar elements, using Mindlin plate elements and Kirchoff plate elements. This project was carried out as a course project at UCLA for CEE 232 Theory of Plates and Shells Course

# Instructions to run the codes

The main routine for running the analysis is the `Program.m` file. All other relevant files contain the sub-routines for appropriate functionalities such as generating the element stiffness matrix, transforming into global level and post-processing. Depending on the parameters such as mesh discretisations, boundary conditions etc, the deflected shapes of the plate will be generated. Also, the impact of shear-locking can be controlled using a reduced order numerical integration rule using the `Int` variable. 

Below are some deflected shapes generated for visualisation. 
