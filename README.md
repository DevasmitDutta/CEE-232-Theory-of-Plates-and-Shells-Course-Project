## This repository contains the finite element simulation routines for generating the deflected shapes of plate using planar elements and using Mindlin plate elements and Kirchoff plate elements. This project was carried out as a course project at UCLA for CEE 232 Theory of Plates and Shells Course in Fall'23

## Instructions to run the codes

The main routine for running the analysis is the `Program.m` file. All other relevant files contain the sub-routines for appropriate functionalities such as generating the element stiffness matrix, transforming into global level and post-processing. Depending on the parameters such as mesh discretisations, boundary conditions etc, the deflected shapes of the plate will be generated. Also, the impact of shear-locking can be controlled using a reduced order numerical integration rule using the `Int` variable. 

Below are some deflected shapes generated for visualisation. 

![untitled](https://github.com/DevasmitDutta/CEE-232-Theory-of-Plates-and-Shells-Course-Project/assets/76597282/31c68774-ac97-4048-a9d5-13b98990bcf8)
![untitled-1](https://github.com/DevasmitDutta/CEE-232-Theory-of-Plates-and-Shells-Course-Project/assets/76597282/e1405d03-f99f-4021-81e3-8764e0d965f4)
![untitled-2](https://github.com/DevasmitDutta/CEE-232-Theory-of-Plates-and-Shells-Course-Project/assets/76597282/06f28a2c-1619-4940-a0a8-ca92713d2303)
![untitled-3](https://github.com/DevasmitDutta/CEE-232-Theory-of-Plates-and-Shells-Course-Project/assets/76597282/acd2214e-96f1-4ebe-89a7-e60b2f34f0bb)
