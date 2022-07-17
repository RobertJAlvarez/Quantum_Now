## Quantum_Now

#Purpose/goal
Quantum Now is an inexpensive based computational learning software aim from high-school to graduate students and K-12 teachers. This would present an affordable quantum simulation software as assistance for early science education.

#Restrictions
The code was written with some restrictions: only addition, subtraction, multiplication, conditional, and iterations were allowed. The main implications of this are that division, sin/cosine, and other basic functions had to be implemented, no calculus was allowed, and no packages for matrix operations were used.

#Current simulations
Current implemented quantum simulations include the particle in a box, harmonic oscillators, hydrogen stark effect, and spin on a ring.

#Structure
The software is divided in 4 main files located at Fortran_2008_Code and Python_Code.
File                                       | Content
-------------------------------------------|-----------------------------------------------------
1: FortranFunctions.f08/PythonFunctions.py | modulus (A%B), absolute value, division, sine, cosine and square root.
2: ArrayFunctions.f08/.py                  | Inverse of an NxN matrix, Jacobian 2x2 diagonalization, and an NxN matrix diagonalization.
3: Applications.f08/.py                    | Current simulations mention above.
4: Main.f08/.py                            | Control the flow of the program.
5: TestingFunctions.f08/.py                | Tests for every function in any of the first three file listed.

#Currently developing
The original code located at Old_Fortran_Code was written in Fortran 77 during a special topics graduate level physics class given by Dr. Pederson at The University of Texas at El Paso during Fall 2021. While the class was given I, Robert Alvarez, was translating the code to Fortran 2008 and giving a more structure way to divide the code as shown in Fortran_2008_Code folder. Now, we are translating the code into Python and creating a [LATEX file](https://www.overleaf.com/read/wrxjwxgrzbtm) guide that explain everything on this code like the mathematics use to develop the algorithms, pseudocode, code, and a Makefile explanation for the Fortran 2008 version.

