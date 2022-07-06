# Quantum_Now

Purpose/goal:
Quantum Now is an inexpensive based computational learning software aim from high-school to graduate students and K-12 teachers. This would present an affordable quantum simulation software as assitance for early sciecne education.

Restrictions:
The code was written with some restrictions: only addition, subtraction, multiplication, conditional, and iterations were allowed. The main implications of this are that division, sin/cosine, and other basic functions had to be implemented, no calculus was allowed, and no packages for matrix operations were used.

Current simulations:
Current implemented quantum simulations include the particle in a box, harmonic oscillators, hydrogen stark effect, and spin on a ring.

Structure:
The software is divided in 4 main files located at Fortran_2008_Code and Python_Code.
1. FortranFunctions.f08/PythonFunctions.py contains the following basic functions: modulus (A%B), absolute value, division, sine, cosine and square root.
2. ArrayFunctions.f08/.py contains matrix operations such as inverse of a n by n matrix, Jacobian 2 by 2 diagonalization, and an n by n matrix diagonalization.
3. Applications.f08/.py contains the 4 different simulations mention above.
4. Main.f08/.py takes care of printing the available options, read input from user, and perform the respective task.

Currently developing:
The original code located at Old_Fortran_Code was written in fortran 77 during a special topics graduate level physics class given by Dr. Pederson at The University of Texas at El Paso during Fall 2021. While the class was given I, Robert Alvarez, was translating the code to Fortran 2008 and giving a more structure way to divide the code as shown in Fortran_2008_Code folder. Now, we are translating the code into Python to generate an easier to understand code and creating a LATEX file guide that explain everything on this code like the mathematics use to develop the algorithms, pseudocode, code, and for the fortran version a Makefile explanation.
