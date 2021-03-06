# cubesolver

This cube solver is a Cube Solver that solves any valid colouring of a Cube.
It is derrived from Herbert Kociemba two-phase algorithm
https://github.com/muodov/kociemba

It constains one .cpp file cubesolver.cpp and header files.

cubesolver.cpp has two main functions:
```
   void initarrays()
   int* solvecube(char*, long, int)
```
## Changes

Moved logic to one .cpp file

Used Initialised arrays for pruning and move tables.
The only arrays that are initialise by initarrays() method are
URFtoDLF_Move and URtoDF_Move.

Parameters to  `int* solvecube(char*, long, int)` function is now char[54] array (values from 0 to 5 corresponding to 6 faces colours) instead of char sequence.
Second parameter is timeout value and third is maximum depth.
This function returns pointer to array with rotations represented by int elements.

The rotations correspond to following int values:
U : 6,
R : 3,
F : 9,
D : -4,
L : -1,
B : -7,
U' : -6,
R' : -3,
F' : -9,
D' : 4,
L' : 1,
B' : 7
For U2, R2, F2, D2, L2, B2 moves, they now become two moves i.e. U2 becomes U,U.

Used casts (short) and (int) to resolve compiler errors.

## Usage

There is no main() function in this project.

To use it there should be created a wrapper (i.e. main function that first calls initarrays() and then solvecube(char*, long, int))


