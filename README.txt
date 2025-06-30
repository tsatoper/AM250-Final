Project Description:
Use OpenMPI to simulate Game of Life. See output in result.out



Compiling:
On Hummingbird, do the following to compile and run my code:

2D decomp:
To change the grid size, edit mpiGOL.f90 and change bigN 
To change number of processors, edit slurm.cmd and change 1234

mpif90 mpiGOL.f90 -o mpiGOL.exe
sbatch slurm.cmd or mpirun –np 1234 mpiGOL.exe
The output data will be in test.out
To look at the data, type more test.out

Column decomp:

To change the grid size, edit column.f90 and change bigM and bigN for a MxN grid 
To change number of processors, edit zcolumn.cmd and change 1234

mpif90 column.f90 -o column.exe
sbatch zcolumn.cmd or mpirun –np 1234 column.exe
The output data will be in testc.out
To look at the data, type more testc.out