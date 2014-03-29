# Parallelize Particle Simulation

We were given a program, in C++, that simulated particle acceleration. The program ran in O(n^2) time where n is the number of particles. 
The reason for this was simple; each particle checked every other particle, too see if they affect each other given their current intermediate distance, 
at each time step. Our main task was to improve this program so that it runs in linear time, O(n). Secondly we were to implement 3 
different parallelized programs using each of the 3 libraries; Pthreads, OpenMP and MPI. All of these 3 should run, as close as possible, 
in time T/p where T is the time taken for our improved linear program and p is the amount of processors used.

## Build
In order to buld our four implementation type the following:

```make serial``` to build the serial implementation.

```make openmp``` to build the openmp implementation.

```make pthreads``` to build the pthreads implementation.

```make mpi``` to build the mpi implementation.

## Run
Run the programs by the following commans:

```./serial``` to run the serial implementation.

```./openmp``` to run the openmp implementation.

```./pthreads``` to run the pthreads implementation.

```mpirun ./mpi``` to run the mpi implementation.

To get the argument list for each program, give the argument ```-h```, for example ```./openmp -h```.

To run an openmp parallel simulation with 4 threads and 40000 particles: ```./openmp -p 4 -n 40000```.

To run an mpi parallel simulation with 4 threads and 40000 particles: ```mpirun -np 4 ./mpi -n 40000```.
