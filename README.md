# mpi-afis
Parallel framework for fingerprint identification.

The project can be executed by executing:

`make`

Requirements:
- MPI (tested with OpenMPI 1.8.7), ideally a thread-safe compilation.
- OpenMP

Example executions:

```
./DPDDFF -h
mpirun -np 2 ./genericMatching -a mcc -k 10 -s m -t template_files.dat -i input_files.dat -N 8 -C LSSR
mpirun -np 2 ./DPDDFF -a jiang,mcc -k 10 -s r -r 10 -t template_files_F1.dat,template_files_F2.dat -i input_files_F1.dat,input_files_F2.dat
```
