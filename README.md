# mpi-afis
Parallel framework for fingerprint identification.
Source code for the framework described in the following papers:
- D. Peralta, I. Triguero, R. Sanchez-Reillo, F. Herrera, J.M. Benítez. **Fast Fingerprint Identification for Large Databases**. *Pattern Recognition* 47:2 (2014) 588–602. doi: 10.1016/j.patcog.2013.08.002
- D. Peralta, I. Triguero, S. García, F. Herrera, J.M. Benítez. **DPD-DFF: A Dual Phase Distributed Scheme with Double Fingerprint Fusion for Fast and Accurate Identification in Large Databases**. Submitted to *Information Fusion*.

The project can be compiled by executing:

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
