# pragma once

#include <mpi.h>
#include <dmumps_c.h>

void solveLinearSystemMUMPS(double* matrix_values, int* row_indices, int* col_indices,
                            double* rhs_vector, double* solution, int matrix_size,
                            int num_nonzeros) {
    // Initialize MPI
    int argc = 0;
    char** argv = NULL;
    MPI_Init(&argc, &argv);
    
    // Define MUMPS data structures
    DMUMPS_STRUC_C mumps_data;
    int INFO = 0;
    
    // Initialize MUMPS structure
    //dmumps_c(&mumps_data);
    mumps_data.job = -1;
    mumps_data.par = 1;   // Set parallel flag (par=1)
    
    // Set matrix size and nonzero structure
    mumps_data.n = matrix_size;
    mumps_data.nz = num_nonzeros;
    mumps_data.irn = row_indices;
    mumps_data.jcn = col_indices;
    mumps_data.a = matrix_values;
    
    // Set RHS vector
    mumps_data.rhs = rhs_vector;
    
    // Perform analysis
    mumps_data.job = 1;
    dmumps_c(&mumps_data);
    
    // Perform factorization
    mumps_data.job = 2;
    dmumps_c(&mumps_data);
    
    // Solve linear system
    mumps_data.job = 3;
    dmumps_c(&mumps_data);
    
    // Access solution vector
    double* sol = mumps_data.a;
    for (int i = 0; i < matrix_size; ++i) {
        solution[i] = sol[i];
    }
    
    // Finalize MUMPS
    mumps_data.job = -2;
    dmumps_c(&mumps_data);
    
    // Finalize MPI
    MPI_Finalize();
}
