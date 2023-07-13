// #include "thesis/sparse-linear-solver.hpp"

#include <mpi.h>
#include <dmumps_c.h>

int main() {

    // Initialize MPI
    int argc = 0;
    char** argv = NULL;
    MPI_Init(&argc, &argv);

    // Example usage
    const int matrix_size = 3;
    int num_nonzeros = 6;
    double matrix_values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    int row_indices[] = {0, 0, 1, 1, 2, 2};
    int col_indices[] = {0, 1, 0, 1, 0, 1};
    double rhs_vector[] = {6.0, 15.0, 24.0};
    double solution[matrix_size];

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

    // // Call MUMPS solver method
    // solveLinearSystemMUMPS(matrix_values, row_indices, col_indices,
    //                        rhs_vector, solution, matrix_size, num_nonzeros);
    
    // // Print solution
    // for (int i = 0; i <matrix_size; ++i) {
    //     std::cout << "x[" << i << "] = " << solution[i] << std::endl;
    // }
    
    // return 0;
}
