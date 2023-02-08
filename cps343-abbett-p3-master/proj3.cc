/* $Smake: mpic++ -Wall -O2 -o %F %f -lhdf5
 * Computes dominant eigenvalue of matrix in HDF5 file via Power Method
 * 
 * This program reads data for a square NxN matrix from an HDF5 data file
 * and uses the Power Method to compute an estimate of the dominant eigenvalue.
 * It assumes the matrix is stored in row-major order (C/C++/Python)
 * 
 * This version computes the linalg operations using for-loops
 * Written by Ben Abbett
 */

//#include <string>
#include <iostream>
//#include "wtime.c"
//#include <stdlib.h>
//#include <unistd.h>
#include <mpi.h>
#include "readMatrixMPI.cc" // The code responsible for reading the matrix
#include <cmath>

using namespace std;

#define IDX(i,j,stride) ((i)*(stride)+(j)) // row-major

/*
 * Display string showing how to run program from command line
 *
 * Parameters:
 *   char* program_name (in)  name of executable
 */
void usage( char* program_name)
{
    fprintf( stderr, "Usage: %s input-file \n", program_name);
}

// Multiply an MxN matrix by a vector (length n) and store the result in y
void mat_vect_mul(double* y, double* A, double* x, int m, int n) {
    for (int ii = 0; ii < m; ii++) {
        y[ii] = 0.0;
        for (int jj = 0; jj < n; jj++) {
            y[ii] += A[IDX(ii,jj,n)] * x[jj];
        }
    }
}

// Compute the dot product of two vectors of size n
double dot_prod(double* x, double* y, int n) {
    double dot_product = 0.0;
    for (int ii = 0; ii < n; ii++) {
        dot_product += x[ii] * y[ii];
    }
    return dot_product;
}

// Divides each element of a vector x of size n by some factor
void scale(double* x, int n, double factor) {
    for (int ii = 0; ii < n; ii++) {
        x[ii] /= factor;
    }
}

int main(int argc, char* argv[]) {
    char* in_name;
    double* x; // the eigenvector
    
    int rank;       // the current processes specific number
    int num_blocks; // the total amount of processes (and thus the amount of matrix blocks)
    double* a_r;    // the matrix block we are computing
    int rows;       // the number of rows in A
    int cols;       // the number of cols in A
    int f_r;        // the first row we are working with in A
    int l_r;        // the last row we are working with in A
    int n_r;        // the number of rows in a_r
    int* F;         // array of offsets to first rows in each process
    int* N;         // array of number of rows assigned to each process

    // Default parameters
    float tolerance = 0.000001; 
    int max_iter = 1000;
    
    // Process command line
    if ( argc != 2 ) {
        usage( argv[0] );
        return EXIT_FAILURE;
    }
    in_name = argv[1];

    // Initialize MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &num_blocks);

    // Read matrix from HDF5 file
    double t1 = MPI_Wtime();   
    readMatrixMPI(in_name, "/A/value", &a_r, &rows, &cols, &f_r, &l_r, rank, num_blocks, MPI_COMM_WORLD );  
    double t2 = MPI_Wtime();
    double read_time = t2 - t1;

    n_r = l_r - f_r + 1;
    
    F = new int[num_blocks];
    F[rank] = f_r; // Set our part of the array
    MPI_Allgather( &f_r,        // The element we are sending
                   1,           // The number of elements we are sending
                   MPI_INT,     // The datatype we send
                   F,           // The element we are receiving
                   1,           // The number of elements we receive from each process
                   MPI_INT,     // The datatype we receive
                   MPI_COMM_WORLD );


    N = new int[num_blocks];
    N[rank] = n_r;  
    MPI_Allgather( &n_r,        // The element we are sending
                   1,           // The number of elements we are sending
                   MPI_INT,     // The datatype we send
                   N,           // The element we are receiving
                   1,           // The number of elements we receive from each process
                   MPI_INT,     // The datatype we receive
                   MPI_COMM_WORLD );

    // Initialize estimate or normalized eigenvector (an array of size n)
    x = new double[rows];
    double nsqrt = sqrt(rows);
    for (int ii = 0; ii < rows; ii++) {
        x[ii] = 1/nsqrt;
    }

    // Main power method loop
    double lambda_new = 0.0;
    double lambda_old = lambda_new + 2 * tolerance;
    double lambda_r = 0.0; // This rank's partial sum of lambda
    double delta = abs(lambda_new - lambda_old);
    int num_iter = 0; // Initialize loop counter
    
    double alpha_r = 0.0; // Result of dot product of y_r with itself
    double alpha = 0.0;   // Agglomerated value that we use to normalize y_r

    double* y_r; // holds the next eigenvector estimate
    y_r = new double[n_r];
    for (int ii = 0; ii < n_r; ii++) {
        y_r[ii] = 0;
    }
    
    t1 = MPI_Wtime();
    while (delta >= tolerance && num_iter <= max_iter) {
        num_iter ++; // Update counter

        // Compute next eigenvector estimate
        mat_vect_mul(y_r, a_r, x, n_r, cols );

        // Save previous eigenvalue estimate
        lambda_old = lambda_new;

        // Compute new estimate  
        lambda_r = dot_prod(&x[f_r], y_r, n_r);

        alpha_r = dot_prod(y_r, y_r, n_r);

        // Reduce all partial lambda sums to one shared sum
        MPI_Allreduce( &lambda_r,       // Value we are sending
                       &lambda_new,     // What we are receiving
                       1,               // Number of elements
                       MPI_DOUBLE,      // Type of element
                       MPI_SUM,         // Operation
                       MPI_COMM_WORLD );

        // Reduce all alpha_r values to one shared alpha value
        MPI_Allreduce( &alpha_r,        // Value we are sending
                       &alpha,          // What we are receiving
                       1,               // Number of elements
                       MPI_DOUBLE,      // Type of element
                       MPI_SUM,         // Operation
                       MPI_COMM_WORLD );

        // Normalize y_r with alpha
        scale(y_r, n_r, sqrt(alpha));

        // Gather y_r to x
        MPI_Allgatherv( y_r,
                        n_r,
                        MPI_DOUBLE,
                        x,
                        N,
                        F,
                        MPI_DOUBLE,
                        MPI_COMM_WORLD );

        delta = abs(lambda_new - lambda_old);        
    }

    t2 = MPI_Wtime();
    double compute_time = t2 - t1;
    double average_time = 0.0;
    // Reduce all compute time values to one value in rank 0 process to be averaged
    MPI_Reduce( &compute_time,        // Value we are sending
                &average_time,        // What we are receiving
                1,                    // Number of elements
                MPI_DOUBLE,           // Type of element
                MPI_SUM,              // Operation
                0,
                MPI_COMM_WORLD );
    average_time /= num_blocks;
    
    if (rank == 0) 
    {
        if (num_iter > max_iter) { 
            cout <<"*** WARNING ****: maximum number of iterations exceeded" << endl; 
        }
        cout << "(eigenvalue = " << lambda_new << "', 'found in " << num_iter << " iterations')" << endl;
        cout << "rank " << rank << ", elapsed read time    =   " << read_time << " seconds" << endl;
        cout << "elapsed compute time =   " << average_time << " seconds" << endl;
    }

    MPI_Finalize();
}
