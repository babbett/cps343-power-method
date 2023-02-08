/* $Smake: mpic++ -Wall -O2 -o %F %f -fopenmp 
 * Computes dominant eigenvalue of matrix in HDF5 file via Power Method
 * 
 * This program reads data for a square NxN matrix from an HDF5 data file
 * and uses the Power Method to compute an estimate of the dominant eigenvalue.
 * It assumes the matrix is stored in row-major order (C/C++/Python)
 * 
 * This version computes the linalg operations using for-loops
 * Written by Ben Abbett
 */

#include <omp.h>
#include <iostream>
#include "wtime.c"
// #include <hdf5.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include "readMatrixMPI.cc" // The code responsible for reading the matrix
#include <cmath>
// extern "C" {
//     #include <cblas.h>
// }

using namespace std;

/* Check return values from HDF5 routines */
#define CHKERR(status,name) if ( status ) \
     fprintf( stderr, "Warning: nonzero status (%d) in %s\n", status, name )
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

void readMatrix(const char* fname, const char* name, double** a, int* n)
{
    hid_t   file_id, dataset_id, file_dataspace_id, dataspace_id;
    herr_t  status;
    hsize_t* dims;
    int rank;
    int ndims;
    hsize_t num_elem;

    /* Open existing HDF5 file */
    file_id = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT );

    /* Open existing first dataset */
    dataset_id = H5Dopen( file_id, name, H5P_DEFAULT );

    /* Determine dataset parameters */
    file_dataspace_id = H5Dget_space( dataset_id );
    rank = H5Sget_simple_extent_ndims( file_dataspace_id );
    //dims = (hsize_t*) malloc( rank * sizeof( hsize_t ) );
    dims = new hsize_t[rank];
    ndims = H5Sget_simple_extent_dims( file_dataspace_id, dims, NULL );
    if ( ndims != rank )
    {
        fprintf( stderr, "Warning: expected dataspace to be dimension " );
        fprintf( stderr, "%d but appears to be %d\n", rank, ndims );
    }

    /* Allocate matrix */
    num_elem = H5Sget_simple_extent_npoints( file_dataspace_id );
    //*a = (double*) malloc( num_elem * sizeof( double ));
    *a = new double[num_elem];
    *n = dims[0]; /* NxN is given, so only need 1 dimension*

    /* Create dataspace */
    dataspace_id = H5Screate_simple( rank, dims, NULL );

    /* Read matrix data from file */
    status = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, dataspace_id,
                      file_dataspace_id, H5P_DEFAULT, *a );
    CHKERR( status, "H5Dread()" );

    /* Close resources */
    status = H5Sclose( dataspace_id ); CHKERR( status, "H5Sclose()" );
    status = H5Sclose( file_dataspace_id ); CHKERR( status, "H5Sclose()" );
    status = H5Dclose( dataset_id ); CHKERR( status, "H5Dclose()" );
    status = H5Fclose( file_id ); CHKERR( status, "H5Fclose()" );
    free( dims );
}

void mat_vect_mul(double* y, double* A, double* x, int n) {
    #pragma omp parallel for
    for (int ii = 0; ii < n; ii++) {
        y[ii] = 0.0;
        for (int jj = 0; jj < n; jj++) {
            y[ii] += A[IDX(ii,jj,n)] * x[jj];
        }
    }
}

double dot_prod(double* x, double* y, int n) {
    double dot_product = 0.0;
    for (int ii = 0; ii < n; ii++) {
        dot_product += x[ii] * y[ii];
    }
    return dot_product;
}

double two_norm(double* x, int n) {
    double norm = 0.0;
    for ( int ii = 0; ii < n; ii++) {
        norm += x[ii] * x[ii];
    }
    norm = sqrt(norm);
    return norm;
}

// Divides each element of a vector x of size n by some value
void scale(double* x, int n, double factor) {
    for (int ii = 0; ii < n; ii++) {
        x[ii] /= factor;
    }
}

void copy(double* x, double* y, int n) {
    for (int ii = 0; ii < n; ii++) {
        x[ii] = y[ii];
    }
}

int main(int argc, char* argv[]) {
    char* in_name;
    double* a; // the matrix we are computing
    int     n; // the N of the NxN matrix
    double* x; // the eigenvector
    
    int rank;  // the current processes specific number
    int num_blocks; // the total amount of processes (and thus the amount of matrix blocks)
    const int tag = 42; // the tag
    double* a_r; // the matrix block we are computing
    int rows_r; // the number of rows in a_r
    int cols_r; // the number of cols in a_r
    int f_r;    // the first row we are working with in A
    int l_r;    // the last row we are working with in A
    int n_r;    // the number of rows in a_r (should just = rows_r)? DELETE THIS AFTER TESTING PROBABLY 
    MPI_Status status;

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
    readMatrixMPI(in_name, "/A/value", &a_r, &rows_r, &cols_r, &f_r, &l_r, rank, num_blocks, MPI_COMM_WORLD );
    // readMatrix(in_name, "/A/value", &a, &n);
    double t2 = MPI_Wtime();
    double read_time = t2 - t1;

    n_r = l_r - f_r + 1;
    printf("Rank %d (of %d) has %d rows, %d cols, %d first row and %d last row, for total size %d\n", rank, num_blocks, rows_r, cols_r, f_r, l_r, n_r);
    /*
    // Initialize estimate or normalized eigenvector (an array of size n)
    //x = (double*) malloc( n * sizeof( double ));
    x = new double[n];
    double nsqrt = sqrt(n);
    for (int ii = 0; ii < n; ii++) {
        x[ii] = 1/nsqrt;
    }

    // Main power method loop
    double lambda_new = 0;
    double lambda_old = lambda_new + 2 * tolerance;
    double delta = abs(lambda_new - lambda_old);
    int num_iter = 0; // Initialize loop counter
    
    double* y; // holds the next eigenvector estimate
    //y = (double*) malloc( n * sizeof( double ));
    y = new double[n];
    for (int ii = 0; ii < n; ii++) {
        y[ii] = 0;
    }

    t1 = wtime();
    while (delta >= tolerance && num_iter <= max_iter) {
        num_iter ++; // Update counter

        // Compute next eigenvector estimate
        mat_vect_mul(y, a, x, n );
       
        // Save previous eigenvalue estimate
        lambda_old = lambda_new;

        // Compute new estimate
        lambda_new = dot_prod(x, y, n);
        
        // Update eigenvector estimate
        copy(x, y, n);

        // Normalize eigenvector estimate
        double norm = two_norm(x, n);

        scale(x, n, norm);

        delta = abs(lambda_old - lambda_new);
    }

    t2 = wtime();
    double compute_time = t2 - t1;
    if (num_iter > max_iter) { 
        cout <<"*** WARNING ****: maximum number of iterations exceeded" << endl; 
    }
    cout << "(eigenvalue = " << lambda_new << "', 'found in " << num_iter << " iterations')" << endl;
    */
    cout << << "rank " << rank << ", elapsed read time    =   " << read_time << " seconds" << endl;
    //cout << "elapsed compute time =   " << compute_time << " seconds" << endl;
}
