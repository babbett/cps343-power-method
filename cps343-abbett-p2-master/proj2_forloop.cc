/* $Smake: g++ -O2 -Wall -o %F %f -lhdf5 -fopenmp
 * 
 * Computes dominant eigenvalue of matrix in HDF5 file via Power Method
 * 
 * This program reads data for a square NxN matrix from an HDF5 data file
 * and uses the Power Method to compute an estimate of the dominant eigenvalue.
 * It assumes the matrix is stored in row-major order (C/C++/Python)
 * 
 * This is a parallelized version of the original (project 1) for-loop program.
 * Written by Ben Abbett
 */

#include <omp.h>
#include <iostream>
#include "wtime.c"
#include <hdf5.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
extern "C" {
    #include <cblas.h>
}

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

void mat_vect_mul(double* y, double* A, double* x, int n, int m) {
    for (int ii = 0 ; ii < m; ii++) {
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
    //norm = sqrt(norm);
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

/*
* Computes the starting and ending displacements for the ith
* subinterval in an n - element array given that there are m
* subintervals of approximately equal size .
*
* Input : int n - length of array ( array indexed [0]..[ n -1])
* int m - number of subintervals
* int i - subinterval number
*
* Output : int * s - location to store subinterval starting index
* int * e - location to store subinterval ending index
*
* Suppose we want to partition a 100 - element array into 3 subintervals
* of roughly the same size . The following three pairs of calls find
* the starting and ending indices of each subinterval :
* decompose1d ( 100 , 3 , 0 , &s , & e ); ( now s = 0 , e = 33)
* decompose1d ( 100 , 3 , 1 , &s , & e ); ( now s = 34 , e = 66)
* decompose1d ( 100 , 3 , 2 , &s , & e ); ( now s = 67 , e = 99)
*
* The subinterval length can be computed with e - s + 1.
*
* Based on the FORTRAN subroutine MPE_DECOMP1D in the file
* UsingMPI / intermediate / decomp . f supplied with the book " Using MPI "
* by Gropp et al . It has been adapted to use 0 - based indexing .
*/
void decompose1d ( int n , int m , int i , int * s , int * e ) {
    const int length = n / m ;
    const int deficit = n % m ;
    * s = i * length + ( i < deficit ? i : deficit );
    * e = * s + length - ( i < deficit ? 0 : 1 );
    if ( ( * e >= n ) || ( i == m - 1 ) ) * e = n - 1;
}

int main(int argc, char* argv[]) {
    char* in_name;
    double* a; // the matrix we are computing
    int     n; // the N of the NxN matrix
    double* x; // the eigenvector

    // Default parameters (will eventually be read through command line)
    float tolerance = 0.000001;
    int max_iter = 1000;

    // Process command line
    if ( argc != 2 ) {
        usage( argv[0] );
        return EXIT_FAILURE;
    }
    in_name = argv[1];

    // Read matrix from HDF5 file
    double t1 = wtime();
    readMatrix(in_name, "/A/value", &a, &n);
    double t2 = wtime();
    double read_time = t2 - t1;

    // Initialize estimate of normalized eigenvector (an array of size n)
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

    t1 = wtime();
    double norm; // the shared variable for accumulation

    // Parallelize the power method
    #pragma omp parallel
    { 
        // Get this thread's number and assign it rows of the matrix 
        const int num_threads = omp_get_max_threads();
        const int rank = omp_get_thread_num();
        int thread_num_rows; // The number of rows that this thread is calculating
        int first_mat_elem;
        int last_mat_elem; 

        // Distribute rows of the matrix to threads approximately evenly
        decompose1d(n, num_threads, rank, &first_mat_elem, &last_mat_elem);

        first_mat_elem *= n;
        last_mat_elem = last_mat_elem * n + n  -1;
        thread_num_rows = (last_mat_elem - first_mat_elem + 1 ) / n;

        const int first_vec_elem = first_mat_elem / n;

        double* y_chunk = new double[thread_num_rows]; // The chunk of y used by this thread
        double* x_chunk = &x[first_vec_elem];          // the chunk of x corresponding to y, used for computing partial lambda
        double* a_chunk = &a[first_mat_elem];          // should be the chunk of a we need
        double lambda_new_chunk;                       // the partial sum of lambda for the chunk we are computing
        double norm_chunk;                             // the part of 'z', accumulated and used for scaling y_chunk

        while (delta >= tolerance && num_iter <= max_iter) {
            #pragma omp single
            {
                num_iter ++; // Update counter
                lambda_old = lambda_new; // Save previous eigenvalue estimate
                lambda_new = 0;
            }

            // Compute next eigenvector estimate
            // We pass pointers to the correct slice of the matrix and vector so that only these parts get multiplied
            mat_vect_mul(y_chunk, a_chunk, x, n, thread_num_rows);

            // Compute new eigenvalue estimate
            lambda_new_chunk = dot_prod(x_chunk, y_chunk, thread_num_rows);

            #pragma omp critical
            lambda_new += lambda_new_chunk;


            // Normalize eigenvector estimate
            norm_chunk = two_norm(y_chunk, thread_num_rows);

	        #pragma omp single
            norm = 0;

            #pragma omp critical
            norm += norm_chunk;

            #pragma omp barrier
            scale(y_chunk, thread_num_rows, sqrt(norm));

            // Update normalized eigenvector estimate
            #pragma omp critical
            copy(&x[first_vec_elem], y_chunk, thread_num_rows);

            delta = abs(lambda_old - lambda_new);
        }
    }
    t2 = wtime();
    double compute_time = t2 - t1;
    if (num_iter > max_iter) {
        cout <<"*** WARNING ****: maximum number of iterations exceeded" << endl;
    }
    cout << "(eigenvalue = " << lambda_new << "', 'found in " << num_iter << " iterations')" << endl;
    cout << "elapsed read time    =   " << read_time << " seconds" << endl;
    cout << "elapsed compute time =   " << compute_time << " seconds" << endl;
}
