#include <cstdio>
#include <cstdlib>
#include <hdf5.h>

/* Check return values from HDF5 routines */
#define CHKERR(status,name) if ( status < 0 ) \
     fprintf( stderr, "Error: nonzero status (%d) in %s\n", status, name )

//----------------------------------------------------------------------------

/*
 * Computes the starting and ending displacements for the ith
 * subinterval in an n-element array given that there are m
 * subintervals of approximately equal size.
 *
 * Input:
 *    int n    - length of array (array indexed [0]..[n-1])
 *    int m    - number of subintervals
 *    int i    - subinterval number
 *
 * Output:
 *    int* s   - location to store subinterval starting index
 *    int* e   - location to store subinterval ending index
 *
 * Suppose we want to partition a 100-element array into 3
 * subintervals of roughly the same size.  The following three
 * pairs of calls find the starting and ending indices of each
 * subinterval:
 *   decompose1d( 100, 3, 0, &s, &e );  (now s =  0, e = 33)
 *   decompose1d( 100, 3, 1, &s, &e );  (now s = 34, e = 66)
 *   decompose1d( 100, 3, 2, &s, &e );  (now s = 67, e = 99)
 *
 * The subinterval length can be computed with e - s + 1.
 *
 * Based on the FORTRAN subroutine MPE_DECOMP1D in the file
 * UsingMPI/intermediate/decomp.f supplied with the book
 * "Using MPI" by Gropp et al.  It has been adapted to use
 * 0-based indexing.
 */
static void decompose1d( int n, int m, int i, int* s, int* e )
{
    const int length  = n / m;
    const int deficit = n % m;
    *s =  i * length + ( i < deficit ? i : deficit );
    *e = *s + length - ( i < deficit ? 0 : 1 );
    if ( ( *e >= n ) || ( i == m - 1 ) ) *e = n - 1;
}

//----------------------------------------------------------------------------
// Opens an HDF5 file, determines the matrix dimensions and what block of
// the matrix this process should load, allocates memory for the matrix
// block, and reads the block from the file.  Size and starting and ending
// row information about the matrix is returned in parameters.

void readMatrixMPI(
    const char* filename, // in  - name of HDF5 file
    const char* path,     // in  - path within HDF5 file to matrix data
    double** a,           // out - location to store pointer to matrix data
    int* rows,            // out - location to store number of rows in matrix
    int* cols,            // out - location to store number of columns in matrix
    int* start,           // out - location to store first row index
    int* end,             // out - location to store last row index
    int rank,             // in  - rank of calling process
    int num_blocks,       // in  - number of blocks of rows
    MPI_Comm comm         // in  - MPI communicator
    )
{
    hid_t plist_id;        // HDF5 id for property list
    hid_t file_id;         // HDF5 id for file
    hid_t dataspace_id;    // HDF5 id for dataspace in file
    hid_t dataset_id;      // HDF5 id for dataset in file
    hid_t memspace_id;     // HDF5 id for dataset in memory
    hsize_t* dims;         // matrix dimensions; converted to block dims
    hsize_t* offset;       // offsets in each dimension to start of block
    herr_t status;         // HDF5 return code
    int ndim;              // number of dimensions in HDF5 dataset
    int partition_dim = 0; // 0 for blocking in groups of rows (1 for columns)

    // Create property list from MPI communicator
    plist_id = H5Pcreate( H5P_FILE_ACCESS );
    if ( plist_id < 0 ) exit( EXIT_FAILURE );
    status = H5Pset_fapl_mpio( plist_id, comm, MPI_INFO_NULL );
    CHKERR( status, "H5Pset_fapl_mpio()" );

    // Open existing shared file using property list
    file_id = H5Fopen( filename, H5F_ACC_RDONLY, plist_id );
    if ( file_id < 0 ) exit( EXIT_FAILURE );
    H5Pclose( plist_id );

    // Open dataset in file
    dataset_id = H5Dopen( file_id, path, H5P_DEFAULT );
    if ( dataset_id < 0 ) exit( EXIT_FAILURE );

    // Determine dataset parameters
    dataspace_id = H5Dget_space( dataset_id );
    ndim = H5Sget_simple_extent_ndims( dataspace_id );
    dims   = new hsize_t [ndim];
    offset = new hsize_t [ndim];

    // Get dimensions for dataset
    ndim = H5Sget_simple_extent_dims( dataspace_id, dims, NULL );
    if ( ndim != 2 )
    {
        if ( rank == 0 )
        {
            fprintf( stderr, "Expected dataspace to be 2-dimensional " );
            fprintf( stderr, "but it appears to be %d-dimensional\n", ndim );
        }
        exit( EXIT_FAILURE );
    }

    // Store matrix dimensions in parameter locations
    *rows = dims[0];
    *cols = dims[1];

    // Determine the location of my rows within the matrix
    decompose1d( dims[partition_dim], num_blocks, rank, start, end );
    dims[partition_dim] = *end - *start + 1;

    // Create memory dataspace
    memspace_id = H5Screate_simple( ndim, dims, NULL );
    if ( memspace_id < 0 ) exit( EXIT_FAILURE );

    // Define hyperslab in file dataspace to make to memory dataspace
    offset[0] = offset[1] = 0;
    offset[partition_dim] = *start;
    status = H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, offset,
                                  NULL, dims, NULL );
    CHKERR( status, "H5Sselect_hyperslab()" );

    // Set transfer mode to be collective
    plist_id = H5Pcreate( H5P_DATASET_XFER );
    if ( plist_id < 0 ) exit( EXIT_FAILURE );

    status = H5Pset_dxpl_mpio( plist_id, H5FD_MPIO_COLLECTIVE );
    CHKERR( status, "H5Pset_dxpl_mpio()" );

    // Allocate memory for portion of matrix and read data from file
    *a = new double [dims[0] * dims[1]];
    status = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, memspace_id,
                      dataspace_id, plist_id, *a );
    CHKERR( status, "H5Dread()" );

    // Close all remaining HDF5 objects
    CHKERR( H5Pclose( plist_id ), "H5Pclose()" );
    CHKERR( H5Sclose( memspace_id ), "H5Sclose()" );
    CHKERR( H5Dclose( dataset_id ), "H5Dclose()" );
    CHKERR( H5Sclose( dataspace_id ), "H5Sclose()" );
    CHKERR( H5Fclose( file_id ), "H5Fclose()" );

    // Clean up
    delete [] dims;
    delete [] offset;
}