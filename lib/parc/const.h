/* const.h:  general constants (including machine specifics)
*/
#ifndef CONST_H
#define CONST_H
#endif
#define TRUE 1
#define FALSE 0
#define SZNAM 32
#define SZLIN 256
#define MAXDIM 7

#define CAPTURE
#define bool int
#define uchar unsigned char
#ifdef DOUBLE
#define uword unsigned long long int
#define float double
#define lint long long int
#define mpi_uword MPI_UNSIGNED_LONG_LONG
#else
#define uword unsigned int
#define lint int
#define mpi_uword MPI_UNSIGNED
#endif
/* Buffer sizes */
#define DAM_MAX_BLOCK 512
#define DAM_READ 101
#define DAM_WRITE 102
