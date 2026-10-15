/* parc.h:  header for PARC utility system ;
*/
#ifndef PARC_H
#define PARC_H
#include "const.h"

/* Constants:
*/
#define MAXPE 256
#define PARC_INT_TYPE 1
#define PARC_FLOAT_TYPE 2
#define PARC_STRING_TYPE 3
#define PARC_BOOL_TYPE 4
#define MAXSTK 31

/* Types:
*/
typedef struct CacheType *CachePtr;
typedef struct CacheType {
/* Is this structure in use */
   bool used;
/* number of PEs */
   int npe;
/* local node limits */
   int lo;
   int me;
   int hi;
/* global node limits */
   int nnpe;
   int nlev;
/* element tree */
   int ntree;
   int *tree;
/* element linked-list */
   int nelem;
   int *elem_far;
   int *elem_near;
/* pointer list */
   int nptr;
   int **ptr_adr;
   int *ptr_far;
   int parent;
/* array details */
   int type;
   int ndim;
   int dims[MAXDIM];
   void *array;
/* tide marks */
   int tider;
   int tidew;
} CacheStruct;

/* Templates:
*/
void parc_args(int, char*[]);
void parc_cerror(char[], int*);
void parc_cwarn(char[]);
void parc_register__(int*, int*, int*);
uword* parc_xdataw(int, int[], uword[], int[], int*);

#endif
