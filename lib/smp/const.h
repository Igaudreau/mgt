/* const.h:  general constants (including machine specifics)
*/
#ifndef CONST_H
#define CONST_H


#define TRUE 1
#define FALSE 0
//#define SZNAM 32
#define SZNAM 1024
#define SZLIN 256
#define MAXDIM 7

#define CAPTURE
#define bool int
#ifdef DOUBLE
#define uword unsigned long long int
#define float double
#else
#define uword unsigned int
#endif
#define uchar unsigned char


#undef MYMALLOC
#endif
