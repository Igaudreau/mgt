/* dr3mat.c:  3-dimensional interface ;  matrix utilities 

  This adopt the convention A[row][col], and have rank MAXD4 unless
  qualified by "3".

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#include "csg.h"

void dr3_zerom(float a[MAXD4][MAXD4]);
void dr3_mult(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4], 
  float c[MAXD4][MAXD4]);
void dr3_copy(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4]);
void dr3_vmult(float a[MAXD4][MAXD4], float x[MAXD4],float y[MAXD4]);
void dr3_vmult3(float a[MAXD4][MAXD4], float x[MAXD4], float y[MAXD]);
void dr3_printmatrix(float a[MAXD4][MAXD4]);


/* dr3_zerom - produce "unit" transformation matrix */
void dr3_zerom(float a[MAXD4][MAXD4])
{
   int row,col;

   for (row=0; row < MAXD4; row++) {
      for (col=0; col < MAXD4; col++) {
         a[row][col] = 0.0;
      }
      a[row][row] = 1.0;
   }
}

/* dr3_mult - multiply a*b, result in c */
void dr3_mult(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4], 
  float c[MAXD4][MAXD4])
{
   int row,col,k;
   float d;

   for (row=0; row < MAXD4; row++) {
      for (col=0; col < MAXD4; col++) {
         d = 0.0;
         for (k=0; k < MAXD4; k++) {
            d = d + a[row][k]*b[k][col];
         }
         c[row][col] = d;
      }
   }
}

/* dr3_copy - copy "a" to "b" */
void dr3_copy(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4])
{
   int row,col;

   for (row=0; row < MAXD4; row++) {
      for (col=0; col < MAXD4; col++) {
         b[row][col] = a[row][col];
      }
   }
}

/* dr3_vmult - y=a*x */
void dr3_vmult(float a[MAXD4][MAXD4], float x[MAXD4],float y[MAXD4])
{
   float d;
   int row,col;

   for (row=0; row < MAXD4; row++) {
      d = 0.0;
      for (col=0; col < MAXD4; col++) {
         d = d + a[row][col]*x[col];
      }
      y[row] = d;
   }
}

/* dr3_vmult3 = y=a*x (3-vector) */
void dr3_vmult3(float a[MAXD4][MAXD4], float x[MAXD4], float y[MAXD])
{
   float d;
   int row,col;

   for (row=0; row < MAXD; row++) {
      d = 0.0;
      for (col=0; col < MAXD; col++) {
         d = d + a[row][col]*x[col];
      }
      y[row] = d;
   }
}

/* dr3_printmatrix - print "a" */
void dr3_printmatrix(float a[MAXD4][MAXD4])
{
   int row,col;

   printf("matrix:-\n");
   for (row=0; row < MAXD4; row++) {
      for (col=0; col < MAXD4; col++) {
         printf("%13.6e ", a[row][col]);
      }
      printf("\n");
   }
}
