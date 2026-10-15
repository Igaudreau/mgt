/* fchr.c:  string conversion, Fortran-C and C-Fortran 

     These procedures are callable from C only.
*/
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"


/* EXPLICIT DECLARATIONS ADDED BY BM */
void fkopy(char a[], int lena, int lenb, char b[]);
void cfopy(char a[], int lenb, char b[]);


/* fkopy - Fortran CHARACTER string to C char string */
// MIR fkopy(a, lena, lenb, b)
void fkopy(a, lena, lenb, b)
   char a[],b[];
   int lena,lenb;
{
   int i,n;

   n = lena;
   for (i=lena-1; i >= 0; i--) {
      if (a[i] != ' ') {
         break;
      }
      n = i;
   }
   if (n > (lenb-1)) {
      n = lenb-1;
   }
   for (i=0; i < n; i++) {
      b[i] = a[i];
   }
   b[n] = '\0';
}

/* cfopy - C char string to Fortran CHARACTER string */
// MIR cfopy(a, lenb, b)
void cfopy(a, lenb, b)
   char a[],b[];
   int lenb;
{
   int i;
   bool ended;

   ended = FALSE;
   for (i=0; i < lenb; i++) {
      if (a[i] == '\0') {
         ended = TRUE;
      }
      if (ended) {
         b[i] = ' ';
      } else {
         b[i] = a[i];
      }
   }
}

