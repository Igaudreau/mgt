/* nam.c:  name table utilities */
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"

void nam_alloc(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);
void nam_find(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);
void nam_free(char name[], int maxnames, bool used[], 
  char names[][SZNAM]);

/* EXPLICIT DECLARATIONS ADDED BY BM */
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);



void nam_alloc(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset)
{
   int i;

   *offset = -1;
   for (i=0; i < maxnames; i++) {
      if (!used[i]) {
         ccopy(name, 0, SZNAM, names[i]);
         used[i] = TRUE;
         *offset = i;
         break;
      }
   }
}

void nam_find(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset)
{
   int i;

   *offset = -1;
   for (i=0; i < maxnames; i++) {
      if (used[i]) {
         if (csame(name, names[i])) {
            *offset = i;
            break;
         }
      }
   }
}

void nam_free(char name[], int maxnames, bool used[], 
  char names[][SZNAM])
{
   int i;

   for (i=0; i < maxnames; i++) {
      if (used[i]) {
         if (csame(name, names[i])) {
            used[i] = FALSE;
            break;
         }
      }
   }
}
