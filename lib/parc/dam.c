/* dam.c:  Data Memory Buffer module (STDC,PARC version)

   Copyright (C) 1999-2002 Mantis Numerics Limited.

   All rights reserved.  This is UNPUBLISHED PROPRIETARY SOURCE CODE
   of Mantis Numerics Limited;  the contents of this file may not be
   disclosed to third parties, copied or duplicated in any form,
   in whole or in part, without the prior written permission of
   Mantis Numerics Limited.

   Mantis Numerics Limited, 46 The Calls, Leeds LS2 7EY, UK.
*/

/* Includes :
*/
#include <stdio.h>
#include <stdlib.h>
#include "dam.h"

#ifdef DOUBLE
#define lint long long int
#define float double		// BM
#else
#define lint int
#endif
/* Constants:
*/
#define BLOCK_GROWTH 1000
/*
   Memory Management Functions
*/
/* dam_alloc:  allocate a memory buffer */
bool dam_alloc(int mm, DamPtr s, bool *error) {
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
/* No need? */
   if (s->mm == mm) {
      s->mp = 0;
      s->mn = 0;
   } else {
/*    Free up any existing block */
      if (s->mb != (uword*)NULL) {
         dam_free(s, error);
      }
/*    Grab raw memory block */
      if ((s->mb=(uword*)malloc((size_t)(mm*sizeof(uword))))==(uword*)NULL) {
         fprintf(stdout, "error:  allocating memory block (dam_alloc)\n");
         *error = TRUE;
      } else {
         s->ib = (lint*)s->mb;
         s->fb = (float*)s->mb;
         s->mm = mm;
         s->mp = 0;
         s->mn = 0;
      }
   }
   return(!(*error));
}
/* dam_free:  free a memory buffer */
bool dam_free(DamPtr s, bool *error) {
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
/* Free raw memory block */
   if (s->mb != (uword*)NULL) {
      free(s->mb);
   }
/* Update structure */
   s->mb = (uword*)NULL;
   s->ib = (lint*)NULL;
   s->fb = (float*)NULL;
   s->mm = 0;
   s->mp = 0;
   s->mn = 0;
   return(!(*error));
}
/* dam_init:  initialize a memory buffer */
bool dam_init(DamPtr s, bool *error) {
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
   s->mb = (uword*)NULL;
   s->ib = (lint*)s->mb;
   s->fb = (float*)s->mb;
   s->mm = 0;
   s->mp = 0;
   s->mn = 0;
   s->st = -1;
   dam_alloc(DAM_MAX_BLOCK, s, error);
   return(!(*error));
}
/* dam_grow:  grow a memory buffer */
bool dam_grow(int mm, DamPtr s, bool *error) {
   uword *newb;
   int newm;
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
/* Grab new raw memory block */
   newm = s->mm + mm;
   if ((newb=(uword*)realloc((void*)s->mb,(size_t)(newm*sizeof(uword))))==(uword*)NULL) {
      fprintf(stdout, "error:  allocating memory block (dam_grow)\n");
      *error = TRUE;
   } else {
      s->mb = newb;
      s->ib = (lint*)s->mb;
      s->fb = (float*)s->mb;
      s->mm = newm;
   }
   return(!(*error));
}
/* dam_space:  check space a memory  (grow) */
bool dam_space(int nw, DamPtr s, bool *error) {
   int need;
/* Return in error state */
   if (*error) {
      return(!(*error));
   }
   if (s->mp + nw > s->mm) {
      if (nw < BLOCK_GROWTH) {
         need = BLOCK_GROWTH;
      } else {
         need = nw;
      }
      dam_grow(need, s, error);
   }
   return(!(*error));
}
