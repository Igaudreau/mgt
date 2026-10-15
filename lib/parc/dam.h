/* dam.h:  header for DAM utility system (POSIX); */
#include "const.h"
/* Types: */
typedef struct DamType *DamPtr;
typedef struct DamType {
   uword *mb;
   lint *ib;
   float *fb;
   int mm;
   int mp;
   int mn;
   int st;
} DamStruct;
/* Templates:
*/
bool dam_init(DamPtr, bool*);
bool dam_alloc(int, DamPtr, bool*);
bool dam_free(DamPtr, bool*);
bool dam_grow(int, DamPtr, bool*);
bool dam_space(int, DamPtr, bool*);
bool dam_rewind(DamPtr, int, bool*);
