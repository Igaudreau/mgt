/* machine.h:  Machine specifics SGI(GL,X11)
*/
#undef SGI
#undef DEC
#undef VAX
#undef SUN
#undef FX
/* Backend specifics:
*/
#define PS
#undef GL
#undef GLX
#define X11
/* Etc :
*/
#define RANCON (float)(RAND_MAX)
#define cfree(ptr,nword,wordsize) (free(ptr))
