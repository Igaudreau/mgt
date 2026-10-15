/* smp.h:  SMP constants 
*/
#ifndef SMP_H
#define SMP_H

#include "const.h"

/* Templates :
*/
void smp_clock(float*);
bool smp_send(int, int, uword[]);
bool smp_recv(int, int, uword[]);
bool smp_register(int, char*[], int*, int*);
void smp_exit();
void smp_abort();
void smp_barrier();
#endif
