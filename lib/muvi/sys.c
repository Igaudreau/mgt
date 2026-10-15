/* sys.c:  various system-related functions */
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"


void sys_date_string(int str_len, char str[]);
void sys_time_string(int str_len, char str[]);
void sys_date(int str_len, char str[]);
extern void ccopy(char a[], int offset, int n, char b[]);

void sys_date_string(int str_len, char str[])
{
   ccopy("<unknown>", 0, str_len, str);
}

void sys_time_string(int str_len, char str[])
{
   ccopy("<unknown>", 0, str_len, str);
}

void sys_date(int str_len, char str[])
{
   ccopy("<unknown>", 0, str_len, str);
}
