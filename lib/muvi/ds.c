/* ds.c:  embryonic data system */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#define MAXPATH 4
#define MAXFILE 12
#define BINEXT "CSB"
#define TXTEXT "CST"
#define REAL "real"
#define INT "int"
#define STRING "string"
#define CHAR "char"
#define BOOL "bool"

static int npath,nvalb,level[MAXFILE],counts[MAXFILE];
static char paths[MAXPATH][SZLIN],typeb[SZNAM],nameb[SZNAM],
  files[MAXFILE][SZLIN],method[SZNAM];
static bool used[MAXFILE],binary[MAXFILE],writen[MAXFILE];
static FILE *fp[MAXFILE];

/* EXPLICIT DECLARATIONS ADDED BY BM */
int ds_prepare();
int ds_flush();
int ds_path(char path[], int *ios);
int ds_method(char met[], int *ios);
int ds_create(char file[], int *unit, int *ios);
int ds_open(char file[], int *unit, int *ios);
int ds_close(int unit, int *ios);
int ds_index(char name[], int no, int lencomb, char comb[]);
int ds_put_begin(int unit, char name[], int *ios);
int ds_put_end(int unit, char name[], int *ios);
int ds_get_begin(int unit, char name[], int *ios);
int ds_get_end(int unit, char name[], int *ios);
int ds_put0i(int unit, char name[], int value, int *ios);
int ds_put1i(int unit, char name[], int nvalue, int value[], int *ios);
int ds_put0r(int unit, char name[], float value, int *ios);
int ds_put1r(int unit, char name[], int nvalue, float value[], int *ios);
int ds_put0b(int unit, char name[], bool value, int *ios);
int ds_put1b(int unit, char name[], int nvalue, bool value[], int *ios);
int ds_put0s(int unit, char name[], char value[], int *ios);
int ds_put0c(int unit, char name[], char value, int *ios);
int ds_put1c(int unit, char name[], int nvalue, char value[], int *ios);
bool ds_get(int unit, char name[], char type[], int *ios);
int ds_get0i(int unit, char name[], int *value, int *ios);
int ds_get1i(int unit, char name[], int nvalue, int value[], int *ios);
int ds_get0r(int unit, char name[], float *value, int *ios);
int ds_get1r(int unit, char name[], int nvalue, float value[], int *ios);
int ds_get0b(int unit, char name[], bool *value, int *ios);
int ds_get1b(int unit, char name[], int nvalue, bool value[], int *ios);
int ds_get0s(int unit, char name[], char value[], int *ios);
int ds_get0c(int unit, char name[], char *value, int *ios);
int ds_get1c(int unit, char name[], int nvalue, char value[], int *ios);
int merror(char text[]);
int ds_filex(char path[], char name[], char ext[], int lenfile, char file[]);
int ds_check_unit(int unit, int *ios);
int ds_put_header(int unit, char name[], char type[], int nvalue, int *ios);
int ds_get_prelude(int unit, char name[], char type[], int nvalue, int *ios);
int ds_get_header(int unit, int *ios);
int ds_check_header(char name[], char type[], int nvalue, int *ios);
int ds_indent(int unit, int indent, int *ios);
int ds_read_error(int unit, char proc[], int *ios);
int ds_write_error(int unit, char proc[], int *ios);
int ds_put1c_binary(int unit, int nvalue, char value[], int *ios);
int ds_get1c_binary(int unit, int *nvalue, char value[], int *ios);
int ds_put1i_binary(int unit, int nvalue, int value[], int *ios);
int ds_get1i_binary(int unit, int nvalue, int value[], int *ios);
int ds_put1r_binary(int unit, int nvalue, float value[], int *ios);
int ds_get1r_binary(int unit, int nvalue, float value[], int *ios);
int ds_put1b_binary(int unit, int nvalue, bool value[], int *ios);
int ds_get1b_binary(int unit, int nvalue, bool value[], int *ios);
int ds_put0s_binary(int unit, char value[], int *ios);
int ds_get0s_binary(int unit, int lenvalue, char value[], int *ios);

extern void ccat(char a[], int lenb, char b[]);
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);
extern void slower(char a[]);
extern int slen(char s[]);


int ds_prepare(){
   int i;

   for (i=0; i < MAXFILE; i++) {
   }
   ccopy("binary", 0, SZNAM, method);
   npath = 1;
   ccopy("", 0, SZLIN, paths[0]);
   return 0;
}

int ds_flush() {
   int i,ios;

   for (i=0; i < MAXFILE; i++) {
      if (used[i]) {
         ios = OK;
         ds_close(i,&ios);
      }
   }
   return 0;
}

int ds_path(path, ios)
   char path[];
   int *ios;
{
   if (!*ios) {
      npath = 0;
      if (slen(path) > 0) {
         npath = 1;
         ccopy(path, 0, SZLIN, paths[npath-1]);
      }
      npath = npath + 1;
      ccopy(" ", 0, SZLIN, paths[npath-1]);
   }
   return 0;
}

int ds_method(met, ios)
   char met[];
   int *ios;
{
   if (!*ios) {
      if (csame("binary",met)) {
         ccopy(met, 0, SZNAM, method);
      } else if (csame("text",met)) {
         ccopy(met, 0, SZNAM, method);
      } else {
         printf("error:  unknown data write method '%s' (ds_method)\n",
           met);
         *ios = ERR;
      }
   }
   return 0;
}

int ds_create(file, unit, ios)
   char file[];
   int *unit,*ios;
{
   int i;
   char filen[SZLIN],ext[SZNAM];
   bool found;

   if (!*ios) {
      found = FALSE;
      for (i=0; i < MAXFILE; i++) {
         if (!used[i]) {
            found = TRUE;
            *unit = i;
            break;
         }
      }
      if (!found) {
         printf("error:  too many file before creating '%s' ",file);
         printf("(ds_create)\n");
         *ios = ERR;
      } else {
         if (csame("binary", method)) {
            ccopy(BINEXT, 0, SZNAM, ext);
         } else if (csame("text", method)) {
            ccopy(TXTEXT, 0, SZNAM, ext);
         } else {
            printf("warning:  illegal method '%s' (ds_create)\n",
              method);
            ccopy(BINEXT, 0, SZNAM, ext);
         }
         ds_filex(paths[0], file, ext, SZLIN, filen);
         if ((fp[*unit] = fopen(filen,"w")) == NULL) {
            printf("error:  creating file '%s' (ds_create)\n", filen);
            *ios = ERR;
         } else {
            used[*unit] = TRUE;
            ccopy(filen, 0, SZLIN, files[*unit]);
            binary[*unit] = csame("binary", method);
            writen[*unit] = TRUE;
            level[*unit] = 0;
         }
      }
   }
   return 0;
}

int ds_open(file, unit, ios)
   char file[];
   int *unit,*ios;
{
   int i;
   char filen[SZLIN];
   bool found;

   if (!*ios) {
      found = FALSE;
      for (i=0; i < MAXFILE; i++) {
         if (!used[i]) {
            found = TRUE;
            *unit = i;
            break;
         }
      }
      if (!found) {
         printf("error:  too many files before opening '%s' ", file);
         printf("(ds_open)\n");
         *ios = ERR;
      } else {
         found = FALSE;
         for (i=0; i < npath; i++) {
            ds_filex(paths[i], file, BINEXT, SZLIN, filen);
            if ((fp[*unit] = fopen(filen,"r")) != NULL) {
               found = TRUE;
               binary[*unit] = TRUE;
               break;
            } else {
               ds_filex(paths[i], file, TXTEXT, SZLIN, filen);
               if ((fp[*unit] = fopen(filen, "r")) != NULL) {
                  found = TRUE;
                  binary[*unit] = FALSE;
                  break;
               }
            }
         }
         if (!found) {
            printf("error:  unable to find file '%s' (ds_open)\n",
              file);
            *ios = ERR;
         } else {
            used[*unit] = TRUE;
            writen[*unit] = FALSE;
            ccopy(filen, 0, SZLIN, files[*unit]);
            level[*unit] = 0;
         }
      }
   }
   return 0;
}

int ds_close(unit, ios)
   int unit,*ios;
{
   ds_check_unit(unit, ios);
   if (!*ios) {
      if (level[unit] > 0) {
         printf("warning:  record nesting wrong on closing file '%s' ",
           files[unit]);
         printf("(ds_close)\n");
      }
      if (writen[unit]) {
         ds_put_header(unit, "eof", "eof", 0, ios);
         if (*ios) {
            printf("error:  writing eof on file '%s' (ds_close)\n",
              files[unit]);
         }
      }
      fclose(fp[unit]);
      used[unit] = FALSE;
   }
   return 0;
}

int ds_index(name, no, lencomb, comb)
   char name[],comb[];
   int no,lencomb;
{
   char text[SZNAM];

   ccopy(name, 0, lencomb, comb);
   sprintf(text, "%d", no);
   ccat(text, lencomb, comb);
   return 0;
}

/* RECORD ROUTINES */

int ds_put_begin(unit, name, ios)
   int unit,*ios;
   char name[];
{
   if (!*ios) {
      ds_put_header(unit, "begin", name, 0, ios);
   }
   return 0;
}

int ds_put_end(unit, name, ios)
   int unit,*ios;
   char name[];
{
   if (!*ios) {
      ds_put_header(unit, "end", name, 0, ios);
   }
   return 0;
}

int ds_get_begin(unit, name, ios)
   int unit,*ios;
   char name[];
{
   if (!*ios) {
      ds_get_prelude(unit, "begin", name, 0, ios);
   }
   return 0;
}

int ds_get_end(unit, name, ios)
   int unit,*ios;
   char name[];
{
   if (!*ios) {
      ds_get_prelude(unit, "end", name, 0, ios);
   }
   return 0;
}

/* PUT ROUTINES */

int ds_put0i(unit, name, value, ios)
   int unit,*ios;
   char name[];
   int value;
{
   int values[1];

   if (!*ios) {
      values[0] = value;
      ds_put1i(unit, name, 1, values, ios);
   }
   return 0;
}

int ds_put1i(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   int value[];
{
   int i,count,indent;

   ds_check_unit(unit, ios);
   if (!*ios) {
      ds_put_header(unit, name, INT, nvalue, ios);
      if (!*ios) {
         if (binary[unit]) {
            ds_put1i_binary(unit, nvalue, value, ios);
         } else {
            indent = level[unit]+1;
            if (nvalue > 1) {
               counts[unit] = 100;
            }
            for (i=0; i < nvalue; i++) {
               if (counts[unit] > 60) {
                  if (fprintf(fp[unit], "\n") < 1) {
                     *ios = ERR;
                     break;
                  } else {
                     counts[unit] = 0;
                  }
                  ds_indent(unit, indent, ios);
                  if (*ios) {
                     break;
                  }
               }
               if ((count = fprintf(fp[unit], " %d ", value[i])) < 1) {
                  *ios = ERR;
                  break;
               } else {
                  counts[unit] = counts[unit] + count;
               }
            }
            fprintf(fp[unit], "\n");
         }
         if (*ios) {
            ds_write_error(unit, "ds_put1i", ios);
         }
      }
   }
   return 0;
}

int ds_put0r(unit, name, value, ios)
   int unit,*ios;
   char name[];
   float value;
{
   float values[1];

   if (!*ios) {
      values[0] = value;
      ds_put1r(unit, name, 1, values, ios);
   }
   return 0;
}


int ds_put1r(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   float value[];
{
   int i,count,indent;

   ds_check_unit(unit, ios);
   if (!*ios) {
      ds_put_header(unit, name, REAL, nvalue, ios);
      if (!*ios) {
         if (binary[unit]) {
            ds_put1r_binary(unit, nvalue, value, ios);
         } else {
            indent = level[unit]+1;
            if (nvalue > 1) {
               counts[unit] = 100;
            }
            for (i=0; i < nvalue; i++) {
               if (counts[unit] > 60) {
                  if (fprintf(fp[unit], "\n") < 1) {
                     *ios = ERR;
                     break;
                  } else {
                     counts[unit] = 0;
                  }
                  ds_indent(unit, indent, ios);
                  if (*ios) {
                     break;
                  }
               }
               if ((count = fprintf(fp[unit],"%13.6E ", value[i])) < 1){
                  *ios = ERR;
                  break;
               } else {
                  counts[unit] = counts[unit] + count;
               }
            }
            fprintf(fp[unit], "\n");
         }
         if (*ios) {
            ds_write_error(unit, "ds_put1r", ios);
         }
      }
   }
   return 0;
}

int ds_put0b(unit, name, value, ios)
   int unit,*ios;
   char name[];
   bool value;
{
   bool values[1];

   if (!*ios) {
      values[0] = value;
      ds_put1b(unit, name, 1, values, ios);
   }
   return 0;
}

int ds_put1b(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   bool value[];
{
   int i,count,indent;

   ds_check_unit(unit, ios);
   if (!*ios) {
      ds_put_header(unit, name, BOOL, nvalue, ios);
      if (!*ios) {
         if (binary[unit]) {
            ds_put1b_binary(unit, nvalue, value, ios);
         } else {
            indent = level[unit]+1;
            if (nvalue > 1) {
               counts[unit] = 100;
            }
            for (i=0; i < nvalue; i++) {
               if (counts[unit] > 60) {
                  if (fprintf(fp[unit], "\n") < 1) {
                     *ios = ERR;
                     break;
                  } else {
                     counts[unit] = 0;
                  }
                  ds_indent(unit, indent, ios);
                  if (*ios) {
                     break;
                  }
               }
               if ((count = fprintf(fp[unit], "%d ", value[i])) < 1) {
                  *ios = ERR;
                  break;
               } else {
                  counts[unit] = counts[unit] + count;
               }
            }
            fprintf(fp[unit], "\n");
         }
         if (*ios) {
            ds_write_error(unit, "ds_put1l", ios);
         }
      }
   }
   return 0;
}

int ds_put0s(unit, name, value, ios)
   int unit,*ios;
   char name[];
   char value[];
{
   int n;

   if (!*ios) {
      n = slen(value);
/*
      printf("ds_put0s:  n=%d\n", n);
      printf("..name=%s\n", name);
      printf("..value=%s\n", value);
*/
      ds_put1c(unit, name, n, value, ios);
   }
   return 0;
}

int ds_put0c(unit, name, value, ios)
   int unit,*ios;
   char name[],value;
{
   char values[1];

   if (!*ios) {
      values[0] = value;
      ds_put1c(unit, name, 1, values, ios);
   }
   return 0;
}

int ds_put1c(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   char value[];
{
   int i;

   ds_check_unit(unit, ios);
   if (!*ios) {
      ds_put_header(unit, name, CHAR, nvalue, ios);
      if (!*ios) {
         if (binary[unit]) {
            ds_put1c_binary(unit, nvalue, value, ios);
         } else {
            for (i=0; i < nvalue; i++) {
               if (fprintf(fp[unit], "%c", value[i]) < 1) {
                  *ios = ERR;
                  break;
               }
            }
            //fprintf(fp[unit], "\n", value);
            fprintf(fp[unit], "\n");			// just print line break BM
         }
         if (*ios) {
            ds_write_error(unit, "ds_put1c", ios);
         }
      }
   }
   return 0;
}


/* GET ROUTINES */

bool ds_get(unit, name, type, ios)
   int unit,*ios;
   char name[],type[];
{
   bool found;

   found = FALSE;
   if (!*ios) {
      ds_get_header(unit, ios);
      if (!*ios) {
         ccopy(nameb, 0, SZNAM, name);
         ccopy(typeb, 0, SZNAM, type);
         if (csame("end",nameb)) {
            found = FALSE;
         } else if (csame("eof",name)) {
            found = FALSE;
         } else {
            found = TRUE;
         }
      }
   }
   return(found);
}

int ds_get0i(unit, name, value, ios)
   int unit,*ios;
   char name[];
   int *value;
{
   int values[1];

   ds_get1i(unit, name, 1, values, ios);
   if (!*ios) {
      *value = values[0];
   }
   return 0;
}

int ds_get1i(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   int value[];
{
   int i;

   ds_get_prelude(unit, name, INT, nvalue, ios);
   if (!*ios) {
      if (binary[unit]) {
         ds_get1i_binary(unit, nvalue, value, ios);
      } else {
         for (i=0; i < nvalue; i++) {
            if (fscanf(fp[unit], " %d ", &(value[i])) < 1) {
               *ios = ERR;
               break;
            }
         }
      }
      if (*ios) {
         ds_read_error(unit, "ds_get1i", ios);
      }
   }
   return 0;
}

int ds_get0r(unit, name, value, ios)
   int unit,*ios;
   char name[];
   float *value;
{
   float values[1];

   ds_get1r(unit, name, 1, values, ios);
   if (!*ios) {
      *value = values[0];
   }
   return 0;
}

int ds_get1r(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   float value[];
{
   int i;

   ds_get_prelude(unit, name, REAL, nvalue, ios);
   if (!*ios) {
      if (binary[unit]) {
         ds_get1r_binary(unit, nvalue, value, ios);
      } else {
         for (i=0; i < nvalue; i++) {
           // Fix for double vs float BM
           #ifdef DOUBLE     
              if (fscanf(fp[unit], " %lf ", &(value[i])) < 1) {
                *ios = ERR;
                break;
              }
           #else
              if (fscanf(fp[unit], " %f ", &(value[i])) < 1) {
                 *ios = ERR;
                 break;
              }
           #endif
         }
      }
      if (*ios) {
         ds_read_error(unit, "ds_get1r", ios);
      }
   }
   return 0;
}

int ds_get0b(unit, name, value, ios)
   int unit,*ios;
   char name[];
   bool *value;
{
   bool values[1];

   ds_get1b(unit, name, 1, values, ios);
   if (!*ios) {
      *value = values[0];
   }
   return 0;
}

int ds_get1b(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   bool value[];
{
   int i;

   ds_get_prelude(unit, name, BOOL, nvalue, ios);
   if (!*ios) {
      if (binary[unit]) {
         ds_get1b_binary(unit, nvalue, value, ios);
      } else {
         for (i=0; i < nvalue; i++) {
            if (fscanf(fp[unit], " %d ", &(value[i])) < 1) {
               *ios = ERR;
               break;
            }
         }
      }
      if (*ios) {
         ds_read_error(unit, "ds_get1b", ios);
      }
   }
   return 0;
}

int ds_get0s(unit, name, value, ios)
   int unit,*ios;
   char name[];
   char value[];
{
   ds_get_header(unit, ios);
   if (!*ios) {
      if (!csame(name, nameb)) {
         printf("error:  expected name '%s' (ds_get0s)\n", name);
         ds_read_error(unit, "ds_get0s", ios);
      } else if (!csame(CHAR, typeb)) {
         printf("error:  expected type '%s' (ds_get0s)\n", CHAR);
         ds_read_error(unit, "ds_get0s", ios);
      } else {
         ds_get1c(unit, "#", nvalb, value, ios);
         if (!*ios) {
            value[nvalb] = '\0';
         }
      }
   }
   return 0;
}

int ds_get0c(unit, name, value, ios)
   int unit,*ios;
   char name[],*value;
{
   char values[1];

   ds_get1c(unit, name, 1, values, ios);
   if (!*ios) {
      *value = values[0];
   }
   return 0;
}

int ds_get1c(unit, name, nvalue, value, ios)
   int unit,nvalue,*ios;
   char name[];
   char value[];
{
   int i;

   ds_get_prelude(unit, name, CHAR, nvalue, ios);
   if (!*ios) {
      if (binary[unit]) {
         ds_get1c_binary(unit, &nvalue, value, ios);
      } else {
         for (i=0; i < nvalue; i++) {
            if (fscanf(fp[unit], "%c", &(value[i])) < 1) {
               ds_read_error(unit, "ds_get1c", ios);
               break;
            }
         }
      }
   }
   return 0;
}

/* ERROR ROUTINE */

int merror(text)
   char text[];
{
   printf("error:  %s\n", text);
   return 0;
}

/* INTERNAL ROUTINES */

int ds_filex(path, name, ext,lenfile, file)
   char path[],name[],ext[],file[];
   int lenfile;
{
   ccopy("", 0, lenfile, file);
   if (slen(path) > 0) {
      ccat(path, lenfile, file);
      ccat("/", lenfile, file);
   }
   ccat(name, lenfile, file);
   if (slen(ext) > 0) {
      ccat(".", lenfile, file);
      ccat(ext, lenfile, file);
   }
   slower(file);
   return 0;
}

int ds_check_unit(unit, ios)
   int unit,*ios;
{
   if (!*ios) {
      if ((unit < 0) || (unit >= MAXFILE)) {
         printf("error:  illegal unit %d (ds_check_unit)\n", unit);
         *ios = ERR;
      } else if (!used[unit]) {
         printf("error:  unit %d out of use (ds_check_unit)\n", unit);
         *ios = ERR;
      }
   }
   return 0;
}

int ds_put_header(unit, name, type, nvalue, ios)
   int unit,nvalue,*ios;
   char name[],type[];
{
   int indent,count,intv[1];

/*
   printf("ds_put_header:  name=%s,nvalue=%d,type=%s\n", name, nvalue,
     type);
*/
   ds_check_unit(unit, ios);
   if (*ios) {
   } else if (!writen[unit]) {
      printf("error:  unable to write on file '%s' (ds_put_header)\n",
        files[unit]);
      ds_write_error(unit, "ds_put_header", ios);
      *ios = ERR;
   } else {
      ccopy(name, 0, SZNAM, nameb);
      ccopy(type, 0, SZNAM, typeb);
      nvalb = nvalue;
      if (csame("begin",name)) {
         indent = level[unit];
         level[unit] = level[unit] + 1;
      } else if (csame("end",name)) {
         level[unit] = level[unit] - 1;
         indent = level[unit];
         if (level[unit] < 0) {
            printf("warning:  too many end-records on file '%s' ",
              files[unit]);
            printf("(ds_put_header)\n");
         }
      } else {
         indent = level[unit];
      }
      if (*ios) {
      } else if (binary[unit]) {
/*
         printf("..write name=%s\n", name);
*/
         ds_put0s_binary(unit, name, ios);
         if (*ios) {
         } else if (csame("begin",name)) {
         } else if (csame("end",name)) {
         } else if (csame("eof",name)) {
         } else {
            intv[0] = nvalue;
/*
            printf("..write nvalue=%d\n",nvalue);
*/
            ds_put1i_binary(unit, 1, intv, ios);
         }
         if (*ios) {
         } else if (csame("eof",name)) {
         } else {
/*
            printf("..write type=%s\n", type);
*/
            ds_put0s_binary(unit, type, ios);
         }
      } else {
         counts[unit] = 0;
         if (indent > 0) {
            ds_indent(unit, indent, ios);
         }
         if ((count = fprintf(fp[unit], "%s ", name)) < 1) {
            *ios = ERR;
         } else {
            counts[unit] = counts[unit] + count;
         }
         if (*ios) {
         } else if (csame("begin",name)) {
         } else if (csame("end",name)) {
         } else if (csame("eof",name)) {
         } else if ((count = fprintf(fp[unit], "%d ", nvalue)) < 1) {
            *ios = ERR;
         } else {
            counts[unit] = counts[unit] + count;
         }
         if (*ios) {
         } else if (csame("eof",name)) {
         } else if ((count = fprintf(fp[unit], "%s ", type)) < 1) {
            *ios = ERR;
         } else {
            counts[unit] = counts[unit] + count;
         }
         if (*ios) {
         } else if (nvalue < 1) {
            fprintf(fp[unit], "\n");
         }
      }
      if (*ios) {
         ds_write_error(unit, "ds_put_header", ios);
      }
   }
   return 0;
}

int ds_get_prelude(unit, name, type, nvalue, ios)
   int unit,nvalue,*ios;
   char name[],type[];
{
   char namew[SZNAM];

   if (!*ios) {
      if (csame("#",name)) {
         ds_check_unit(unit, ios);
         ccopy(nameb, 0, SZNAM, namew);
      } else {
         ds_get_header(unit, ios);
         ccopy(name, 0, SZNAM, namew);
      }
      ds_check_header(namew, type, nvalue, ios);
   }
   return 0;
}

int ds_get_header(unit, ios)
   int unit,*ios;
{
   int intv[1];

   if (!*ios) {
      ds_check_unit(unit, ios);
      if (*ios) {
      } else if (writen[unit]) {
         printf("error:  cannot read from file '%s' while writing ",
           files[unit]);
         printf("(ds_get_header)\n");
         *ios = ERR;
      } else {
         if (binary[unit]) {
/*
            printf("ds_get_header:  read nameb(binary)\n");
*/
            ds_get0s_binary(unit, SZNAM, nameb, ios);
/*
            printf("..got nameb='%s'\n", nameb);
*/
            if (*ios) {
            } else if (csame("begin",nameb)) {
               nvalb = 0;
            } else if (csame("end", nameb)) {
               nvalb = 0;
            } else if (csame("eof", nameb)) {
               nvalb = 0;
            } else {
/*
               printf("..read nvalb\n");
*/
               ds_get1i_binary(unit, 1, intv, ios);
               nvalb = intv[0];
/*
               printf("..got nvalb=%d\n", nvalb);
*/
            }
            if (*ios) {
            } else if (csame("eof",nameb)) {
               ccopy("eof", 0, SZNAM, typeb);
            } else {
/*
               printf("..read typeb\n");
*/
               ds_get0s_binary(unit, SZNAM, typeb, ios);
/*
               printf("..got typeb='%s'\n", typeb);
*/
            }
         } else {
            if (fscanf(fp[unit], "%s ", nameb) < 1) {
               printf("error:  reading name from file '%s' ",
                 files[unit]);
               printf("(ds_get_header)\n");
               *ios = ERR;
            }
            if (*ios) {
            } else if (csame("begin",nameb)) {
               nvalb = 0;
            } else if (csame("end", nameb)) {
               nvalb = 0;
            } else if (csame("eof", nameb)) {
               nvalb = 0;
            } else if (fscanf(fp[unit], "%d ", &nvalb) < 1) {
               printf("error:  rading size from file '%s' ",
                  files[unit]);
               printf("(ds_get_header)\n");
               *ios = ERR;
            }
            if (*ios) {
            } else if (csame("eof",nameb)) {
               ccopy("eof", 0, SZNAM, typeb);
            } else if (fscanf(fp[unit], "%s ", typeb) < 1) {
               printf("error:  reading type from file '%s' ",
                 files[unit]);
               printf("(ds_get_header)\n");
               *ios = ERR;
            }
         }
         if (*ios) {
         } else if (csame("begin", nameb)) {
            level[unit] = level[unit] + 1;
         } else if (csame("end", nameb)) {
            level[unit] = level[unit] - 1;
            if (level[unit] < 0) {
               printf("warning:  too many end-records in file '%s' ",
                 files[unit]);
               printf("ds_get_header)\n");
            }
         }
      }
   }
   return 0;
}

int ds_check_header(name, type, nvalue, ios)
   char name[],type[];
   int nvalue,*ios;
{
   if (!*ios) {
      if (!csame(name,nameb)) {
         printf("error:  expect name '%s' found '%s' (ds_check_header)",
           name, nameb);
         *ios = ERR;
      } else if (!csame(type,typeb)) {
         printf("error:  expect type '%s' found '%s' (ds_check_header)",
           type, typeb);
         *ios = ERR;
      } else if (nvalue != nvalb) {
         printf("error:  expect size %d found %d (ds_check_header)",
           nvalue, nvalb);
         *ios = ERR;
     }
   }
   return 0;
}

int ds_indent(unit, indent, ios)
   int unit,indent,*ios;
{
   int i,count;

   for (i=0; i < indent; i++) {
      if ((count = fprintf(fp[unit], "  ")) < 1) {
         *ios = ERR;
         break;
      } else {
         counts[unit] = counts[unit] + count;
      }
   }
   return 0;
}

int ds_read_error(unit, proc, ios)
   int unit,*ios;
   char proc[];
{
   printf("error:  reading '%s[%d]%s' from file '%s' (%s)\n",
     nameb, nvalb, typeb, files[unit], proc);
   *ios = ERR;
   return 0;
}

int ds_write_error(unit, proc, ios)
   int unit,*ios;
   char proc[];
{
   printf("error:  writing '%s[%d]%s' to file '%s' (%s)\n",
     nameb, nvalb, typeb, files[unit], proc);
   *ios = ERR;
   return 0;
}

/* SIMPLIFIED BINARY INPUT-OUTPUT ROUTINES */

int ds_put1c_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   char value[];
{
   int int_values[1],nc;

   nc = nvalue+1;
   int_values[0] = nc;
/*
   printf("ds_put1c_binary:  nvalue=%d\n", nvalue);
*/
   if (fwrite(int_values, sizeof(int), 1, fp[unit]) < 1) {
      *ios = ERR;
   } else {
      if (fwrite(value, sizeof(char), nc, fp[unit]) < nc) {
         *ios = ERR;
      }
   }
   return 0;
}

int ds_get1c_binary(unit, nvalue, value, ios)
   int unit,*nvalue,*ios;
   char value[];
{
   int nc,int_values[1];

   if (fread(int_values, sizeof(int), 1, fp[unit]) < 1) {
      *ios = ERR;
   } else {
      nc = int_values[0];
      *nvalue = nc;
      if (fread(value, sizeof(char), nc, fp[unit]) < nc) {
         *ios = ERR;
/*
      } else {
         printf("ds_get1c_binary:  nvalue=%d\n",*nvalue);
*/
      }
   }
   return 0;
}

int ds_put1i_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   int value[];
{

/*
   printf("ds_put1i_binary:  nvalue=%d,value=%d\n", nvalue,value[0]);
*/
   if (fwrite(value, sizeof(int), nvalue, fp[unit]) < nvalue) {
      *ios = ERR;
   }
   return 0;
}

int ds_get1i_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   int value[];
{

   if (fread(value, sizeof(int), nvalue, fp[unit]) < nvalue) {
      *ios = ERR;
/*
   } else {
      printf("ds_get1i_binary:  nvalue=%d,value=%d\n", nvalue,value[0]);
*/
   }
   return 0;
}

int ds_put1r_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   float value[];
{

   if (fwrite(value, sizeof(float), nvalue, fp[unit]) < nvalue) {
      *ios = ERR;
   }
   return 0;
}

int ds_get1r_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   float value[];
{

   if (fread(value, sizeof(float), nvalue, fp[unit]) < nvalue) {
      *ios = ERR;
   }
   return 0;
}

int ds_put1b_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   bool value[];
{

   if (fwrite(value, sizeof(bool), nvalue, fp[unit]) < nvalue) {
      *ios = ERR;
   }
   return 0;
}

int ds_get1b_binary(unit, nvalue, value, ios)
   int unit,nvalue,*ios;
   bool value[];
{

   if (fread(value, sizeof(bool), nvalue, fp[unit]) < nvalue) {
      *ios = ERR;
   }
   return 0;
}

int ds_put0s_binary(unit, value, ios)
   int unit,*ios;
   char value[];
{
   int int_values[1],nc;

   nc = slen(value)+1;
/*
   printf("ds_put0s_binary:  value=%s\n", value);
*/
   int_values[0] = nc;
   if (fwrite(int_values, sizeof(int), 1, fp[unit]) < 1) {
      *ios = ERR;
   } else {
      if (fwrite(value, sizeof(char), nc, fp[unit]) < nc) {
         *ios = ERR;
      }
   }
   return 0;
}

int ds_get0s_binary(unit, lenvalue, value, ios)
   int unit,lenvalue,*ios;
   char value[];
{
   int nc,int_values[1];

   if (fread(int_values, sizeof(int), 1, fp[unit]) < 1) {
      *ios = ERR;
   } else {
      nc = int_values[0];
      if (nc > lenvalue) {
         printf("error:  string overflow\n");
         printf("..nc(in)=%d,lenvalue=%d\n", nc, lenvalue);
         *ios = ERR;
      } else if (fread(value, sizeof(char), nc, fp[unit]) < nc) {
         *ios = ERR;
/*
      } else {
         printf("ds_get0s_binary:  nc=%d,value=%s\n", nc,value);
*/
      }
   }
   return 0;
}

