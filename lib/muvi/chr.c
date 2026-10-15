#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"

void ccopy(char a[], int offset, int n, char b[]);
int slen(char s[]);
int rem(int a,int b);
/* EXPLICIT DECLARATIONS ADDED BY BM */
int cindex(char a[], char c);
void ccat(char a[], int lenb, char b[]);
void cbuild(char delim, char a[], int lenb, char b[]);
void cscan(int lena, char a[], int offset, int lenb, char b[], int *newoff);
bool csame(char a[], char b[]);
void slower(char a[]);
bool cabrev(char a[], char b[], int *offset);
void srmblk(char a[]);
void cvector(int ndiv, float v0, float dv, int *ipower, float *scale, int *maxfield, int *maxdp);
void cfloat(float v, int field, int ndp, int lena, char a[]);


/* slen() - return string length (exclude EOS) */
int slen(s)
   char s[];
{
   int i;

   i = 0;
   while (s[i] != '\0') {
      i++;
   }
   return(i);
}

int cindex(a, c)
   char a[],c;
{
   int i,n;

   n = slen(a);
   for (i=0; i < n; i++) {
      if (a[i] == c) {
         return(i);
      }
   }
   return(n);
}


/* ccat - add string "a" onto the end of "b" */
// MIR ccat(a, lenb, b)
void ccat(a, lenb, b)
   char a[],b[];
   int lenb;
{
   int offset;

   offset = slen(b);
   ccopy(a, offset, lenb, b);
}

/* ccopy() - copy string a[] into string b[], starting at offset */
// MIR ccopy(a, offset, n, b)
void ccopy(a, offset, n, b)
   int offset,n;
   char a[],b[];
{
   int i,m;

   m = slen(a);
   if (m+offset < n) {
      for (i=0; i<m+1; i++) {
         b[i+offset] = a[i];
      }
   }
}

/* cbuild - build a composite string (assume starts delim) */
//MIR cbuild(delim, a, lenb, b)
void cbuild(delim, a, lenb, b)
   char delim,a[],b[];
   int lenb;
{
   char delims[2];

   delims[0] = delim;
   delims[1] = '\0';
   ccat(a, lenb, b);
   ccat(delims, lenb, b);
}

/* cscan - scan a composite (starts delim) */
// MIR cscan(lena, a, offset, lenb, b, newoff)
void cscan(lena, a, offset, lenb, b, newoff)
   char a[],b[];
   int lena,offset,lenb,*newoff;
{
   char delim;
   int i,k,n;

   delim = a[0];
   if (offset == 0) {
      offset = 1;
   }
   n = 0;
   for (k = offset; k < lena; k++) {
      i = k - offset;
      n = i;
      if (a[k] == '\0') {
         break;
      } else if (a[k] == delim) {
         break;
      }
   }
   *newoff = offset+n+1;
   if (n > (lenb-1)) {
      n = lenb-1;
   }
   for (i=0; i < n; i++) {
      b[i] = a[offset+i];
   }
   b[n] = '\0';
}

/* csame() - compare strings, return whether same (ignore case) */
bool csame(a, b)
   char a[],b[];
{
   char ca,cb;
   int i,let;

   i = 0;
   while (TRUE) {
      ca = a[i];
      let = ca - 'A';
      if (let < 0) {
      } else if (let < 26) {
         ca = 'a' + let;
      }
      cb = b[i++];
      let = cb - 'A';
      if (let < 0) {
      } else if (let < 26) {
         cb = 'a' + let;
      }
      if (ca != cb) {
         break;
      } else if (ca == '\0') {
         return(TRUE);
      }
   }
   return(FALSE);
}

/* clower - convert a string in situe to lower case */
// MIR slower(a)
void slower(a)
   char a[];
{
   int n,i,let;

   n = slen(a);
   for (i=0; i < n; i++) {
      let = a[i] - 'A';
      if (let < 0) {
      } else if (let < 26) {
         a[i] = 'a' + let;
      }
   }
}

/* cabrev - see if string "b" starts "a/" (ignode case)*/
bool cabrev(a, b, offset)
   char a[],b[];
   int *offset;
{
   int i,lena,lenb,count,let;
   char ca,cb;

   lena = slen(a);
   lenb = slen(b);
   if (lena+1 > lenb) {
      return(FALSE);
   } else {
      count = 0;
      for (i=0; i < lena; i++) {
         ca = a[i];
         let = ca - 'A';
         if (let < 0) {
         } else if (let < 26) {
            ca = 'a' + let;
         }
         cb = b[i];
         let = cb - 'A';
         if (let < 0) {
         } else if (let < 26) {
            cb = 'a' + let;
         }
         if (ca != cb) {
            break;
         } else {
            count = count + 1;
         }
      }
      if (count < lena) {
         return(FALSE);
      } else if (b[lena] != '/') {
         return(FALSE);
      } else {
         *offset = lena+1;
         return(TRUE);
      }
   }
}

/* srmblk - remove leading blanks */
// MIR srmblk(a)
void srmblk(a)
   char a[];
{
   int i,n,k,j;

   n = slen(a);
   k = 0;
   for (i=0; i < n; i++) {
      if (a[i] != ' ') {
         k = i+1;
         break;
      }
   }
   if (k > 0) {
      j = 0;
      for (i=k-1; i < (n+1); i++) {
         a[j] = a[i];
         j = j + 1;
      }
   }
}

// MIR cvector(ndiv, v0, dv, ipower, scale, maxfield, maxdp)
void cvector(ndiv, v0, dv, ipower, scale, maxfield, maxdp)
   int ndiv,*ipower,*maxfield,*maxdp;
   float v0,dv,*scale;
{
   char str[SZLIN];
   int i,maxip,n,j,k,ndec,maxdec,maxfd;
   float v,s,p,vmax;

/*
   printf("cvector:  ndiv=%d,v0=%f,dv=%f\n", ndiv, v0, dv);
*/
/* find max absolute value */
   vmax = -BIG;
   for (i=0; i < (ndiv+1); i++) {
      v = (float)i*dv + v0;
/*
      printf("..i=%d,v=%f\n", i, v);
*/
      if (v < 0.0) {
         v = -v;
      }
      if (v > vmax) {
         vmax = v;
      }
   }
/* find corresponding exponent */
   if (vmax <= 0.0) {
      maxip = 0;
   } else {
/*
      glog10(vmax, &p);
*/
      p = log10(vmax);
      maxip = (int)p;
   }
/*
   printf("..maxip=%d\n", maxip);
*/
/* choose to ignore small exponents */
   if ((maxip > -2) && (maxip < 3)) {
      *ipower = 0;
      s = 1.0;
   } else {
      *ipower = maxip;
      v = (float)maxip;
/*
      gpow(10.0, v, &s);
*/
      s = pow(10.0, v);
   }
   *scale = s;
/*
   printf("..ipower=%d,scale=%f\n", *ipower, *scale);
*/
/* format numbers and determine ndp and so on */
   maxdec = 1;
   maxfd = 0;
   for (i=0; i < (ndiv+1); i++) {
      v = ((float)i*dv + v0)/s;
      sprintf(str, "%f", v);
      n = slen(str);
      j = cindex(str, '.');
      ndec = (n-1)-j;
      for (k=(n-1); k > (j+1); k--) {
         if (str[k] != '0') {
            break;
         }
         ndec = ndec - 1;
      }
      if (ndec > maxdec) {
         maxdec = ndec;
      }
      if (j+1 > maxfd) {
         maxfd = j+1;
      }
   }
   *maxfield = maxfd+maxdec;
   *maxdp = maxdec;
/*
   printf("cvector:  maxfd=%d,maxdec=%d,maxfield=%d\n",maxfd, maxdec,
     *maxfield);
*/
}

//MIR cfloat(v, field, ndp, lena, a)
void cfloat(v, field, ndp, lena, a)
   float v;
   int field,ndp,lena;
   char a[];
{
   char f[SZLIN];

   sprintf(f, "%%%d.%df", field, ndp);
   sprintf(a, f, v);
}

int rem(a,b)
   int a,b;
{
   return ( a - ( (a/b)*b ) );
}

