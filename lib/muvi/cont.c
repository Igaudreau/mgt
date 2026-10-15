/* cont.c:  pauls contour routine, isolated */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#define bit char
#define ACROSS 1
#define UP 2
#define ACROSS_DOWN 3
#define UP_LEFT 4

/* private copies of called arguments */

static int contour_nx,contour_ny;
static float *contour_xa,*contour_ya,*contour_qa;
static int *contour_sa;
static char *contour_lab;

/* internal scratch data */

static bit *bitmap,*wall1,*wall2;

/* contour label control */

static float contour_x,contour_y,contour_dl,contour_dl0,contour_left;
static bool contour_first,contour_shown;

/* Declarations BM */
void setbitmap();
void locate(int line, float level);


/* External Declarations BM */
extern void vis_bp(void);
extern void vis_ap(float h, float v);
extern void vis_text(float h, float v, char str[], char just[]);
extern void vis_ep(void);
extern int slen(char s[]);
extern int rem(int a,int b);

/* draw a single contour level */
//MIR contour_map(nx, ny, xa, ya, sa, qa, qlevel, qlab, ios)
void contour_map(int nx, int ny, float xa[], float ya[], int sa[], float qa[], float qlevel, char qlab[], int *ios);

void contour_map(nx, ny, xa, ya, sa, qa, qlevel, qlab, ios)
   int nx,ny,*ios;
   float xa[],ya[],qa[],qlevel;
   int sa[];
   char qlab[];
{
int k,tot,nbit;
float dx,dy;

/* check for memory overflow */
   nbit = 2*nx*ny;
   bitmap = (bit*)calloc(nbit,sizeof(bit));
   wall1 = (bit*)calloc(nbit,sizeof(bit));
   wall2 = (bit*)calloc(nbit,sizeof(bit));
   if (bitmap == NULL) {
      printf("Error:  no memory for contours (bitmap)\n");
      *ios = ERR;
   } else if (wall1 == NULL) {
      printf("Error:  no memory for contours (wall1)\n");
      *ios = ERR;
   } else if (wall2 == NULL) {
      printf("Error:  no memory for contours (wall2)\n");
      *ios = ERR;
   } else {
/*    take addresses */
      contour_nx = nx;
      contour_ny = ny;
      contour_xa = xa;
      contour_ya = ya;
      contour_qa = qa;
      contour_sa = sa;
      contour_lab = qlab;

/*    suitable contour label spacing */
      dx = xa[nx-1] - xa[0];
      dy = ya[ny-1] - ya[0];
      contour_dl = 0.15*sqrt(dx*dx + dy*dy);
      contour_dl0 = 0.5*contour_dl;
/*    lets go */
      setbitmap();
      tot= ( (nx-1)*ny ) +( (ny-1)*nx );
      for (k=1 ; k<=tot ;k++ ) {
         if ( (wall1[k]==TRUE) && (bitmap[k]==TRUE) ) {
            locate(k,qlevel);
         }
      }
      for (k=1 ; k<=tot ; k++ ) {
         if ( bitmap[k]==TRUE ) {
            locate(k,qlevel);
         }
      }
   }
/* return any non-NULL memory */
   if (bitmap != NULL) {
      cfree(bitmap, nbit, sizeof(char));
   }
   if (wall1 != NULL) {
      cfree(wall1, nbit, sizeof(char));
   }
   if (wall2 != NULL) {
      cfree(wall2, nbit, sizeof(char));
   }
}


// MIR contour_begin() {
void contour_begin() {
   contour_first = TRUE;
   vis_bp();
}

// MIR contour_append(ax, ay)
void contour_append(ax, ay)
   float ax,ay;
{
   float x,y,dx,dy,dl;
   int ix,iy;

   ix = ax;
   if (ix < 0) {
      x = contour_xa[0];
   } else if (ix >= contour_nx-1) {
      x = contour_xa[contour_nx-1];
   } else {
      dx = ax - (float)ix;
      x = (contour_xa[ix+1]-contour_xa[ix])*dx + contour_xa[ix];
   }
   iy = ay;
   if (iy < 0) {
      y = contour_ya[0];
   } else if (iy >= contour_ny-1) {
      y = contour_ya[contour_ny-1];
   } else {
      dy = ay - (float)iy;
      y = (contour_ya[iy+1]-contour_ya[iy])*dy + contour_ya[iy];
   }
   vis_ap(x, y);
   if (contour_first) {
      contour_x = x;
      contour_y = y;
      contour_first = FALSE;
      contour_shown = FALSE;
      contour_left = contour_dl;
   } else {
      dx = x - contour_x;
      dy = y - contour_y;
      contour_x = x;
      contour_y = y;
      dl = sqrt(dx*dx + dy*dy);
      contour_left = contour_left - dl;
      if (contour_left < 0.0) {
         if (slen(contour_lab) > 0) {
            vis_text(x, y, contour_lab, "mm");
         }
         contour_left = contour_dl;
         contour_shown = TRUE;
      }
   }
}

// MIR contour_end() {
void contour_end() {
   vis_ep();
   if (!contour_first) {
      if (!contour_shown) {
         if (contour_left < 0.8*contour_dl) {
            if (slen(contour_lab) > 0) {
               vis_text(contour_x, contour_y, contour_lab, "mm");
            }
         }
      }
   }
}

/* BELOW HERE IS A 3RD GENERATION RATS NETS, ENTER AT YOUR PERIL */

int getdirection(line,lines_across)
   int line,lines_across;
{
   if (line <= lines_across ){
      return ACROSS ;
   }  else {return UP;}
}
bool intersect(a,b,c)
   float a,b,c;
{
   if ( (a-c)*(b-c) <= 0 ) {
      return TRUE;
   } else {
      return FALSE;
   }
}


float lin_interp(y1,y2,level)
   float y1,y2,level;
{
   if (y1==y2) {
      return 0;
   } else {
      return ( (level-y1)/(y2-y1) );
   }
}
// MIR getij(line,lines_across,i,j)
void getij(line,lines_across,i,j)
   int line,lines_across,*i,*j;
{
int temp;
   temp=line-1;
  if ( line<=lines_across ) {
      *j=temp/(contour_nx-1);
      *i=rem(temp, (contour_nx-1) );
   } else {
      temp=temp-lines_across;
      *i=temp/(contour_ny-1);
      *j=rem(temp, (contour_ny-1) );
   }
}

//MIR linetf(line,in)
void linetf(line,in)
   int line,in;
{
   if (bitmap[line]==TRUE) {
      if ( in==1 ) {
         bitmap[line]=FALSE;
      } else {
        if ( wall1[line]==TRUE ) {
           bitmap[line]=FALSE;
        } else {
           wall1[line]=TRUE;
           if ( in==3 ) { wall2[line]=TRUE ; }
        }
      }
   }
}

// MIR setbitmap()
void setbitmap()
{
int left,down,ix,iy,k,tot,lines_across;
   lines_across=(contour_nx-1)*contour_ny;
   tot=lines_across + ( (contour_ny-1)*contour_nx );
   for ( k=1 ; k<=tot ;k++ ) {
      bitmap[k]=TRUE;
      wall1[k]=FALSE;
      wall2[k]=FALSE;
   }
   for ( k=1 ;k<contour_nx; k++ ) {
      wall1[k]=TRUE;
      wall2[k]=FALSE;
      wall1[lines_across- (k-1)]=TRUE;
      wall2[lines_across- (k-1)]=TRUE;
   }
   for ( k=1 ; k<contour_ny ;k++ ) {
      wall1[lines_across+k]=TRUE;
      wall2[lines_across+k]=FALSE;
      wall1[tot-(k-1)]=TRUE;
      wall2[tot- (k-1)]=TRUE;
   }
   for (iy=0; iy<contour_ny; iy++ ) {
      for (ix=0; ix<contour_nx; ix++ ) {
         k=iy*contour_nx + ix;
         if ( contour_sa[k]!=FLUID ) {
            left= (contour_nx-1)*iy + ix;
            down= (lines_across + iy) + ( (contour_ny-1)*ix );
            if (ix>0) {
               linetf( left,1 );
               wall1[ left ]=FALSE; 
            }
            if (ix<(contour_nx-1)) {
               linetf( left+1,1);
               wall1[left+1]=FALSE;
            }
            if (iy>0) {
               linetf(down,1);
               wall1[down]=FALSE;
            }
            if (iy<(contour_ny-1)) {
               linetf(down+1,1);
               wall1[down+1]=FALSE;
            }
            if ( (ix>0)&&(ix<(contour_nx-1))&&(iy>0)&&(iy<(contour_ny-1)) ) {
               linetf( (down-(contour_ny-2)) , 3);
               linetf( (down-(contour_ny-1)) , 3);
               linetf( (left-(contour_nx-2)) , 3);
               linetf( (left-(contour_nx-1)) , 3);
               linetf( (down+(contour_ny-1)) , 2);
               linetf( (down+contour_ny)   , 2);
               linetf( (left+contour_nx)   , 2);
               linetf( (left+(contour_nx-1)) , 2);
            }
         }
      }
   }
}

//MIR  getab(i,j,direction,a,b)
void getab(i,j,direction,a,b)
   int i,j,direction;
   float *a,*b;
{
   int k;
   k=(j*contour_nx)+i;
   *a=contour_qa[k];
   if ( (direction==ACROSS) || (direction==ACROSS_DOWN) ) {
      k=k+1;
   } else {
      k=k+contour_nx;
   }
   *b=contour_qa[k];
}
// MIR getcd(i,j,direction,c,d)
void getcd(i,j,direction,c,d)
   int i,j,direction;
   float *c,*d;
{
int k;
   if ( (direction==ACROSS) || (direction==ACROSS_DOWN) ) {
      k=(j+1)*contour_nx + i;
      if ( direction==ACROSS_DOWN ) { k=k-(2*contour_nx); }
      *d=contour_qa[k];
      k=k+1;
   } else {
      k=(j*contour_nx) +(i+1);
      if ( direction==UP_LEFT ) { k=k-2; }
      *d=contour_qa[k];
      k=k+contour_nx;
   }
   *c=contour_qa[k];
}

// MIR verify(a,b,c,d,level,line,i,j,direction_in,
void verify(a,b,c,d,level,line,i,j,direction_in,
       start,lineout,direction_out,located)
   int line,i,j,direction_in,*lineout,*direction_out;
   bool *located;
   float a,b,c,d,level,*start;
{
int intersectcount,lineintersected,line_to_blank,store,
    count,flag,flag2,k;
int blank[5];
bool path;
float square[6];

   *located=FALSE;
   path=FALSE;
   square[0]=d;
   square[1]=a;
   square[2]=b;
   square[3]=c;
   square[4]=d;
   square[5]=a;
   
   k=(contour_nx-1)*contour_ny;
   if ( (direction_in==ACROSS) || (direction_in==ACROSS_DOWN)) {
      store=k + ( i*(contour_ny-1) ) + (j+1);
      store=store - ( direction_in/2 );
      blank[4]=store;
      blank[3]=line + ( (contour_nx-1)*(2-direction_in) );
      blank[2]=store + (contour_ny-1);
   } else {
      store= j*(contour_nx-1) + (i+1);
      store=store-(direction_in/3);
      blank[4]=store;
      blank[3]=line + ( (contour_ny-1)*(3-direction_in) );
      blank[2]=store+ (contour_nx-1);
   }
   count=0;
   for (k=1; k<5; k++) {
      if ( square[k]==level ) {
         count=count+1;
         flag2=k;
      } else {
         flag=k;
      } 
   }
   if ( count==4 ) {
      *located=FALSE;
   } else if (count==3) {
     *located=TRUE;
     if ( (flag==3) || (flag==4) ) {
        *lineout=blank[3];
        *direction_out=2;
        bitmap[ blank[2] ]=FALSE;
        bitmap[ blank[4] ]=FALSE;
        if ( flag==3 ) { path=TRUE; }
      } else {
        if (flag==2) {
           *lineout=blank[2];
           *direction_out=1;
        } else {
           *lineout=blank[4];
           *direction_out=3;
        }
        bitmap[ blank[flag+1] ]=FALSE;
        bitmap[ blank[flag+2] ]=FALSE;
      }
   } else if (count==2) {
      if ( flag2==2 ) {
         *lineout=blank[2];
         *direction_out=1;
         bitmap[ blank[3] ]=FALSE;
         bitmap[ blank[4] ]=FALSE;
      } else {
         *direction_out=2;
         *lineout=blank[3];
         bitmap[ blank[2] ]=FALSE;
         bitmap[ blank[4] ]=FALSE;
      }
      *located=TRUE;
   } else {
      intersectcount=0;
      lineintersected=0;
      for (k=2; k<5; k++) {
         line_to_blank= blank[k];
         if ( bitmap[line_to_blank]==TRUE ) {
            if ( intersect( square[k],square[k+1],level)==TRUE ){
               intersectcount=intersectcount+1;
               lineintersected=k;
               *located=TRUE;
            } else {
               bitmap[ line_to_blank ]=FALSE;
            }
         }
      }
      if ( *located==TRUE) {
         *lineout= blank[ lineintersected ];
         *direction_out= lineintersected-1;
      }
      if ( count==1 ) {
        if ( ((flag2==1) || (flag2==2)) && (intersectcount==3)) {
           *located=FALSE;
        }
        if ( (flag2==1) && (intersectcount>1) ) {
           *direction_out=1;
        }
        if ( (flag2==4) ) { 
           *direction_out=2;
           *lineout=blank[3];
        }
/*
        if ( (flag2==1) && (lineintersected==4) ) {
           *located=FALSE; 
        }
        if ( (flag2==2) && (lineintersected==2) ) {
           *located=FALSE;
        }
*/
      }
   }
   if ( path==FALSE ) {
   *start=lin_interp(a,b,level);
   } else { *start=1; }
}


void locate(line,level)
   float level;
   int line;
{
float a,b,c,d,xpos,ypos,start;
int   i,j,currentline,lineout,startline,old_direction,
      direction,direction_out,lines_across;
bool located;
   lines_across=(contour_nx-1)*contour_ny;
   direction=getdirection(line,lines_across);
   getij(line,lines_across,&i,&j);
   getab(i,j,direction,&a,&b);
   if ( wall2[line]==TRUE ) {
      direction=direction+2;
   }
   if ( intersect(a,b,level)==TRUE ) {
      getcd(i,j,direction,&c,&d);
      verify(a,b,c,d,level,line,i,j,direction,&start,
            &lineout,&direction_out,&located);
      startline=line;
      if ( located==TRUE ) {
         if ( (direction==ACROSS) || (direction==ACROSS_DOWN) ) {
            xpos=(float)i + start;
            ypos=(float)j;
         } else {
            xpos=(float)i;
            ypos=(float)j + start;
         }
         contour_begin();
         contour_append(xpos,ypos);
         while ( (bitmap[line]==TRUE) &&  (located==TRUE ) ) {
            currentline=lineout;
            getij(currentline,lines_across,&i,&j);
            old_direction=direction;
            direction=getdirection(currentline,lines_across);
            if ( (old_direction==3) && (direction==1) ) {
               direction=3;
            } else if ( (old_direction==4)&&(direction==2)) {
               direction=4;
            } else if ( (direction_out)==3 ) {
               direction=direction+2;
            }
            getab(i,j,direction,&a,&b);
            if ( wall1[currentline]==FALSE ) {
               getcd(i,j,direction,&c,&d);
               verify(a,b,c,d,level,currentline,i,j,direction,
                       &start,&lineout,&direction_out,&located);
               bitmap[ currentline ]=FALSE;
            } else {
              start=lin_interp(a,b,level);
              located=FALSE;
              bitmap[ currentline ]=FALSE;
            }
/*
*/
            if ( (direction==ACROSS) || (direction==ACROSS_DOWN) ) {
               xpos=(float)i + start;
               ypos=(float)j;
            } else {
               xpos=(float)i;
               ypos=(float)j + start;
            }
            contour_append(xpos,ypos);
/*

*/
         }   
         contour_end();

      }
      bitmap[ startline ]=FALSE;
   } else {
      bitmap[ line ]=FALSE;
   }
}
