/* dr3solid.c:  a collection of solid primitives */
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"

#undef DEBUG

#ifdef STL
  void dr3_solid(char shape[], float par[], int npar);
#else
  void dr3_solid(char shape[], float par[]);
#endif
void dr3_solid_cylinder(float radius, float height);
void dr3_solid_rectangle(float width, float height);
void dr3_solid_cone(float radius, float height);
void dr3_solid_block(float width, float height, float depth);
void dr3_solid_sphere(float radius);
void dr3_solid_rect(float width, float height);

/* EXPLICIT DECLARATIONS ADDED BY BM */
extern bool csame(char a[], char b[]);
extern int dr3_begin_polygon();
extern void dr3_polygon_vertex(float x,float y,float z);
extern void dr3_polygon_normal(float nx,float ny,float nz);
extern void dr3_end_polygon();




#ifdef STL
  void dr3_solid(char shape[], float par[], int npar)
#else
  void dr3_solid(char shape[], float par[])
#endif
{
   if (csame("cylinder",shape)) {
      dr3_solid_cylinder(par[0], par[1]);
   } else if (csame("cone",shape)) {
      dr3_solid_cone(par[0], par[1]);
   } else if (csame("block",shape)) {
      dr3_solid_block(par[0], par[1], par[2]);
   } else if (csame("sphere",shape)) {
      dr3_solid_sphere(par[0]);
   } else if (csame("rect",shape)) {
      dr3_solid_rect(par[0],par[1]);
   } else {
   }
}


void dr3_solid_cylinder(float radius, float height)
{
   int nr,ir,nz,iz,ntheta,itheta;
   float dr,rl,rr, dz,zl,zr, dtheta,thetar,
     sintl,sintr,costl,costr, xl,xr,yl,yr;

   ntheta = 24;
   dtheta = 2.0*PI/(float)ntheta;
   nz = (int)((float)ntheta*height/(2.0*PI*radius));
   if (nz < 1) {
      nz = 1;
   }
   nr = ntheta*(int)(0.5/PI);
   if (nr < 1) {
      nr = 1;
   }
   dr = radius/(float)nr;
   dz = height/(float)nz;
   thetar = 0.0;
   sintr = sin(thetar);
   costr = cos(thetar);
   for (itheta=0; itheta < ntheta; itheta++) {
      sintl = sintr;
      costl = costr;
      thetar = thetar + dtheta;
      sintr = sin(thetar);
      costr = cos(thetar);
      rr = 0.0;
      for (ir=0; ir < nr; ir++) {
         rl = rr;
         rr = rr + dr;
         dr3_begin_polygon();
         dr3_polygon_normal(0.0, 0.0, -1.0);
         if (ir > 0) {
            dr3_polygon_vertex(rl*costl, rl*sintl, 0.0);
         }
         dr3_polygon_vertex(rl*costr, rl*sintr, 0.0);
         dr3_polygon_vertex(rr*costr, rr*sintr, 0.0);
         dr3_polygon_vertex(rr*costl, rr*sintl, 0.0);
         dr3_end_polygon();
         dr3_begin_polygon();
         dr3_polygon_normal(0.0, 0.0, 1.0);
         if (ir > 0) {
            dr3_polygon_vertex(rl*costl, rl*sintl, height);
         }
         dr3_polygon_vertex(rl*costr, rl*sintr, height);
         dr3_polygon_vertex(rr*costr, rr*sintr, height);
         dr3_polygon_vertex(rr*costl, rr*sintl, height);
         dr3_end_polygon();
      }
      xl = radius*costl;
      xr = radius*costr;
      yl = radius*sintl;
      yr = radius*sintr;
      zr = 0.0;
      for (iz=0; iz < nz; iz++) {
         zl = zr;
         zr = zr + dz;
         dr3_begin_polygon();
         dr3_polygon_normal(costl, sintl, 0.0);
         dr3_polygon_vertex(xl, yl, zl);
         dr3_polygon_normal(costr, sintr, 0.0);
         dr3_polygon_vertex(xr, yr, zl);
         dr3_polygon_normal(costr, sintr, 0.0);
         dr3_polygon_vertex(xr, yr, zr);
         dr3_polygon_normal(costl, sintl, 0.0);
         dr3_polygon_vertex(xl, yl, zr);
         dr3_end_polygon();
      }
   }
}

void dr3_solid_rectangle(float width, float height)
{
   dr3_begin_polygon();
   dr3_polygon_normal(0.0, 0.0, 1.0);
   dr3_polygon_vertex(0.0, 0.0, 0.0);
   dr3_polygon_vertex(width, 0.0, 0.0);
   dr3_polygon_vertex(width, height, 0.0);
   dr3_polygon_vertex(0.0, height, 0.0);
   dr3_end_polygon();
}

void dr3_solid_cone(float radius, float height)
{
   int itheta,ntheta;
   float dtheta,thetar, sintl,sintr,costl,costr, sintm,costm,
     sinp,cosp,side;

/* azimuthal angles */
   ntheta = 24;
   dtheta = 2.0*PI/(float)ntheta;
/* conical angle */
   side = sqrt(height*height + radius*radius);
   cosp = height/side;
   sinp = radius/side;
/* work around azimuthal segments */
   thetar = 0.0;
   sintr = sin(thetar);
   costr = cos(thetar);
   for (itheta=0; itheta < ntheta; itheta++) {
      sintl = sintr;
      costl = costr;
      thetar = thetar + dtheta;
      sintr = sin(thetar);
      costr = cos(thetar);
      sintm = 0.5*(costl+costr);
      costm = 0.5*(sintl+sintr);
/*    base segment */
      dr3_begin_polygon();
      dr3_polygon_normal(0.0, 0.0, -1.0);
      dr3_polygon_vertex(0.0, 0.0, 0.0);
      dr3_polygon_vertex(radius*costl, radius*sintl, 0.0);
      dr3_polygon_vertex(radius*costr, radius*sintr, 0.0);
      dr3_end_polygon();
/*    top segment */
      dr3_begin_polygon();
      dr3_polygon_normal(costm*cosp, sintm*cosp, sinp);
      dr3_polygon_vertex(0.0, 0.0, height);
      dr3_polygon_normal(costr*cosp, sintr*cosp, sinp);
      dr3_polygon_vertex(radius*costr, radius*sintr, 0.0);
      dr3_polygon_normal(costl*cosp, sintl*cosp, sinp);
      dr3_polygon_vertex(radius*costl, radius*sintl, 0.0);
      dr3_end_polygon();
   }
}

void dr3_solid_block(float width, float height, float depth)
{
/* front */
   dr3_begin_polygon();
   dr3_polygon_normal(0.0, 0.0, -1.0);
   dr3_polygon_vertex(0.0, 0.0, 0.0);
   dr3_polygon_vertex(width, 0.0, 0.0);
   dr3_polygon_vertex(width, height, 0.0);
   dr3_polygon_vertex(0.0, height, 0.0);
   dr3_end_polygon();
/* back */
   dr3_begin_polygon();
   dr3_polygon_normal(0.0, 0.0, 1.0);
   dr3_polygon_vertex(0.0, 0.0, depth);
   dr3_polygon_vertex(width, 0.0, depth);
   dr3_polygon_vertex(width, height, depth);
   dr3_polygon_vertex(0.0, height, depth);
   dr3_end_polygon();
/* left side */
   dr3_begin_polygon();
   dr3_polygon_normal(-1.0, 0.0, 0.0);
   dr3_polygon_vertex(0.0, 0.0, 0.0);
   dr3_polygon_vertex(0.0, height, 0.0);
   dr3_polygon_vertex(0.0, height, depth);
   dr3_polygon_vertex(0.0, 0.0, depth);
   dr3_end_polygon();
/* right side */
   dr3_begin_polygon();
   dr3_polygon_normal(1.0, 0.0, 0.0);
   dr3_polygon_vertex(width, 0.0, 0.0);
   dr3_polygon_vertex(width, height, 0.0);
   dr3_polygon_vertex(width, height, depth);
   dr3_polygon_vertex(width, 0.0, depth);
   dr3_end_polygon();
/* bottom side */
   dr3_begin_polygon();
   dr3_polygon_normal(0.0, -1.0, 0.0);
   dr3_polygon_vertex(0.0, 0.0, 0.0);
   dr3_polygon_vertex(width, 0.0, 0.0);
   dr3_polygon_vertex(width, 0.0, depth);
   dr3_polygon_vertex(0.0, 0.0, depth);
   dr3_end_polygon();
/* top side */
   dr3_begin_polygon();
   dr3_polygon_normal(0.0, 1.0, 0.0);
   dr3_polygon_vertex(0.0, height, 0.0);
   dr3_polygon_vertex(width, height, 0.0);
   dr3_polygon_vertex(width, height, depth);
   dr3_polygon_vertex(0.0, height, depth);
   dr3_end_polygon();
}

void dr3_solid_sphere(float radius)
{
   int ntheta,itheta, nphi,iphi;
   float dtheta,thetar, sintl,sintr,costl,costr, x,y,z,divr,
     dphi,phir,sinpl,sinpr,cospl,cospr;

   divr = 1.0/radius;
   ntheta = 24;
   dtheta = 2.0*PI/(float)ntheta;
   nphi = 12;
   dphi = PI/(float)nphi;
   thetar = 0.0;
   sintr = sin(thetar);
   costr = cos(thetar);
   for (itheta=0; itheta < ntheta; itheta++) {
      sintl = sintr;
      costl = costr;
      thetar = thetar + dtheta;
      sintr = sin(thetar);
      costr = cos(thetar);
      phir = -PI/2.0;
      sinpr = sin(phir);
      cospr = cos(phir);
      for (iphi=0; iphi < nphi; iphi++) {
         sinpl = sinpr;
         cospl = cospr;
         phir = phir + dphi;
         sinpr = sin(phir);
         cospr = cos(phir);
         dr3_begin_polygon();
         if (iphi == 0) {
            x = 0.0;
            y = 0.0;
            z = -radius;
            dr3_polygon_normal(x*divr, y*divr, z*divr);
            dr3_polygon_vertex(x, y, z);
         } else {
            x = radius*cospl*costl;
            y = radius*cospl*sintl;
            z = radius*sinpl;
            dr3_polygon_normal(x*divr, y*divr, z*divr);
            dr3_polygon_vertex(x, y, z);
            x = radius*cospl*costr;
            y = radius*cospl*sintr;
            z = radius*sinpl;
            dr3_polygon_normal(x*divr, y*divr, z*divr);
            dr3_polygon_vertex(x, y, z);
         }
         if (iphi < (nphi-1)) {
            x = radius*cospr*costr;
            y = radius*cospr*sintr;
            z = radius*sinpr;
            dr3_polygon_normal(x*divr, y*divr, z*divr);
            dr3_polygon_vertex(x, y, z);
            x = radius*cospr*costl;
            y = radius*cospr*sintl;
            z = radius*sinpr;
            dr3_polygon_normal(x*divr, y*divr, z*divr);
            dr3_polygon_vertex(x, y, z);
         } else {
            x = 0.0;
            y = 0.0;
            z = radius;
            dr3_polygon_normal(x*divr, y*divr, z*divr);
            dr3_polygon_vertex(x, y, z);
         }
         dr3_end_polygon();
      }
   }
}

void dr3_solid_rect(float width, float height)
{
/* front */
   dr3_begin_polygon();
   dr3_polygon_normal(0.0, 0.0, 1.0);
   dr3_polygon_vertex(0.0, 0.0, 0.0);
   dr3_polygon_vertex(width, 0.0, 0.0);
   dr3_polygon_vertex(width, height, 0.0);
   dr3_polygon_vertex(0.0, height, 0.0);
   dr3_end_polygon();
}

