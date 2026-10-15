/* dr3wire.c:  a collection of wire primitives */
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"

/* EXPLICIT DECLARATIONS ADDED BY BM */
extern void dr3_begin_polyline();
extern void dr3_polyline_vertex(float x,float y,float z);
extern void dr3_end_polyline();
void dr3_wire_block(float xmin, float xmax, float ymin, float ymax,
  float zmin, float zmax);

void dr3_wire_block(float xmin, float xmax, float ymin, float ymax,
  float zmin, float zmax)
{
/* front */
   dr3_begin_polyline();
   dr3_polyline_vertex(xmin, ymin, zmin);
   dr3_polyline_vertex(xmax, ymin, zmin);
   dr3_polyline_vertex(xmax, ymax, zmin);
   dr3_polyline_vertex(xmin, ymax, zmin);
   dr3_polyline_vertex(xmin, ymin, zmin);
   dr3_end_polyline();
/* back */
   dr3_begin_polyline();
   dr3_polyline_vertex(xmin, ymin, zmax);
   dr3_polyline_vertex(xmax, ymin, zmax);
   dr3_polyline_vertex(xmax, ymax, zmax);
   dr3_polyline_vertex(xmin, ymax, zmax);
   dr3_polyline_vertex(xmin, ymin, zmax);
   dr3_end_polyline();
/* left side */
   dr3_begin_polyline();
   dr3_polyline_vertex(xmin, ymin, zmin);
   dr3_polyline_vertex(xmin, ymax, zmin);
   dr3_polyline_vertex(xmin, ymax, zmax);
   dr3_polyline_vertex(xmin, ymin, zmax);
   dr3_polyline_vertex(xmin, ymin, zmin);
   dr3_end_polyline();
/* right side */
   dr3_begin_polyline();
   dr3_polyline_vertex(xmax, ymin, zmin);
   dr3_polyline_vertex(xmax, ymax, zmin);
   dr3_polyline_vertex(xmax, ymax, zmax);
   dr3_polyline_vertex(xmax, ymin, zmax);
   dr3_polyline_vertex(xmax, ymin, zmin);
   dr3_end_polyline();
/* bottom side */
   dr3_begin_polyline();
   dr3_polyline_vertex(xmin, ymin, zmin);
   dr3_polyline_vertex(xmax, ymin, zmin);
   dr3_polyline_vertex(xmax, ymin, zmax);
   dr3_polyline_vertex(xmin, ymin, zmax);
   dr3_polyline_vertex(xmin, ymin, zmin);
   dr3_end_polyline();
/* top side */
   dr3_begin_polyline();
   dr3_polyline_vertex(xmin, ymax, zmin);
   dr3_polyline_vertex(xmax, ymax, zmin);
   dr3_polyline_vertex(xmax, ymax, zmax);
   dr3_polyline_vertex(xmin, ymax, zmax);
   dr3_polyline_vertex(xmin, ymax, zmin);
   dr3_end_polyline();
}
