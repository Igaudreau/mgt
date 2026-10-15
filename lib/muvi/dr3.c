/* dr3.c:  3-dimensional drawing module 

   This module provides facilities for constructing 3-dimensional
   pictures consisting of text, markers, lines, surfaces and
   solids.  It includes various shading models, and an image
   stack which supports alpha-blending.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#include "drp.h"
#include "csg.h"
#include "dr3.h"

#undef DEBUG /* general debug */
#undef DEBUG_SOLID /* debug solid */
#undef DEBUG_POL /* debug polygons */
#undef DEBUG_SCENE /* debug scene */
#undef DEBUG_VIEW /* debug view */
#undef DEBUG_PROJECT /* debug projection */
#undef DEBUG_IMAGES /* debug images */
#undef DEBUG_SHADE /* debug shade */
#undef DEBUG_LIGHT /* debug light */
#undef DEBUG_DATA
#undef DEBUG_DRAW

#define MAXOP 100 /* number of nib-elements along the line */
#define MAXV 2048 /* vertical scan lines */
#define MAXBUF 128 /* vertices per polygon */
#define MAXMATRIX 8 /* matrices in stack */

/* COLOUR MODES */

static bool vector_bitmap,colour_bitmap,rgb_bitmap,rgbz_bitmap;
static int dev_class;

/* DATA FIELDS */

static bool no_s=TRUE;
static int nd_s,nc_s[MAXD],*sa_s;
static float x0_s[MAXD],dxdj_s[MAXD],djdx_s[MAXD],qmin_s,qmax_s,
  q0_s,dqdp_s,dpdq_s,*qa_s, *xe_s[MAXD];

/* SCENE */

static bool noscene=TRUE;
static int scene_index,scene_count=0;

/* VIEWING AND PROJECTION */

static float view_azim,view_altd,view_twist,view_distance,
  view_focus[MAXD4],view_point[MAXD4], view_matrix[MAXD4][MAXD4],
  project_anglev,project_near,project_far,project_hmin,project_hmax,
  project_vmin,project_vmax, dhdxz,dvdyz,dddz,h0,v0,d0,
  view_world[MAXD4][MAXD4], view_nx[MAXD],view_ny[MAXD],view_nz[MAXD],
  view_w0[MAXD];
static bool noview=TRUE,noproject=TRUE;

/* MATRIX STACK */

static float vm[MAXD4][MAXD4],vms[MAXMATRIX][MAXD4][MAXD4],
  vw[MAXD4][MAXD4],vws[MAXMATRIX][MAXD4][MAXD4];
static int nmatrix;

/* IMAGE BITMAP */

static int ih0,iv0,nh,nv,*cim=NULL;
static float dhidh,dvidv,hi0,vi0,*dim=NULL,*rgbim=NULL;
static bool noimage=TRUE,image_owned;

/* LIGHTING */

static float normal[MAXD],scalar;
static bool normal_defined=FALSE,scalar_defined=FALSE;

static float light_azim,light_altd,light_distance,light_point[MAXD4],
  light_ambient,light_candle,light_view[MAXD4];
static bool nolight=TRUE;

/* MATERIALS */

static int nmaterial,material_index;
static float material_reds[MAXMAT],material_greens[MAXMAT],
  material_blues[MAXMAT],material_roughs[MAXMAT],material_glows[MAXMAT],
  material_red,material_green,material_blue,material_rough,
  material_glow;
static bool nomaterial=TRUE,material_used[MAXMAT];
static char material_names[MAXMAT][SZNAM],material_name[SZNAM];

/* SHADING MODEL */

static float sky_red=0.9,sky_green=1.0,sky_blue=1.0,shade_red[MAXCM],
  shade_green[MAXCM],shade_blue[MAXCM],shade_didc[MAXCD],
  ranf=1.0,divran=(1.0/RANCON),shade_mat_r[MAXMAT],shade_mat_g[MAXMAT],
  shade_mat_b[MAXMAT];
static int sky_index,nd,nc[MAXCD],ncm,shade_index,shade_offset,
  ncolour=8,nshade=8,ncd,matin=0;
static bool random_shade=TRUE,hard_surface=TRUE,lucid_surface=FALSE,
  random_colour=FALSE,noshade=TRUE,wire_frame=FALSE;
static char shade_name[SZNAM],shade_mat[MAXMAT][SZNAM];

/* RANDOM SEQUENCE */

static bool norand=TRUE;
#define MAXRAND 32765 /* Random sequence size */
static float randaf[MAXRAND];
static int rand_count;

/* SOLID RENDERING */

static bool solid_mode=FALSE,surface_mode=TRUE;
static float solid_op_sign=1.0,nib_buffer_overflow=FALSE;
static int nib_count=0,nib_size,nib_ivmin,nib_ivmax,
  nib_ihmin[MAXV],nib_ihmax[MAXV], nibop_ivmin,nibop_ivmax,
  nibop_ihmin[MAXV],nibop_ihmax[MAXV];
typedef struct nib_type {
   struct nib_type *prev;
   float d;
   int op;
   int c;
   float rgb[3];
} nib_struct,*nib_ptr;
static nib_struct *nib_buffer;
static nib_ptr *nib,*nibop;

/* POLYGON BUFFER */

static float hbuf[MAXBUF],vbuf[MAXBUF],dbuf[MAXBUF],cbuf[MAXBUF][MAXCD];
static int nbuf;
static int ivmin,ivmax,ihmin[MAXV],ihmax[MAXV];
static float dmin[MAXV],dmax[MAXV],cmin[MAXV][MAXCD],cmax[MAXV][MAXCD];
static float hi,vi,di;
static bool first;

/* POLYLINE DATA */

static int np,ip0[MAXSEG],nps[MAXSEG],ns;
static float  hp[MAXPOL],vp[MAXPOL],dp[MAXPOL];
static bool polyline_buffer_fault=FALSE;


/* EXPLICIT DECLARATIONS ADDED BY BM */
int dr3_setup_rand();
void dr3_create_scene(float anglev, float azim, float altd, float twist, float distance, float f[MAXD], float near, float far, float hmin, float hmax, float vmin, float vmax);
void dr3_fresh_scene();
void dr3_zap_scene();
void dr3_view(float azim, float altd, float twist, float distance, float f[MAXD]);
void dr3_project(float anglev, float near, float far, float hmin, float hmax, float vmin, float vmax);
void dr3_put_scene(scene_struct *scp);
void dr3_get_scene(scene_struct *scp);
void dr3_pop_image(scene_struct *scp);
void dr3_push_image(scene_struct *scp);
void dr3_create_image(int *ios);
void dr3_alloc_image(int *ios);
void dr3_free_image();   
void dr3_clear_image();
void dr3_shade(char name[], int *ios);
void dr3_get_nl(int *n);
void dr3_get_nd(int *n);
void dr3_sky(float red, float green, float blue);
void dr3_lights(float azim, float altd, float distance, float fusion);
void dr3_zap_materials();   
void dr3_define_material(char material[], float red, float green, float blue, float rough, float glow);
void dr3_select_material(char material[]);
void dr3_begin_solid();
int dr3_end_solid();
#ifdef STL
  void dr3_cut_solid(char shape[], float par[], int npar);
#else
  void dr3_cut_solid(char shape[], float par[]);
#endif
#ifdef STL
  void dr3_add_solid(char shape[], float par[], int npar);
#else
  void dr3_add_solid(char shape[], float par[]);
#endif
void dr3_alloc_nibs();
void dr3_add_nibs();
void dr3_cut_nibs();
void dr3_fix_nib(int n, float d[], int op[], int c[], int *nnn);
void dr3_render_nibs();
void dr3_sort_nib(int n, float d[], int op[], int c[]);
void dr3_pop_nib(nib_struct *nib, int *n, float d[], int op[], int c[]); 
void dr3_push_nib(int n, float d[], int op[], int c[], nib_struct *(*nib));
void dr3_free_nibs();   
int dr3_begin_polygon();
void dr3_polygon_vertex(float x,float y,float z);
void dr3_polygon_normal(float nx,float ny,float nz);
void dr3_polygon_scalar(float s);
void dr3_end_polygon();
void dr3_zap_polylines();
void dr3_display_polylines();
void dr3_hide_polylines();
void dr3_begin_polyline();
void dr3_polyline_vertex(float x,float y,float z);
void dr3_end_polyline();
void dr3_push_matrix();
void dr3_pop_matrix();
void dr3_get_matrix(float a[MAXD4][MAXD4]);
void dr3_put_matrix(float a[MAXD4][MAXD4]);
void dr3_zap_matrix();
void dr3_rotran_matrix(float x, float y, float z, float azim, float altd, float twist);
void dr3_rot_matrix(int id, float angled);
void dr3_tran_matrix(float tx, float ty, float tz);
void dr3_study_polygon();
void dr3_begin_surface();
void dr3_end_surface();
void dr3_render_surface();
void dr3_render_solid();
void dr3_shade_cm(float h, float v, float d, float c[], float *dimk, int *cimk);
void dr3_shade_rgb(float h, float v, float d, float c[], float *dimk, float rgbimk[3]);
void dr3_shade_rgbz(float xv[MAXD4], float rgb[3]);
void dr3_light(float xv[MAXD4], float normal[MAXD], float *light);
void dr3_wtod(float x,float y,float z, float *h, float *v, float *d);
void dr3_get_class(int *class);
void dr3_rain(float phi, float *red, float *green, float *blue);
void dr3_scalar(float xw[MAXD4], float *s, bool *found);
void dr3_vtod(float xv[MAXD4], float *h, float *v, float *d);
void dr3_wtov(float xo[MAXD4], float xv[MAXD4]);
void dr3_dtov(float h,float v,float d, float xv[MAXD4]);
void dr3_vtow(float xv[MAXD4], float xw[MAXD4]);
void dr3_clock(float *cp);
void dr3_loads(int nc[], float *xe[MAXD4], int sa[], float qa[], float qmin, float qmax, int *ios);
void dr3_zaps();
void dr3_data_plane(float distance, int *ios);
void dr3_data_shade(int *ios);
void dr3_coerce_wire(int *ios);
void dr3_coerce_ld(int *ios);
void dr3_coerce_lm(int *ios);
void dr3_coerce_d(int *ios);
void dr3_seta(char name[]);
void dr3_set0i(char name[], int i);
void dr3_get_shade(scene_struct *scp);
void dr3_get_image_type(scene_struct *scp);
void dr3_default_limits(int nh,int nv, int *ivmin,int *ivmax,int ihmin[], int ihmax[]);
void dr3_extend_limits(int mh, int mv, int jvmin, int jvmax, int jhmin[], int jhmax[],
  int *ivmin, int *ivmax, int ihmin[], int ihmax[]);
void dr3_dump_scene(int unit, int *ios);
void dr3_load_scene(int unit, int *ios);
void dr3_wire_icon(float size, float x, float y, float z);
void dr3_solid_icon(float s, float size, float x, float y, float z);
void dr3_begin();
void dr3_end();

extern void dr3_zerom(float a[MAXD4][MAXD4]);
extern void dr3_mult(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4], 
  float c[MAXD4][MAXD4]);
extern void dr3_copy(float a[MAXD4][MAXD4], float b[MAXD4][MAXD4]);
extern void dr3_vmult(float a[MAXD4][MAXD4], float x[MAXD4],float y[MAXD4]);
extern void dr3_vmult3(float a[MAXD4][MAXD4], float x[MAXD4], float y[MAXD]);
extern void dr3_printmatrix(float a[MAXD4][MAXD4]);
extern void ccopy(char a[], int offset, int n, char b[]);
extern int cindex(char a[], char c);
extern void ccat(char a[], int lenb, char b[]);
extern void cbuild(char delim, char a[], int lenb, char b[]);
extern void cscan(int lena, char a[], int offset, int lenb, char b[], int *newoff);
extern bool csame(char a[], char b[]);
extern void slower(char a[]);
extern bool cabrev(char a[], char b[], int *offset);
extern void srmblk(char a[]);
extern void cvector(int ndiv, float v0, float dv, int *ipower, float *scale, int *maxfield, int *maxdp);
extern void cfloat(float v, int field, int ndp, int lena, char a[]);

extern void nam_alloc(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);
extern void nam_find(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);
extern void nam_free(char name[], int maxnames, bool used[], 
  char names[][SZNAM]);
#ifdef STL
  extern void dr3_solid(char shape[], float par[], int npar);
#else
  extern void dr3_solid(char shape[], float par[]);
#endif
extern void dr3_solid_cylinder(float radius, float height);
extern void dr3_solid_rectangle(float width, float height);
extern void dr3_solid_cone(float radius, float height);
extern void dr3_solid_block(float width, float height, float depth);
extern void dr3_solid_sphere(float radius);
extern void dr3_solid_rect(float width, float height);
extern int drp_prepare();
extern int drp_flush();
extern int drp_print();
extern int drp_avail();
extern int drp_pos_r(float rh, float rv);
extern int drp_pos_a(float ah, float av);
extern int drp_size_r(float rh, float rv);
extern int drp_size_a(float ah, float av);
extern int drp_open(char *name, char *type, int *ios);
extern int drp_select(char name[], int *ios);
extern int drp_close(char name[], int *ios);
extern int drp_get_count(int *count);
extern int drp_active(bool *active);
extern int drp_get_name(char name[]);
extern int drp_get_type(char type[]);
extern int drp_set_hv(float h1, float h2, float v1, float v2);
extern int drp_zap_hv();
extern int drp_get_hv(float *h1,float *h2,float *v1,float *v2);
extern int drp_setup_hv();
extern int drp_setup_clip();
extern int drp_erase();
extern int drp_rub(float h1, float h2, float v1, float v2);
extern int drp_avail_rub(bool *avail);
extern int drp_line(float h1, float v1, float h2, float v2);
extern int drp_polyline(int na, float ha[], float va[]);
extern int drp_set_line_width(float line_width);
extern int drp_get_line_width(float *line_width);
extern int drp_set_line_style(int style);
extern int drp_get_line_styles(int *n);
extern int drp_get_line_style(int *style);
extern int drp_bp();
extern int drp_ap(float h, float v);
extern int drp_ep();
extern int drp_clip(float *hh1, float *vv1, bool *clip1, float *hh2, float *vv2, bool *clip2, bool *empty);
extern int drp_mode(float cmin, float cmax, float c, int *ic);
extern int drp_rect(float h1, float h2, float v1, float v2);
extern int drp_avail_rect(bool *avail);
extern int drp_polygon(int na, float ha[], float va[]);
extern int drp_avail_polygon(bool *avail);
extern int drp_bg();
extern int drp_ag(float h, float v);
extern int drp_eg();
extern int drp_text(float h, float v, char str[], char just[2]);
extern int drp_text_box(char str[], float *dh, float *dv);
extern int drp_set_font_style(char fontstyle[]);
extern int drp_get_font_style(char fontstyle[]);
extern int drp_set_font_size(float fontsize);
extern int drp_get_font_size(float *fontsize);
extern int drp_avail_fontscale(bool *avail);
extern int drp_set_font();
extern int drp_get_char(float *height,float *width,float *depth);
extern int drp_avail_cursor(bool *avail);
extern int drp_cursor(float *h, float *v, int *button);
extern int drp_avail_colour(bool *avail);
extern int drp_set_colours(int offset, int size, float reds[], float greens[], float blues[]);
extern int drp_set_colour(int colour_index);
extern int drp_change_colour(bool *change);
extern int drp_get_colour_shape(int *offset, int *size);
extern int drp_get_colour(int *colour_index);
extern int drp_get_colours(int offset, int size, float reds[], float greens[], float blues[]);
extern int drp_avail_rgb(bool *avail);
extern int drp_set_rgb(float red, float green, float blue);
extern int drp_set_rgb_background(float red, float green, float blue);
extern int drp_get_rgb(float *red, float *green, float *blue);
extern int drp_get_rgb_background(float *red_background,float *green_background,float *blue_background);
extern int drp_avail_rgbz(bool *avail);
extern int rgbz_setup(float anglev, float near, float far, float hmin, float hmax, float vmin, float vmax);
extern int rgbz_begin();
extern int rgbz_end();
extern int rgbz_clear();
extern int rgbz_begin_polygon();
extern int rgbz_polygon_vertex(float x, float y, float z);
extern int rgbz_end_polygon();
extern int rgbz_begin_polyline();
extern int rgbz_polyline_vertex(float x, float y, float z);
extern int rgbz_end_polyline();
extern int rgbz_put_bitmap(int ih0,int iv0,int nh,int nv,float rgb_image[],float d_image[]);
extern int rgbz_get_bitmap(int ih0,int iv0,int nh,int nv,float rgb_image[],float d_image[]);
extern int drp_set_marker_scale(float scale);
extern int drp_set_marker_size(float marker_size);
extern int drp_get_marker_scale(float *scale);
extern int drp_get_marker_size(float *marker_size);
extern int drp_marker(int marker, float h, float v);
extern int drp_avail_bitmap(bool *avail);
extern int drp_bitmap(int *nh,int *nv,float *dh,float *dv,float *h0,float *v0);
extern int drp_get_bitmap(int *ih0,int *iv0,int *nh,int *nv,float *dhidh,float *dvidv,float *hi0,float *vi0);
extern int drp_colour_bitmap(int ih0,int iv0,int nh,int nv,int image[]);
extern int drp_rgb_bitmap(int ih0,int iv0,int nh,int nv,float rgb_image[]);
extern int drp_zap_colour();
extern int drp_zap_ramps();
extern int drp_print_ramps();
extern int drp_put_ramp(char name[], int nd, int nc[], int *offset);
extern int drp_get_ramp(char name[], int *nd, int nc[], int *offset);
extern int drp_find_ramp(int nd, char name[]);
extern int drp_zap_paints();
extern int drp_put_paint(char name[], int *offset);
extern int drp_get_paint(char name[], int *offset);
extern int drp_set_paint(int offset, float red, float green, float blue);
extern int drp_print_paints();
extern int drp_get_colour_space(int *space);
extern int crap();
extern int drp_gflush();
extern int ds_prepare();
extern int ds_flush();
extern int ds_path(char path[], int *ios);
extern int ds_method(char met[], int *ios);
extern int ds_create(char file[], int *unit, int *ios);
extern int ds_open(char file[], int *unit, int *ios);
extern int ds_close(int unit, int *ios);
extern int ds_index(char name[], int no, int lencomb, char comb[]);
extern int ds_put_begin(int unit, char name[], int *ios);
extern int ds_put_end(int unit, char name[], int *ios);
extern int ds_get_begin(int unit, char name[], int *ios);
extern int ds_get_end(int unit, char name[], int *ios);
extern int ds_put0i(int unit, char name[], int value, int *ios);
extern int ds_put1i(int unit, char name[], int nvalue, int value[], int *ios);
extern int ds_put0r(int unit, char name[], float value, int *ios);
extern int ds_put1r(int unit, char name[], int nvalue, float value[], int *ios);
extern int ds_put0b(int unit, char name[], bool value, int *ios);
extern int ds_put1b(int unit, char name[], int nvalue, bool value[], int *ios);
extern int ds_put0s(int unit, char name[], char value[], int *ios);
extern int ds_put0c(int unit, char name[], char value, int *ios);
extern int ds_put1c(int unit, char name[], int nvalue, char value[], int *ios);
extern bool ds_get(int unit, char name[], char type[], int *ios);
extern int ds_get0i(int unit, char name[], int *value, int *ios);
extern int ds_get1i(int unit, char name[], int nvalue, int value[], int *ios);
extern int ds_get0r(int unit, char name[], float *value, int *ios);
extern int ds_get1r(int unit, char name[], int nvalue, float value[], int *ios);
extern int ds_get0b(int unit, char name[], bool *value, int *ios);
extern int ds_get1b(int unit, char name[], int nvalue, bool value[], int *ios);
extern int ds_get0s(int unit, char name[], char value[], int *ios);
extern int ds_get0c(int unit, char name[], char *value, int *ios);
extern int ds_get1c(int unit, char name[], int nvalue, char value[], int *ios);
extern int merror(char text[]);
extern int ds_filex(char path[], char name[], char ext[], int lenfile, char file[]);
extern int ds_check_unit(int unit, int *ios);
extern int ds_put_header(int unit, char name[], char type[], int nvalue, int *ios);
extern int ds_get_prelude(int unit, char name[], char type[], int nvalue, int *ios);
extern int ds_get_header(int unit, int *ios);
extern int ds_check_header(char name[], char type[], int nvalue, int *ios);
extern int ds_indent(int unit, int indent, int *ios);
extern int ds_read_error(int unit, char proc[], int *ios);
extern int ds_write_error(int unit, char proc[], int *ios);
extern int ds_put1c_binary(int unit, int nvalue, char value[], int *ios);
extern int ds_get1c_binary(int unit, int *nvalue, char value[], int *ios);
extern int ds_put1i_binary(int unit, int nvalue, int value[], int *ios);
extern int ds_get1i_binary(int unit, int nvalue, int value[], int *ios);
extern int ds_put1r_binary(int unit, int nvalue, float value[], int *ios);
extern int ds_get1r_binary(int unit, int nvalue, float value[], int *ios);
extern int ds_put1b_binary(int unit, int nvalue, bool value[], int *ios);
extern int ds_get1b_binary(int unit, int nvalue, bool value[], int *ios);
extern int ds_put0s_binary(int unit, char value[], int *ios);
extern int ds_get0s_binary(int unit, int lenvalue, char value[], int *ios);
extern void dr3_wire_block(float xmin, float xmax, float ymin, float ymax,
  float zmin, float zmax);





/* SETUP */

int dr3_setup_rand() {
   int i;

   for (i=0; i < MAXRAND; i++) {
      randaf[i] = ranf*((float)rand()*divran-0.5);
   }
   rand_count = 0;
   norand = FALSE;
   return 0;
}

/* SCENE PRIMITIVES */

/* dr3_create_scene - setup viewing and projection transforms */
void dr3_create_scene(anglev, azim, altd, twist, distance, f, near, far,
  hmin, hmax, vmin, vmax)
   float anglev, azim,altd,twist, distance,f[MAXD],near,far,
     hmin,hmax,vmin,vmax;
{
   float hoff,voff,dh,dv,r,g,b;
   int mh,mv,ios;

/* Create random numbers */
   if (norand) {
      dr3_setup_rand();
   }
/* Find out what kind of device we have */
   dr3_get_class(&dev_class);
/* Setup boolean switchs */
   rgbz_bitmap = dev_class == DEV_RGBZ;
   rgb_bitmap = dev_class == DEV_RGB;
   colour_bitmap = dev_class == DEV_CM;
   vector_bitmap = dev_class == DEV_VEC;
#ifdef DEBUG_SCENE
   printf("dr3_create_scene:  dev_class=%d\n", dev_class);
#endif
/* View */
   dr3_view(azim, altd, twist, distance, f);
/* Projection */
   dr3_project(anglev, near, far, hmin, hmax, vmin, vmax);
/* lights (default) */
   if (nolight) {
      dr3_lights(30.0, 45.0, 10.0, 0.2);
   }
/* Get an image bitmap for current 2D viewport */
/*
   drp_set_hv(hmin, hmax, vmin, vmax);
   drp_get_bitmap(&ih0,&iv0,&nh,&nv,&dhidh,&dvidv,&hi0,&vi0);
*/
   drp_bitmap(&mh,&mv,&dh,&dv,&hoff,&voff);
   ih0 = (hmin-hoff)/dh;
   iv0 = (vmin-voff)/dv;
   nh = (hmax-hmin)/dh;
   nv = (vmax-vmin)/dv;
   dhidh = 1.0/dh;
   dvidv = 1.0/dv;
   hi0 = -hmin/dh;
   vi0 = -vmin/dv;
#ifdef DEBUG_SCENE
   printf("..ih0=%d,iv0=%d,nh=%d,nv=%d\n", ih0, iv0, nh, nv);
   printf("..dhidh=%f,dvidv=%f,hi0=%f,vi0=%f\n", dhidh,dvidv,hi0,vi0);
#endif
/* prepare polyline data structure */
   dr3_zap_polylines();
/* zap materials */
   dr3_zap_materials();
/* we have noshading model as yet */
   noshade = TRUE;
/* we have a scene! */
   noscene = FALSE;
   scene_count++;
   scene_index = scene_count;
/* RGBZ is special:  must allocate bitmaps here rather than in create*/
   if (rgbz_bitmap) {
      dr3_alloc_image(&ios);
      if (ios) {
         printf("error:  failed to allocate rgbz image\n");
      } else {
         dr3_shade("WIRE", &ios);
         if (ios) {
            printf("error:  failed to setup prelim rgbz wire\n");
         } else {
            drp_get_rgb(&r, &g, &b);
            drp_set_rgb(sky_red, sky_green, sky_blue);
            drp_rect(hmin, hmax, vmin, vmax);
            drp_set_rgb(r, g, b);
         }
      }
   }
}

/* dr3_fresh_scene - fresh scene (by displaying it) */
void dr3_fresh_scene()
{
   int ipaint;

#ifdef DEBUG_SCENE
   printf("dr3_fresh_scene:  called\n");
#endif
/* draw resultant bitmap image */
   if (noimage) {
   } else if (colour_bitmap) {
      drp_colour_bitmap(ih0, iv0, nh, nv, cim);
   } else if (rgb_bitmap) {
      drp_rgb_bitmap(ih0, iv0, nh, nv, rgbim);
   } else if (rgbz_bitmap) {
/*
      printf("dr3_fresh_scene:  not calling rgbz_put_bitmap\n");
      rgbz_put_bitmap(ih0, iv0, nh, nv, rgbim, dim);
*/
   }
/* refresh polylines if necessary */
   if (!rgbz_bitmap) {
      if (colour_bitmap) {
         drp_get_paint("BLACK", &ipaint);
         drp_set_colour(ipaint);
      }
      drp_set_hv(project_hmin,project_hmax,project_vmin,project_vmax);
      if (noimage) {
         dr3_display_polylines();
      } else {
         dr3_hide_polylines();
      }
   }
}

/* dr3_zap_scene - zap the entire scene */
void dr3_zap_scene() {
#ifdef DEBUG_SCENE
   printf("zap_scene:  called\n");
#endif
/* current image */
   dr3_free_image();
/* data fields */
   no_s = TRUE;
/* view,project,light,shade */
   noview = TRUE;
   noproject = TRUE;
   nolight = TRUE;
   noshade = TRUE;
/* materials */
   dr3_zap_materials();
/* matrix stack */
   nmatrix = TRUE;
/* polyline buffer */
   np = 0;
   ns = 0;
/* shading model */
   noshade = TRUE;
/* no scene! */
   noscene = TRUE;
}

/* dr3_view - setup viewing transformation */
void dr3_view(azim, altd, twist, distance, f)
   float azim,altd,twist, distance,f[MAXD];
{
   float sint,cost,sinp,cosp,nx[MAXD4],ny[MAXD4],nz[MAXD4],vf[MAXD4];
   int col,id;

#ifdef DEBUG_VIEW
   printf("dr3_view:  azim=%f,altd=%f,twist=%f,distance=%f\n",
     azim, altd, twist, distance);
   printf("..f=[%f,%f,%f]\n", f[0], f[1], f[2]);
#endif
/* copy parameters */
   view_azim = azim;
   view_altd = altd;
   view_twist = twist;
   view_distance = distance;
   for (id=0; id < MAXD; id++) {
      view_focus[id] = f[id];
   }
   view_focus[3] = 1.0;
/* Construct a matrix which translates and rotates world coordinates
   (x,y,z) so that they are positioned in the viewing frustrum
   frame (X,Y,Z)
*/
/* X-normal vector */
   sint = sin(DEGRAD*(90.0+azim));
   cost = cos(DEGRAD*(90.0+azim));
   sinp = 0.0;
   cosp = 1.0;
   nx[0] = cost*cosp;
   nx[1] = sint*cosp;
   nx[2] = sinp;
/* Y-normal vector */
   sint = sin(DEGRAD*(180.0+azim));
   cost = cos(DEGRAD*(180.0+azim));
   sinp = sin(DEGRAD*(90.0-altd));
   cosp = cos(DEGRAD*(90.0-altd));
   ny[0] = cost*cosp;
   ny[1] = sint*cosp;
   ny[2] = sinp;
/* Z-normal vector */
   nz[0] = cos(azim*DEGRAD)*cos(altd*DEGRAD);
   nz[1] = sin(azim*DEGRAD)*cos(altd*DEGRAD);
   nz[2] = sin(altd*DEGRAD);
/* viewpoint */
   for (col=0; col < MAXD; col++) {
      vf[col] = distance*nz[col];
      view_point[col] = f[col] + vf[col];
   }
   view_point[3] = 1.0;
#ifdef DEBUG_VIEW
   printf("..vf=[%f,%f,%f]\n",vf[0],vf[1],vf[2]);
   printf("..view_point=[%f,%f,%f]\n",view_point[0],view_point[1],
     view_point[2]);
#endif
/* construct viewing matrix */
   dr3_zerom(view_matrix);
   for (col=0; col < MAXD; col++) {
      view_matrix[0][col] = nx[col];
      view_matrix[1][col] = ny[col];
      view_matrix[2][col] = nz[col];
      view_matrix[0][3] = view_matrix[0][3] - nx[col]*view_focus[col];
      view_matrix[1][3] = view_matrix[1][3] - ny[col]*view_focus[col];
      view_matrix[2][3] = view_matrix[2][3] - nz[col]*view_focus[col];
   }
   view_matrix[2][3] = view_matrix[2][3] - view_distance;
   view_matrix[3][3] = 1.0;
   noview = FALSE;
   dr3_zerom(view_world);
   dr3_zap_matrix();
/* take a copy of the [nx,ny,nz] */
   for (id=0; id < MAXD; id++) {
      view_nx[id] = nx[id];
      view_ny[id] = ny[id];
      view_nz[id] = nz[id];
      view_w0[id] = view_focus[id] + view_distance*nz[id];
   }
#ifdef DEBUG_VIEW
   printf("nx=%f %f %f\n", nx[0],nx[1],nx[2]);
   printf("ny=%f %f %f\n", ny[0],ny[1],ny[2]);
   printf("nz=%f %f %f\n", nz[0],nz[1],nz[2]);
   printf("w0=vf=%f %f %f\n", view_w0[0],view_w0[1],view_w0[2]);
#endif
}

void dr3_project(anglev, near, far, hmin, hmax, vmin, vmax)
   float anglev,near,far,hmin,hmax,vmin,vmax;
{
   float tant;

#ifdef DEBUG_PROJECT
   printf("dr3_project:  anglev=%f,near=%f,far=%f\n", anglev,near,far);
   printf("..hmin=%f,hmax=%f,vmin=%f,vmax=%f\n", hmin, hmax,
     vmin, vmax);
#endif
/* copy parameters */
   project_anglev = anglev;
   project_near = near;
   project_far = far;
   project_hmin = hmin;
   project_hmax = hmax;
   project_vmin = vmin;
   project_vmax = vmax;
/* projecting [X,Y,Z] onto [h,v,d] */
   tant = tan(anglev*DEGRAD);
   dvdyz = (-0.5)*(hmax-hmin)/tant;
   v0 = 0.5*(vmin+vmax);
   dhdxz = dvdyz;
   h0 = 0.5*(hmin+hmax);
   dddz = 1.0/(near-far);
   d0 = near*dddz;
#ifdef DEBUG_PROJECT
   printf("..dhdxz=%f,h0=%f,dvdyz=%f,v0=%f\n", dhdxz, h0, dvdyz, v0);
   printf("..dddz=%f,d0=%f\n", dddz, d0);
   printf("..near=%f,far=%f\n", near, far);
#endif
/* rgbz-mode:  impose [X,Y,Z] */
   if (rgbz_bitmap) {
      rgbz_setup(anglev, near, far, hmin, hmax, vmin, vmax);
      rgbz_clear();
   }
   noproject = FALSE;
}

/* dr3_put_scene - caller puts new scene into dr3 */
void dr3_put_scene(scp)
   scene_struct *scp;
{
   int id,jd;

#ifdef DEBUG_SCENE
   printf("dr3_put_scene:  called\n");
#endif
/* view parameters */
   noview = scp->noview;
   if (!noview) {
      view_azim = scp->view_azim;
      view_altd = scp->view_altd;
      view_twist = scp->view_twist;
      view_distance = scp->view_distance;
      for (id=0; id < MAXD4; id++) {
         view_focus[id] = scp->view_focus[id];
         view_point[id] = scp->view_point[id];
         for (jd=0; jd < MAXD4; jd++) {
            view_matrix[id][jd] = scp->view_matrix[id][jd];
            view_world[id][jd] = scp->view_world[id][jd];
         }
      }
      for (id=0; id < MAXD; id++) {
         view_nx[id] = scp->view_nx[id];
         view_ny[id] = scp->view_ny[id];
         view_nz[id] = scp->view_nz[id];
         view_w0[id] = scp->view_w0[id];
      }
      dr3_zap_matrix();
   }
/* projection parameters */
   noproject = scp->noproject;
   if (!noproject) {
      project_anglev = scp->project_anglev;
      project_near = scp->project_near;
      project_far = scp->project_far;
      project_hmin = scp->project_hmin;
      project_hmax = scp->project_hmax;
      project_vmin = scp->project_vmin;
      project_vmax = scp->project_vmax;
      dhdxz = scp->project_dhdxz;
      h0 = scp->project_h0;
      dvdyz = scp->project_dvdyz;
      v0 = scp->project_v0;
      dddz = scp->project_dddz;
      d0 = scp->project_d0;
   }
/* light parameters */
   nolight = scp->nolight;
   if (!nolight) {
      light_azim = scp->light_azim;
      light_altd = scp->light_altd;
      light_distance = scp->light_distance;
      for (id=0; id < MAXD4; id++) {
         light_point[id] = scp->light_point[id];
         light_view[id] = scp->light_view[id];
      }
      light_ambient = scp->light_ambient;
      light_candle = scp->light_candle;
   }
/* scene */
   noscene = FALSE;
   scene_index = scp->scene_index;
}

/* dr3_get_scene - get scene from dr3 */
void dr3_get_scene(scp)
   scene_struct *scp;
{
   int id,jd;

#ifdef DEBUG_SCENE
   printf("dr3_get_scene:  called\n");
#endif
/* view parameters */
   scp->noview = noview;
   if (!noview) {
      scp->view_azim = view_azim;
      scp->view_altd = view_altd;
      scp->view_twist = view_twist;
      scp->view_distance = view_distance;
      for (id=0; id < MAXD4; id++) {
         scp->view_focus[id] = view_focus[id];
         scp->view_point[id] = view_point[id];
         for (jd=0; jd < MAXD4; jd++) {
            scp->view_matrix[id][jd] = view_matrix[id][jd];
            scp->view_world[id][jd] = view_world[id][jd];
         }
      }
      for (id=0; id < MAXD; id++) {
         scp->view_nx[id] = view_nx[id];
         scp->view_ny[id] = view_ny[id];
         scp->view_nz[id] = view_nz[id];
         scp->view_w0[id] = view_w0[id];
      }
   }
/* projection parameters */
   scp->noproject = noproject;
   if (!noproject) {
      scp->project_anglev = project_anglev;
      scp->project_near = project_near;
      scp->project_far = project_far;
      scp->project_hmin = project_hmin;
      scp->project_hmax = project_hmax;
      scp->project_vmin = project_vmin;
      scp->project_vmax = project_vmax;
      scp->project_dhdxz = dhdxz;
      scp->project_h0 = h0;
      scp->project_dvdyz = dvdyz;
      scp->project_v0 = v0;
      scp->project_dddz = dddz;
      scp->project_d0 = d0;
   }
/* light parameters */
   scp->nolight = nolight;
   if (!nolight) {
      scp->light_azim = light_azim;
      scp->light_altd = light_altd;
      scp->light_distance = light_distance;
      for (id=0; id < MAXD4; id++) {
         scp->light_point[id] = light_point[id];
         scp->light_view[id] = light_view[id];
      }
      scp->light_ambient = light_ambient;
      scp->light_candle = light_candle;
   }
   scp->scene_index = scene_index;
}

/* IMAGE PRIMITIVES */

/* dr3_pop_image - pop image from dr3 to caller */
void dr3_pop_image(scp)
   scene_struct *scp;
{
   int i,is,id;

#ifdef DEBUG_IMAGE
   printf("dr3_pop_image:  called\n");
#endif
/* copy image (if any) and its shading model */
   scp->noimage = noimage;
   if (noimage) {
      scp->noshade = TRUE;
   } else {
/*    Image part */
      scp->dhidh = dhidh;
      scp->hi0 = hi0;
      scp->dvidv = dvidv;
      scp->vi0 = vi0;
      scp->ih0 = ih0;
      scp->iv0 = iv0;
      scp->nh = nh;
      scp->nv = nv;
      scp->dev_class = dev_class;
      scp->vector_bitmap = vector_bitmap;
      scp->colour_bitmap = colour_bitmap;
      scp->rgb_bitmap = rgb_bitmap;
      scp->rgbz_bitmap = rgbz_bitmap;
      if (!image_owned) {
         printf("crazy:  attempt to pop image not owned by dr3\n");
      }
      scp->image_owned = TRUE;
      image_owned = FALSE;
      scp->dim = dim;
#ifdef DEBUG
      printf("dr3_get_scene:  dim=%d\n", (int)dim);
#endif
      dim = NULL;
      if (cim != NULL) {
         scp->cim = cim;
#ifdef DEBUG
         printf("dr3_get_scene:  cim=%d\n", (int)cim);
#endif
         cim = NULL;
      }
      if (rgbim != NULL) {
         scp->rgbim = rgbim;
#ifdef DEBUG
         printf("dr3_get_scene:  rgbim=%d\n", (int)rgbim);
#endif
         rgbim = NULL;
      }
/*    Shading model part */
      scp->noshade = noshade;
      if (!noshade) {
         scp->shade_index = shade_index;
         ccopy(shade_name, 0, SZNAM, scp->shade_name);
         scp->shade_nd = nd;
         for (id=0; id < nd; id++) {
            scp->shade_nc[id] = nc[id];
            scp->shade_didc[id] = (float)(nc[id]-1);
         }
         if (colour_bitmap) {
            scp->shade_ncm = ncm;
            for (i=0; i < ncm; i++) {
               scp->shade_red[i] = shade_red[i];
               scp->shade_green[i] = shade_green[i];
               scp->shade_blue[i] = shade_blue[i];
            }
            scp->shade_offset = shade_offset;
         }
         if (shade_index == SHADE_LM) {
            for (i=0; i < nc[1]; i++) {
               ccopy(shade_mat[i],0,SZNAM,scp->shade_mat[i]);
               scp->shade_mat_r[i] = shade_mat_r[i];
               scp->shade_mat_g[i] = shade_mat_g[i];
               scp->shade_mat_b[i] = shade_mat_b[i];
            }
         }
         scp->sky_red = sky_red;
         scp->sky_green = sky_green;
         scp->sky_blue = sky_blue;
      }
   }
/* Polyline data */
   scp->ns = ns;
   for (is=0; is < ns; is++) {
      scp->ip0[is] = ip0[is];
      scp->nps[is] = nps[is];
   }
   scp->np = np;
   for (i=0; i < np; i++) {
      scp->hp[i] = hp[i];
      scp->vp[i] = vp[i];
      scp->dp[i] = dp[i];
   }
/* zap image flags */
   noimage = TRUE;
   noshade = TRUE;
   np = 0;
   ns = 0;
}

/* dr3_push_image - push image from caller into dr3 */
void dr3_push_image(scp)
   scene_struct *scp;
{
   int is,i,id;

#ifdef debug_SCENE
   printf("dr3_push_image:  called\n");
#endif
/* copy image (if any) with its shading model */
   noimage = scp->noimage;
   if (noimage) {
      noshade = TRUE;
   } else {
/*    The image part */
      dhidh = scp->dhidh;
      hi0 = scp->hi0;
      dvidv = scp->dvidv;
      vi0 = scp->vi0;
      ih0 = scp->ih0;
      iv0 = scp->iv0;
      nh = scp->nh;
      nv = scp->nv;
      dev_class = scp->dev_class;
      vector_bitmap = scp->vector_bitmap;
      colour_bitmap = scp->colour_bitmap;
      rgb_bitmap = scp->rgb_bitmap;
      rgbz_bitmap = scp->rgbz_bitmap;
      if (!scp->image_owned) {
         printf("crazy:  attempt to push image not owned by scp\n");
      }
      image_owned = TRUE;
      scp->image_owned = FALSE;
      dim = scp->dim;
      scp->dim = NULL;
#ifdef DEBUG_SCENE
      printf("dr3_put_scene:  dim=%d\n", (int)dim);
#endif
      if (scp->cim != NULL) {
         cim = scp->cim;
         scp->cim = NULL;
#ifdef DEBUG_SCENE
         printf("dr3_put_scene:  cim=%d\n", (int)cim);
#endif
      }
      if (scp->rgbim != NULL) {
         rgbim = scp->rgbim;
         scp->rgbim = NULL;
#ifdef DEBUG_SCENE
         printf("dr3_put_scene:  rgbim=%d\n", (int)rgbim);
#endif
      }
/*    Shading model part */
      noshade = scp->noshade;
      if (!noshade) {
         shade_index = scp->shade_index;
         ccopy(scp->shade_name, 0, SZNAM, shade_name);
         nd = scp->shade_nd;
         for (id=0; id < nd; id++) {
            nc[id] = scp->shade_nc[id];
            shade_didc[id] = (float)(nc[id]-1);
         }
         if (colour_bitmap) {
            ncm = scp->shade_ncm;
            for (i=0; i < ncm; i++) {
               shade_red[i] = scp->shade_red[i];
               shade_green[i] = scp->shade_green[i];
               shade_blue[i] = scp->shade_blue[i];
            }
            shade_offset = scp->shade_offset;
         }
         if (shade_index == SHADE_LM) {
            for (i=0; i < nc[1]; i++) {
               ccopy(scp->shade_mat[i], 0, SZNAM, shade_mat[i]);
               shade_mat_r[i]  = scp->shade_mat_r[i];
               shade_mat_g[i]  = scp->shade_mat_g[i];
               shade_mat_b[i]  = scp->shade_mat_b[i];
            }
         }
         sky_red = scp->sky_red;
         sky_green = scp->sky_green;
         sky_blue = scp->sky_blue;
      }
   }
/* Polyline data */
   ns = scp->ns;
   for (is=0; is < ns; is++) {
      ip0[is] = scp->ip0[is];
      nps[is] = scp->nps[is];
   }
   np = scp->np;
   for (i=0; i < np; i++) {
      hp[i] = scp->hp[i];
      vp[i] = scp->vp[i];
      dp[i] = scp->dp[i];
   }
}

/* dr3_create_image - create a new  image */
void dr3_create_image(ios)
   int *ios;
{
/* if there is an image, free its memory and report */
   if (!rgbz_bitmap) {
      dr3_free_image();
/*    we do need a shading model */
      if (noshade) {
         printf("error:  need shading model to begin image\n");
         *ios = ERR;
      } else {
/*       allocate image memory */
         dr3_alloc_image(ios);
         if (*ios) {
            printf("unable to allocate memory for image\n");
         } else {
/*          clear image depth and colour */
            dr3_clear_image();
         }
      }
   }
}

/* dr3_alloc_image - allocate memory for image */
void dr3_alloc_image(ios) 
   int *ios;
{

/* free any pre-existing arrays */
   dr3_free_image();
/* Create depth-array dim[] */
   *ios = OK;
   dim = (float*)malloc((unsigned)(nh*nv*sizeof(float)));
   if (dim == NULL) {
      printf("error:  allocating memory for depth-array\n");
      *ios = ERR;
   }
/* Create cmap-arrays: cim[] */
   if (colour_bitmap) {
      cim = (int*)malloc((unsigned)(nh*nv*sizeof(int)));
      if (cim == NULL) {
         printf("error:  allocating memory for colour-array\n");
         *ios = ERR;
      }
/* Create rgb-array: rgbim[] */
   } else if (rgb_bitmap || rgbz_bitmap) {
      rgbim = (float*)malloc((unsigned)(3*nh*nv*sizeof(float)));
      if (rgbim == NULL) {
         printf("error:  allocating memory for rgb-array\n");
         *ios = ERR;
      }
   }
   if (*ios) {
      printf("error:  allocating image memory\n");
      dr3_free_image();
   } else {
      noimage = FALSE;
      image_owned = TRUE;
   }
}

/* dr3_free_image - free image memory */
void dr3_free_image() {

/* Free depth-array dim[] */
   if (dim != NULL) {
      free(dim);
      dim = NULL;
   }
   if (cim != NULL) {
      free(cim);
      cim = NULL;
   }
   if (rgbim != NULL) {
      free(rgbim);
      rgbim = NULL;
   }
   noimage = TRUE;
}

/* dr3_clear_image - clear image depth and colour */
void dr3_clear_image() {
   int k,n,m;

#ifdef DEBUG_IMAGES
   printf("dr3_clear_image:  called\n");
#endif
   if (noimage) {
      printf("error:  no image to clear\n");
   } else {
      n = nh*nv;
      for (k=0; k < n; k++) {
         dim[k] = BIG;
      }
      if (colour_bitmap) {
         for (k=0; k < n; k++) {
            cim[k] = shade_offset;
         }
      } else if (rgb_bitmap || rgbz_bitmap) {
         m = 0;
         for (k=0; k < n; k++) {
            rgbim[m++] = sky_red;
            rgbim[m++] = sky_green;
            rgbim[m++] = sky_blue;
         }
      }
/*    bitmap pixels affected
      for (iv=0; iv < nv; iv++) {
         bitmap_ihmin[iv] = nh;
         bitmap_ihmax[iv] = (-1);
      }
      bitmap_ivmin = nv;
      bitmap_ivmax = (-1);
*/
   }
}

/* SHADE PRIMITIVES */

/* dr3_shade - setup shading model */
void dr3_shade(name, ios)
   char name[];
   int *ios;
{
   int n,i,j,k,offset,nl,id;
   bool write_ramp;
   float red,green,blue,phi,dcdi[MAXCD];
   for (i=0; i < MAXCD; i++)
     dcdi[i] = 0.0;			// BM

#ifdef DEBUG_SHADE
   printf("dr3_shade:  name=%s\n", name);
#endif
/* set some hard parameters */
   *ios = OK;
   noshade = TRUE;
   write_ramp = FALSE;
   if (csame("LM",name)) {
      shade_index = SHADE_LM;
      ncd = 1;
      if (colour_bitmap) {
         drp_get_ramp(name, &nd, nc, &offset);
         if (offset < 0) {
            nd = 2;
            dr3_get_nl(&nl);
#ifdef DEBUG_SHADE
            printf("dr3_shade:  nl=%d\n", nl);
#endif
            drp_get_colour_space(&n);
            if (nl > 0) {
               nc[0] = nl;
               if (nl*nmaterial > n) {
                  nc[1] = 1;
               } else {
                  nc[1] = nmaterial;
               }
            } else if ((nshade+1)*nmaterial > n) {
               if (nmaterial*8 > n) {
                  if ((nshade+1) > n) {
                     nc[0] = n;
                  } else {
                     nc[0] = (nshade+1);
                  }
                  nc[1] = 1;
               } else {
                  nc[0] = n/nmaterial;
                  nc[1] = nmaterial;
               }
            } else {
               nc[0] = (nshade+1);
               nc[1] = nmaterial;
            }
            ncm = nc[0]*nc[1];
            if (ncm > n) {
               printf("error:  no colour space for LM image\n");
               *ios = ERR;
            } else if (nc[0] < 3) {
               printf("error:  nshade reduced too much for LM image\n");
               *ios = ERR;
            } else {
               if (nc[0] < (nshade+1)) {
                  printf("warning:  nshade reduced for LM image\n");
               }
               if (nc[1] < nmaterial) {
                  printf("warning:  nmaterial reduced for LM image\n");
               }
               dcdi[0] = 1.0/(float)(nc[0]-2);
               shade_didc[0] = (float)(nc[0]-1);
               k = 0;
/*             Material loop */
               for (i=0; i < nc[1]; i++) {
                  ccopy(material_names[i], 0, SZNAM, shade_mat[i]);
                  red = material_reds[i];
                  green = material_greens[i];
                  blue = material_blues[i];
                  if (i == 0) {
                     shade_red[k] = sky_red;
                     shade_green[k] = sky_green;
                     shade_blue[k] = sky_blue;
                     k++;
                  } else {
                     shade_red[k] = red;
                     shade_green[k] = green;
                     shade_blue[k] = blue;
                     k++;
                  }
                  for (j=1; j < nc[0]; j++) {
                     phi = dcdi[0]*(float)(j-1);
                     phi = sqrt(phi);
                     shade_red[k] = red*phi;
                     shade_green[k] = green*phi;
                     shade_blue[k] = blue*phi;
                     k++;
                  }
               }
               write_ramp = TRUE;
            }
         } else if (nd != 2) {
            printf("error:  LM ramp has wrong shape\n");
            *ios = ERR;
         } else {
            printf("note:  use existing LM ramp\n");
            shade_offset = offset;
            shade_didc[0] = (float)(nc[0]-1);
         }
      } else if (rgb_bitmap || rgbz_bitmap) {
         nd = 2;
         nc[0] = (nshade+1);
         nc[1] = nmaterial;
         for (id=0; id < nd; id++) {
            shade_didc[id] = (float)(nc[id]-1);
         }
         for (i=0; i < nmaterial; i++) {
            shade_mat_r[i] = material_reds[i];
            shade_mat_g[i] = material_greens[i];
            shade_mat_b[i] = material_blues[i];
            ccopy(material_names[i], 0, SZNAM, shade_mat[i]);
         }
      }
   } else if (csame("LD",name)) {
      shade_index = SHADE_LD;
      ncd = 2;
      if (colour_bitmap) {
         drp_get_ramp(name, &nd, nc, &offset);
         if (offset < 0) {
            nd = 2;
            dr3_get_nl(&nl);
#ifdef DEBUG_SHADE
            printf("dr3_shade:  nl=%d\n", nl);
#endif
            drp_get_colour_space(&n);
            if (nl > 0) {
               nc[0] = nl;
               if (nl*(ncolour+1) > n) {
                  nc[1] = n/nl;
               } else {
                  nc[1] = (ncolour+1);
               }
            } else if ((nshade+1)*(ncolour+1) > n) {
               nc[0] = (int)sqrt((float)n);
               nc[1] = n/nc[0];
            } else {
               nc[0] = (nshade+1);
               nc[1] = (ncolour+1);
            }
            ncm = nc[0]*nc[1];
            if (ncm > n) {
               printf("error:  no colour space for LD image\n");
               *ios = ERR;
            } else if (nc[0] < 3) {
               printf("error:  nshade reduced too much for LD image\n");
               *ios = ERR;
            } else if (nc[1] < 3) {
               printf("error:  ncolour reduced too much for LD image\n");
               *ios = ERR;
            } else {
               if (nc[0] < (nshade+1)) {
                  printf("warning:  nshade reduced for LD image\n");
               }
               if (nc[1] < (ncolour+1)) {
                  printf("warning:  ncolour reduced for LD image\n");
               }
               dcdi[0] = 1.0/(float)(nc[0]-2);
               shade_didc[0] = (float)(nc[0]-1);
               dcdi[1] = 1.0/(float)(nc[1]-2);
               shade_didc[1] = (float)(nc[1]-1);
               k = 0;
/*             Data loop */
               for (i=0; i < nc[1]; i++) {
/*                No data:  L only */
                  if (i == 0) {
/*                   Sky light */
                     shade_red[k] = sky_red;
                     shade_green[k] = sky_green;
                     shade_blue[k] = sky_blue;
                     k++;
/*                   Grey light loop */
                     for (j=1; j < nc[0]; j++) {
                        phi = dcdi[0]*(float)(j-1);
                        phi = sqrt(phi);
                        shade_red[k] = phi;
                        shade_green[k] = phi;
                        shade_blue[k] = phi;
                        k++;
                     }
/*                Data:  LD or D */
                  } else {
                     phi = (float)(i-1)*dcdi[1];
                     dr3_rain(phi, &red, &green, &blue);
/*                   Data only */
                     shade_red[k] = red;
                     shade_green[k] = green;
                     shade_blue[k] = blue;
                     k++;
/*                   Light loop */
                     for (j=1; j < nc[0]; j++) {
                        phi = dcdi[0]*(float)(j-1);
                        phi = sqrt(phi);
                        shade_red[k] = red*phi;
                        shade_green[k] = green*phi;
                        shade_blue[k] = blue*phi;
                        k++;
                     }
                  }
               }
               write_ramp = TRUE;
             }
         } else if (nd != 2) {
            printf("error:  LD ramp has wrong shape\n");
            *ios = ERR;
         } else {
            printf("note:  use existing LD ramp\n");
            shade_offset = offset;
            shade_didc[0] = (float)(nc[0]-1);
            shade_didc[1] = (float)(nc[1]-1);
         }
      } else if (rgb_bitmap || rgbz_bitmap) {
         nd = 2;
         nc[0] = (nshade+1);
         nc[1] = (ncolour+1);
         for (id=0; id < nd; id++) {
            shade_didc[id] = (float)(nc[id]-1);
         }
      }
   } else if (csame("L",name)) {
      shade_index = SHADE_L;
      ncd = 1;
      if (colour_bitmap) {
         drp_get_ramp(name, &nd, nc, &offset);
         if (offset < 0) {
            nd = 1;
            dr3_get_nl(&nl);
#ifdef DEBUG_SHADE
            printf("dr3_shade:  nl=%d\n", nl);
#endif
            drp_get_colour_space(&n);
            if (nl > n) {
               printf("error:  no colour space for L image\n");
               *ios = ERR;
            } else {
               if ((nshade+1) > n) {
                  nc[0] = n;
               } else {
                  nc[0] = (nshade+1);
               }
               if (nc[0] < 3) {
                  printf("error:  nshade reduced too much for L image\n");
                  *ios = ERR;
               } else {
                  if (nc[0] < (nshade+1)) {
                     printf("warning:  nshade reduced for L image\n");
                  }
                  ncm = nc[0];
                  dcdi[0] = 1.0/(float)(nc[0]-2);
                  shade_didc[0] = (float)(nc[0]-1);
                  shade_red[0] = sky_red;
                  shade_green[0] = sky_green;
                  shade_blue[0] = sky_blue;
                  for (i=1; i < nc[0]; i++) {
                     phi = dcdi[0]*(float)(i-1);
                     phi = sqrt(phi);
                     shade_red[i] = phi;
                     shade_green[i] = phi;
                     shade_blue[i] = phi;
                  }
                  write_ramp = TRUE;
               }
            }
         } else if (nd != 1) {
            printf("error:  L ramp wrong shape\n");
            *ios = ERR;
         } else {
            printf("note:  use existing L ramp\n");
            shade_offset = offset;
            shade_didc[0] = (float)(nc[0]-1);
         }
      } else if (rgb_bitmap || rgbz_bitmap) {
         nd = 1;
         nc[0] = (nshade+1);
         for (id=0; id < nd; id++) {
            shade_didc[id] = (float)(nc[id]-1);
         }
      }
   } else if (csame("WIRE",name)) {
      shade_index = SHADE_WIRE;
      nd = 0;
   } else if (csame("D",name)) {
      shade_index = SHADE_D;
      ncd = 1;
      if (colour_bitmap) {
         drp_get_ramp(name, &nd, nc, &offset);
         if (offset < 0) {
            nd = 2;
            dr3_get_nd(&nl);
            drp_get_colour_space(&n);
            if (nl > 0) {
               nc[0] = nl;
            } else {
               nc[0] = (ncolour+1);
            }
            ncm = nc[0];
            if (ncm > n) {
               printf("error:  no colour space for D image\n");
               *ios = ERR;
            } else if (nc[0] < 3) {
               printf("error:  ncolour reduced too much for D image\n");
               *ios = ERR;
            } else {
               if (nc[0] < (ncolour+1)) {
                  printf("warning:  ncolour reduced for D image\n");
               }
               dcdi[0] = 1.0/(float)(nc[0]-2);
               shade_didc[0] = (float)(nc[0]-1);
/*             Sky light */
               k = 0;
               shade_red[k] = sky_red;
               shade_green[k] = sky_green;
               shade_blue[k] = sky_blue;
               k++;
/*             Data loop */
               for (i=1; i < nc[0]; i++) {
                  phi = (float)(i-0)*dcdi[1];
                  dr3_rain(phi, &red, &green, &blue);
                  shade_red[k] = red;
                  shade_green[k] = green;
                  shade_blue[k] = blue;
                  k++;
               }
               write_ramp = TRUE;
             }
         } else if (nd != 1) {
            printf("error:  D ramp has wrong shape\n");
            *ios = ERR;
         } else {
            printf("note:  use existing D ramp\n");
            shade_offset = offset;
            shade_didc[0] = (float)(nc[0]-1);
         }
      } else if (rgb_bitmap || rgbz_bitmap) {
         nd = 1;
         nc[0] = (ncolour+1);
         for (id=0; id < nd; id++) {
            shade_didc[id] = (float)(nc[id]-1);
         }
      }
   } else {
      printf("error:  unknown shading model\n");
      *ios = ERR;
   }
/* if all is well, create colour ramp */
   if (!*ios) {
      ccopy(name, 0, SZNAM, shade_name);
      if (shade_index != SHADE_WIRE) {
         if (colour_bitmap) {
            if (write_ramp) {
               drp_put_ramp(shade_name, nd, nc, &shade_offset);
               if (shade_offset < 0) {
                  printf("error:  putting shading model ramp\n");
                  *ios = ERR;
               } else {
                  drp_set_colours(shade_offset, ncm, shade_red, 
                    shade_green, shade_blue);
               }
            }
         }
      }
   }
   noshade = (*ios);
}

/* dr3_get_nl - get L-size for any existing Lx shade */
void dr3_get_nl(n)
   int *n;
{
   int ndi,nci[MAXCD],offset;

/* Try for LM */
   drp_get_ramp("LM", &ndi, nci, &offset);
   if (offset >= 0) {
      *n = nci[0];
/* Otherwise try for L */
   } else {
      drp_get_ramp("L", &ndi, nci, &offset);
      if (offset >= 0) {
         *n = nc[0];
      } else {
         *n = 0;
      }
   }
}

/* dr3_get_nd - get D-size for any existing xD shade */
void dr3_get_nd(n)
   int *n;
{
   int ndi,nci[MAXCD],offset;

/* Try for LD */
   drp_get_ramp("LD", &ndi, nci, &offset);
   if (offset >= 0) {
      *n = nci[1];
/* Otherwise try for D */
   } else {
      drp_get_ramp("D", &ndi, nci, &offset);
      if (offset >= 0) {
         *n = nc[0];
      } else {
         *n = 0;
      }
   }
}

/* dr3_sky - define sky colour */
void dr3_sky(red, green, blue)
   float red,green,blue;
{
   float reds[1],greens[1],blues[1];

   sky_red = red;
   sky_green = green;
   sky_blue = blue;
#ifdef DEBUG_SHADE
   printf("dr3_sky:  rgb=[%f,%f,%f]\n", red, green, blue);
#endif
   if (colour_bitmap) {
      drp_put_paint("SKY", &sky_index);
      if (sky_index < 0) {
         printf("error:  no paint room for sky\n");
      } else {
         reds[0] = red;
         greens[0] = green;
         blues[0] = blue;
         drp_set_colours(sky_index, 1, reds, greens, blues);
      }
   }
}

/* LIGHTING PRIMITIVES */

/* dr3_lights - define lighting model */
void dr3_lights(azim, altd, distance, fusion)
   float azim,altd,distance,fusion;
{

#ifdef DEBUG_LIGHT
   printf("dr3_lights:  azim=%f,altd=%f,distance=%f,fusion=%f\n",
     azim, altd, distance, fusion);
#endif
/* define light position relative to view-focus */
   light_point[0] = distance*cos(azim*DEGRAD)*cos(altd*DEGRAD) + 
     view_focus[0];
   light_point[1] = distance*sin(azim*DEGRAD)*cos(altd*DEGRAD) + 
     view_focus[1];
   light_point[2] = distance*sin(altd*DEGRAD) + view_focus[2];
   light_point[3] = 1.0;
#ifdef DEBUG_LIGHT
   printf("dr3_lights:  light_point=%f,%f,%f\n", light_point[0],
     light_point[1], light_point[2]);
#endif
/* copy other values */
   light_azim = azim;
   light_altd = altd;
   light_distance = distance;
   if (fusion > 0.1) {
      light_ambient = fusion;
   } else {
      light_ambient = 0.1;
   }
   light_candle = 1.0 - fusion;
#ifdef DEBUG_LIGHT
   printf("..light_ambient=%f,light_candle=%f\n", light_ambient,
     light_candle);
#endif
/* transform to viewing coordinates */
   if (noview) {
      light_view[0] = light_point[0];
      light_view[1] = light_point[1];
      light_view[2] = light_point[2];
      light_view[3] = 1.0;
   } else {
      dr3_wtov(light_point, light_view);
   }
   nolight = FALSE;
}

/* MATERIAL PRIMITIVES */

void dr3_zap_materials() {
   nmaterial = 1;
   ccopy("DEFAULT", 0, SZNAM, material_names[0]);
   material_used[0] = TRUE;
   material_reds[0] = 1.0;
   material_greens[0] = 1.0;
   material_blues[0] = 1.0;
   material_roughs[0] = 0.5;
   material_glows[0] = 0.0;
   ccopy("DEFAULT", 0, SZNAM, material_name);
   material_index = 0;
   nomaterial = FALSE;
   dr3_select_material("DEFAULT");
}

void dr3_define_material(material, red, green, blue, rough, glow)
   char material[];
   float red,green,blue,rough,glow;
{
   int i;

   nam_find(material, nmaterial, material_used,
     material_names, &i);
   if (i >= 0) {
      material_reds[i] = red;
      material_greens[i] = green;
      material_blues[i] = blue;
      material_roughs[i] = rough;
      material_glows[i] = glow;
      dr3_select_material(material);
   } else if (nmaterial > MAXMAT) {
      printf("error:  too many materials\n");
      ccopy(material, 0, SZNAM, material_name);
   } else {
      ccopy(material, 0, SZNAM, material_names[nmaterial]);
      material_used[nmaterial] = TRUE;
      material_reds[nmaterial] = red;
      material_greens[nmaterial] = green;
      material_blues[nmaterial] = blue;
      material_roughs[nmaterial] = rough;
      material_glows[nmaterial] = glow;
      ccopy(material, 0, SZNAM, material_name);
      nmaterial++;
      dr3_select_material(material);
   }
}

void dr3_select_material(material)
   char material[];
{
   int i;

   if (nmaterial < 1) {
      printf("error:  no materials\n");
   } else {
      nam_find(material, nmaterial, material_used,
        material_names, &i);
      if (i < 0) {
         printf("error:  material=%s not found\n", material);
         material_index = 0;
         material_red = material_reds[0];
         material_green = material_greens[0];
         material_blue = material_blues[0];
         material_rough = material_roughs[0];
         material_glow = material_glows[0];
      } else {
         ccopy(material_names[i], 0, SZNAM, material_name);
         material_red = material_reds[i];
         material_green = material_greens[i];
         material_blue = material_blues[i];
         material_rough = material_roughs[i];
         material_glow = material_glows[i];
         material_index = i;
         nomaterial = FALSE;
      }
   }
}

/* SOLID PRIMITIVES */

void dr3_begin_solid() {
/* Select solid mode */
   if (!rgbz_bitmap) {
      solid_mode = TRUE;
      surface_mode = FALSE;
   }
/* select-check material index (assume LM) */
   if (shade_index == SHADE_LM) {
      if (material_index < nc[1]) {
         matin = material_index;
      } else {
         matin = 0;
      }
   }
/* Create a null "nib" images */
   if (solid_mode) {
      dr3_alloc_nibs();
   }
}

int dr3_end_solid() {
/* Only if in solid mode */
   if (solid_mode) {
/*    Transfer the "nib" image to the current image */
      dr3_render_nibs();
/*    Free the nib image */
      dr3_free_nibs();
      solid_mode = FALSE;
      surface_mode = TRUE;
   }
/* update screen? */
   if (colour_bitmap) {
      drp_colour_bitmap(ih0, iv0, nh, nv, cim);
   } else if (rgb_bitmap) {
      drp_rgb_bitmap(ih0, iv0, nh, nv, rgbim);
   }
   return 0;
}

#ifdef STL
  void dr3_cut_solid(shape, par, npar)
   int npar;
#else
  void dr3_cut_solid(shape, par)
#endif
   char shape[];
   float par[];
{
   if (solid_mode) {
/*    put it into cut-mode */
      solid_op_sign = (-1.0);
#ifdef STL
      dr3_solid(shape, par, npar);
#else
      dr3_solid(shape, par);
#endif
/*    apply nibop to current nib array */
      dr3_cut_nibs();
/*    revert back to add-mode */
      solid_op_sign = 1.0;
   }
}

#ifdef STL
  void dr3_add_solid(shape, par, npar)
   int npar;
#else
  void dr3_add_solid(shape, par)
#endif
   char shape[];
   float par[];
{
#ifdef STL
  #ifdef DEBUG
    printf("dr3_add_solid:  shape=%s,npar=%d\n", shape, npar);
  #endif
  dr3_solid(shape, par, npar);
#else
  dr3_solid(shape, par);
#endif
/* apply nibop to nib array */
   if (solid_mode) {
      dr3_add_nibs();
   }
}

/* dr3_alloc_nibs - allocate temporary (current,op) nib arrays */
void dr3_alloc_nibs() {

#ifdef DEBUG_SOLID
   printf("dr3_alloc_nibs:  called\n");
#endif
/* nib buffer */
   nib_buffer_overflow = FALSE;
   nib_count = 0;
/* clear (via limits) */
   dr3_default_limits(nh, nv, &nib_ivmin, &nib_ivmax, 
     nib_ihmin, nib_ihmax);
   dr3_default_limits(nh, nv, &nibop_ivmin, &nibop_ivmax, 
     nibop_ihmin, nibop_ihmax);
/* free any existing arrays (to be sure) */
   dr3_free_nibs();
/* allocate nib array */
   nib = (nib_ptr*)malloc((unsigned)(nh*nv*sizeof(nib_ptr)));
   nibop = (nib_ptr*)malloc((unsigned)(nh*nv*sizeof(nib_ptr)));
   nib_size = 4*nh*nv;
   nib_buffer = (nib_struct*)malloc((unsigned)(nib_size*sizeof(nib_struct)));
   if ((nib == NULL) || (nibop == NULL) || (nib_buffer == NULL)) {
      printf("error:  not enough memory for solid cuts\n");
      solid_mode = FALSE;
      surface_mode = TRUE;
#ifdef DEBUG_SOLID
   } else {
      printf("dr3_alloc_nibs:  #nib[op]=%d\n", (int)sizeof(nib_ptr));
      printf("dr3_alloc_nibs:  #nib_=%d\n", (int)sizeof(nib_struct));
      printf("dr3_alloc_nibs:  nib(%d)=%d\n", (nh*nv), (int)nib);
      printf("dr3_alloc_nibops:  nibop(%d)=%d\n", (nh*nv), (int)nibop);
      printf("dr3_alloc_nib_s:  nib_(%d)=%d\n",nib_size,(int)nib_buffer);
#endif
   }
}

/* dr3_add_nib - add nibop into nib array */
void dr3_add_nibs() {
   int i,k,n1,nn,op1[MAXOP],op[MAXOP],c1[MAXOP],c[MAXOP],ih,iv,nm,im;
   float d1[MAXOP],d[MAXOP];
   bool out;

#ifdef DEBUG_SOLID
   printf("dr3_add_nibs:  called\n");
#endif
/* extend initialisation of nib-array to include nibop-array */
   dr3_extend_limits(nh, nv, nibop_ivmin,nibop_ivmax,
     nibop_ihmin, nibop_ihmax, 
     &nib_ivmin, &nib_ivmax, nib_ihmin, nib_ihmax);
/* apply nibop to nib, then zap nibop limits*/
   for (iv=nibop_ivmin; iv <= nibop_ivmax; iv++) {
      k = nh*iv + nibop_ihmin[iv];
      for (ih=nibop_ihmin[iv]; ih <= nibop_ihmax[iv]; ih++) {
         if (nibop[k] != NULL) {
            n1 = 0;
            dr3_pop_nib(nib[k], &n1, d1, op1, c1);
#ifdef    DEBUG_SOLID
            printf("..[%d,%d],n1(nib)=%d", ih, iv, n1);
#endif
            dr3_pop_nib(nibop[k], &n1, d1, op1, c1);
#ifdef    DEBUG_SOLID
            printf(",n1(nib+nibop)=%d\n", n1);
#endif
            dr3_sort_nib(n1,d1,op1,c1);
#ifdef    DEBUG_SOLID
            printf("..(sorted):-\n");
            for (i=0; i < n1; i++) {
               printf("..i=%d,d1=%8.6f,op1=%d,c=%d\n",i,d1[i],op1[i],c1[i]);
            }
#endif
/*          new list is first +, last -, first + etc. */
            i = 0;
            nn = 0;
            out = TRUE;
            while (i < n1) {
               if (out) {
                  while(i < n1) {
                     if (op1[i] > 0) {
                        d[nn] = d1[i];
                        op[nn] = op1[i];
                        c[nn] = c1[i];
                        i++;
                        nn++;
                        out = FALSE;
                        break;
                     }
                     i++;
                  }
               } else {
                  nm = 0;
                  while (i < n1) {
                     if (op1[i] < 0) {
                        im = i;
                        nm++;
                     } else if (nm > 0) {
/*
                        i++;
*/
                        break;
                     }
                     i++;
                  }
                  if (nm > 0) {
                     d[nn] = d1[im];
                     op[nn] = op1[im];
                     c[nn] = c1[im];
                     nn++;
                     out = TRUE;
                  }
               }
            }
            if (nn < 1) {
#ifdef BUG
               nib[k] == NULL;
#else
               nib[k] = NULL;
#endif
            } else {
               if (op[nn-1] > 0) {
                  nn--;
/*
                  printf("add_nibs:  [%d,%d], op ends +; nn->%d\n", 
                    ih,iv,nn);
*/
               }
               if (nn < 1) {
                  nib[k] = NULL;
               } else {
                  dr3_push_nib(nn,d,op,c, &(nib[k]));
               }
            }
#ifdef    DEBUG_SOLID
            printf("..(out):-\n");
            for (i=0; i < nn; i++) {
               printf("..i=%d,d=%8.6f,op=%d,c=%d\n",i,d[i],op[i],c[i]);
            }
            printf("..nibk=%d\n", (int)nib[k]);
#endif
         }
         k++;
      }
      nibop_ihmin[iv] = nh;
      nibop_ihmax[iv] = (-1);
   }
   nibop_ivmin = nv;
   nibop_ivmax = (-1);
}

/* dr3_cut_nib - cut nibop into nib array */
void dr3_cut_nibs() {
   int i,k,n1,nn,op1[MAXOP],op[MAXOP],c1[MAXOP],c[MAXOP],ih,iv,
     n2,op2[MAXOP],c2[MAXOP],c1f,c2b,j;
   float d1[MAXOP],d[MAXOP],d2[MAXOP],d1f,d1b,d2f,d2b;
   bool pend;

#ifdef DEBUG_ETC
   printf("dr3_cut_nib:  called\n");
#endif
/* extend initialisation of nib-array to include nibop-array */
   dr3_extend_limits(nh, nv, nibop_ivmin,nibop_ivmax,
     nibop_ihmin, nibop_ihmax, 
     &nib_ivmin, &nib_ivmax, nib_ihmin, nib_ihmax);
/* apply nibop to nib */
   for (iv=nibop_ivmin; iv <= nibop_ivmax; iv++) {
      k = nh*iv + nibop_ihmin[iv];
      for (ih=nibop_ihmin[iv]; ih <= nibop_ihmax[iv]; ih++) {
         if ((nib[k] != NULL) && (nibop[k] != NULL)) {
            n1 = 0;
            dr3_pop_nib(nib[k], &n1, d1, op1, c1);
            n2 = 0;
            dr3_pop_nib(nibop[k], &n2, d2, op2, c2);
#ifdef    DEBUG
            printf("..[%d,%d],n1(nib),n2(nibop)=%d", ih, iv, n1, n2);
#endif
            dr3_sort_nib(n1,d1,op1,c1);
            dr3_sort_nib(n2,d2,op2,c2);
#ifdef    DEBUG
            printf("..(nib:sorted):-\n");
            for (i=0; i < n1; i++) {
               printf("..i=%d,d1=%8.6f,op1=%d,c1=%d\n",
                 i, d1[i], op1[i], c1[i]);
            }
            printf("..(nibop:sorted):-\n");
            for (i=0; i < n2; i++) {
               printf("..i=%d,d2=%8.6f,op2=%d,c2=%d\n",
                 i, d2[i], op2[i], c2[i]);
            }
#endif
            dr3_fix_nib(n2, d2, op2, c2, &n2);
#ifdef DEBUG_ETC
            printf("..(nibop:fixed):-\n");
            for (i=0; i < n2; i++) {
               printf("..i=%d,d2=%8.6f,op2=%d,c2=%d\n", i, d2[i],
                 op2[i], c2[i]);
            }
#endif
/*          For each (nib) solid section */
            nn = 0;
            i = 1;
            while (i < n1) {
/*             extract copy */
               c1f = c1[i-1];
               d1f = d1[i-1];
               d1b = d1[i];
#ifdef DEBUG_ETC
               printf("..solid,i=%d,c1f=%d,d1f=%8.6f,d1b=%8.6f\n",
                 i, c1f, d1f, d1b);
#endif
               pend = TRUE;
/*             For each (nibop) cut section */
               j = 1;
               while (j < n2) {
                  c2b = c2[j];
                  d2f = d2[j-1];
                  d2b = d2[j];
#ifdef DEBUG_ETC
                  printf("..cut,j=%d,d2f=%8.6f,d2b=%8.6f,c2b=%d\n",
                    j, d2f, d2b, c2b);
#endif
/*                Cut totally infront:  next */
                  if (d2b < d1f) {
                     j = j + 2;
/*                Cut totally behind:  break */
                  } else if (d2f > d1b) {
                     break;
/*                Cuts front of solid */
                  } else if (d2f < d1f) {
/*                   Cut front off solid */
                     if (d2b <= d1b) {
                        d1f = d2b;
                        c1f = c2b;
                        j = j + 2;
/*                   Cut across entire solid */
                     } else {
                        pend = FALSE;
                        break;
                     }
/*                Cut rear off solid */
                  } else if (d2b >= d1b) {
                     d1b = d2f;
                     j = j + 2;
/*                Cut hole in solid */
                  } else {
/*                   Push (c1f,d1f,d2f) */
                     d[nn] = d1f;
                     op[nn] = 1;
                     c[nn] = c1f;
                     nn++;
                     d[nn] = d2f;
                     op[nn] = (-1);
                     c[nn] = 0;
                     nn++;
                     d1f = d2b;
                     c1f = c2b;
                     j = j + 2;
                  }
               }
/*             Flush any pending solid */
               if (pend) {
/*                Push (c1f,d1f,d1b) */
                  d[nn] = d1f;
                  op[nn] = 1;
                  c[nn] = c1f;
                  nn++;
                  d[nn] = d1b;
                  op[nn] = (-1);
                  c[nn] = 0;
                  nn++;
               }
               i = i + 2;
            }
            if (nn < 1) {
               nib[k] = NULL;
            } else {
               dr3_push_nib(nn,d,op,c, &(nib[k]));
            }
#ifdef    DEBUG
            printf("..(out):-\n");
            for (i=0; i < nn; i++) {
               printf("..i=%d,d=%8.6f,op=%d,c=%d\n", i, d[i], op[i], c[i]);
            }
            printf("..nibk=%d\n", (int)nib[k]);
#endif
         }
         k++;
      }
      nibop_ihmin[iv] = nh;
      nibop_ihmax[iv] = (-1);
   }
   nibop_ivmin = nv;
   nibop_ivmax = (-1);
}

void dr3_fix_nib(n, d, op, c, nnn)
   int n,op[],c[],*nnn;
   float d[];
{
   int nn,i,np,ip;
   bool out;

/* new list is last +, first -, last + etc. */
   nn = 0;
   i = 0;
   out = FALSE;
   while (i < n) {
      if (out) {
         np = 0;
         while (i < n) {
            if (op[i] > 0) {
               np++;
               ip = i;
            } else if (np > 0) {
               break;
            }
            i++;
         }
         if (np > 0) {
            d[nn] = d[ip];
            op[nn] = op[ip];
            c[nn] = c[ip];
            nn++;
            out = FALSE;
         }
      } else {
         while (i < n) {
            if (op[i] < 0) {
               d[nn] = d[i];
               op[nn] = op[i];
               c[nn] = c[i];
               nn++;
               i++;
               out = TRUE;
               break;
            }
            i++;
         }
      }
   }
   if (nn > 0) {
      if (op[nn-1] < 0) {
         nn--;
      }
   }
   *nnn = nn;
}

/* dr3_render_nibs - transfer nibs to bitmap */
void dr3_render_nibs() {
   int i,iprev,k,nn,op[MAXOP],c[MAXOP],ih,iv,kk;
   float d[MAXOP],light;

#ifdef DEBUG_SOLID
   printf("dr3_render_nibs:  called\n");
#endif
   for (iv=nib_ivmin; iv <= nib_ivmax; iv++) {
      k = nh*iv + nib_ihmin[iv];
      for (ih=nib_ihmin[iv]; ih <= nib_ihmax[iv]; ih++) {
         if (nib[k] != NULL) {
/*          Pop nibs into local array */
            nn = 0;
            dr3_pop_nib(nib[k], &nn, d, op, c);
#ifdef    DEBUG_SOLID
            printf("..[%d,%d],nn=%d\n", ih, iv, nn);
            for (i=0; i < nn; i++) {
               printf("..i=%d,d=%5.2f,op=%d,c=%d\n",i,d[i],op[i],c[i]);
            }
#endif
/*          Find first front-facing nib and maybe use it */
            iprev = (-1);
            for (i=0; i < nn; i++) {
               if (op[i] < 0) {
               } else if (iprev < 0) {
                  iprev = i;
               } else if (d[i] < d[iprev]) {
                  iprev = i;
               }
            }
/*          z-buffer it onto bitmap */
#ifdef    DEBUG_SOLID
            printf("..iprev=%d\n", iprev);
#endif
            if (iprev < 0) {
            } else if (d[iprev] < dim[k]) {
               dim[k] = d[iprev];
               if (colour_bitmap) {
                  cim[k] = c[iprev];
               } else if (rgb_bitmap) {
                  kk = k*3;
                  light = (float)c[iprev]*0.01;
                  rgbim[kk] = light*material_red;
                  rgbim[kk+1] = light*material_green;
                  rgbim[kk+2] = light*material_blue;
               }
            }
         }
         k++;
      }
      nib_ihmin[iv] = nh;
      nib_ihmax[iv] = (-1);
   }
   nib_ivmin = nv;
   nib_ivmax = (-1);
}

/* dr3_sort_nib - sort nib array into order */
void dr3_sort_nib(n, d, op, c) 
   int n,op[],c[];
   float d[];
{
   int k,i,opk,ck;
   float dk;

   for (k=0; k < (n-1); k++) {
      for (i=k+1; i < n; i++) {
         if (d[i] < d[k]) {
            dk = d[k];
            opk = op[k];
            ck = c[k];
            d[k] = d[i];
            op[k] = op[i];
            c[k] = c[i];
            d[i] = dk;
            op[i] = opk;
            c[i] = ck;
         }
      }
   }
}

/* dr3_pop_nib - pop nibs from structure list */
void dr3_pop_nib(nib, n, d, op, c) 
   nib_struct *nib;
   int *n,op[],c[];
   float d[];
{
   nib_struct *nibc;

   nibc = nib;
   while (nibc != NULL) {
      d[*n] = nibc->d;
      c[*n] = nibc->c;
      op[*n] = nibc->op;
      *n = (*n) + 1;
      nibc = nibc->prev;
   }
}

/* dr3_push_nib - push nibs to structure list */
void dr3_push_nib(n, d, op, c, nib)
   nib_struct *(*nib);
   int n,op[],c[];
   float d[];
{
   nib_struct *nibc,*prev;
   int i;

   nibc = NULL;
   for (i=(n-1); i >= 0; i--) {
      if (nib_buffer_overflow) {
         break;
      } else if (nib_count < nib_size) {
         prev = nibc;
         nibc = &(nib_buffer[nib_count++]);
         nibc->prev = prev;
         nibc->d = d[i];
         nibc->op = op[i];
         nibc->c = c[i];
      } else {
         printf("error:  nib buffer overflow\n");
         nib_buffer_overflow = TRUE;
         break;
      }
   }
   *nib = nibc;
}

/* dr3_free_nibs - free temporary (curremt,op) nib arrays */
void dr3_free_nibs() {
   if (nib != NULL) {
      free(nib);
      nib = NULL;
   }
   if (nibop != NULL) {
      free(nibop);
      nibop = NULL;
   }
   if (nib_buffer != NULL) {
      free(nib_buffer);
      nib_buffer = NULL;
   }
}

/* POLYGON PRIMITIVES */

/* dr3_begin_polygon - start polygon construction */
int dr3_begin_polygon() {
   normal_defined = FALSE;
   scalar_defined = FALSE;
/* rgbz-mode:  pass it on */
   if (rgbz_bitmap) {
      rgbz_begin_polygon();
   }
   nbuf = 0;
   first=TRUE;
/* vector-mode means we deal in polylines */
   if (vector_bitmap || wire_frame) {
      if (polyline_buffer_fault) {
      } else if (ns >= MAXSEG) {
         printf("error:  3D polyline segment buffer overflow\n");
         polyline_buffer_fault = TRUE;
      } else {
         ip0[ns] = np;
      }
   }
   return 0;
}

/* dr3_polygon_vertex - add a vertex to current polygon */
void dr3_polygon_vertex(x,y,z)
   float x,y,z;
{
   float xo[MAXD4],xv[MAXD4],h,v,d,light,rgb[3];

/* plant object coords in vector */
   xo[0] = x;
   xo[1] = y;
   xo[2] = z;
   xo[3] = 1.0;
/* [x,y,z]    to [X,Y,Z] */
   dr3_vmult(vm, xo, xv);
/* rgbz-mode:  pass it on */
   if (rgbz_bitmap) {
      dr3_shade_rgbz(xv, rgb);
      drp_set_rgb(rgb[0], rgb[1], rgb[2]);
      rgbz_polygon_vertex(xv[0], xv[1], xv[2]);
      hbuf[nbuf] = xv[0];
      vbuf[nbuf] = xv[1];
      dbuf[nbuf] = xv[2];
      nbuf++;
      if (first) {
         hi = xv[0];
         vi = xv[1];
         di = xv[2];
         first = FALSE;
      }
/* otherwise: */
   } else {
/*    [X,Y,Z] to [h,v,d] */
      dr3_vtod(xv, &h, &v, &d);
/*    append to outline polyline buffer */
      if (surface_mode && (vector_bitmap || wire_frame)) {
         if (polyline_buffer_fault) {
         } else if (np < MAXPOL) {
            hp[np] = h;
            vp[np] = v;
            dp[np] = d;
            np++;
         } else {
            printf("error:  3D polyline data buffer overflow\n");
            polyline_buffer_fault = TRUE;
         }
      }
/*    append [h,v,d] to buffer */
      hbuf[nbuf] = h;
      vbuf[nbuf] = v;
      dbuf[nbuf] = d;
/*    append shade-model data to buffer */
      dr3_light(xv, normal, &light);
      if (surface_mode) {
         cbuf[nbuf][0] = light;
      } else if (normal[2] > 0.0) {
         cbuf[nbuf][0] = light;
      } else {
         cbuf[nbuf][0] = (-light);
      }
      if (shade_index == SHADE_LD) {
         if (scalar_defined) {
            cbuf[nbuf][1] = scalar;
         } else {
            cbuf[nbuf][1] = 0.0;
         }
      }
#ifdef DEBUG_SOLID
      printf("..polygon_vertex:  normal=[%f,%f,%f]\n", normal[0],
        normal[1], normal[2]);
      printf("..light=%f,cbuf=%f\n", light, cbuf[nbuf][0]);
#endif
/*    increment polygon vertex counter */
      nbuf = nbuf + 1;
      if (first) {
         hi = h;
         vi = v;
         di = d;
         first = FALSE;
      }
   }
}

/* dr3_polygon_normal - define surface normal for subsequent vertices */
void dr3_polygon_normal(nx,ny,nz)
   float nx,ny,nz;
{
   float a[MAXD4];

   normal_defined = TRUE;
/* compose world coordinate 3-vec (factor by sign) */
   a[0] = solid_op_sign*nx;
   a[1] = solid_op_sign*ny;
   a[2] = solid_op_sign*nz;
/* convert it to viewing coordinates [X,Y,Z] */
   dr3_vmult3(vm, a, normal);
#ifdef DEBUG_SOLID
   printf("dr3_polygon_normal:  nx=%f,ny=%f,nz=%f\n", nx, ny, nz);
   printf("..solid_op_sign=%f,normal=[%f,%f,%f]\n", solid_op_sign,
     normal[0],normal[1],normal[2]);
#endif
}

/* dr3_polygon_scalar - define surface scalar for subsequent vertices */
void dr3_polygon_scalar(s)
   float s;
{
   scalar_defined = TRUE;
   scalar = s;
}

/* dr3_end_polygon - complete polygon construction */
void dr3_end_polygon() {
   int i;

/* rgbz-mode:  pass it on */
   if (rgbz_bitmap) {
      rgbz_end_polygon();
      drp_set_rgb(0.0, 0.0, 0.0);
      hbuf[nbuf] = hi;
      vbuf[nbuf] = vi;
      dbuf[nbuf] = di;
      nbuf++;
      if (wire_frame) {
         rgbz_begin_polyline();
         for (i=0; i < nbuf; i++) {
            rgbz_polyline_vertex(hbuf[i], vbuf[i], dbuf[i]);
         }
         rgbz_end_polyline();
      }
/* otherwise:  private z-buffer */
   } else {
/*    construct patch */
      dr3_study_polygon();
/*    render patch */
      if (solid_mode) {
         dr3_render_solid();
      } else {
         dr3_render_surface();
      }
/*    add perimeter to polyline buffer */
      if (surface_mode && (vector_bitmap || wire_frame)) {
         if (polyline_buffer_fault) {
         } else if (ns < MAXSEG) {
            if (np >= MAXPOL) {
               printf("error:  3D polyline data buffer overflow\n");
               polyline_buffer_fault = TRUE;
            } else {
               hp[np] = hi;
               vp[np] = vi;
               dp[np] = di;
               np++;
            }
            nps[ns] = np;
            ns++;
         } else {
            printf("error:  3D polyline segment buffer overflow\n");
            polyline_buffer_fault = TRUE;
         }
      }
   }
   normal_defined = FALSE;
   scalar_defined = FALSE;
}

/* POLYLINE PRIMITIVES */

/* dr3_zap_polylines - begin set of polylines to be "z-buffered" */
void dr3_zap_polylines() {
   if (!rgbz_bitmap) {
      np = 0;
      ns = 0;
      polyline_buffer_fault = FALSE;
   }
}

/* dr3_display_polylines - end set of polylines */
void dr3_display_polylines() {
   int is,i;

   for (is=0; is < ns; is++) {
      drp_bp();
      for (i=ip0[is]; i < nps[is]; i++) {
         drp_ap(hp[i], vp[i]);
      }
      drp_ep();
   }
}

/* dr3_hide_polylines - end set of polylines, z-buffer them */
void dr3_hide_polylines() {
   float h1,h2,v1,v2,dh,dv,dd,dl,dddl,d1,d2,h,v,d,dix;
   int is,i,n,j,ih,iv,lh,lv,k, ihl,ihr,ivl,ivr;
   bool empty,clip1,clip2,in;
   dddl = 0.0;

   for (is=0; is < ns; is++) {
      in = FALSE;
      for (i=(ip0[is]+1); i < nps[is]; i++) {
/*       clip segment */
         h1 = hp[i-1];
         v1 = vp[i-1];
         h2 = hp[i];
         v2 = vp[i];
         drp_clip(&h1, &v1, &clip1, &h2, &v2, &clip2, &empty);
/*       empty segment */
         if (empty) {
            if (in) {
               drp_ep();
               in = FALSE;
            }
/*       non-empty segment */
         } else {
/*          any clipping means we must interpolate */
            if (clip1 || clip2) {
               dh = (hp[i]-hp[i-1]);
               dv = (vp[i]-vp[i-1]);
               dd = (dp[i]-dp[i-1]);
               dl = dh*dh + dv*dv;
               if (dl > 0.0) {
                  dddl = dd/sqrt(dl);
               } else {
                  dddl = 0.0;
               }
            }
/*          cliped at start? */
            if (clip1) {
               if (in) {
                  drp_ep();
                  in = FALSE;
               }
               dh = h1 - hp[i-1];
               dv = v1 - vp[i-1];
               d1 = sqrt(dh*dh+dv*dv)*dddl + dp[i-1];
            } else {
               d1 = dp[i-1];
            }
/*          clipped at end? */
            if (clip2) {
               dh = h2 - hp[i-1];
               dv = v2 - vp[i-1];
               d2 = sqrt(dh*dh+dv*dv)*dddl + dp[i-1];
            } else {
               d2 = dp[i];
            }
/*          breakup line into a grid */
            lh = (int)(fabs(h2-h1)*dhidh) + 1;
            lv = (int)(fabs(v2-v1)*dvidv) + 1;
            n = 2;
            if (lh > n) {
               n = lh;
            }
            if (lv > n) {
               n = lv;
            }
            dh = (h2-h1)/(float)(n-1);
            dv = (v2-v1)/(float)(n-1);
            dd = (d2-d1)/(float)(n-1);
/*          work through the line */
            h = h1;
            v = v1;
            d = d1;
            for (j=0; j < n; j++) {
               if ((d < 0.0) || (d > 1.0)) {
                  if (in) {
                     drp_ep();
                     in = FALSE;
                  }
               } else {
                  ih = (int)(dhidh*h+hi0);
                  if (ih < 1) {
                    ihl = 0;
                    ihr = 1;
                  } else if (ih > (nh-1)) {
                    ihl = (nh-2);
                    ihr = (nh-1);
                  } else {
                     ihl = ih-1;
                     ihr = ih+1;
                  }
                  iv = (int)(dvidv*v+vi0);
                  if (iv < 1) {
                     ivl = 0;
                     ivr = 1;
                  } else if (iv > (nv-1)) {
                     ivl = (nv-2);
                     ivr = (nv-1);
                  } else {
                     ivl = iv-1;
                     ivr = iv+1;
                  }
                  dix = 0.0;
                  for (iv=ivl; iv <=ivr; iv++) {
                     for (ih=ihl; ih <=ihr; ih++) {
                        k = nh*iv + ih;
                        if (dim[k] > dix) {
                           dix = dim[k];
                        }
                     }
                  }
                  if (d > dix) {
                     if (in) {
                        drp_ep();
                        in = FALSE;
                     }
                  } else if (in) {
                     drp_ap(h, v);
                  } else {
                     drp_bp();
                     in = TRUE;
                     drp_ap(h, v);
                  }
               }
               h = h + dh;
               v = v + dv;
               d = d + dd;
            }
         }
      }
      if (in) {
         drp_ep();
      }
   }
}

/* dr3_begin_polyline - begin polyline construction */
void dr3_begin_polyline() {
   if (rgbz_bitmap) {
      rgbz_begin_polyline();
   } else if (polyline_buffer_fault) {
   } else if (ns >= MAXSEG) {
      printf("error:  3D polyline segment buffer overflow\n");
      polyline_buffer_fault = TRUE;
   } else {
      ip0[ns] = np;
   }
}

/* dr3_polyline_vertex - add a vertex to current polyline */
void dr3_polyline_vertex(x,y,z)
   float x,y,z;
{
   float h,v,d, xo[MAXD4],xv[MAXD4];

   if (rgbz_bitmap) {
      xo[0] = x;
      xo[1] = y;
      xo[2] = z;
      xo[3] = 1.0;
      dr3_wtov(xo, xv);
      rgbz_polyline_vertex(xv[0], xv[1], xv[2]);
   } else if (polyline_buffer_fault) {
   } else if (np >= MAXPOL) {
      printf("error:  3D polyline data buffer overflow\n");
      polyline_buffer_fault = TRUE;
   } else {
      dr3_wtod(x, y, z, &h, &v, &d);
      hp[np] = h;
      vp[np] = v;
      dp[np] = d;
      np++;
   }
}

/* dr3_end_polyline - complete polyline construction */
void dr3_end_polyline() {
   if (rgbz_bitmap) {
      rgbz_end_polyline();
   } else if (polyline_buffer_fault) {
   } else if (ns < MAXSEG) {
      nps[ns] = np;
      ns++;
   }
}

/* TRANSLATION-ROTATION PRIMITIVES */

/* dr3_push_matrix - push current matrix onto stack */
void dr3_push_matrix() {
   if (nmatrix >= MAXMATRIX) {
      printf("error:  matrix stack full\n");
   } else {
      dr3_copy(vm, vms[nmatrix]);
      dr3_copy(vw, vws[nmatrix]);
      nmatrix++;
   }
}

/* dr3_pop_matrix - pop matrix from stack */
void dr3_pop_matrix() {
   if (nmatrix == 0) {
      printf("error:  matrix stack empty\n");
   } else {
      nmatrix--;
      dr3_copy(vms[nmatrix], vm);
      dr3_copy(vws[nmatrix], vw);
   }
}

/* dr3_get_matrix - get current view-matrix */
void dr3_get_matrix(a)
   float a[MAXD4][MAXD4];
{
   dr3_copy(vm, a);
}

/* dr3_put_matrix - put new view-matrix */
void dr3_put_matrix(a)
   float a[MAXD4][MAXD4];
{
   dr3_copy(a, vm);
}

/* dr3_zap_matrix - zap matrix stack */
void dr3_zap_matrix() {
   nmatrix = 0;
   if (!noview) {
      dr3_copy(view_matrix, vm);
      dr3_copy(view_world, vw);
   }
}

/* dr3_rotran_matrix - perform complete rot-tran on matrix */
void dr3_rotran_matrix(x, y, z, azim, altd, twist)
   float x,y,z, azim,altd,twist;
{
   dr3_tran_matrix(x, y, z);
/*
   dr3_rot_matrix(2, twist);
*/
   dr3_rot_matrix(2, -azim);
   dr3_rot_matrix(1, -altd);
}

/* dr3_rot_matrix -  rotate by "angled" (degrees) about "id" */
void dr3_rot_matrix(id, angled)
   int id;
   float angled;
{
   float r[MAXD4][MAXD4],b[MAXD4][MAXD4],cosa,sina;
   int row,col;

/* construct r-matrix */
   for (row=0; row < MAXD4; row++) {
      for (col=0; col < MAXD4; col++) {
         r[row][col] = 0.0;
      }
      r[row][row] = 1.0;
   }
   cosa = cos(angled*DEGRAD);
   sina = sin(angled*DEGRAD);
   if (id == 0) {
      r[1][1] = cosa;
      r[1][2] = sina;
      r[2][1] = -sina;
      r[2][2] = cosa;
   } else if (id == 1) {
      r[0][0] = cosa;
      r[0][2] = -sina;
      r[2][0] = sina;
      r[2][2] = cosa;
   } else {
      r[0][0] = cosa;
      r[0][1] = sina;
      r[1][0] = -sina;
      r[1][1] = cosa;
   }
/* apply to VM-matrix,result - B-matrix */
   dr3_mult(vm, r, b);
/* copy back into VM */
   dr3_copy(b, vm);
   dr3_mult(vw, r, b);
   dr3_copy(b, vw);
}

/* dr3_tran_matrix - modify "a" by translation "[tx,ty,tz]" */
void dr3_tran_matrix(tx, ty, tz)
   float tx,ty,tz;
{
   float t[MAXD4][MAXD4],b[MAXD4][MAXD4];
   int row,col;

/* construct T-matrix */
   for (row=0; row < MAXD4; row++) {
      for (col=0; col < MAXD4; col++)  {
         t[row][col] = 0.0;
      }
      t[row][row] = 1.0;
   }
   col = MAXD;
   t[0][col] = tx;
   t[1][col] = ty;
   t[2][col] = tz;
/* apply to VM-matrix */
   dr3_mult(vm, t, b);
/* copy back into VM-matrix */
   dr3_copy(b, vm);
   dr3_mult(vw, t, b);
   dr3_copy(b, vw);
}

/* PRIVATE FUNCTIONS */

/* dr3_study_polygon - deduce image share for polygon */
void dr3_study_polygon()
{
   float denom,hmu,dmu,cmu[MAXCD];
   int i,iprev,iv,ivprev,ih,ihprev,jv,jh,jvmin,jvmax,id;

/* set safe limits for each line */
   for (iv=0; iv < nv; iv++) {
      ihmin[iv] = nh;
      ihmax[iv] = (-1);
   }
/* set safe line limits */
   ivmin = nv;
   ivmax = (-1);
/* map end point onto array */
   iprev = (nbuf-1);
   ivprev = (int)(dvidv*vbuf[iprev]+vi0);
   ihprev = (int)(dhidh*hbuf[iprev]+hi0);
/* work through segments */
   for (i=0; i < nbuf; i++) {
/*    map current point onto array */
      iv = (int)(dvidv*vbuf[i]+vi0);
      ih = (int)(dhidh*hbuf[i]+hi0);
/*    only consider inclined segments */
      if (iv != ivprev) {
/*       find line range (if any) */
         if (iv < ivprev) {
            jvmin = iv;
            jvmax = ivprev;
         } else {
            jvmin = ivprev;
            jvmax = iv;
         }
/*       clip to array */
         if (jvmin < 0) {
            jvmin = 0;
         }
         if (jvmax > (nv-1)) {
            jvmax = (nv-1);
         }
/*       avoid degenerate or other bad cases */
         if (jvmin < jvmax) {
/*          affect vertical extents */
            if (jvmin < ivmin) {
               ivmin = jvmin;
            }
            if (jvmax > ivmax) {
               ivmax = jvmax;
            }
/*          interpolation factors */
            denom = 1.0/(float)(iv-ivprev);
            hmu = (float)(ih-ihprev)*denom;
            dmu = (dbuf[i]-dbuf[iprev])*denom;
            for (id=0; id < ncd; id++) {
               cmu[id] = (cbuf[i][id]-cbuf[iprev][id])*denom;
            }
/*          work through lines in vertical range */
            for (jv=jvmin; jv <= jvmax; jv++) {
/*             find "ih" by interpolation */
               jh = hmu*(jv-ivprev) + ihprev;
/*             affect both extents */
               if (ihmin[jv] > ihmax[jv]) {
                  ihmin[jv] = jh;
                  dmin[jv] = dmu*(jv-ivprev) + dbuf[iprev];
                  for (id=0; id < ncd; id++) {
                     cmin[jv][id] = cmu[id]*(jv-ivprev)+cbuf[iprev][id];
                     cmax[jv][id] = cmin[jv][id];
                  }
                  ihmax[jv] = jh;
                  dmax[jv] = dmin[jv];
/*             affect minimum h-extent only */
               } else if (jh < ihmin[jv]) {
                  ihmin[jv] = jh;
                  dmin[jv] = dmu*(jv-ivprev) + dbuf[iprev];
                  for (id=0; id < ncd; id++) {
                     cmin[jv][id] = cmu[id]*(jv-ivprev)+cbuf[iprev][id];
                  }
/*             affect maximum h-extent only */
               } else if (jh > ihmax[jv]) {
                  ihmax[jv] = jh;
                  dmax[jv] = dmu*(jv-ivprev) + dbuf[iprev];
                  for (id=0; id < ncd; id++) {
                     cmax[jv][id] = cmu[id]*(jv-ivprev)+cbuf[iprev][id];
                  }
               }
            }
         }
      }
/*    update previous povoid */
      iprev = i;
      ivprev = iv;
      ihprev = ih;
   }
}

/* SURFACES */

void dr3_begin_surface() {
/* Select surface mode */
   solid_mode = FALSE;
   surface_mode = TRUE;
/* select-check material index (assume LM) */
   if (shade_index == SHADE_LM) {
      if (material_index < nc[1]) {
         matin = material_index;
      } else {
         matin = 0;
      }
   }
}

void dr3_end_surface() {
/* update screen? */
   if (colour_bitmap) {
      drp_colour_bitmap(ih0, iv0, nh, nv, cim);
   } else if (rgb_bitmap) {
#ifndef NO_UPDATE_RGB
      drp_rgb_bitmap(ih0, iv0, nh, nv, rgbim);
#endif
   }
}

/* dr3_render_surface - render surface patch onto current image */
void dr3_render_surface() {
   float denom,dmu,d0,cmu[MAXCD],c0[MAXCD],d,c[MAXCD],h,v;
   int ih,iv,ihmn,ihmx,k,id;

/* only for some vertical extent */
   if (ivmin <= ivmax) {
/*    random  fix */
/*    for each (already clipped) line in vertical */
      for (iv=ivmin; iv <= ivmax; iv++) {
/*       only if there is some (unclipped) horizontal extent */
         if (ihmin[iv] < ihmax[iv]) {
/*          clip to array */
            if (ihmin[iv] < 0) {
               ihmn = 0;
            } else {
               ihmn = ihmin[iv];
            }
            if (ihmax[iv] > (nh-1)) {
               ihmx = (nh-1);
            } else {
               ihmx = ihmax[iv];
            }
/*          only for some overlap */
            if (ihmn <= ihmx) {
/*             interpolates for depth and colour */
               denom = 1.0/(float)(ihmax[iv]-ihmin[iv]);
               dmu = (dmax[iv]-dmin[iv])*denom;
               d0 = dmin[iv];
               for (id=0; id < ncd; id++) {
                  cmu[id] = (cmax[iv][id]-cmin[iv][id])*denom;
                  c0[id] = cmin[iv][id];
               }
/*             cut across valid pixels */
               k = nh*iv + ihmn;
               v = ((float)iv-vi0)/dvidv;
               for (ih=ihmn; ih <= ihmx; ih++) {
                  d = dmu*(float)(ih-ihmin[iv]) + d0;
                  if ((d < 0.0) || (d > 1.0)) {
                  } else if (d < dim[k]) {
                     h = ((float)ih-hi0)/dhidh;
                     for (id=0; id < ncd; id++) {
                        c[id] = cmu[id]*(float)(ih-ihmin[iv]) + c0[id];
                     }
                     if (colour_bitmap) {
                        dr3_shade_cm(h,v,d, c, &(dim[k]), &(cim[k]));
                     } else if (rgb_bitmap) {
                        dr3_shade_rgb(h,v,d,c,&(dim[k]),&(rgbim[k*3]));
                     } else {
                        dim[k] = d;
                     }
                  }
                  k++;
               }
            }
         }
      }
   }
}

/* dr3_render_solid - render solid patch onto current image */
void dr3_render_solid() {
   float denom,dmu,d0,cmu[MAXCD],c0[MAXCD],d,c[MAXCD],light;
   int ih,iv,ihmn,ihmx,k,id,ic;
   nib_struct *newc,*prevc;

/* only for some vertical extent */
   if (ivmin <= ivmax) {
#ifdef DEBUG_SOLID
      printf("dr3_render_solid:  ivmin=%d,ivmax=%d\n", ivmin, ivmax);
#endif
/*    Ensure that the nibop region includes this vertical extent */
      if (nibop_ivmin > nibop_ivmax) {
#ifdef DEBUG_SOLID
         printf("..(new)nibop\n");
#endif
         for (iv=ivmin; iv <= ivmax; iv++) {
            nibop_ihmin[iv] = nh;
            nibop_ihmax[iv] = (-1);
         }
         nibop_ivmin = ivmin;
         nibop_ivmax = ivmax;
      } else {
#ifdef DEBUG_SOLID
         printf("..(adjust)nibop\n");
#endif
         if (ivmin < nibop_ivmin) {
#ifdef DEBUG_SOLID
            printf("..reduce nibop_ivmin\n");
#endif
            for (iv=ivmin; iv < nibop_ivmin; iv++) {
               nibop_ihmin[iv] = nh;
               nibop_ihmax[iv] = (-1);
            }
            nibop_ivmin = ivmin;
         }
         if (nibop_ivmax < ivmax) {
#ifdef DEBUG_SOLID
            printf("..increase nibop_ivmax\n");
#endif
            for (iv=ivmax; iv > nibop_ivmax; iv--) {
               nibop_ihmin[iv] = nh;
               nibop_ihmax[iv] = (-1);
            }
            nibop_ivmax = ivmax;
         }
      }
/*    for each (already clipped) line in vertical */
      for (iv=ivmin; iv <= ivmax; iv++) {
/*       only if there is some (unclipped) horizontal extent */
         if (ihmin[iv] < ihmax[iv]) {
/*          clip to array */
            if (ihmin[iv] < 0) {
               ihmn = 0;
            } else {
               ihmn = ihmin[iv];
            }
            if (ihmax[iv] > (nh-1)) {
               ihmx = (nh-1);
            } else {
               ihmx = ihmax[iv];
            }
/*          only for some overlap */
            if (ihmn <= ihmx) {
/*             Ensure that the nibop region includes this */
               if (nibop_ihmin[iv] > nibop_ihmax[iv]) {
                  k = nh*iv + ihmn;
                  for (ih=ihmn; ih <= ihmx; ih++) {
                     nibop[k++] = NULL;
                  }
                  nibop_ihmin[iv] = ihmn;
                  nibop_ihmax[iv] = ihmx;
               } else {
                  if (ihmn < nibop_ihmin[iv]) {
                     k = nh*iv + ihmn;
                     for (ih=ihmn; ih < nibop_ihmin[iv]; ih++) {
                        nibop[k++] = NULL;
                     }
                     nibop_ihmin[iv] = ihmn;
                  }
                  if (ihmx > nibop_ihmax[iv]) {
                     k = nh*iv + ihmx;
                     for (ih=ihmx; ih > nibop_ihmax[iv]; ih--) {
                        nibop[k--] = NULL;
                     }
                     nibop_ihmax[iv] = ihmx;
                  }
               }
/*             interpolates for depth and colour */
               denom = 1.0/(float)(ihmax[iv]-ihmin[iv]);
               dmu = (dmax[iv]-dmin[iv])*denom;
               d0 = dmin[iv];
               for (id=0; id < ncd; id++) {
                  cmu[id] = (cmax[iv][id]-cmin[iv][id])*denom;
                  c0[id] = cmin[iv][id];
               }
/*             cut across valid pixels */
               k = nh*iv + ihmn;
/*
               v = ((float)iv-vi0)/dvidv;
*/
               for (ih=ihmn; ih <= ihmx; ih++) {
                  d = dmu*(float)(ih-ihmin[iv]) + d0;
/*
                  h = ((float)ih-hi0)/dhidh;
*/
                  for (id=0; id < ncd; id++) {
                     c[id] = cmu[id]*(float)(ih-ihmin[iv]) + c0[id];
                  }
/*                ColourMap mode */
                  if (colour_bitmap | rgb_bitmap) {
/*                   Append nib to array */
                     if (nib_buffer_overflow) {
                     } else if (nib_count < nib_size) {
                        newc = &nib_buffer[nib_count++];
                        prevc = nibop[k];
                        if (prevc == NULL) {
                           nibop[k] = newc;
                        } else {
                           nibop[k] = newc;
                        }
                        newc->prev = prevc;
                        if (c[0] < 0.0) {
                           c[0] = (-c[0]);
                           newc->op = (-1);
                        } else {
                           newc->op = 1;
                        }
                        newc->d = d;
                        if (colour_bitmap) {
                           if (random_shade) {
                              light = randaf[rand_count++] + c[0]*
                                shade_didc[0];
                              if (rand_count == MAXRAND) {
                                 rand_count = 0;
                              }
                           } else {
                              light = c[0]*shade_didc[0];
                           }
                           ic = (int)light;
                           if (ic < 0) {
                              ic = 0;
                           } else if (ic > (nc[0]-2)) {
                              ic = nc[0]-2;
                           }
                           newc->c = ((ic+1)+nc[0]*matin)+shade_offset;
                        } else {
                           newc->c = (int)(c[0]*100.0);
                        }
                     } else {
                        printf("error:  nib buffer overflow\n");
                        nib_buffer_overflow = TRUE;
                     }
                  } else {
                     dim[k] = d;
                  }
                  k++;
               }
            }
         }
      }
   }
}

/* dr3_shade_cm - construct colours for povoid */
void dr3_shade_cm(h, v, d, c, dimk, cimk)
   float h,v,d,c[],*dimk;
   int *cimk;
{
   float light,s;
   int ic,cam[MAXCD],ilight,iother;
   //bool found;

   switch (shade_index) {
      case SHADE_LM:
         *dimk = d;
         light = c[0]*shade_didc[0];
         if (random_shade) {
/*
            light = light + ranf*((float)rand()*divran - 0.5);
*/
            light = light + randaf[rand_count++];
            if (rand_count == MAXRAND) {
               rand_count = 0;
            }
         }
         ic = (int)light;
         if (ic < 0) {
            ic = 0;
         } else if (ic > (nc[0]-2)) {
            ic = nc[0]-2;
         }
         cam[0] = ic + 1;
         if (material_index < nc[1]) {
            ic = material_index;
         } else {
            ic = 0;
         }
         cam[1] = ic;
         *cimk = (cam[0] + nc[0]*cam[1]) + shade_offset;
         break;
      case SHADE_LD:
         if (hard_surface) {
            *dimk = d;
         }
         if (lucid_surface) {
            ilight = (*cimk - shade_offset);
            iother = ilight/nc[0];
            cam[0] = ilight - iother*nc[0];
         } else {
            light = c[0]*shade_didc[0];
            if (random_shade) {
/*
               light = light + ranf*((float)rand()*divran - 0.5);
*/
               light = light + randaf[rand_count++];
               if (rand_count == MAXRAND) {
                  rand_count = 0;
               }
            }
            ic = (int)light;
            if (ic < 0) {
               ic = 0;
            } else if (ic > (nc[0]-2)) {
               ic = nc[0]-2;
            }
            cam[0] = ic + 1;
         }
/*
         dr3_dtov(h,v,d, xv);
         dr3_vtow(xv, xw);
*/
         if (scalar_defined) {
            s = c[1]*shade_didc[1];
            if (random_colour) {
               s = s + randaf[rand_count++];
               if (rand_count == MAXRAND) {
                  rand_count = 0;
               }
            }
            ic = (int)s;
            if (ic < 0) {
               ic = 0;
            } else if (ic > (nc[1]-2)) {
               ic = nc[1]-2;
            }
            ic = ic+1;
            cam[1] = ic;
         } else {
            cam[1] = 0;
         }
         *cimk = (cam[0] + nc[0]*cam[1]) + shade_offset;
         break;
      case SHADE_L:
         *dimk = d;
         light = c[0]*shade_didc[0];
         if (random_shade) {
/*
            light = light + ranf*((float)rand()*divran - 0.5);
*/
            light = light + randaf[rand_count++];
            if (rand_count == MAXRAND) {
               rand_count = 0;
            }
         }
         ic = (int)light;
         if (ic < 0) {
            ic = 0;
         } else if (ic > (nc[0]-1)) {
            ic = nc[0]-1;
         }
         cam[0] = ic;
         *cimk = cam[0] + shade_offset;
         break;
      default:
         break;
   }
}

/* dr3_shade_rgb - construct colours for point */
void dr3_shade_rgb(h, v, d, c, dimk, rgbimk)
   float h,v,d,c[],*dimk,rgbimk[3];
{
   float light,red,green,blue;
   //bool found;

   switch (shade_index) {
      case SHADE_LM:
         *dimk = d;
         light = c[0];
         rgbimk[0] = light*material_red;
         rgbimk[1] = light*material_green;
         rgbimk[2] = light*material_blue;
         break;
      case SHADE_LD:
         *dimk = d;
         light = c[0];
         light = 1.0;
/*
         dr3_dtov(h,v,d, xv);
         dr3_vtow(xv, xw);
         dr3_scalar(xw, &s, &found);
*/
         if (scalar_defined) {
            dr3_rain(c[1], &red, &green, &blue);
            rgbimk[0] = light*red;
            rgbimk[1] = light*green;
            rgbimk[2] = light*blue;
         } else {
            rgbimk[0] = light;
            rgbimk[1] = light;
            rgbimk[2] = light;
         }
         break;
/*    Light:  c=[light],cam=[l] */
      case SHADE_L:
         *dimk = d;
         light = c[0];
         rgbimk[0] = light;
         rgbimk[1] = light;
         rgbimk[2] = light;
         break;
/*    Other:  do nothing */
      default:
         break;
   }
}

/* dr3_shade_rgbz - construct colours for point */
void dr3_shade_rgbz(xv, rgb)
   float xv[MAXD4],rgb[3];
{
   float light,s,xw[MAXD4];
   bool found;

   switch (shade_index) {
      case SHADE_LM:
         dr3_light(xv, normal, &light);
         rgb[0] = light*material_red;
         rgb[1] = light*material_green;
         rgb[2] = light*material_blue;
         break;
      case SHADE_LD:
         dr3_light(xv, normal, &light);
         dr3_vtow(xv, xw);
         dr3_scalar(xw, &s, &found);
         if (found) {
            dr3_rain(s, &rgb[0], &rgb[1], &rgb[2]);
            rgb[0] = light*rgb[0];
            rgb[1] = light*rgb[1];
            rgb[2] = light*rgb[2];
         } else {
            rgb[0] = light;
            rgb[1] = light;
            rgb[2] = light;
         }
         break;
      case SHADE_L:
         dr3_light(xv, normal, &light);
         rgb[0] = light;
         rgb[1] = light;
         rgb[2] = light;
         break;
      default:
         rgb[0] = sky_red;
         rgb[1] = sky_green;
         rgb[2] = sky_blue;
         break;
   }
}

/* dr3_light - define light from model */
void dr3_light(xv, normal, light)
   float xv[MAXD4],normal[MAXD],*light;
{
   float cost,dist2,lp[MAXD4];
   int id;

/* ambient component */
   *light = light_ambient;
/* normal means can look at directed light */
   if (normal_defined) {
/*    vector to light-point */
      dist2 = 0.0;
      for (id=0; id < MAXD; id++) {
         lp[id] = light_view[id] - xv[id];
         dist2 = dist2 + lp[id]*lp[id];
      }
/*    only do point source distance from vertex */
      if (dist2 > 0.0) {
         cost = 0.0;
         for (id=0; id < MAXD; id++) {
            cost = cost + normal[id]*lp[id];
         }
         cost = cost/sqrt(dist2);
/*       only surfaces that face the light are illuminated */
         if (cost > 0.0) {
            *light = *light + light_candle*cost;
         }
      }
      if (material_glow > 0.0) {
         if (normal[2] > 0.0) {
            *light = (*light) + normal[2]*material_glow*light_candle;
          }
      }
   }
   *light = (*light)*(*light);
}

/* dr3_wtod - transform "world" to "device" coordinates */
void dr3_wtod(x,y,z, h, v, d)
   float x,y,z, *h,*v,*d;
{
   float a[MAXD4],b[MAXD4];

/* compose world coordinate 4-vector */
   a[0] = x;
   a[1] = y;
   a[2] = z;
   a[3] = 1.0;
/* convert it to viewing coordinates [X,Y,Z] */
   dr3_vmult(vm, a, b);
/* apply projection transformation */
   *h = dhdxz*(b[0]/b[2]) + h0;
   *v = dvdyz*(b[1]/b[2]) + v0;
   *d = dddz*b[2] + d0;
}

/* dr3_get_class - classify underlying device */
void dr3_get_class(class) 
   int *class;
{
   bool got_rgbz,got_rgb,got_colour,got_change;
   int i,n;

   drp_avail_rgbz(&got_rgbz);
   drp_avail_rgb(&got_rgb);
   drp_avail_colour(&got_colour);
   drp_change_colour(&got_change);
   drp_get_colour_shape(&i,&n);
   if (got_rgbz) {
      *class = DEV_RGBZ;
   } else if (got_rgb) {
      *class = DEV_RGB;
   } else if (got_colour && got_change && (n > 16)) {
      *class = DEV_CM;
   } else {
      *class = DEV_VEC;
   }
}

void dr3_rain(phi, red, green, blue)
   float phi,*red,*green,*blue;
{
   if (phi < 0.5) {
      *red = 0.0;
      *green = 2.0*phi;
      *blue = 2.0*(0.5-phi);
   } else {
      *red = 2.0*(phi-0.5);
      *green = 2.0*(1.0-phi);
      *blue = 0.0;
   }
/*
   b = 2.0*(phi-1.0);
   *red = exp(-b*b);
   b = 2.0*(phi-0.5);
   *green = exp(-b*b);
   b = 2.0*phi;
   *blue = exp(-b*b);
*/
}

void dr3_scalar(xw, s, found)
   float xw[MAXD4],*s;
   bool *found;
{
   int k,id,j,con[MAXD],conl[MAXD],conr[MAXD],n,ix,iy,iz;
   float x0,xn,dx,wx,dy,wy,dz,wz,w,wsum,qsum,xoff[MAXD];

/* predefined scalar? */
   if (scalar_defined) {
      *found = TRUE;
      *s = scalar;
      return;			// was return(0) BM
   }
/* quick safe exist */
   if (no_s) {
      *found = FALSE;
      return;			// was return(0) BM
   }
/* find whether we are inside or not */
   *found = TRUE;
   for (id=0; id < MAXD; id++) {
/*    very quickly exclude points outside data box */
      x0 = xe_s[id][0];
      n = nc_s[id];
      xn = xe_s[id][n];
      if ((xw[id] < x0) || (xw[id] > xn)) {
         *found = FALSE;
         break;
      } else {
/*       assume uniform grid (even for boundary cells) */
         j = (int)(djdx_s[id]*(xw[id]-x0));
         con[id] = j;
         conl[id] = j-1;
         if (conl[id] < 0) {
            conl[id] = 0;
         }
         conr[id] = j+1;
         if (conr[id] > (n-1)) {
            conr[id] = n-1;
         }
         if (conl[id] > conr[id]) {
            *found = FALSE;
            break;
         }
         xoff[id] = djdx_s[id]*(x0-xw[id]) + 0.5;
      }
   }
   if (*found) {
      wsum = 0.0;
      qsum = 0.0;
      for (iz=conl[2]; iz <= conr[2]; iz++) {
         dz = (float)iz + xoff[2];
         wz = 1.0/(1.0 + dz*dz);
         for (iy=conl[1]; iy <= conr[1]; iy++) {
            dy = (float)iy + xoff[1];
            wy = 1.0/(1.0 + dy*dy);
            for (ix=conl[0]; ix <= conr[0]; ix++) {
               k = ix + nc_s[0]*(iy + nc_s[1]*iz);
               if (sa_s[k] == FLUID) {
                  dx = (float)ix + xoff[0];
                  wx = 1.0/(1.0 + dx*dx);
                  w = wx*wy*wz;
                  wsum = wsum + w;
                  qsum = qsum + w*qa_s[k];
               }
            }
         }
      }
      if (wsum > 0.0) {
         *s = (qsum/wsum - qmin_s)*dpdq_s;
      } else {
         *found = FALSE;
      }
   }
   (void)con; /* Silences the unused-but-set-variable warning safely BM */
}

/* dr3_vtod - transform "viewing" to "device" coordinates */
void dr3_vtod(xv, h, v, d)
   float xv[MAXD4], *h,*v,*d;
{
   *h = dhdxz*(xv[0]/xv[2]) + h0;
   *v = dvdyz*(xv[1]/xv[2]) + v0;
   *d = dddz*xv[2] + d0;
}

/* dr3_wtov - transform "world" to "viewing" coordinates */
void dr3_wtov(xo, xv)
   float xo[MAXD4], xv[MAXD4];
{

/* convert it to viewing coordinates [X,Y,Z] */
   dr3_vmult(vm, xo, xv);
}

/* dr3_dtov - transform "device" to "viewing" coordinates */
void dr3_dtov(h,v,d, xv)
   float xv[MAXD4], h,v,d;
{
   xv[2] = (d-d0)/dddz;
   xv[0] = (h-h0)*xv[2]/dhdxz;
   xv[1] = (v-v0)*xv[2]/dvdyz;
   xv[3] = 1.0;
}

/* dr3_vtow - transform "viewing" to "world" coordinates */
void dr3_vtow(xv, xw)
   float xv[MAXD4],xw[MAXD4];
{
   int id;

   for (id=0; id < MAXD; id++) {
      xw[id] = xv[0]*view_nx[id] + xv[1]*view_ny[id] + 
        xv[2]*view_nz[id] + view_w0[id];
   }
   xw[3] = 1.0;
}

void dr3_clock(cp)
   float *cp;
{
/*
   cpop_(cp);
*/
}

/* LOADING FIELD DATA */

/* dr3_loads - load scalar field (+grid,-vector) */
void dr3_loads(nc, xe, sa, qa, qmin, qmax, ios)
   int nc[],sa[],*ios;
   float *xe[MAXD],qa[],qmin,qmax;
{
   int id,n;

#ifdef DEBUG_DATA
   printf("dr3_loads:  called\n");
#endif
/* zap any others */
/* private copy of coordinate mesh */
   nd_s = nd;
   for (id=0; id < MAXD; id++) {
      n = nc[id];
      nc_s[id] = n;
      x0_s[id] = xe[id][0];
      dxdj_s[id] = (xe[id][n] - xe[id][0])/(float)n;
      djdx_s[id] = 1.0/dxdj_s[id];
      xe_s[id] = xe[id];
   }
/* scalar mapping */
   qmin_s = qmin;
   qmax_s = qmax;
   q0_s = qmin;
   dqdp_s = (qmax-qmin);
   dpdq_s = 1.0/dqdp_s;
/* pointers to sa[] and qa[] */
   sa_s = sa;
   qa_s = qa;
   no_s = FALSE;
   *ios = OK;
}

/* dr3_zaps - zap scalar field (free any memory) */
void dr3_zaps() {
   int id;

#ifdef DEBUG_DATA
   printf("dr3_zaps:  called\n");
#endif
/* if no scalar field report! */
   if (no_s) {
      printf("warning:  no scalar field to zap (dr3_zaps)\n");
   } else {
      for (id=0; id < MAXD; id++) {
         xe_s[id] = NULL;
      }
      sa_s = NULL;
      qa_s = NULL;
      no_s = TRUE;
   }
}

/* dr3_data_plane - data shade existing Lx image */
void dr3_data_plane(distance, ios)
   float distance;
   int *ios;
{
   int k,ih,iv,ilight,is,iother;
   float s,h,v,xv[MAXD4],xw[MAXD4],dcut;
   bool found;

/* determine "d" for this distance */
   if (distance < project_near) {
      printf("error:  data plane too near, use distance=%f\n",
         project_near);
      dcut = 0.0;
   } else if (distance > project_far) {
      printf("error:  data plane too far, use distance=%f\n",
        project_far);
      dcut = 1.0;
   } else {
      dcut = (-dddz*distance + d0);
   }
/* arrange an LD image */
   dr3_coerce_ld(ios);
   if (*ios) {
      printf("error:  unable to coerce LD image\n");
   } else if (!colour_bitmap) {
      printf("error:  need colour-bitmap image\n");
      *ios = ERR;
   } else {
      k=0;
      for (iv=0; iv < nv; iv++) {
         v = ((float)iv-vi0)/dvidv;
         for (ih=0; ih < nh; ih++) {
            h = ((float)ih-hi0)/dhidh;
            ilight = (cim[k] - shade_offset);
            iother = ilight/nc[0];
            ilight = ilight - iother*nc[0];
            dr3_dtov(h,v,dcut, xv);
            dr3_vtow(xv, xw);
            dr3_scalar(xw, &s, &found);
            if (found) {
               s = s*shade_didc[1];
               if (random_colour) {
                  s = s + randaf[rand_count++];
                  if (rand_count == MAXRAND) {
                     rand_count = 0;
                  }
               }
               is = (int)s;
               if (is < 0) {
                  is = 0;
               } else if (is > (nc[1]-2)) {
                  is = nc[1]-2;
               }
               is = is+1;
               if (dcut < dim[k]) {
                  cim[k] = (ilight+nc[0]*is) + shade_offset;
               } else {
                  cim[k] = ilight + shade_offset;
               }
            } else {
               cim[k] = ilight + shade_offset;
            }
            k++;
         }
      }
      drp_colour_bitmap(ih0,iv0,nh,nv,cim);
   }
}

/* dr3_data_shade - data shade existing Lx image */
void dr3_data_shade(ios)
   int *ios;
{
   int k,ih,iv,ilight,is,iother,kk;
   float s,h,v,xv[MAXD4],xw[MAXD4],red,green,blue;
   bool found;

   dr3_coerce_ld(ios);
   if (*ios) {
      printf("error:  unable to coerce LD image\n");
   } else if (colour_bitmap) {
      k=0;
      for (iv=0; iv < nv; iv++) {
         v = ((float)iv-vi0)/dvidv;
         for (ih=0; ih < nh; ih++) {
            ilight = (cim[k] - shade_offset);
            iother = ilight/nc[0];
            ilight = ilight - iother*nc[0];
            if (dim[k] < 1.01) {
               h = ((float)ih-hi0)/dhidh;
               dr3_dtov(h,v,dim[k], xv);
               dr3_vtow(xv, xw);
               dr3_scalar(xw, &s, &found);
               if (found) {
                  s = s*shade_didc[1];
                  if (random_colour) {
                     s = s + randaf[rand_count++];
                     if (rand_count == MAXRAND) {
                        rand_count = 0;
                     }
                  }
                  is = (int)s;
                  if (is < 0) {
                     is = 0;
                  } else if (is > (nc[1]-2)) {
                     is = nc[1]-2;
                  }
                  is = is+1;
                  cim[k] = (ilight+nc[0]*is) + shade_offset;
               } else {
                  cim[k] = ilight + shade_offset;
               }
            } else {
               cim[k] = ilight + shade_offset;
            }
            k++;
         }
      }
      drp_colour_bitmap(ih0,iv0,nh,nv,cim);
   } else if (rgb_bitmap) {
      k=0;
      for (iv=0; iv < nv; iv++) {
         v = ((float)iv-vi0)/dvidv;
         for (ih=0; ih < nh; ih++) {
            if (dim[k] < 1.01) {
               h = ((float)ih-hi0)/dhidh;
               dr3_dtov(h,v,dim[k], xv);
               dr3_vtow(xv, xw);
               dr3_scalar(xw, &s, &found);
               if (found) {
                  dr3_rain(s, &red, &green, &blue);
                  kk = k*3;
                  rgbim[kk] = rgbim[kk]*red;
                  rgbim[kk+1] = rgbim[kk+1]*green;
                  rgbim[kk+2] = rgbim[kk+2]*blue;
               }
            }
            k++;
         }
      }
      drp_rgb_bitmap(ih0,iv0,nh,nv,rgbim);
   } else {
      printf("error:  need (colour,rgb)-bitmap\n");
      *ios = ERR;
   }
}

/* dr3_coerce_wire - coerce existing shading to WIRE */
void dr3_coerce_wire(ios)
   int *ios;
{

/* is there one already? */
   if (noimage) {
      dr3_shade("WIRE", ios);
      if (!*ios) {
         dr3_create_image(ios);
         if (*ios) {
            printf("error:  creating WIRE image\n");
         }
      }
/* if so always can use it */
   } else {
      *ios = OK;
   }
}

/* dr3_coerce_ld - create or coerce existing Lx to LD */
void dr3_coerce_ld(ios)
   int *ios;
{
   int mc,id,k,ih,iv,old_offset,ilight,iother,kk;
   float grey;

/* is there one already? */
   if (noimage) {
      dr3_shade("LD", ios);
      if (!*ios) {
         dr3_create_image(ios);
         if (*ios) {
            printf("error:  creating LD image\n");
         }
      }
/* rgbz always ok */
   } else if (rgbz_bitmap) {
      dr3_shade("LD", ios);
      if (*ios) {
         printf("error:  doing rgbz LD shade\n");
      }
/* if so check it out against LD */
   } else if (shade_index == SHADE_LD) {
/* otherwise try to coerce Lx */
   } else if (shade_name[0] == 'L') {
      if (colour_bitmap) {
         mc = 1;
         for (id=0; id < (nd-1); id++) {
            mc = mc*nc[id];
         }
         old_offset = shade_offset;
         dr3_shade("LD", ios);
         if (!*ios) {
            k=0;
            for (iv=0; iv < nv; iv++) {
               for (ih=0; ih < nh; ih++) {
                  ilight = (cim[k] - old_offset);
                  iother = ilight/mc;
                  ilight = ilight - iother*mc;
                  cim[k] = ilight + shade_offset;
                  k++;
               }
            }
            drp_colour_bitmap(ih0,iv0,nh,nv,cim);
         }
      } else {
         dr3_shade("LD", ios);
         if (!*ios) {
            k = 0;
            for (iv=0; iv < nv; iv++) {
               for (ih=0; ih < nh; ih++) {
                  if (dim[k] < 1.01) {
                     kk = k*3;
                     grey = 0.333*(rgbim[kk]+rgbim[kk+1]+rgbim[kk+2]);
                     rgbim[kk] = grey;
                     rgbim[kk+1] = grey;
                     rgbim[kk+2] = grey;
                  }
                  k++;
               }
            }
         }
      }
   }
}

/* dr3_coerce_lm - coerce existing Lx to LM */
void dr3_coerce_lm(ios)
   int *ios;
{
   int mc,id,k,ih,iv,old_offset,ilight,iother;

/* is there one already? */
   if (noimage) {
      dr3_shade("LM", ios);
      if (*ios) {
         printf("..problem creating shade LM\n");
      } else {
         dr3_create_image(ios);
         if (*ios) {
            printf("error:  creating LM image\n");
         }
      }
/* if so check it out against LM */
   } else if (shade_index == SHADE_LM) {
/* rgbz is always ok */
   } else if (rgbz_bitmap) {
      dr3_shade("LM", ios);
      if (*ios) {
         printf("error:  creating LM for rgb(coerce_lm)\n");
      }
/* otherwise try to coerce Lx */
   } else if (shade_name[0] == 'L') {
      if (colour_bitmap) {
         mc = 1;
         for (id=0; id < (nd-1); id++) {
            mc = mc*nc[id];
         }
         old_offset = shade_offset;
         dr3_shade("LM", ios);
         if (!*ios) {
            k=0;
            for (iv=0; iv < nv; iv++) {
               for (ih=0; ih < nh; ih++) {
                  ilight = (cim[k] - old_offset);
                  iother = ilight/mc;
                  ilight = ilight - iother*mc;
                  cim[k] = ilight + shade_offset;
                  k++;
               }
            }
            drp_colour_bitmap(ih0,iv0,nh,nv,cim);
         }
      } else {
         dr3_shade("LD", ios);
      }
/* otherwise we have failed! */
   } else {
      printf("error:  cannot coerce %s image to LM\n", shade_name);
      *ios = ERR;
   }
}

void dr3_coerce_d(ios)
   int *ios;
{

/* is there one already? */
   if (noimage) {
      dr3_shade("D", ios);
      if (!*ios) {
         dr3_create_image(ios);
         if (*ios) {
            printf("error:  creating D image\n");
         }
      }
/* rgbz always OK */
   } else if (rgbz_bitmap) {
      dr3_shade("D", ios);
      if (*ios) {
         printf("error:  doing rgbz D shade\n");
      }
/* if so check it out against D or LD */
   } else if (shade_index == SHADE_D) {
      *ios = OK;
   } else if (shade_index == SHADE_LD) {
      *ios = OK;
/* otherwise we have failed */
   } else {
      printf("error:  cannot coerce %s image to xD\n", shade_name);
      *ios = ERR;
   }
}

/* ATTRIBUTE FUNCTIONS */

void dr3_seta(name)
   char name[];
{
   if (csame("random_shade",name)) {
      random_shade = TRUE;
   } else if (csame("norandom_shade",name)) {
      random_shade = FALSE;
   } else if (csame("random_colour",name)) {
      random_colour = TRUE;
   } else if (csame("norandom_colour",name)) {
      random_colour = FALSE;
   } else if (csame("wire_frame",name)) {
      wire_frame = TRUE;
   } else if (csame("nowire_frame",name)) {
      wire_frame = FALSE;
   } else if (csame("hard_surface",name)) {
      hard_surface = TRUE;
   } else if (csame("nohard_surface",name)) {
      hard_surface = FALSE;
   } else if (csame("lucid_surface",name)) {
      lucid_surface = TRUE;
   } else if (csame("nolucid_surface",name)) {
      lucid_surface = FALSE;
   } else {
      printf("error:  unknown attribute=%s (dr3_seta)\n", name);
   }
}

void dr3_set0i(name, i)
   char name[];
   int i;
{
   if (csame("ncolour",name)) {
      if (i < 3) {
         printf("error:  ncolour=%d not accepted (dr3_set0i)\n", i);
      } else {
         ncolour = i;
      }
   } else if (csame("nshade",name)) {
      if (i < 3) {
         printf("error:  nshade=%d not accepted (dr3_set0i)\n", i);
      } else {
         nshade = i;
      }
   } else {
      printf("error:  unknown attribute=%s (dr3_set0i)\n", name);
   }
}

/* SOME SPECIALIST ENQUIRY FUNCTIONS */

/* dr3_get_shade - get shading model to caller */
void dr3_get_shade(scp)
   scene_struct *scp;
{
   int id,i;

   scp->dev_class = dev_class;
   scp->vector_bitmap = vector_bitmap;
   scp->colour_bitmap = colour_bitmap;
   scp->rgb_bitmap = rgb_bitmap;
   scp->rgbz_bitmap = rgbz_bitmap;
   scp->noshade = noshade;
   if (!noshade) {
      scp->shade_index = shade_index;
      ccopy(shade_name, 0, SZNAM, scp->shade_name);
      scp->shade_nd = nd;
      for (id=0; id < nd; id++) {
         scp->shade_nc[id] = nc[id];
         scp->shade_didc[id] = (float)(nc[id]-1);
      }
      if (colour_bitmap) {
         scp->shade_ncm = ncm;
         for (i=0; i < ncm; i++) {
            scp->shade_red[i] = shade_red[i];
            scp->shade_green[i] = shade_green[i];
            scp->shade_blue[i] = shade_blue[i];
         }
         scp->shade_offset = shade_offset;
      }
      if (shade_index == SHADE_LM) {
         for (i=0; i < nc[1]; i++) {
            ccopy(shade_mat[i], 0, SZNAM, scp->shade_mat[i]);
            scp->shade_mat_r[i] = shade_mat_r[i];
            scp->shade_mat_g[i] = shade_mat_g[i];
            scp->shade_mat_b[i] = shade_mat_b[i];
         }
      }
      scp->sky_red = sky_red;
      scp->sky_green = sky_green;
      scp->sky_blue = sky_blue;
   }
}

/* dr3_get_image_type - pop scene to caller */
void dr3_get_image_type(scp)
   scene_struct *scp;
{

   scp->noimage = noimage;
   if (!noimage) {
      scp->dev_class = dev_class;
      scp->vector_bitmap = vector_bitmap;
      scp->colour_bitmap = colour_bitmap;
      scp->rgb_bitmap = rgb_bitmap;
      scp->rgbz_bitmap = rgbz_bitmap;
   }
}

void dr3_default_limits(nh,nv, ivmin,ivmax,ihmin,ihmax)
   int nh,nv, *ivmin,*ivmax,ihmin[],ihmax[];
{
   int iv;

   *ivmin = nv;
   *ivmax = (-1);
   for (iv=0; iv < nv; iv++) {
      ihmin[iv] = nh;
      ihmax[iv] = (-1);
   }
}

void dr3_extend_limits(mh, mv, jvmin, jvmax, jhmin, jhmax,
  ivmin, ivmax, ihmin, ihmax)
  int mh,mv,jvmin,jvmax,jhmin[],jhmax[], *ivmin,*ivmax,
    ihmin[],ihmax[];
{
   int jh,jv,k;

/* new profiles */
   if ((*ivmin) > (*ivmax)) {
      for (jv=jvmin; jv <= jvmax; jv++) {
         ihmin[jv] = jhmin[jv];
         ihmax[jv] = jhmax[jv];
         k = mh*jv + jhmin[jv];
         for (jh=jhmin[jv]; jh <= jhmax[jv]; jh++) {
            nib[k++] = NULL;
         }
      }
      *ivmin = jvmin;
      *ivmax = jvmax;
/* adjust old ones */
   } else {
/*    Bottom */
      if (jvmin < (*ivmin)) {
         for (jv=jvmin; jv < (*ivmin); jv++) {
            k = mh*jv + jhmin[jv];
            for (jh=jhmin[jv]; jh <= jhmax[jv]; jh++) {
               nib[k++] = NULL;
            }
            ihmin[jv] = jhmin[jv];
            ihmax[jv] = jhmax[jv];
         }
         *ivmin = jvmin;
      }
/*    Top */
      if ((*ivmax) < jvmax) {
         for (jv=jvmax; jv > (*ivmax); jv--) {
            k = mh*jv + jhmin[jv];
            for (jh=jhmin[jv]; jh <= jhmax[jv]; jh++) {
               nib[k++] = NULL;
            }
            ihmin[jv] = jhmin[jv];
            ihmax[jv] = jhmax[jv];
         }
         *ivmax = jvmax;
      }
/*    Interior */
      for (jv=jvmin; jv <= jvmax; jv++) {
/*       Left */
         if (jhmin[jv] < ihmin[jv]) {
            k = mh*jv + jhmin[jv];
            for (jh=jhmin[jv]; jh < ihmin[jv]; jh++) {
               nib[k++] = NULL;
            }
            ihmin[jv] = jhmin[jv];
         }
/*       Right */
         if (jhmax[jv] > ihmax[jv]) {
            k = mh*jv + jhmax[jv];
            for (jh=jhmax[jv]; jh > ihmax[jv]; jh--) {
               nib[k--] = NULL;
            }
            ihmax[jv] = jhmax[jv];
         }
      }
   }
}

/* DUMP AND LOAD */

/* dr3_dump_scene - dump current dr3 scene to unit */
void dr3_dump_scene(unit, ios)
   int unit,*ios;
{
   int n,i;

#ifdef DEBUG_SCENE
   printf("dr3_dump_scene:  called\n");
#endif
/* dr3scene */
   ds_put_begin(unit, "dr3scene", ios);
   if (*ios) {
      printf("error:  begin(dr3scene)\n");
   } else {
/*    view */
      if (!noview) {
         ds_put_begin(unit, "view", ios);
         ds_put0r(unit, "azim", view_azim, ios);
         ds_put0r(unit, "altd", view_altd, ios);
         ds_put0r(unit, "twist", view_twist, ios);
         ds_put0r(unit, "distance", view_distance, ios);
         ds_put1r(unit, "focus", MAXD4, view_focus, ios);
         ds_put1r(unit, "point", MAXD4, view_point, ios);
         //ds_put1r(unit, "matrix", MAXD4*MAXD4, view_matrix, ios);
         //ds_put1r(unit, "world", MAXD4*MAXD4, view_world, ios);
         ds_put1r(unit, "matrix", MAXD4*MAXD4, &view_matrix[0][0], ios);	// BM
         ds_put1r(unit, "world", MAXD4*MAXD4, &view_world[0][0], ios);		// BM
         ds_put1r(unit, "nx", MAXD, view_nx, ios);
         ds_put1r(unit, "ny", MAXD, view_ny, ios);
         ds_put1r(unit, "nz", MAXD, view_nz, ios);
         ds_put1r(unit, "w0", MAXD, view_w0, ios);
         ds_put_end(unit, "view", ios);
      }
/*    projection */
      if (!noproject) {
         ds_put_begin(unit, "project", ios);
         ds_put0r(unit, "anglev", project_anglev, ios);
         ds_put0r(unit, "near", project_near, ios);
         ds_put0r(unit, "far", project_far, ios);
         ds_put0r(unit, "hmin", project_hmin, ios);
         ds_put0r(unit, "hmax", project_hmax, ios);
         ds_put0r(unit, "vmin", project_vmin, ios);
         ds_put0r(unit, "vmax", project_vmax, ios);
         ds_put0r(unit, "dhdxz", dhdxz, ios);
         ds_put0r(unit, "h0", h0, ios);
         ds_put0r(unit, "dvdyz", dvdyz, ios);
         ds_put0r(unit, "v0", v0, ios);
         ds_put0r(unit, "dddz", dddz, ios);
         ds_put0r(unit, "d0", d0, ios);
         ds_put_end(unit, "project", ios);
      }
/*    light */
      if (!nolight) {
         ds_put_begin(unit, "light", ios);
         ds_put0r(unit, "azim", light_azim, ios);
         ds_put0r(unit, "altd", light_altd, ios);
         ds_put0r(unit, "distance", light_distance, ios);
         ds_put1r(unit, "point", MAXD4, light_point, ios);
         ds_put1r(unit, "view", MAXD4, light_view, ios);
         ds_put0r(unit, "ambient", light_ambient, ios);
         ds_put0r(unit, "candle", light_candle, ios);
         ds_put_end(unit, "light", ios);
      }
/*    image */
      if (!noimage) {
         ds_put_begin(unit, "image", ios);
         ds_put0r(unit, "dhidh", dhidh, ios);
         ds_put0r(unit, "hi0", hi0, ios);
         ds_put0r(unit, "dvidv", dvidv, ios);
         ds_put0r(unit, "vi0", vi0, ios);
         ds_put0i(unit, "ih0", ih0, ios);
         ds_put0i(unit, "iv0", iv0, ios);
         ds_put0i(unit, "nh", nh, ios);
         ds_put0i(unit, "nv", nv, ios);
/*       Device classification */
         ds_put0i(unit, "dev_class", dev_class, ios);
         ds_put0b(unit, "vector_bitmap", vector_bitmap, ios);
         ds_put0b(unit, "colour_bitmap", colour_bitmap, ios);
         ds_put0b(unit, "rgb_bitmap", rgb_bitmap, ios);
         ds_put0b(unit, "rgbz_bitmap", rgbz_bitmap, ios);
/*       Actual image */
         if (image_owned) {
            n = nh*nv;
            if (dim != NULL) {
               ds_put1r(unit, "dim", n, dim, ios);
            }
            if (colour_bitmap) {
               if (cim != NULL) {
                  ds_put1i(unit, "cim", n, cim, ios);
               }
            } else if (rgb_bitmap || rgbz_bitmap) {
               if (rgbim != NULL) {
                  ds_put1r(unit, "rgbim", n*3, rgbim, ios);
               }
            }
         }
/*       Shade */
         if (!noshade) {
            ds_put_begin(unit, "shade", ios);
            ds_put0i(unit, "index", shade_index, ios);
            ds_put0s(unit, "name", shade_name, ios);
            ds_put0i(unit, "nd", nd, ios);
            ds_put1i(unit, "nc", nd, nc, ios);
            if (colour_bitmap) {
               ds_put0i(unit, "ncm", ncm, ios);
               ds_put1r(unit, "red", ncm, shade_red, ios);
               ds_put1r(unit, "green", ncm, shade_green, ios);
               ds_put1r(unit, "blue", ncm, shade_blue, ios);
               ds_put0i(unit, "offset", shade_offset, ios);
            }
            ds_put0r(unit, "sky_red", sky_red, ios);
            ds_put0r(unit, "sky_green", sky_green, ios);
            ds_put0r(unit, "sky_blue", sky_blue, ios);
            if (shade_index == SHADE_LM) {
               ds_put_begin(unit, "materials", ios);
               for (i=0; i < nc[1]; i++) {
                  ds_put0s(unit, "material", shade_mat[i], ios);
               }
               ds_put_end(unit, "materials", ios);
            }
            ds_put1r(unit, "shade_mat_r", nc[1], shade_mat_r, ios);
            ds_put1r(unit, "shade_mat_g", nc[1], shade_mat_g, ios);
            ds_put1r(unit, "shade_mat_b", nc[1], shade_mat_b, ios);
            ds_put_end(unit, "shade", ios);
         }
         ds_put_end(unit, "image", ios);
      }
/*    polylines */
      if (ns > 0) {
         ds_put_begin(unit, "polylines", ios);
         ds_put0i(unit, "ns", ns, ios);
         ds_put0i(unit, "np", np, ios);
         ds_put1i(unit, "ip0", ns, ip0, ios);
         ds_put1i(unit, "nps", ns, nps, ios);
         ds_put1r(unit, "hp", np, hp, ios);
         ds_put1r(unit, "vp", np, vp, ios);
         ds_put1r(unit, "dp", np, dp, ios);
         ds_put_end(unit, "polylines", ios);
      }
/*    end dr3scene */
      ds_put_end(unit, "dr3scene", ios);
   }
}

/* dr3_load_scene - dump current dr3 scene to unit */
void dr3_load_scene(unit, ios)
   int unit,*ios;
{
   char name[SZNAM],type[SZNAM];
   int n,alt_ih0,alt_iv0,alt_dev_class,alt_shade_offset,mh,mv,
     ndi,nci[MAXCD],id,i,icoff,k;
   float alt_project_hmin,alt_project_hmax,alt_project_vmin,
     alt_project_vmax,alt_h0,alt_v0,alt_hi0,alt_vi0,hgmin,hgmax,
     vgmin,vgmax,delh,delv,hglen,vglen,hlen,vlen,hoff,voff,dh,dv,
     reds[1],greens[1],blues[1];

#ifdef DEBUG_SCENE
   printf("dr3_load_scene:  called\n");
#endif
/* setup defaults */
   noview = TRUE;
   noproject = TRUE;
   nolight = TRUE;
   noimage = TRUE;
   noshade = TRUE;
   ns = 0;
   np = 0;
   noscene = TRUE;
/* assume that the caller has found begin(dr3scene)? */
   ds_get_begin(unit, "dr3scene", ios);
   while (ds_get(unit, name, type, ios)) {
/*    view */
      if (!csame("begin",name)) {
         printf("error:  %s unexpected in dr3-record\n", name);
         *ios = ERR;
      } else if (csame("view",type)) {
         noview = FALSE;
         ds_get0r(unit, "azim", &view_azim, ios);
         ds_get0r(unit, "altd", &view_altd, ios);
         ds_get0r(unit, "twist", &view_twist, ios);
         ds_get0r(unit, "distance", &view_distance, ios);
         ds_get1r(unit, "focus", MAXD4, view_focus, ios);
         ds_get1r(unit, "point", MAXD4, view_point, ios);
         //ds_get1r(unit, "matrix", MAXD4*MAXD4, view_matrix, ios);
         //ds_get1r(unit, "world", MAXD4*MAXD4, view_world, ios);
         ds_get1r(unit, "matrix", MAXD4*MAXD4, &view_matrix[0][0], ios);	// BM
         ds_get1r(unit, "world", MAXD4*MAXD4, &view_world[0][0], ios);		// BM
         ds_get1r(unit, "nx", MAXD, view_nx, ios);
         ds_get1r(unit, "ny", MAXD, view_ny, ios);
         ds_get1r(unit, "nz", MAXD, view_nz, ios);
         ds_get1r(unit, "w0", MAXD, view_w0, ios);
         ds_get_end(unit, "view", ios);
         dr3_zap_matrix();
/*    projection */
      } else if (csame("project",type)) {
         noproject = FALSE;
         ds_get0r(unit, "anglev", &project_anglev, ios);
         ds_get0r(unit, "near", &project_near, ios);
         ds_get0r(unit, "far", &project_far, ios);
         ds_get0r(unit, "hmin", &alt_project_hmin, ios);
         ds_get0r(unit, "hmax", &alt_project_hmax, ios);
         ds_get0r(unit, "vmin", &alt_project_vmin, ios);
         ds_get0r(unit, "vmax", &alt_project_vmax, ios);
         ds_get0r(unit, "dhdxz", &dhdxz, ios);
         ds_get0r(unit, "h0", &alt_h0, ios);
         ds_get0r(unit, "dvdyz", &dvdyz, ios);
         ds_get0r(unit, "v0", &alt_v0, ios);
         ds_get0r(unit, "dddz", &dddz, ios);
         ds_get0r(unit, "d0", &d0, ios);
         ds_get_end(unit, "project", ios);
/*    light */
      } else if (csame("light",type)) {
         nolight = FALSE;
         ds_get0r(unit, "azim", &light_azim, ios);
         ds_get0r(unit, "altd", &light_altd, ios);
         ds_get0r(unit, "distance", &light_distance, ios);
         ds_get1r(unit, "point", MAXD4, light_point, ios);
         ds_get1r(unit, "view", MAXD4, light_view, ios);
         ds_get0r(unit, "ambient", &light_ambient, ios);
         ds_get0r(unit, "candle", &light_candle, ios);
         ds_get_end(unit, "light", ios);
/*    image */
      } else if (csame("image",type)) {
         noimage = FALSE;
         ds_get0r(unit, "dhidh", &dhidh, ios);
         ds_get0r(unit, "hi0", &alt_hi0, ios);
         ds_get0r(unit, "dvidv", &dvidv, ios);
         ds_get0r(unit, "vi0", &alt_vi0, ios);
         ds_get0i(unit, "ih0", &alt_ih0, ios);
         ds_get0i(unit, "iv0", &alt_iv0, ios);
         ds_get0i(unit, "nh", &nh, ios);
         ds_get0i(unit, "nv", &nv, ios);
/*       Device classification */
         ds_get0i(unit, "dev_class", &alt_dev_class, ios);
         ds_get0b(unit, "vector_bitmap", &vector_bitmap, ios);
         ds_get0b(unit, "colour_bitmap", &colour_bitmap, ios);
         ds_get0b(unit, "rgb_bitmap", &rgb_bitmap, ios);
         ds_get0b(unit, "rgbz_bitmap", &rgbz_bitmap, ios);
/*       Actual image */
         if (!*ios) {
            dr3_alloc_image(ios);
            if (!*ios) {
               n = nh*nv;
               ds_get1r(unit, "dim", n, dim, ios);
               if (colour_bitmap) {
                  ds_get1i(unit, "cim", n, cim, ios);
               } else if (rgb_bitmap || rgbz_bitmap) {
                  ds_get1r(unit, "rgbim", n*3, rgbim, ios);
               }
            }
         }
/*       Shade */
         noshade = FALSE;
         ds_get_begin(unit, "shade", ios);
         ds_get0i(unit, "index", &shade_index, ios);
         ds_get0s(unit, "name", shade_name, ios);
         ds_get0i(unit, "nd", &nd, ios);
         ds_get1i(unit, "nc", nd, nc, ios);
         for (id=0; id < nd; id++) {
            shade_didc[id] = (float)(nc[id]-1);
         }
         if (colour_bitmap) {
            ds_get0i(unit, "ncm", &ncm, ios);
            ds_get1r(unit, "red", ncm, shade_red, ios);
            ds_get1r(unit, "green", ncm, shade_green, ios);
            ds_get1r(unit, "blue", ncm, shade_blue, ios);
            ds_get0i(unit, "offset", &alt_shade_offset, ios);
         }
/*       Options in shade */
         sky_red = 0.8;
         sky_green = 0.9;
         sky_blue = 1.0;
         if (shade_index == SHADE_LM) {
            for (i=0; i < nc[1]; i++) {
               ccopy(" ", 0, SZNAM, shade_mat[i]);
               shade_mat_r[i] = 0.5;
               shade_mat_g[i] = 0.5;
               shade_mat_b[i] = 0.5;
            }
         }
         while (ds_get(unit, name, type, ios)) {
            if (csame("begin",name)) {
               if (csame("materials",type)) {
                  for (i=0; i < nc[1]; i++) {
                     ds_get0s(unit, "material", shade_mat[i], ios);
                  }
                  ds_get_end(unit, "materials", ios);
               } else {
                  printf("error:  %s unexpected in materials\n", type);
                  *ios = ERR;
               }
            } else if (csame("sky_red",name)) {
               ds_get0r(unit, "#", &sky_red, ios);
            } else if (csame("sky_green",name)) {
               ds_get0r(unit, "#", &sky_green, ios);
            } else if (csame("sky_blue",name)) {
               ds_get0r(unit, "#", &sky_blue, ios);
            } else if (csame("shade_mat_r",name)) {
               if (shade_index == SHADE_LM) {
                  ds_get1r(unit, "#", nc[1], shade_mat_r, ios);
               }
            } else if (csame("shade_mat_g",name)) {
               if (shade_index == SHADE_LM) {
                  ds_get1r(unit, "#", nc[1], shade_mat_g, ios);
               }
            } else if (csame("shade_mat_b",name)) {
               if (shade_index == SHADE_LM) {
                  ds_get1r(unit, "#", nc[1], shade_mat_b, ios);
               }
            } else {
               printf("error:  %s unexpected in shade-record\n", name);
               *ios = ERR;
            }
         }
/*       End image */
         ds_get_end(unit, "image", ios);
/*    polylines */
      } else if (csame("polylines",type)) {
         ds_get0i(unit, "ns", &ns, ios);
         ds_get0i(unit, "np", &np, ios);
         ds_get1i(unit, "ip0", ns, ip0, ios);
         ds_get1i(unit, "nps", ns, nps, ios);
         ds_get1r(unit, "hp", np, hp, ios);
         ds_get1r(unit, "vp", np, vp, ios);
         ds_get1r(unit, "dp", np, dp, ios);
         ds_get_end(unit, "polylines", ios);
      } else {
         printf("error:  %s unexpected in dr3scene-record\n", name);
         *ios = ERR;
      }
   }
/* Get the actual device class */
   dr3_get_class(&dev_class);
/* Check for fails */
   if (*ios) {
      dr3_zap_scene();
/* Test device class */
   } else if (alt_dev_class != dev_class) {
      printf("error:  wrong class of device\n");
      *ios = ERR;
/* Try to display */
   } else {
/*    Map the image onto available display surface */
      drp_get_hv(&hgmin, &hgmax, &vgmin, &vgmax);
      hglen = hgmax - hgmin;
      vglen = vgmax - vgmin;
      hlen = alt_project_hmax - alt_project_hmin;
      vlen = alt_project_vmax - alt_project_vmin;
      drp_bitmap(&mh,&mv,&dh,&dv,&h0,&v0);
      hoff = 0.5*(hglen - hlen);
      voff = 0.5*(vglen - vlen);
      ih0 = (hgmin + hoff)/dh;
      iv0 = (vgmin + voff)/dv;
/*    What we must add to [h,v] coordinates for new placement */
      delh = (float)(ih0-alt_ih0)*dh;
      delv = (float)(iv0-alt_iv0)*dv;
/*    data window */
      project_hmin = alt_project_hmin + delh;
      project_hmax = alt_project_hmax + delh;
      project_vmin = alt_project_vmin + delv;
      project_vmax = alt_project_vmax + delv;
/*    fixup hv-projection offset */
      h0 = 0.5*(project_hmin + project_hmax);
      v0 = 0.5*(project_vmin + project_vmax);
/*    fixup pixel<>coordinate offset */
      hi0 = alt_hi0 - delh/dh;
      vi0 = alt_vi0 - delv/dv;
/*    What about the matrix stack??? */
/*    In colour_bitmap mode we must address colour-ramp */
      if (shade_index == SHADE_WIRE) {
      } else if (colour_bitmap) {
/*       Find a suitable ramp to support shade model? */
         drp_get_ramp(shade_name, &ndi, nci, &shade_offset);
         if (shade_offset < 0) {
            drp_put_ramp(shade_name, nd, nc, &shade_offset);
            if (shade_offset < 0) {
               printf("error:  not enough colour space\n");
               *ios = ERR;
            } else {
               drp_set_colours(shade_offset, ncm, shade_red, 
                 shade_green, shade_blue);
               drp_put_paint("SKY", &sky_index);
               if (sky_index < 0) {
                  printf("error:  no paint room for sky\n");
               } else {
                  reds[0] = sky_red;
                  greens[0] = sky_green;
                  blues[0] = sky_blue;
                  drp_set_colours(sky_index, 1, reds, greens, blues);
               }
            }
         } else if (ndi != nd) {
            printf("error:  %s ramp mismatch\n", shade_name);
            *ios = ERR;
         } else {
            for (id=0; id < nd; id++) {
               if (nci[id] != nc[id]) {
                  printf("error:  %s ramp mismatch\n", shade_name);
                  *ios = ERR;
                  break;
               }
            }
         }
/*       If all is well, adjust bitmap */
         if (!*ios) {
            icoff = shade_offset - alt_shade_offset;
            n = nh*nv;
            for (k=0; k < n; k++) {
               cim[k] = cim[k] + icoff;
            }
         }
      }
/*    Polylines get adjustment */
      for (i=0; i < np; i++) {
         hp[i] = hp[i] + delh;
         vp[i] = vp[i] + delv;
      }
   }
/* If at the end of all this we are OK, refresh display */
   if (!*ios) {
      dr3_fresh_scene();
   } else {
      dr3_zap_scene();
   }
}

void dr3_wire_icon(size, x, y, z)
   float size,x,y,z;
{
   dr3_wire_block(x-size, x+size, y-size, y+size, z-size, z+size);
}

void dr3_solid_icon(s, size, x, y, z)
   float s,size,x,y,z;
{
/* front */
   dr3_begin_polygon();
   dr3_polygon_scalar(s);
   dr3_polygon_normal(0.0, 0.0, -1.0);
   dr3_polygon_vertex(x-size,y-size, z-size);
   dr3_polygon_vertex(x+size,y-size, z-size);
   dr3_polygon_vertex(x+size, y+size, z-size);
   dr3_polygon_vertex(x-size, y+size, z-size);
   dr3_end_polygon();
/* back */
   dr3_begin_polygon();
   dr3_polygon_scalar(s);
   dr3_polygon_normal(0.0, 0.0, 1.0);
   dr3_polygon_vertex(x-size,y-size, z+size);
   dr3_polygon_vertex(x+size,y-size, z+size);
   dr3_polygon_vertex(x+size, y+size, z+size);
   dr3_polygon_vertex(x-size, y+size, z+size);
   dr3_end_polygon();
/* left side */
   dr3_begin_polygon();
   dr3_polygon_scalar(s);
   dr3_polygon_normal(-1.0, 0.0, 0.0);
   dr3_polygon_vertex(x-size,y-size, z-size);
   dr3_polygon_vertex(x-size, y+size, z-size);
   dr3_polygon_vertex(x-size, y+size, z+size);
   dr3_polygon_vertex(x-size,y-size, z+size);
   dr3_end_polygon();
/* right side */
   dr3_begin_polygon();
   dr3_polygon_scalar(s);
   dr3_polygon_normal(1.0, 0.0, 0.0);
   dr3_polygon_vertex(x+size,y-size, z-size);
   dr3_polygon_vertex(x+size, y+size, z-size);
   dr3_polygon_vertex(x+size, y+size, z+size);
   dr3_polygon_vertex(x+size,y-size, z+size);
   dr3_end_polygon();
/* bottom side */
   dr3_begin_polygon();
   dr3_polygon_scalar(s);
   dr3_polygon_normal(0.0, -1.0, 0.0);
   dr3_polygon_vertex(x-size,y-size, z-size);
   dr3_polygon_vertex(x+size,y-size, z-size);
   dr3_polygon_vertex(x+size,y-size, z+size);
   dr3_polygon_vertex(x-size,y-size, z+size);
   dr3_end_polygon();
/* top side */
   dr3_begin_polygon();
   dr3_polygon_scalar(s);
   dr3_polygon_normal(0.0, 1.0, 0.0);
   dr3_polygon_vertex(x-size, y+size, z-size);
   dr3_polygon_vertex(x+size, y+size, z-size);
   dr3_polygon_vertex(x+size, y+size, z+size);
   dr3_polygon_vertex(x-size, y+size, z+size);
   dr3_end_polygon();
}

void dr3_begin() {
   drp_set_hv(project_hmin, project_hmax, project_vmin, project_vmax);
   if (rgbz_bitmap) {
      rgbz_begin();
/*
      printf("dr3_begin:  not calling rgbz_put_bitmap\n");
      rgbz_put_bitmap(ih0, iv0, nh, nv, rgbim, dim);
*/
   }
}

void dr3_end() {
   if (rgbz_bitmap) {
      rgbz_end();
/*
      printf("dr3_end:  not calling rgbz_get_bitmap\n");
      rgbz_get_bitmap(ih0, iv0, nh, nv, rgbim, dim);
*/
   }
}
