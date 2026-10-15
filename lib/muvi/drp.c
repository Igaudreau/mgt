/* drp.c:  2D drawing package (generic) */
#include <stdio.h>
#include "drp.h"
#include "machine.h"
#define MAXDEV 12
#define MAXLINESTYLE 1
#define MAXP 4096
#define MAXMARKER 4

struct drp_surface {
   char name[SZNAM],type[SZNAM];
   bool used,clear;
   struct dev *p;
   float marker_size,line_size; float hmin,hmax,vmin,vmax;
   int ihmin,ihmax,ivmin,ivmax;
};

static struct drp_surface surface[MAXDEV], *wp;
static struct dev *p;
static int wc=(-1);


static char line_styles[MAXLINESTYLE][SZNAM];
static float line_scale=1.0,marker_scale=0.2;
static int line_style=0;

static bool drp_pos_rel=TRUE,drp_size_rel=TRUE;
static float drp_pos_rh=1.0,drp_pos_rv=1.0,drp_size_rh=0.5,
  drp_size_rv=0.5, drp_pos_ah,drp_pos_av,drp_size_ah,drp_size_av;

static float hp[MAXP],vp[MAXP],hg[MAXP],vg[MAXP];
static int np=0,ng=0;
static int is0[]={0,2,4,5},is1[]={1,3,4,5},ip0[]={0,2,4,6,8,13},
  ip1[]={1,3,5,7,12,17},
  dh[]={-1,1,0,0, -1,1,-1,1, -1,1,1,-1,-1, 0,1,0,-1,0},
  dv[]={0,0,-1,1, 1,-1,-1,1, -1,-1,1,1,-1, -1,0,1,0,-1};

static float drp_hc,drp_vc,drp_dh,drp_dv, drp_hmin,drp_hmax,
  drp_vmin,drp_vmax;



/* EXPLICIT DECLARATIONS ADDED BY BM */
int drp_prepare();
int drp_flush();
int drp_print();
int drp_avail();
int drp_pos_r(float rh, float rv);
int drp_pos_a(float ah, float av);
int drp_size_r(float rh, float rv);
int drp_size_a(float ah, float av);
int drp_open(char *name, char *type, int *ios);
int drp_select(char name[], int *ios);
int drp_close(char name[], int *ios);
int drp_get_count(int *count);
int drp_active(bool *active);
int drp_get_name(char name[]);
int drp_get_type(char type[]);
int drp_set_hv(float h1, float h2, float v1, float v2);
int drp_zap_hv();
int drp_get_hv(float *h1,float *h2,float *v1,float *v2);
int drp_setup_hv();
int drp_setup_clip();
int drp_erase();
int drp_rub(float h1, float h2, float v1, float v2);
int drp_avail_rub(bool *avail);
int drp_line(float h1, float v1, float h2, float v2);
int drp_polyline(int na, float ha[], float va[]);
int drp_set_line_width(float line_width);
int drp_get_line_width(float *line_width);
int drp_set_line_style(int style);
int drp_get_line_styles(int *n);
int drp_get_line_style(int *style);
int drp_bp();
int drp_ap(float h, float v);
int drp_ep();
int drp_clip(float *hh1, float *vv1, bool *clip1, float *hh2, float *vv2, bool *clip2, bool *empty);
int drp_mode(float cmin, float cmax, float c, int *ic);
int drp_rect(float h1, float h2, float v1, float v2);
int drp_avail_rect(bool *avail);
int drp_polygon(int na, float ha[], float va[]);
int drp_avail_polygon(bool *avail);
int drp_bg();
int drp_ag(float h, float v);
int drp_eg();
int drp_text(float h, float v, char str[], char just[2]);
int drp_text_box(char str[], float *dh, float *dv);
int drp_set_font_style(char fontstyle[]);
int drp_get_font_style(char fontstyle[]);
int drp_set_font_size(float fontsize);
int drp_get_font_size(float *fontsize);
int drp_avail_fontscale(bool *avail);
int drp_set_font();
int drp_get_char(float *height,float *width,float *depth);
int drp_avail_cursor(bool *avail);
int drp_cursor(float *h, float *v, int *button);
int drp_avail_colour(bool *avail);
int drp_set_colours(int offset, int size, float reds[], float greens[], float blues[]);
int drp_set_colour(int colour_index);
int drp_change_colour(bool *change);
int drp_get_colour_shape(int *offset, int *size);
int drp_get_colour(int *colour_index);
int drp_get_colours(int offset, int size, float reds[], float greens[], float blues[]);
int drp_avail_rgb(bool *avail);
int drp_set_rgb(float red, float green, float blue);
int drp_set_rgb_background(float red, float green, float blue);
int drp_get_rgb(float *red, float *green, float *blue);
int drp_get_rgb_background(float *red_background,float *green_background,float *blue_background);
int drp_avail_rgbz(bool *avail);
int rgbz_setup(float anglev, float near, float far, float hmin, float hmax, float vmin, float vmax);
int rgbz_begin();
int rgbz_end();
int rgbz_clear();
int rgbz_begin_polygon();
int rgbz_polygon_vertex(float x, float y, float z);
int rgbz_end_polygon();
int rgbz_begin_polyline();
int rgbz_polyline_vertex(float x, float y, float z);
int rgbz_end_polyline();
int rgbz_put_bitmap(int ih0,int iv0,int nh,int nv,float rgb_image[],float d_image[]);
int rgbz_get_bitmap(int ih0,int iv0,int nh,int nv,float rgb_image[],float d_image[]);
int drp_set_marker_scale(float scale);
int drp_set_marker_size(float marker_size);
int drp_get_marker_scale(float *scale);
int drp_get_marker_size(float *marker_size);
int drp_marker(int marker, float h, float v);
int drp_avail_bitmap(bool *avail);
int drp_bitmap(int *nh,int *nv,float *dh,float *dv,float *h0,float *v0);
int drp_get_bitmap(int *ih0,int *iv0,int *nh,int *nv,float *dhidh,float *dvidv,float *hi0,float *vi0);
int drp_colour_bitmap(int ih0,int iv0,int nh,int nv,int image[]);
int drp_rgb_bitmap(int ih0,int iv0,int nh,int nv,float rgb_image[]);
int drp_zap_colour();
int drp_zap_ramps();
int drp_print_ramps();
int drp_put_ramp(char name[], int nd, int nc[], int *offset);
int drp_get_ramp(char name[], int *nd, int nc[], int *offset);
int drp_find_ramp(int nd, char name[]);
int drp_zap_paints();
int drp_put_paint(char name[], int *offset);
int drp_get_paint(char name[], int *offset);
int drp_set_paint(int offset, float red, float green, float blue);
int drp_print_paints();
int drp_get_colour_space(int *space);
int crap();
int drp_gflush();

extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);
extern int ps_size_r(float rh, float rv);
extern int ps_size_a(float ah, float av);
extern int ps_open(char *name, char *type, struct dev *(*p), int *ios);
extern int ps_setup_window(float hlen,float vlen, float *hsize,float *vsize);
extern int ps_close(char name[], int *ios);
extern int ps_select(char name[], int *ios);
extern int ps_erase();
extern int ps_rub(float h1, float h2, float v1, float v2); 
extern int ps_rect(float h1, float h2, float v1, float v2);
extern int ps_polyline(int na, float ha[], float va[]);
extern int ps_polygon(int na, float ha[], float va[]);
extern int ps_text(float h, float v, char str[], char just[2]);
extern int ps_text_box(char str[], float *dh, float *dv);
extern int ps_cursor(float *h, float *v, int *button);
extern int ps_sample(float *h, float *v, int *button);
extern int ps_set_font();
extern int ps_set_colours(int offset, int size, float reds[], float greens[], float blues[]);
extern int ps_set_colour(int index);
extern int ps_set_rgb(float red, float green, float blue);
extern int ps_set_rgb_background(float red, float green, float blue);
extern int ps_set_line_width(float width);
extern int ps_colour_bitmap(int ih0,int iv0, int nh,int nv, int image[]);
extern int ps_rgb_bitmap(int ih0,int iv0, int nh,int nv, float rgbim[]);
extern int  ps_flush();
extern int ps_prepare();
extern int ps_setup_colour();
extern int ps_set_paint(int offset, float red, float green, float blue);
extern int x11_prepare();
extern int x11_open(char *name, char *type, struct dev *(*p), int *ios);
extern int x11_pos_r(float rh, float rv);
extern int x11_pos_a(float ah, float av);
extern int x11_size_r(float rh, float rv);
extern int x11_size_a(float ah, float av);
extern void nam_find(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);




/* SETUP */

int drp_prepare() {
   int i;

   for (i=0; i < MAXDEV; i++) {
      surface[i].used = FALSE;
   }
   ccopy("SOLID", 0, SZNAM, line_styles[0]);
#ifdef TEK
   tek_prepare();
#endif
#ifdef MX
   mx_prepare();
#endif
#ifdef SUNVIEW
   sun_prepare();
#endif
#ifdef PS
   ps_prepare();
#endif
#ifdef GL
   sgi_prepare();
#endif
#ifdef GLX
   glx_prepare();
#endif
#ifdef X11
   x11_prepare();
#endif
   return 0;
}

int drp_flush() {
   int i,ios;

   for (i=0; i < MAXDEV; i++) {
      if (surface[i].used) {
         drp_close(surface[i].name, &ios);
      }
   }
   return 0;
}

int drp_print() {
   int i;
   char c;

   printf("%2s %12s  %12s\n","No", "Name", "Type");
   for (i=0; i < MAXDEV; i++) {
      if (surface[i].used) {
         if (wc == i) {
            c = '@';
         } else {
            c = ' ';
         };
         printf("%2i %12s%c %12s\n",i,surface[i].name,c,
           surface[i].type);
      }
   }
   return 0;
}

int drp_avail() {
   printf("Available worksurface types:-\n");
#ifdef TEK
   printf("tek:  Tektronix 4010/4014 variants\n");
   printf("mg:   Pericom MG 200/600/700 series (tek)\n");
   printf("ln03: DEC LN03+ laser printer (tek)\n");
#endif
#ifdef MX
   printf("mx:   Pericom MX 7000 series 4107 mode\n");
#endif
#ifdef SUNVIEW
   printf("sun:  Sunview colour window\n");
   printf("sun_mono:  Sunview mono window\n");
   printf("sun_rgb:  Sunview rgb window\n");
#endif
#ifdef PS
   printf("ps:   PostScript\n");
   printf("ps_rgb:   Colour PostScript\n");
#endif
#ifdef GL
   printf("sgi:  Silicon Graphics GL (8bit)\n");
   printf("sgi_rgbz:  Silicon Graphics GL (z-buffer RGB)\n");
   printf("sgi_rgbs:  Silicon Graphics GL (8-bit RGB)\n");
   printf("sgi_4096:  Silicon Graphics GL (12-bit)\n");
#endif
#ifdef GLX
   printf("glx:  OpenGL/X (8bit)\n");
   printf("glx_rgbz:  OpenGL/X (z-buffer RGB)\n");
   printf("glx_rgbs:  OpenGL/X (8-bit RGB)\n");
   printf("glx_4096:  OpenGL/X (12-bit)\n");
#endif
#ifdef GLX
   printf("glx:  X11 (8bit)\n");
#endif
   return 0;
}

/* WORKSURFACE MANAGEMENT */

int drp_pos_r(rh, rv)
   float rh,rv;
{
   drp_pos_rh = rh;
   drp_pos_rv = rv;
   drp_pos_rel = TRUE;
   return 0;
}

int drp_pos_a(ah, av)
   float ah,av;
{
   drp_pos_ah = ah;
   drp_pos_av = av;
   drp_pos_rel = FALSE;
   return 0;
}

int drp_size_r(rh, rv)
   float rh,rv;
{
   drp_size_rh = rh;
   drp_size_rv = rv;
   drp_size_rel = TRUE;
   return 0;
}

int drp_size_a(ah, av)
   float ah,av;
{
   drp_size_ah = ah;
   drp_size_av = av;
   drp_size_rel = FALSE;
   return 0;
}

int drp_open(name, type, ios)
   char *name,*type;
   int *ios;
{
   int i,wd;

   *ios = OK;
   wd = -1;
   for (i=0; i < MAXDEV; i++) {
      if (*ios) {
      } else if (!surface[i].used) {
         if (wd < 0) {
            wd = i;
         };
      } else if (csame(name,surface[i].name)) {
         printf("Error:  duplicate surface named %s\n", name);
         *ios = ERR;
      }
   }
   if (*ios) {
   } else if (wd < 0) {
      printf("error:  no room in device table (drp_open)\n");
      *ios = ERR;
#ifdef TEK
   } else if (csame("tek",type)) {
      tek_open(name, type, &p, ios);
   } else if (csame("mg",type)) {
      tek_open(name, type, &p, ios);
   } else if (csame("ln03",type)) {
      tek_open(name, type, &p, ios);
#endif
#ifdef MX
   } else if (csame("mx",type)) {
      mx_open(name, "mx", &p, ios);
#endif
#ifdef SUNVIEW
   } else if (csame("sun",type) || csame("sun_mono",type) ||
     csame("sun_rgb",type) || csame("sun_rgbz",type)) {
      if (drp_pos_rel) {
         sun_pos_r(drp_pos_rh, drp_pos_rv);
      } else {
         sun_pos_a(drp_pos_ah, drp_pos_av);
      }
      if (drp_size_rel) {
         sun_size_r(drp_size_rh, drp_size_rv);
      } else {
         sun_size_a(drp_size_ah, drp_size_av);
      }
      sun_open(name, type, &p, ios);
#endif
#ifdef PS
   } else if (csame("ps",type) || csame("ps_rgb",type)) {
      if (drp_size_rel) {
         ps_size_r(drp_size_rh, drp_size_rv);
      } else {
         ps_size_a(drp_size_ah, drp_size_av);
      }
      ps_open(name, type, &p, ios);
#endif
#ifdef GL
   } else if (csame("sgi",type)||csame("sgi_rgbz",type) ||
     csame("sgi_4096",type) || csame("sgi_rgbs",type)) {
      if (drp_pos_rel) {
         sgi_pos_r(drp_pos_rh, drp_pos_rv);
      } else {
         sgi_pos_a(drp_pos_ah, drp_pos_av);
      }
      if (drp_size_rel) {
         sgi_size_r(drp_size_rh, drp_size_rv);
      } else {
         sgi_size_a(drp_size_ah, drp_size_av);
      }
      sgi_open(name, type, &p, ios);
#endif
#ifdef GLX
   } else if (csame("glx",type)||csame("glx_rgbz",type) ||
     csame("glx_4096",type) || csame("glx_rgbs",type)) {
      if (drp_pos_rel) {
         glx_pos_r(drp_pos_rh, drp_pos_rv);
      } else {
         glx_pos_a(drp_pos_ah, drp_pos_av);
      }
      if (drp_size_rel) {
         glx_size_r(drp_size_rh, drp_size_rv);
      } else {
         glx_size_a(drp_size_ah, drp_size_av);
      }
      glx_open(name, type, &p, ios);
#endif
#ifdef X11
   } else if (csame("x",type) ||
     csame("x_my_cmap",type) || csame("x_my_hsv",type) ||
     csame("x_my_hv",type) || csame("x_my_rgbs",type) ||
     csame("x_def_cmap",type) || csame("x_def_hsv",type) ||
     csame("x_def_hv",type) || csame("x_def_rgbs",type) ||
     csame("x_clone_cmap",type) || csame("x_clone_hsv",type) ||
#define TRUE_TECRA16
#ifdef  TRUE_TECRA16
     csame("x_tecra",type) ||
#endif
     csame("x_clone_hv",type) || csame("x_clone_rgbs",type)) {
      if (drp_pos_rel) {
         x11_pos_r(drp_pos_rh, drp_pos_rv);
      } else {
         x11_pos_a(drp_pos_ah, drp_pos_av);
      }
      if (drp_size_rel) {
         x11_size_r(drp_size_rh, drp_size_rv);
      } else {
         x11_size_a(drp_size_ah, drp_size_av);
      }
      x11_open(name, type, &p, ios);
#endif
   } else {
      printf("error:  unknown type=%s (drp_open)\n",type);
      *ios = ERR;
   }
   if (*ios) {
      printf("error:  failed to open worksurface (drp_open)\n");
      if (wc >= 0) {
         p = wp->p;
      }
   } else {
      wc = wd;
      wp = &surface[wc];
      wp->used = TRUE;
      ccopy(name, 0, SZNAM, wp->name);
      ccopy(type, 0, SZNAM, wp->type);
      wp->p = p;
      wp->clear = TRUE;
      wp->marker_size = marker_scale*(p->char_height);
      wp->line_size = line_scale*(p->char_height);
/*    keep copy of original coordinate limits */
      wp->hmin = p->hmin;
      wp->hmax = p->hmax;
      wp->vmin = p->vmin;
      wp->vmax = p->vmax;
/*    copy bitmap equivalents too */
      wp->ihmin = p->ihmin;
      wp->ihmax = p->ihmax;
      wp->ivmin = p->ivmin;
      wp->ivmax = p->ivmax;
/*    setup clipping data */
      drp_setup_clip();
/*    zap colour ramps and palets */
      if (p->count == 1) {
         drp_zap_colour();
      }
   }
   return 0;
}

int drp_select(name, ios)
   char name[];
   int *ios;
{
   int i,wd;

   wd = -1;
   for (i=0; i < MAXDEV; i++) {
      if (surface[i].used) {
         if (csame(name, surface[i].name)) {
            wd = i;
         }
      }
   }
   if (wd < 0) {
      printf("error:  unknown worksurface named %s\n", name);
      *ios = ERR;
   } else {
/*    worksurface number */
      wc = wd;
/*    worksurface pointer */
      wp = &surface[wc];
/*    device structure pointer */
      p = wp->p;
/*    restore original (total viewport) coordinate limits */
      p->hmin = wp->hmin;
      p->hmax = wp->hmax;
      p->vmin = wp->vmin;
      p->vmax = wp->vmax;
/*    restore bitmap equivalents too */
      p->ihmin = wp->ihmin;
      p->ihmax = wp->ihmax;
      p->ivmin = wp->ivmin;
      p->ivmax = wp->ivmax;
/*    establish clipping data */
      drp_setup_clip();
#ifdef REPAIR
/*    invoke "gflush" to implement any required window repair */
      drp_gflush();
#endif
      *ios = OK;
      (*(p->select))(name, ios);
   }
   return 0;
}

int drp_close(name, ios)
   char name[];
   int *ios;
{
   int i,wd;

   wd = -1;
   for (i=0; i < MAXDEV; i++) {
      if (surface[i].used) {
         if (csame(name, surface[i].name)) {
            wd = i;
         }
      }
   }
   if (wd < 0) {
      printf("error:  unknown worksurface named %s\n", name);
      *ios = ERR;
   } else {
      (*((surface[wd].p)->close))(name, ios);
      surface[wd].used = FALSE;
      if (wd == wc) {
         wc = -1;
      }
      *ios = OK;
   }
   return 0;
}

int drp_get_count(count)
   int *count;
{
   *count = p->count;
   return 0;
}

int drp_active(active)
   bool *active;
{
   *active = (wc >= 0);
   return 0;
}

int drp_get_name(name)
   char name[];
{
   if (wc < 0) {
      ccopy("?", 0, SZNAM, name);
   } else {
      ccopy((p->name),0,SZNAM,name);
   }
   return 0;
}

int drp_get_type(type)
   char type[];
{
   if (wc < 0) {
      ccopy("?", 0, SZNAM, type);
   } else {
      ccopy((p->type),0,SZNAM,type);
   }
   return 0;
}

/* VIEWPORTS */

int drp_set_hv(h1, h2, v1, v2)
   float h1,h2,v1,v2;
{
/* take copies */
   p->hmin = h1;
   p->hmax = h2;
   p->vmin = v1;
   p->vmax = v2;
/* clip to maximal viewport */
   if (p->hmin < wp->hmin) {
      p->hmin = wp->hmin;
   }
   if (p->hmax > wp->hmax) {
      p->hmax = wp->hmax;
   }
   if (p->vmin < wp->vmin) {
      p->vmin = wp->vmin;
   }
   if (p->vmax > wp->vmax) {
      p->vmax = wp->vmax;
   }
/* avoid null extents (zap) */
   if ((p->hmin >= p->hmax) || (p->vmin >= p->vmax)) {
      printf("error:  null viewport extent\n");
      drp_zap_hv();
   } else {
      drp_setup_hv();
   }
   return 0;
}

int drp_zap_hv() {
   p->hmin = wp->hmin;
   p->hmax = wp->hmax;
   p->vmin = wp->vmin;
   p->vmax = wp->vmax;
   p->ihmin = wp->ihmin;
   p->ihmax = wp->ihmax;
   p->ivmin = wp->ivmin;
   p->ivmax = wp->ivmax;
   drp_setup_clip();
   return 0;
}

int drp_get_hv(h1,h2,v1,v2)
   float *h1,*h2,*v1,*v2;
{
   *h1 = p->hmin;
   *h2 = p->hmax;
   *v1 = p->vmin;
   *v2 = p->vmax;
   return 0;
}

int drp_setup_hv() {
   p->ihmin = (wp->ihmax - wp->ihmin)*(p->hmin - wp->hmin)/
     (wp->hmax - wp->hmin) + wp->ihmin;
   if (p->ihmin < wp->ihmin) {
      p->ihmin = wp->ihmin;
   } else if (p->ihmin > wp->ihmax) {
      p->ihmin = wp->ihmax;
   }
   p->ihmax = (wp->ihmax - wp->ihmin)*(p->hmax - wp->hmin)/
     (wp->hmax - wp->hmin) + wp->ihmin;
   if (p->ihmax < wp->ihmin) {
      p->ihmax = wp->ihmin;
   } else if (p->ihmax > wp->ihmax) {
      p->ihmax = wp->ihmax;
   }
   p->ivmin = (wp->ivmax - wp->ivmin)*(p->vmin - wp->vmin)/
     (wp->vmax - wp->vmin) + wp->ivmin;
   if (p->ivmin < wp->ivmin) {
      p->ivmin = wp->ivmin;
   } else if (p->ivmin > wp->ivmax) {
      p->ivmin = wp->ivmax;
   }
   p->ivmax = (wp->ivmax - wp->ivmin)*(p->vmax - wp->vmin)/
     (wp->vmax - wp->vmin) + wp->ivmin;
   if (p->ivmax < wp->ivmin) {
      p->ivmax = wp->ivmin;
   } else if (p->ivmax > wp->ivmax) {
      p->ivmax = wp->ivmax;
   }
   drp_setup_clip();
   return 0;
}

int drp_setup_clip() {
   drp_hc = 0.5*(p->hmin + p->hmax);
   drp_vc = 0.5*(p->vmin + p->vmax);
   drp_dh = 0.5*(p->hmax - p->hmin);
   drp_dv = 0.5*(p->vmax - p->vmin);
   drp_hmin = p->hmin;
   drp_hmax = p->hmax;
   drp_vmin = p->vmin;
   drp_vmax = p->vmax;
   return 0;
}

/* ERASURE */

int drp_erase(){
   wp->clear = TRUE;
   (*(p->erase))();
   return 0;
}

int drp_rub(h1, h2, v1, v2)
   float h1,h2,v1,v2;
{
   (*(p->rub))(h1,h2,v1,v2);
   return 0;
}

int drp_avail_rub(avail)
   bool *avail;
{
   *avail = p->avail_rub;
   return 0;
}

/* LINE DRAWING */

int drp_line(h1, v1, h2, v2)
   float h1,v1,h2,v2;
{
/*
   (*(p->line))(h1, v1, h2, v2);
*/
   drp_bp();
   drp_ap(h1, v1);
   drp_ap(h2, v2);
   drp_ep();
   return 0;
}

int drp_polyline(na, ha, va)
   int na;
   float ha[],va[];
{
   (*(p->polyline))(na, ha, va);
   return 0;
}

int drp_set_line_width(line_width)
   float line_width;
{
   (*(p->set_line_width))(line_width);
   return 0;
}

int drp_get_line_width(line_width)
   float *line_width;
{
   *line_width = p->line_width;
   return 0;
}

int drp_set_line_style(style)
   int style;
{
   if (style < 0) {
   } else if (style > (MAXLINESTYLE-1)) {
   } else {
      line_style = style;
   }
   return 0;
}

int drp_get_line_styles(n)
   int *n;
{
   *n = MAXLINESTYLE;
   return 0;
}

int drp_get_line_style(style)
   int *style;
{
   *style = line_style;
   return 0;
}

int drp_bp() {
   np = 0;
   return 0;
}

int drp_ap(h, v)
   float h,v;
{
   bool clip1,clip2,empty;
   float h1,h2,v1,v2;

   if (np == 0) {
      hp[np] = h;
      vp[np] = v;
      np++;
   } else {
      if (np == MAXP) {
         drp_polyline(np, hp, vp);
         np = 0;
         hp[np] = hp[MAXP-1];
         vp[np] = vp[MAXP-1];
         np++;
      }
      h1 = hp[np-1];
      v1 = vp[np-1];
      h2 = h;
      v2 = v;
      drp_clip(&h1, &v1, &clip1, &h2, &v2, &clip2, &empty);
      if (empty) {
         np = 0;
      } else {
         if (clip1) {
            np = 0;
            hp[np] = h1;
            vp[np] = v1;
            np++;
         }
         hp[np] = h2;
         vp[np] = v2;
         np++;
         if (clip2) {
            drp_polyline(np, hp, vp);
            np = 0;
            hp[np] = h;
            vp[np] = v;
            np++;
         }
      }
   }
   return 0;
}

int drp_ep() {
   if (np > 1) {
      drp_polyline(np, hp, vp);
   }
   np = 0;
   return 0;
}

int drp_clip(hh1, vv1, clip1, hh2, vv2, clip2, empty)
   float *hh1,*vv1,*hh2,*vv2;
   bool *clip1,*clip2,*empty;
{
   int ih1,ih2,iv1,iv2;
   float h1,h2,v1,v2,hh,vv;

/* make local copies */
   h1 = *hh1;
   v1 = *vv1;
   h2 = *hh2;
   v2 = *vv2;
/* find modes */
   drp_mode(drp_hmin, drp_hmax, h1, &ih1);
   drp_mode(drp_hmin, drp_hmax, h2, &ih2);
   drp_mode(drp_vmin, drp_vmax, v1, &iv1);
   drp_mode(drp_vmin, drp_vmax, v2, &iv2);
/* is line entirely outside viewport? */
   if (((ih1*ih2) == 1) || ((iv1*iv2) == 1)) {
      *empty = TRUE;
/* no... */
   } else {
      *empty = FALSE;
      *clip1 = FALSE;
      if (iv1 != 0) {
         vv = drp_dv*(float)iv1 + drp_vc;
         h1 = h1 + (h2-h1)*(vv-v1)/(v2-v1);
         v1 = vv;
         drp_mode(drp_hmin, drp_hmax, h1, &ih1);
         *clip1 = TRUE;
      }
      if (ih1 != 0) {
         hh = drp_dh*(float)ih1 + drp_hc;
         v1 = v1 + (v2-v1)*(hh-h1)/(h2-h1);
         h1 = hh;
         *clip1 = TRUE;
      }
      *clip2 = FALSE;
      if (iv2 != 0) {
         vv = drp_dv*(float)iv2 + drp_vc;
         h2 = h1 + (h2-h1)*(vv-v1)/(v2-v1);
         v2 = vv;
         drp_mode(drp_hmin, drp_hmax, h2, &ih2);
         *clip2 = TRUE;
      }
      if (ih2 != 0) {
         hh = drp_dh*(float)ih2 + drp_hc;
         v2 = v1 + (v2-v1)*(hh-h1)/(h2-h1);
         h2 = hh;
         *clip2 = TRUE;
      }
/*    copy back to caller */
      *hh1 = h1;
      *vv1 = v1;
      *hh2 = h2;
      *vv2 = v2;
   }
   return 0;
}

int drp_mode(cmin, cmax, c, ic)
   float cmin,cmax,c;
   int *ic;
{
   if (c < cmin) {
      *ic = -1;
   } else if (c <= cmax) {
      *ic = 0;
   } else {
      *ic = 1;
   }
   return 0;
}

/* RECTANGLE FILL */

int drp_rect(h1, h2, v1, v2)
   float h1,h2,v1,v2;
{
   float ha,hb,va,vb;

   if (h1 > drp_hmax) {
   } else if (h2 < drp_hmin) {
   } else if (v1 > drp_vmax) {
   } else if (v2 < drp_vmin) {
   } else {
      if (h1 < drp_hmin) {
         ha = drp_hmin;
      } else {
         ha = h1;
      }
      if (h2 > drp_hmax) {
         hb = drp_hmax;
      } else {
         hb = h2;
      }
      if (v1 < drp_vmin) {
         va = drp_vmin;
      } else {
         va = v1;
      }
      if (v2 > drp_vmax) {
         vb = drp_vmax;
      } else {
         vb = v2;
      }
      (*(p->rect))(ha, hb, va, vb);
   }
   return 0;
}

int drp_avail_rect(avail)
   bool *avail;
{
   *avail = p->avail_rect;
   return 0;
}

/* POLYGON FILL */

int drp_polygon(na, ha, va)
   int na;
   float ha[],va[];
{
   (*(p->polygon))(na, ha, va);
   return 0;
}

int drp_avail_polygon(avail)
   bool *avail;
{
   *avail = p->avail_polygon;
   return 0;
}

int drp_bg() {
   ng = 0;
   return 0;
}

int drp_ag(h, v)
   float h,v;
{
   if (ng >= MAXP) {
      drp_polygon(ng, hg, vg);
      ng = 1;
      hg[ng-1] = hg[MAXP-1];
      vg[np-1] = vg[MAXP-1];
   }
   ng = ng+1;
   if (h < drp_hmin) {
      hg[ng-1] = drp_hmin;
   } else if (h <= drp_hmax) {
      hg[ng-1] = h;
   } else {
      hg[ng-1] = drp_hmax;
   }
   if (v < drp_vmin) {
      vg[ng-1] = drp_vmin;
   } else if (v <= drp_vmax) {
      vg[ng-1] = v;
   } else {
      vg[ng-1] = drp_vmax;
   }
   return 0;
}

int drp_eg() {
   if (ng > 2) {
      drp_polygon(ng, hg, vg);
   }
   ng = 0;
   return 0;
}

/* TEXT */

int drp_text(h, v, str, just)
   float h,v;
   char str[],just[2];
{
   float dh,dv,h0,v0;

   drp_text_box(str, &dh, &dv);
   if (just[0] == 'r') {
      h0 = h - dh;
   } else if (just[0] == 'm') {
      h0 = h - 0.5*dh;
   } else {
      h0 = h;
   }
   if (just[1] == 't') {
      v0 = v - dv;
   } else if (just[1] == 'm') {
      v0 = v - 0.5*dv;
   } else {
      v0 = v;
   }
   if (h0 > drp_hmax) {
   } else if ((h0+dh) < drp_hmin) {
   } else if (h0 < drp_hmin) {
   } else if ((h0+dh) > drp_hmax) {
   } else if (v0 > drp_vmax) {
   } else if ((v0+dv) < drp_vmin) {
   } else if (v0 < drp_vmin) {
   } else if ((v0+dv) > drp_vmax) {
   } else {
      (*(p->text))(h, v, str, just);
   }
   return 0;
}

int drp_text_box(str, dh, dv)
   char str[];
   float *dh,*dv;
{
   (*(p->text_box))(str, dh, dv);
   return 0;
}

int drp_set_font_style(fontstyle)
   char fontstyle[];
{
   ccopy(fontstyle, 0, SZNAM, p->font_style);
   return 0;
}

int drp_get_font_style(fontstyle)
   char fontstyle[];
{
   ccopy(p->font_style, 0, SZNAM, fontstyle);
   return 0;
}

int drp_set_font_size(fontsize)
   float fontsize;
{
   p->font_size = fontsize;
   return 0;
}

int drp_get_font_size(fontsize)
   float *fontsize;
{
   *fontsize = p->font_size;
   return 0;
}

int drp_avail_fontscale(avail)
   bool *avail;
{
   *avail = p->avail_fontscale;
   return 0;
}

int drp_set_font() {
   (*(p->set_font))();
   return 0;
}

int drp_get_char(height,width,depth)
   float *height,*width,*depth;
{
   *height = p->char_height;
   *width = p->char_width;
   *depth = p->char_depth;
   return 0;
}

/* CURSOR */

int drp_avail_cursor(avail)
   bool *avail;
{
   *avail = p->avail_cursor;
   return 0;
}

int drp_cursor(h, v, button)
   float *h,*v;
   int *button;
{
   (*(p->cursor))(h, v, button);
   return 0;
}

/* COLOUR PRIMITIVES */

/* drp_avail_colour - determine whether colourmap mode is available */
int drp_avail_colour(avail)
   bool *avail;
{
   *avail = p->avail_colour;
   return 0;
}

/* drp_set_colours - set RGB for colourmap */
int drp_set_colours(offset, size, reds, greens, blues)
   int offset,size;
   float reds[],greens[],blues[];
{
   (*(p->set_colours))(offset, size, reds, greens, blues);
   return 0;
}

/* drp_set_colour - set colourmap index for drawing */
int drp_set_colour(colour_index)
   int colour_index;
{
   (*(p->set_colour))(colour_index);
   return 0;
}

/* drp_change_colour - can we change colourmap entries? */
int drp_change_colour(change)
   bool *change;
{
   *change = p->change_colour;
   return 0;
}

/* drp_get_colour_shape - get non-special part of colourmap */
int drp_get_colour_shape(offset, size)
   int *offset,*size;
{
   *offset = p->colour_offset;
   *size = p->colour_size;
   return 0;
}

/* drp_get_colour - get colour index for drawing */
int drp_get_colour(colour_index)
   int *colour_index;
{
   *colour_index = p->colour_index;
   return 0;
}

/* drp_get_colours - get RGB for colourmap */
int drp_get_colours(offset, size, reds, greens, blues)
   int offset,size;
   float reds[],greens[],blues[];
{
   int i,k;

   k = 0;
   for (i=offset; i < size; i++) {
      reds[k] = (p->reds)[i];
      greens[k] = (p->greens)[i];
      blues[k] = (p->blues)[i];
      k++;
   }
   return 0;
}

/* RGB PRIMITIVES */

int drp_avail_rgb(avail)
   bool *avail;
{
   *avail = p->avail_rgb;
   return 0;
}

int drp_set_rgb(red, green, blue)
   float red,green,blue;
{
   (*(p->set_rgb))(red, green, blue);
   return 0;
}

int drp_set_rgb_background(red, green, blue)
   float red,green,blue;
{
   (*(p->set_rgb_background))(red, green, blue);
   return 0;
}

int drp_get_rgb(red, green, blue)
   float *red,*green,*blue;
{
  *red = p->red;
  *green = p->green;
  *blue = p->blue;
   return 0;
}

int drp_get_rgb_background(red_background,green_background,blue_background)
   float *red_background,*green_background,*blue_background;
{
  *red_background = p->red_background;
  *green_background = p->green_background;
  *blue_background = p->blue_background;
   return 0;
}

/* RGBZ PRIMITIVES */

int drp_avail_rgbz(avail)
   bool *avail;
{
   *avail = p->avail_rgbz;
   return 0;
}

/* rgbz_setup - establish the 3d projection etc. */
int rgbz_setup(anglev, near, far, hmin, hmax, vmin, vmax)
   float anglev,near,far,hmin,hmax,vmin,vmax;
{
   if (p->avail_rgbz) {
      (*(p->rgbz_setup))(anglev, near, far, hmin, hmax, vmin, vmax);
   }
   return 0;
}

/* rgbz_begin - begin 3d projection */
int rgbz_begin() {
   if (p->avail_rgbz) {
      (*(p->rgbz_begin))();
   }
   return 0;
}

/* rgbz_end - end 3d projection */
int rgbz_end() {
   if (p->avail_rgbz) {
      (*(p->rgbz_end))();
   }
   return 0;
}

/* rgbz_clear - clear the 3d region */
int rgbz_clear() {
   if (p->avail_rgbz) {
      (*(p->rgbz_clear))();
   }
   return 0;
}

int rgbz_begin_polygon() {
   if (p->avail_rgbz) {
      (*(p->rgbz_begin_polygon))();
   }
   return 0;
}

int rgbz_polygon_vertex(x, y, z)
   float x,y,z;
{
   if (p->avail_rgbz) {
      (*(p->rgbz_polygon_vertex))(x, y, z);
   }
   return 0;
}

int rgbz_end_polygon() {
   if (p->avail_rgbz) {
      (*(p->rgbz_end_polygon))();
   }
   return 0;
}

int rgbz_begin_polyline() {
   if (p->avail_rgbz) {
      (*(p->rgbz_begin_polyline))();
   }
   return 0;
}

int rgbz_polyline_vertex(x, y, z)
   float x,y,z;
{
   if (p->avail_rgbz) {
      (*(p->rgbz_polyline_vertex))(x, y, z);
   }
   return 0;
}

int rgbz_end_polyline() {
   if (p->avail_rgbz) {
      (*(p->rgbz_end_polyline))();
   }
   return 0;
}

int rgbz_put_bitmap(ih0,iv0,nh,nv,rgb_image,d_image)
   int ih0,iv0,nh,nv;
   float rgb_image[],d_image[];
{
   (*(p->rgbz_put_bitmap))(ih0,iv0,nh,nv,rgb_image,d_image);
   return 0;
}

int rgbz_get_bitmap(ih0,iv0,nh,nv,rgb_image,d_image)
   int ih0,iv0,nh,nv;
   float rgb_image[],d_image[];
{
   (*(p->rgbz_get_bitmap))(ih0,iv0,nh,nv,rgb_image,d_image);
   return 0;
}

/* MARKERS */

int drp_set_marker_scale(scale)
   float scale;
{
   marker_scale = scale;
   wp->marker_size = marker_scale*(p->char_height);
   return 0;
}

int drp_set_marker_size(marker_size)
   float marker_size;
{
   wp->marker_size = marker_size;
   return 0;
}

int drp_get_marker_scale(scale)
   float *scale;
{
   *scale = marker_scale;
   return 0;
}

int drp_get_marker_size(marker_size)
   float *marker_size;
{
   *marker_size = wp->marker_size;
   return 0;
}

int drp_marker(marker, h, v)
  int marker;
   float h,v;
{
   int is,j,temp_line_style;
   float hm,vm,mksiz;

   if (marker < 1) {
   } else if (marker > MAXMARKER) {
   } else {
      temp_line_style = line_style;
      line_style = 0;
      mksiz = wp->marker_size;
      for (is=is0[marker-1]; is <(is1[marker-1]+1); is++) {
         drp_bp();
         for (j=ip0[is]; j < (ip1[is]+1); j++) {
            hm = h + mksiz*(float)dh[j];
            vm = v + mksiz*(float)dv[j];
            drp_ap(hm, vm);
         }
         drp_ep();
      }
      line_style = temp_line_style;
   }
   return 0;
}

/* BITMAP PRIMITIVES */

int drp_avail_bitmap(avail)
   bool *avail;
{
   *avail = p->avail_bitmap;
   return 0;
}

int drp_bitmap(nh,nv,dh,dv,h0,v0)
   int *nh,*nv;
   float *dh,*dv,*h0,*v0;
{
   *nh = (wp->ihmax - wp->ihmin) + 1;
   *nv = (wp->ivmax - wp->ivmin) + 1;
/*
   *nh = (*nh)/4;
   *nh = (*nh)*4;
*/
   *dh = p->bitmap_dh;
   *dv = p->bitmap_dv;
   *h0 = wp->hmin;
   *v0 = wp->vmin;
   return 0;
}

int drp_get_bitmap(ih0,iv0,nh,nv,dhidh,dvidv,hi0,vi0)
   int *nh,*nv,*ih0,*iv0;
   float *dhidh,*dvidv,*hi0,*vi0;
{
   *ih0 = p->ihmin;
   *iv0 = p->ivmin;
/*
   *nh = (p->ihmax - p->ihmin)/4;
   *nh = *nh*4;
   *nv = (p->ivmax - p->ivmin)/4;
   *nv = *nv*4;
*/
   *nh = p->ihmax - p->ihmin;
   *nv = p->ivmax - p->ivmin;
   *dhidh = 1.0/p->bitmap_dh;
   *dvidv = 1.0/p->bitmap_dv;
   *hi0 = (-(p->hmin/p->bitmap_dh));
   *vi0 = (-(p->vmin/p->bitmap_dv));
   return 0;
}

int drp_colour_bitmap(ih0,iv0,nh,nv,image)
   int ih0,iv0,nh,nv,image[];
{
   (*(p->colour_bitmap))(ih0,iv0,nh,nv,image);
   return 0;
}

int drp_rgb_bitmap(ih0,iv0,nh,nv,rgb_image)
   int ih0,iv0,nh,nv;
   float rgb_image[];
{
   (*(p->rgb_bitmap))(ih0,iv0,nh,nv,rgb_image);
   return 0;
}

/* COLOUR RAMP MANAGER */

int drp_zap_colour() {
   drp_zap_paints();
   drp_zap_ramps();
   return 0;
}

int drp_zap_ramps() {
   int i;

   p->ramp_top = p->colour_offset;
   p->ramp_count = 0;
   for (i=0; i < MAXRAMP; i++) {
      p->ramp_used[i] = FALSE;
   }
   return 0;
}

int drp_print_ramps() {
   int i,id;

   printf("Colour ramps:-\n");
   for (i=0; i < p->ramp_count; i++) {
      printf("  i=%d,name=%s,nd=%d,nc=", i, p->ramp_names[i],
        p->ramp_nd[i]);
      for (id=0; id < p->ramp_nd[i]; id++) {
         printf("%d", p->ramp_nc[i][id]);
         if (id == (p->ramp_nd[i]-1)) {
            printf("],");
         } else {
            printf(",");
         }
      }
      printf("offset=%d,size=%d\n", p->ramp_offset[i], p->ramp_size[i]);
   }
   return 0;
}

int drp_put_ramp(name, nd, nc, offset)
   char name[];
   int nd,nc[],*offset;
{
   int i,id,space,size;

   nam_find(name, p->ramp_count, p->ramp_used, p->ramp_names, &i);
   if (i < 0) {
      space = p->paint_top - p->ramp_top;
      size = 1;
      for (id=0; id < nd; id++) {
         size = size*nc[id];
      }
      if (size > space) {
         *offset = (-1);
      } else if (p->ramp_count == MAXRAMP) {
         *offset = (-1);
      } else {
         i = p->ramp_count;
         p->ramp_size[i] = size;
         *offset = p->ramp_top;
         p->ramp_offset[i] = *offset;
         ccopy(name, 0, SZNAM, p->ramp_names[i]);
         p->ramp_top = p->ramp_top + size;
         p->ramp_nd[i] = nd;
         for (id=0; id < nd; id++) {
            p->ramp_nc[i][id] = nc[id];
         }
         p->ramp_used[i] = TRUE;
         p->ramp_count = p->ramp_count + 1;
      }
   } else {
      *offset = p->ramp_offset[i];
   }
   return 0;
}

int drp_get_ramp(name, nd, nc, offset)
   char name[];
   int *nd,nc[],*offset;
{
   int i,id;

   nam_find(name, p->ramp_count, p->ramp_used, p->ramp_names, &i);
   if (i < 0) {
      *offset = i;
   } else {
      *nd = p->ramp_nd[i];
      *offset = p->ramp_offset[i];
      for (id=0; id <  *nd; id++) {
         nc[id] = p->ramp_nc[i][id];
      }
   }
   return 0;
}

int drp_find_ramp(nd, name)
   int nd;
   char name[];
{
   int i;

   ccopy(" ", 0, SZNAM, name);
   for (i=0; i < p->ramp_count; i++) {
     if (p->ramp_nd[i] == nd) {
         ccopy(name, 0, SZNAM, p->ramp_names[i]);
         break;
      }
   }
   return 0;
}

/* COLOUR PAINT MANAGER */

int drp_zap_paints() {
   int i;

   p->paint_top = p->colour_size + p->colour_offset;
   p->paint_count = 0;
   for (i=0; i < MAXPAINT; i++) {
      p->paint_used[i] = FALSE;
   }
   return 0;
}

int drp_put_paint(name, offset)
   char name[];
   int *offset;
{
   int i;

   nam_find(name, p->paint_count, p->paint_used, p->paint_names, &i);
   if (i < 0) {
      if ((p->paint_top-1) < p->ramp_top) {
         *offset = (-1);
      } else if (p->paint_count == MAXPAINT) {
         *offset = (-1);
      } else {
         *offset = p->paint_top-1;
         i = p->paint_count;
         p->paint_offset[i] = *offset;
         ccopy(name, 0, SZNAM, p->paint_names[i]);
         p->paint_used[i] = TRUE;
         p->paint_top = p->paint_top - 1;
         p->paint_count = p->paint_count + 1;
      }
   } else {
      *offset = p->paint_offset[i];
   }
   return 0;
}

int drp_get_paint(name, offset)
   char name[];
   int *offset;
{
   int i;

   nam_find(name, p->paint_count, p->paint_used, p->paint_names, &i);
   if (i < 0) {
      *offset = i;
   } else {
      *offset = p->paint_offset[i];
   }
   return 0;
}

int drp_set_paint(offset, red, green, blue)
   int offset;
   float red, green, blue;
{
   float reds[1],greens[1],blues[1];

   reds[0] = red;
   greens[0] = green;
   blues[0] = blue;
   drp_set_colours(offset, 1, reds, greens, blues);
   return 0;
}

int drp_print_paints() {
   int i;

   printf("Colour paints:-\n");
   for (i=0; i < p->paint_count; i++) {
      printf("  i=%d,name=%s,", i, p->paint_names[i]);
      printf("offset=%d\n", p->paint_offset[i]);
   }
   return 0;
}

int drp_get_colour_space(space)
   int *space;
{
   *space = p->paint_top - p->ramp_top;
   return 0;
}

int crap() {
   float a,b,c;

   c = 0.0;
   b = 1.0;
   printf("crap out\n");
   a = b/c;
   return a;
}

/* GFLUSH */

int drp_gflush() {
   (*(p->gflush))();
   return 0;
}

