/* x11.c:  Vanilla X11 backend */

/* Includes :
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <time.h>
#include "drp.h"
#include "machine.h"

/* Types :
*/
#define bool int
#define TRUE24
#undef DEBUG_TRUE24
#define TRUE_TECRA16

/* Constants :
*/
#define RANDY
#define DEFAULT_FONT "fixed"
#undef DEBUG
#define WM_DELETE
#define REPAIR
#define MAXWIND 12
#define MAXPIXEL 4096
#define MAXRAND 32765

/* Current and general X data :
*/
static int screen;
static Display *display;
static XFontStruct *font;
static char font_name[SZLIN];
static int npixel;
static unsigned long pixels[MAXPIXEL];
static bool my_cmap,def_cmap,clone_cmap;
#ifdef TRUE24
static bool true_colour;
static unsigned int x11_red_offset,x11_green_offset,x11_blue_offset;
static int x11_background_pixel;
#ifdef TRUE_TECRA16
static float x11_red_fsd,x11_green_fsd,x11_blue_fsd;
#endif
#endif
static Window window;
static Pixmap pix;
static GC gc;
static Visual *visual;
static Colormap cmap,defmap;
static int cmap_size;
static int cmap_offset;
static XGCValues gcv;
static unsigned long gc_mask;
//static XEvent event;			// not used BM
static Window rootwin;
static unsigned int width,height;


static struct dev x11;

static bool x11_avail_rgbs,x11_avail_hsv;
static int r_count,g_count,b_count,rgb_offset,h_count,s_count,
  v_count,hsv_offset,r_mult,g_mult,b_mult,h_mult,s_mult,v_mult,
  grey_offset;
static float r_space,g_space,b_space,h_space,s_space,v_space,
  grey_space;
static bool x11_pos_rel=TRUE,x11_size_rel=TRUE;
static float x11_pos_rh=1.0,x11_pos_rv=1.0,x11_size_rh=0.5,
  x11_size_rv=0.5, x11_pos_ah,x11_pos_av,x11_size_ah,x11_size_av;
static int x11_background_index;
static bool x11_prev_font;

/* Window-specific data :
*/
static int x11_count=0;
static bool x11_used[MAXWIND];
static char x11_names[MAXWIND][SZNAM],x11_type[SZNAM];
static int x11_ihmin[MAXWIND],x11_ihmax[MAXWIND],x11_ivmin[MAXWIND],
  x11_ivmax[MAXWIND];
static float x11_hmin[MAXWIND],x11_hmax[MAXWIND],x11_vmin[MAXWIND],
  x11_vmax[MAXWIND];
static Window x11_windows[MAXWIND];
static GC x11_gcs[MAXWIND];
static Pixmap  x11_pixs[MAXWIND];
static int x11_now;

/* Radom number data :
*/
static int norand=TRUE;
static float ranf=1.0,divran=(1.0/RANCON);
static float randaf[MAXRAND];
static int rand_count;


/* EXPLICIT DECLARATIONS ADDED BY BM */
int x11_pos_r(float rh, float rv);
int x11_pos_a(float ah, float av);
int x11_size_r(float rh, float rv);
int x11_size_a(float ah, float av);
void x11_await_map();
void x11_pend_events();
unsigned long x11_rgb_pixel(float r, float g, float b);
int x11_open(char *name, char *type, struct dev *(*p), int *ios);
int x11_close(char name[], int *ios);
int x11_select(char name[], int *ios);
int x11_erase();
int x11_rub(float h1, float h2, float v1, float v2);
int x11_rect(float h1, float h2, float v1, float v2);
int x11_polyline(int na, float ha[], float va[]);
int x11_polygon(int na, float ha[], float va[]);
int x11_text(float h, float v, char str[], char just[2]);
int x11_text_box(char str[], float *dh, float *dv);
int x11_cursor(float *h, float *v, int *button);
int x11_sample(float *h, float *v, int *button);
int x11_invent_font(char style[], float size, char name[]);
int x11_set_font();
int x11_set_colours(int offset, int size, float reds[], float greens[], float blues[]);
int x11_set_colour(int index);
int x11_set_rgb(float red, float green, float blue);
int x11_set_rgb_background(float red, float green, float blue);
int x11_set_line_width(float width);
int x11_colour_bitmap(int ih0,int iv0, int nh,int nv, int cim[]);
int x11_rgb_bitmap(int ih0,int iv0, int nh,int nv, float rgb_image[]);
int x11_rgbz_put_bitmap(int ih0,int iv0, int nh,int nv, float rgb_image[], float d_image[]);
int x11_rgbz_get_bitmap(int ih0,int iv0, int nh,int nv, float rgb_image[], float d_image[]);
int x11_rgbz_setup(float anglev, float near, float far, float hmin, float hmax, float vmin, float vmax);
int x11_rgbz_begin();
int x11_rgbz_end();
int x11_rgbz_clear();
int x11_rgbz_begin_polygon();
int x11_rgbz_polygon_vertex(float x, float y, float z);
int x11_rgbz_end_polygon();
int x11_rgbz_begin_polyline();
int x11_rgbz_polyline_vertex(float x, float y, float z);
int x11_rgbz_end_polyline();
int x11_flush();
int x11_prepare();
int x11_setup_window(int ixpmax,int iypmax,int ixmmax,int iymmax, float *h0,float *v0, float *hsize,float *vsize);
int x11_dump_colour(int offset, int size);
int x11_set_paint(int offset, float red, float green, float blue);
int x11_setup_rand();
int x11_mrs_videoc(char file[], int *ios);
int x11_wm_hints(int argc, char *argv[], char name[], int h, int v);
int x11_hsv_to_rgb(float hsv[3], float rgb[3]);
int x11_rgb_to_hsv(float rgb[3], float hsv[3]);
int x11_setup_hsv(int nhue, int nsat);
int x11_hsv_index(float hsv[3]);
int x11_setup_rgbs();
int x11_rgb_index(float rgb[3]);

extern void nam_find(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);
extern void nam_alloc(char name[], int maxnames, bool used[], 
  char names[][SZNAM], int *offset);
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);
extern int slen(char s[]);







int x11_pos_r(rh, rv)
   float rh,rv;
{
   x11_pos_rh = rh;
   x11_pos_rv = rv;
   x11_pos_rel = TRUE;
   return 0;
}

int x11_pos_a(ah, av)
   float ah,av;
{
   x11_pos_ah = ah;
   x11_pos_av = av;
   x11_pos_rel = FALSE;
   return 0;
}

int x11_size_r(rh, rv)
   float rh,rv;
{
   x11_size_rh = rh;
   x11_size_rv = rv;
   x11_size_rel = TRUE;
   return 0;
}

int x11_size_a(ah, av)
   float ah,av;
{
   x11_size_ah = ah;
   x11_size_av = av;
   x11_size_rel = FALSE;
   return 0;
}

void x11_await_map() {
   XEvent event;
   bool more=TRUE;

#ifdef DEBUG
   printf("x11_await_map:  called\n");
#endif
   while(more) {
      XNextEvent(display, &event);
      switch (event.type) {
         case Expose:
#ifdef DEBUG
            printf("expose event\n");
#endif
            break;
         case NoExpose:
#ifdef DEBUG
            printf("noexpose event (ignore)\n");
#endif
            break;
         case MapNotify:
#ifdef DEBUG
            printf("mapnotify event (expected)\n");
#endif
            more = FALSE;
            break;
         case ReparentNotify:
#ifdef DEBUG
            printf("reparentnotify event (ignore)\n");
#endif
            break;
         case UnmapNotify:
#ifdef DEBUG
            printf("unmap event (ignore)\n");
#endif
            break;
         case DestroyNotify:
#ifdef DEBUG
            printf("destroy event (ignore)\n");
#endif
            break;
         case ButtonPressMask:
#ifdef DEBUG
            printf("button event\n");
#endif
            break;
         case ConfigureNotify:
#ifdef DEBUG
            printf("configure event\n");
#endif
            break;
         default:
#ifdef DEBUG
            printf("unexpected event=%d\n", (int)event.type);
#else
            ;
#endif
      }
   }
}

void x11_pend_events() {
   XEvent event;
   int i,count,j;
#ifdef WM_DELETE
   Atom protocols_atom,delete_atom;
#endif

#ifdef DEBUG
   printf("x11_pend_events:  called\n");
#endif
   if ((count=XPending(display)) > 0) {
#ifdef DEBUG
      printf("x11_pend_events:  count=%d\n", count);
#endif
      for (i=0; i < count; i++) {
#ifdef DEBUG
         printf("x11_pend_events:  i=%d\n", i);
#endif
         XNextEvent(display, &event);
         switch (event.type) {
            case Expose:
#ifdef DEBUG
               printf("expose event\n");
#endif
#ifdef REPAIR
#ifdef DEBUG
               printf("x11_pend_events:  repairs window\n");
#endif
               for (j=0; j < MAXWIND; j++) {
                  if (x11_used[j]) {
                     if (event.xexpose.window == x11_windows[j]) {
                        XCopyArea(display, x11_pixs[j], x11_windows[j],
                          gc, 0, 0, width, height, 0, 0);
                     }
                  }
               }
#endif
               break;
            case NoExpose:
#ifdef DEBUG
               printf("noexpose event (ignore)\n");
#endif
               break;
#ifdef WM_DELETE
            case ClientMessage:
               protocols_atom = XInternAtom(display, "WM_PROTOCOLS",
                 FALSE);
               delete_atom = XInternAtom(display, "WM_DELETE_WINDOW",
                  FALSE);
               if ((event.xclient.message_type == protocols_atom) &&
                  (event.xclient.data.l[0] == delete_atom) &&
                  (event.xclient.format == 32)) {
                  fprintf(stderr, "warning:  reveived delete atom\n");
               } else {
                  fprintf(stderr, "warning:  weird client message\n");
               }
               break;
#endif
            case MapNotify:
#ifdef DEBUG
               printf("mapnotify event (expected)\n");
#endif
               break;
            case ReparentNotify:
#ifdef DEBUG
               printf("reparentnotify event (ignore)\n");
#endif
               break;
            case UnmapNotify:
#ifdef    DEBUG
               printf("unmap event (ignore)\n");
#endif
               break;
            case DestroyNotify:
#ifdef    DEBUG
               printf("destroy event (ignore)\n");
#endif
               break;
            case ButtonPressMask:
#ifdef DEBUG
               printf("button event\n");
#endif
               break;
            case ConfigureNotify:
#ifdef DEBUG
               printf("configure event\n");
#endif
#ifdef REPAIR
#ifdef DEBUG
               printf("x11_pend_events:  repairs window\n");
#endif
               for (j=0; j < MAXWIND; j++) {
                  if (x11_used[j]) {
                     if (event.xconfigure.window == x11_windows[j]) {
                        XCopyArea(display, x11_pixs[j], x11_windows[j],
                          gc, 0, 0, width, height, 0, 0);
                     }
                  }
               }
#endif
               break;
            default:
#ifdef DEBUG
               printf("unexpected event=%d\n", (int)event.type);
#else
               ;
#endif
         }
      }
   }
}

#ifdef TRUE24
/* x11_rgb_pixel:  calculate pixel from RGB compoents */
unsigned long x11_rgb_pixel(float r, float g, float b) {
   unsigned long pixel;
   int ir,ig,ib;

   r = (r > 1.0) ? 1.0 : r;
   r = (r < 0.0) ? 0.0 : r;
   g = (g > 1.0) ? 1.0 : g;
   g = (g < 0.0) ? 0.0 : g;
   b = (b > 1.0) ? 1.0 : b;
   b = (b < 0.0) ? 0.0 : b;
#ifdef TRUE_TECRA16
   ir = (unsigned long)(r*x11_red_fsd);
   ig = (unsigned long)(g*x11_green_fsd);
   ib = (unsigned long)(b*x11_blue_fsd);
#else
   ir = (unsigned long)(r*255.0);
   ig = (unsigned long)(g*255.0);
   ib = (unsigned long)(b*255.0);
#endif
   pixel = (unsigned long) (ir*x11_red_offset +ig*x11_green_offset  +
     ib*x11_blue_offset);
#ifdef DUMP
   printf("ir=%d,ig=%d,ib=%d,pixel=%6x\n",  ir, ig, ib,
     pixel);
#endif
   return(pixel);
}
#endif

int x11_open(name, type, p, ios)
   char *name,*type;
   struct dev *(*p);
   int *ios;
{
   int argc;
   char *argv[2];
   XSetWindowAttributes attributes;
   unsigned long attribute_mask;
   //XColor c;
   int depth,k,window_h,window_v;
   int ihlen,ivlen,ihoff,ivoff,i,nhue,nsat;
   nhue = 0;
   nsat = 0;
#ifdef WM_DELETE
   Atom protocols_atom,delete_atom;
#endif
   float h0,v0,hlen,vlen;
   //int gdtmp,mmtmp;
   static int ixpmax,iypmax,ixmmax,iymmax;

/* check name and allocate new one if unique */
   nam_find(name, MAXWIND, x11_used, x11_names, &i);
   if (i < 0) {
      nam_alloc(name, MAXWIND, x11_used, x11_names, &i);
      if (i < 0) {
         fprintf(stderr,"error:  too many x11 windows\n");
         *ios = ERR;
//         return;  MIR
         return 0;
      } else {
/*       check for a single type */
         if (x11_count > 0) {
            if (!csame(type,x11_type)) {
               fprintf(stderr,"error:  mixed x11 window types\n");
               *ios = ERR;
//               return;  MIR
               return 0;
            }
         } else {
            ccopy(type, 0, SZNAM, x11_type);
         }
         x11_now = i;
      }
   } else {
      fprintf(stderr,"error:  duplicate x11 window name\n");
      *ios = ERR;
//      return;  MIR
      return 0;
   }
/* connect to X server */
   if (x11_count == 0) {
/*    Determine colour model */
      def_cmap = FALSE;
      my_cmap = FALSE;
      clone_cmap = FALSE;
      x11_avail_rgbs = FALSE;
      x11_avail_hsv = FALSE;
      x11.avail_colour = FALSE;
      x11.change_colour = FALSE;
      x11.avail_rgbz = FALSE;
#ifdef TRUE24
/*    Open X display */
      display = XOpenDisplay((char*)NULL);
      if (display == NULL) {
         fprintf(stderr, "error:  opening X server display\n");
         *ios = ERR;
//         return; MIR
         return 0;
      }
      screen = DefaultScreen(display);
      rootwin = RootWindow(display,screen);
      visual = DefaultVisual(display,screen);
      true_colour = (visual->class==TrueColor);
      if (true_colour) {
         x11.avail_rgb = TRUE;
#ifdef CRAP24
         x11_red_offset = visual->red_mask >> visual->bits_per_rgb;
         x11_green_offset = visual->green_mask >> visual->bits_per_rgb;
         x11_blue_offset = visual->blue_mask >> visual->bits_per_rgb;
#endif
#ifdef DEBUG_TRUE24
         printf("red_mask=%08X\n", visual->red_mask, visual->red_mask);
         printf("green_mask=%08X\n", visual->green_mask, visual->green_mask);
         printf("blue_mask=%08X\n", visual->blue_mask, visual->blue_mask);
         printf("bits_per_rgb=%d\n", visual->bits_per_rgb, visual->bits_per_rgb);
#endif
#ifdef TRUE_TECRA16
         if (csame("x_tecra",type)) {
            x11_red_offset = 2048;
            x11_green_offset =  32;
            x11_blue_offset = 1;
            x11_red_fsd = 31.0;
            x11_green_fsd = 63.0;
            x11_blue_fsd = 31.0;
         } else {
            x11_red_offset=256*256;
            x11_green_offset=256;
            x11_blue_offset=1;
            x11_red_fsd = 255.0;
            x11_green_fsd = 255.0;
            x11_blue_fsd = 255.0;
         }
#else
         x11_red_offset=256*256;
         x11_green_offset=256;
         x11_blue_offset=1;
#endif
#ifdef DEBUG_TRUE24
         printf("x11_red_offset=%d\n", x11_red_offset, x11_red_offset);
         printf("x11_green_offset=%d\n", x11_green_offset, x11_green_offset);
         printf("x11_blue_offset=%d\n", x11_blue_offset, x11_blue_offset);
#ifdef TRUE_TECRA16
         printf("x11_red_fsd=%6.1f\n", x11_red_fsd);
         printf("x11_green_fsd=%6.1f\n", x11_green_fsd);
         printf("x11_blue_fsd=%6.1f\n", x11_blue_fsd);
#endif
#endif
      } else if (csame("x_my_cmap", type) || csame("x_my_hsv",type) ||
#else
      if (csame("x_my_cmap", type) || csame("x_my_hsv",type) ||
#endif
         csame("x_my_hv",type) || csame("x_my_rgbs",type)) {
         my_cmap = TRUE;
         if (csame("x_my_hsv",type)) {
            x11_avail_hsv = TRUE;
            x11.avail_rgb = TRUE;
            nhue = 12;
            nsat = 2;
         } else if  (csame("x_my_hv",type)) {
            x11_avail_hsv = TRUE;
            x11.avail_rgb = TRUE;
            nhue = 24;
            nsat = 1;
         } else if (csame("x_my_rgbs",type)) {
            x11_avail_rgbs = TRUE;
            x11.avail_rgb = TRUE;
         } else {
            x11.avail_colour = TRUE;
            x11.change_colour = TRUE;
         }
      } else if (csame("x_def_cmap", type) || csame("x_def_hsv",type) ||
         csame("x_def_hv",type) || csame("x_def_rgbs",type) ||
         csame("x",type)) {
         def_cmap = TRUE;
         if (csame("x_def_hsv",type)) {
            x11_avail_hsv = TRUE;
            x11.avail_rgb = TRUE;
            nhue = 12;
            nsat = 2;
         } else if  (csame("x_def_hv",type)) {
            x11_avail_hsv = TRUE;
            x11.avail_rgb = TRUE;
            nhue = 24;
            nsat = 1;
         } else if (csame("x_def_rgbs",type)) {
            x11_avail_rgbs = TRUE;
            x11.avail_rgb = TRUE;
         } else {
            x11.avail_colour = TRUE;
            x11.change_colour = TRUE;
         }
      } else if (csame("x_clone_cmap", type) || csame("x_clone_hsv",type) ||
         csame("x_clone_hv",type) || csame("x_clone_rgbs",type)) {
         clone_cmap = TRUE;
         if (csame("x_clone_hsv",type)) {
            x11_avail_hsv = TRUE;
            x11.avail_rgb = TRUE;
            nhue = 12;
            nsat = 2;
         } else if  (csame("x_clone_hv",type)) {
            x11_avail_hsv = TRUE;
            x11.avail_rgb = TRUE;
            nhue = 24;
            nsat = 1;
         } else if (csame("x_clone_rgbs",type)) {
            x11_avail_rgbs = TRUE;
            x11.avail_rgb = TRUE;
         } else {
            x11.avail_colour = TRUE;
            x11.change_colour = TRUE;
         }
        def_cmap = TRUE;
      }
#ifndef TRUE24
/*    Open X display */
      display = XOpenDisplay((char*)NULL);
      if (display == NULL) {
         fprintf(stderr, "error:  opening X server display\n");
         *ios = ERR;
         return;
      }
      screen = DefaultScreen(display);
      rootwin = RootWindow(display,screen);
#endif
      ixpmax = DisplayWidth(display, screen);
      iypmax = DisplayHeight(display, screen);
      ixmmax = DisplayWidthMM(display, screen);
      iymmax = DisplayHeightMM(display, screen);
/*    Acquire colormap and initialise it */
#ifdef TRUE24
      if (true_colour) {
         npixel = 0;
         x11_background_index = npixel;
         x11.reds[npixel] = 1.0;
         x11.greens[npixel] = 1.0;
         x11.blues[npixel] = 1.0;
         x11_background_pixel = x11_rgb_pixel(1.0,1.0,1.0);
         pixels[npixel] = x11_background_pixel;
         npixel++;
         x11.colour_index = npixel;
         x11.reds[npixel] = 0.0;
         x11.greens[npixel] = 0.0;
         x11.blues[npixel] = 0.0;
         pixels[npixel] = x11_rgb_pixel(0.0,0.0,0.0);
         npixel++;
      } else {
#endif
      if (my_cmap) {
         visual = DefaultVisual(display,screen);
         cmap = XCreateColormap(display, rootwin, visual, AllocAll);
         cmap_size = DisplayCells(display,screen);
         cmap_offset = 0;
#ifdef DEBUG
         printf("main:  my_cmap,cmap_size=%d\n", cmap_size);
#endif
         npixel = 0;
         for (k=cmap_offset; k < cmap_size; k++) {
            pixels[npixel++] = k;
         }
#ifdef DEBUG
         printf("main:  my_cmap:  offset=%d,npixel=%d\n", cmap_offset,
           npixel);
#endif
      } else if (clone_cmap) {
         defmap = DefaultColormap(display, screen);
         cmap_size = DisplayCells(display, screen);
#ifdef DEBUG
         printf("main:  clone_cmap,cmap_size=%d\n", cmap_size);
#endif
         npixel = 0;
         while (XAllocColorCells(display,defmap,FALSE,NULL,0,&pixels[npixel],1)) {
            npixel++;
         }
#ifdef DEBUG
         printf("main:  cloned cmap, found npixel=%d colour cells\n",npixel);
#endif
         cmap = XCopyColormapAndFree(display, defmap);
      } else {
         cmap = DefaultColormap(display, screen);
         cmap_size = DisplayCells(display, screen);
#ifdef DEBUG
         printf("main:  default cmap,cmap_size=%d\n", cmap_size);
#endif
         npixel = 0;
         while (XAllocColorCells(display,cmap,FALSE,NULL,0,&pixels[npixel],1)) {
            npixel++;
         }
#ifdef DEBUG
         printf("main:  default cmap, found npixel=%d colour cells\n",npixel);
#endif
      }
/*    Check that we have SOME colours */
      if (npixel < 64) {
         fprintf(stderr, "error:  X server has too few colours (%d) \n", npixel);
         *ios = ERR;
//         return;  MIR
         return 0;
      }
      x11.colour_offset = 0;
      x11.colour_size = npixel;
/*    Zap colourmap with RED */
      for (i=0; i < x11.colour_size; i++) {
         x11.reds[i] = 1.0;
         x11.greens[i] = 0.0;
         x11.blues[i] = 0.0;
      }
/*    Background is first */
      x11_background_index = 0;
      x11.reds[0] = 1.0;
      x11.greens[0] = 1.0;
      x11.blues[0] = 1.0;
/*    Foreground is first */
      x11.colour_index = 1;
      x11.reds[1] = 0.0;
      x11.greens[1] = 0.0;
      x11.blues[1] = 0.0;
/*    Dump those colours */
      x11_dump_colour(x11.colour_offset, x11.colour_size);
/*    Now increment offset to secure background */
      x11.colour_offset++;
      x11.colour_size--;
/*    Create RGBS colour cube */
      if (x11_avail_rgbs) {
         x11_setup_rgbs();
/*    Create HSV or HV colour cube */
      } else if (x11_avail_hsv) {
         x11_setup_hsv(nhue, nsat);
      }
#ifdef TRUE24
      }
#endif
   }
/* return pointer to "x11" structure */
   x11_count++;
   x11.count = x11_count;
   *p = &x11;
/* record name and type */
   ccopy(name,0,SZNAM,x11.name);
   ccopy(type,0,SZNAM,x11.type);
/* setup window size and position */
   x11_setup_window(ixpmax, iypmax, ixmmax, iymmax,
       &h0, &v0, &hlen, &vlen);
/* Realise requested width and height (inches) */
   ihlen = (int)hlen;
   ivlen = (int)vlen;
/* Realise window origin in screen frame */
   ihoff = (int)h0;
   ivoff = (int)v0;
/* open window with desired name, size and position */
   x11.ihmin = ihoff + 1;
   x11.ihmax = ihoff + ihlen;
   x11.ivmin = ivoff + 1;
   x11.ivmax = ivoff + ivlen;
/* Static window attributes */
   attribute_mask = CWBorderPixel | CWBackPixel;
   attributes.background_pixel = pixels[x11_background_index];
   attributes.border_pixel = pixels[x11_background_index];
/* Create window */
   window_h = x11.ihmin;
   window_v = iypmax - x11.ivmin - ivlen;
   width = ihlen;
   height = ivlen;
#ifdef DEBUG
   printf("x11_open:  width=%d,height=%d\n", width, height);
#endif
   window = XCreateWindow(display, rootwin,
     window_h, window_v, width, height, 5, CopyFromParent, InputOutput,
     CopyFromParent, attribute_mask, &attributes);
   if (window == 0) {
      fprintf(stderr, "error:  creating window\n");
      *ios = ERR;
//      return;  MIR
      return 0;
   }
/* Create pixmap for backing restoration */
   depth = DefaultDepth(display, screen);
   pix = XCreatePixmap(display, window, width, height, depth);
   if (pix == 0) {
      fprintf(stderr, "error:  creating backing pixmap\n");
      *ios = ERR;
//      return;  MIR
      return 0;
   }
/* Set colourmap for window */
#ifdef DEBUG
   printf("main:  set my cmap for window\n");
#endif
#ifdef TRUE24
   if (!true_colour) {
#endif
   XSetWindowColormap(display, window, cmap);
#ifdef TRUE24
   }
#endif
#ifdef OLD_FONT
/* Create a font */
   strcpy(font_name, "fixed");
   font = XLoadQueryFont(display, font_name);
   if (font == (XFontStruct*)NULL) {
      fprintf(stderr,"warning:  failed to load font=%s\n", font_name);
      *ios = ERR;
      return;
   }
#endif
/* Create a drawing context */
   gc_mask =  (unsigned long)0;
   gc = XCreateGC(display, window, gc_mask, &gcv);
#ifdef OLD_FONT
/* Set font for gc */
   XSetFont(display, gc, font->fid);
#endif
/* Select event types wanted */
   XSelectInput(display, window, ExposureMask | ButtonPressMask |
     StructureNotifyMask);
/* Set Standard props */
   XSetStandardProperties(display, window, name, name, None,
     0, 0, 0);
/* Tell WM about hints */
   argc=0;
   argv[argc++] = "muvi";
   argv[argc++] = name;
   x11_wm_hints(argc, argv, name, window_h, window_v);
#ifdef WM_DELETE
/* Deal with WM Protocol relating to destruction */
   protocols_atom = XInternAtom(display, "WM_PROTOCOLS", FALSE);
   delete_atom = XInternAtom(display, "WM_DELETE_WINDOW", FALSE);
   if ((delete_atom != (Atom)None) &&
       (protocols_atom != (Atom)None)) {
      XSetWMProtocols(display, window, &delete_atom, 1);
   }
#endif
/* Map window (make it appear) */
   XMapWindow(display, window);
/* Tell WM about colourmap */
#ifdef DEBUG
   printf("main:  set cmap for WM\n");
#endif
#ifdef TRUE24
   if (!true_colour) {
#endif
   XSetWMColormapWindows(display, window, &window, 1);
#ifdef TRUE24
   }
#endif
/* Flush display */
   XFlush(display);
/* Await map notification */
   x11_await_map();
#ifdef DEBUG
   printf("x11_open:  received map notify\n");
#endif
/* get actual coordinate limits*/
   x11.ihmin = 0;
   x11.ihmax = ihlen;
   x11.ivmin = 0;
   x11.ivmax = ivlen;
/* real coordinate limits */
   x11.hmin = (float)x11.ihmin;
   x11.hmax = (float)x11.ihmax;
   x11.vmin = (float)x11.ivmin;
   x11.vmax = (float)x11.ivmax;
   x11_prev_font = FALSE;
#ifdef OLD_FONT
/* Open a suitable font */
   x11.avail_fontscale = FALSE;
   ccopy("Screen", 0, SZNAM, x11.font_style);
   x11.font_size = 12.0;
   x11_set_font();
#endif
/* rectangles, polygons, bitmap */
   x11.avail_rect = TRUE;
   x11.avail_rub = TRUE;
   x11.avail_polygon = TRUE;
   x11.avail_bitmap = TRUE;
   x11.bitmap_dh = 1.0;
   x11.bitmap_dv = 1.0;
/* cursor */
   x11.avail_cursor = FALSE;
   x11.avail_sample = FALSE;
   x11.cursor_buttons = 3;
/* line width */
   x11.line_width = 0.25;
/* copy window-specifics into table */
   x11_windows[x11_now] = window;
   x11_gcs[x11_now] = gc;
   x11_pixs[x11_now] = pix;
   x11_ihmin[x11_now] = x11.ihmin;
   x11_ihmax[x11_now] = x11.ihmax;
   x11_ivmin[x11_now] = x11.ivmin;
   x11_ivmax[x11_now] = x11.ivmax;
   x11_hmin[x11_now] = x11.hmin;
   x11_hmax[x11_now] = x11.hmax;
   x11_vmin[x11_now] = x11.vmin;
   x11_vmax[x11_now] = x11.vmax;
/* clear things */
   x11_erase();
#ifndef OLD_FONT
/* Load fonts if on primary window */
   if (x11_count == 1) {
      ccopy("default", 0, SZNAM, x11.font_style);
      x11.font_size = 13.0;
      x11_set_font();
   } else {
#ifdef DEBUG
      printf("x11_open:  set font for new window\n");
#endif
      XSetFont(display, gc, font->fid);
   }
#endif
/* it worked! */
   *ios = OK;
   return 0;
}

int x11_close(name, ios)
   char name[];
   int *ios;
{
   int i;

#ifdef DEBUG
   printf("x11_close:  closing %s\n", name);
#endif
   nam_find(name, MAXWIND, x11_used, x11_names, &i);
   if (i < 0) {
      fprintf(stderr,"error:  cannot close unknown window\n");
      *ios = ERR;
   } else {
#ifdef BUG
      window = x11_windows[i];
#endif
#ifdef REPAIR
      XFreePixmap(display, x11_pixs[i]);
      x11_pixs[i] = 0;
#endif
      XDestroyWindow(display, x11_windows[i]);
      XFlush(display);
      x11_count--;
      x11.count = x11_count;
      *ios = OK;
      if (i == x11_now) {
         x11_now = -1;
      }
      x11_used[i] = FALSE;
      if (x11_count == 0) {
         if (x11_prev_font) {
            XFreeFont(display, font);
            x11_prev_font = FALSE;
         }
         XCloseDisplay(display);
      }
   }
   return 0;
}

int x11_select(name, ios)
   char name[];
   int *ios;
{
   int i;

   nam_find(name, MAXWIND, x11_used, x11_names, &i);
   if (i < 0) {
      fprintf(stderr,"error:  cannot select unknown x11 window\n");
      *ios = ERR;
   } else {
      window = x11_windows[i];
      gc = x11_gcs[i];
      pix = x11_pixs[i];
      x11.ihmin = x11_ihmin[i];
      x11.ihmax = x11_ihmax[i];
      x11.ivmin = x11_ivmin[i];
      x11.ivmax = x11_ivmax[i];
      x11.hmin = x11_hmin[i];
      x11.hmax = x11_hmax[i];
      x11.vmin = x11_vmin[i];
      x11.vmax = x11_vmax[i];
      x11_flush();
      x11_now = i;
      *ios = OK;
   }
   return 0;
}

int x11_erase() {

   x11_pend_events();
#ifdef DEBUG
   printf("x11_erase:  background_index=%d,pixel=%d\n",
     x11_background_index, pixels[x11_background_index]);
#endif
   XSetForeground(display, gc, pixels[x11_background_index]);
   XFillRectangle(display, window, gc, 0, 0, width, height);
#ifdef REPAIR
   XFillRectangle(display, pix, gc, 0, 0, width, height);
#endif
   XSetForeground(display, gc, pixels[x11.colour_index]);
   x11_flush();
   return 0;
}

int x11_rub(h1, h2, v1, v2) 
   float h1,h2,v1,v2;
{
   int dh,dv,h0,v0;

   dh = (int)(h2-h1)+1;
   dv = (int)(v2-v1)+1;
   h0 = (int)h1;
   v0 = height - (int)v1 - dv;
#ifdef DEBUG
   printf("x11_rub:  h0=%d,v0=%d,dh=%d,dv=%d\n", h0, v0, dh, dv);
#endif
   XSetForeground(display, gc, pixels[x11_background_index]);
   XFillRectangle(display, window, gc, h0, v0, dh, dv);
#ifdef REPAIR
   XFillRectangle(display, pix, gc, h0, v0, dh, dv);
#endif
   XSetForeground(display, gc, pixels[x11.colour_index]);
   x11_flush();
   return 0;
}


int x11_rect(h1, h2, v1, v2)
   float h1, h2, v1, v2;
{
   int dh,dv,h0,v0;

   dh = (int)(h2-h1)+1;
   dv = (int)(v2-v1)+1;
   h0 = (int)h1;
   v0 = height - (int)v1 - dv;
#ifdef DEBUG
   printf("x11_rect:  h1=%f,h2=%f,v1=%f,v2=%f\n", h1, h2, v1, v2);
   printf("x11_rect:  v0=%d\n", v0);
#endif
   XFillRectangle(display, window, gc, h0, v0, dh, dv);
#ifdef REPAIR
   XFillRectangle(display, pix, gc, h0, v0, dh, dv);
#endif
   return 0;
}

int x11_polyline(na, ha, va)
   float ha[], va[];
   int na;
{
   int h1,h2,v1,v2;
   int i;

   for (i=1; i < na; i++){
      h1 = (int)ha[i-1];
      v1 = height - (int)va[i-1];
      h2 = (int)ha[i];
      v2 = height - (int)va[i];
      XDrawLine(display, window, gc, h1, v1, h2, v2);
#ifdef REPAIR
      XDrawLine(display, pix, gc, h1, v1, h2, v2);
#endif
   };
   return 0;
}

int x11_polygon(na, ha, va)
   float ha[], va[];
   int na;
{
   XPoint points[256];
   int i;

   for (i=0; i < na; i++) {
      points[i].x = (int)ha[i];
      points[i].y = height - (int)va[i];
   }
   XFillPolygon(display, window, gc, points, na, Convex, 
     CoordModeOrigin);
#ifdef REPAIR
   XFillPolygon(display, pix, gc, points, na, Convex, 
     CoordModeOrigin);
#endif
   return 0;
}

int x11_text(h, v, str, just)
   float h, v;
   char str[],just[2];
{

   float hlen,vlen;
   int h0,v0;

   x11_text_box(str,&hlen,&vlen);

   if (just[0] == 'r') {
      h0 = (int)(h - hlen);
   } else if (just[0] == 'm') {
      h0 = (int)(h - (hlen/2.0));
   } else {
      h0 = (int)h;
   }
   if (just[1] == 't') {
      v0 = (int)(v - x11.char_height);
   } else if (just[1] == 'm') {
      v0 = (int)(v - x11.char_height/2.0);
   } else {
      v0 = (int)v;;
   }
   v0 = height - v0;
#ifdef DEBUG
   printf("x11_text:  h0=%d,v0=%d,str=%s\n", h0, v0, str);
#endif
   XDrawString(display, window, gc, h0, v0, str, strlen(str));
#ifdef REPAIR
   XDrawString(display, pix, gc, h0, v0, str, strlen(str));
#endif
   return 0;
}

int x11_text_box(str, dh, dv)
   char str[];
   float *dh,*dv;
{
   int nchar;

   nchar = slen(str);
   *dh = XTextWidth(font, str, nchar);
   *dv = x11.char_height+x11.char_depth;
   return 0;
}

int x11_cursor(h, v, button)
   float *h, *v;
   int *button;
{
   return 0;
}

int x11_sample(h, v, button)
   float *h,*v;
   int *button;
{
   return 0;
}


int x11_invent_font(char style[], float size, char name[]) {
   int pt,ptc;
   FILE *fp;
   bool found;
   char line[SZLIN],stylec[SZNAM],file[SZLIN],*ep;

   pt = (int)(size+0.5);
   if ((ep=getenv("MUVI_FONT_TABLE"))==(char*)NULL) {
      ccopy("muvi.fonts", 0, SZLIN, file);
      fprintf(stderr, "warning:  no MUVI_FONT_TABLE undefined\n");
   } else {
      ccopy(ep, 0, SZLIN, file);
   }
   if (csame("default",style)) {
      ccopy(DEFAULT_FONT, 0, SZLIN, name);
   } else if ((fp=fopen(file,"r")) == NULL) {
      fprintf(stderr, "warning:  unable to open MUVI font table %s\n", file);
      if (csame(style,"variable")) {
         ccopy("variable", 0, SZLIN, name);
      } else {
         ccopy(DEFAULT_FONT, 0, SZLIN, name);
      }
   } else {
#ifdef DEBUG
      printf("x11_invent_font:  opened file=%s\n", file);
#endif
      found = FALSE;
      while(fgets(line, SZLIN,fp)!=(char*)NULL) {
#ifdef DEBUG
            printf("..line=%s\n", line);
#endif
         if (sscanf(line, "%s %d %s\n", stylec, &ptc, name)!=3) {
            fprintf(stderr, "error:  bad MUVI font table line\n");
            break;
         } else {
#ifdef DEBUG
            printf("..stylec=%s,ptc=%d\n", stylec, ptc);
#endif
            if (csame(style,stylec) && (pt==ptc)) {
#ifdef DEBUG
               printf("..matches\n");
#endif
               found = TRUE;
               break;
            }
         }
      }
      fclose(fp);
      if (found) {
      } else {
         fprintf(stderr, "error:  did not find font match\n");
         ccopy(DEFAULT_FONT, 0, SZLIN, name);
      }
   }
   return 0;
}

int x11_set_font() {
   char c;
   int i,w,wmax;

   if (x11_prev_font) {
      XFreeFont(display, font);
      x11_prev_font = FALSE;
#ifdef DEBUG
      printf("x11_set_font:  free previous one\n");
#endif
   }
   x11_invent_font(x11.font_style, x11.font_size, font_name);
#ifdef DEBUG
   printf("x11_set_font:  font_style=%s,font_size=%6.2f,font_name=%s\n",
     x11.font_style, x11.font_size, font_name);
#endif
   font = XLoadQueryFont(display, font_name);
   if (font == (XFontStruct*)NULL) {
      fprintf(stderr, "error:  cannot load font %s %6.2f\n",
        x11.font_style, x11.font_size);
      fprintf(stderr, "note:  that was font_name=%s\n", font_name);
      fprintf(stderr, "note:  will revert to default\n");
      ccopy("Default", 0, SZNAM, x11.font_style);
      x11.font_size = 13.0;
      ccopy(DEFAULT_FONT, 0, SZLIN, font_name);
      font = XLoadQueryFont(display, font_name);
      if (font == (XFontStruct*)NULL) {
         fprintf(stderr, "error:  reverting to default font\n");
         fprintf(stderr, "warning:  no characters will appear!\n");
      } else {
         for (i=0; i < MAXWIND; i++) {
            if (x11_used[i]) {
#ifdef DEBUG
               printf("x11_set_font:  setting for window %d\n", i);
#endif
               XSetFont(display, x11_gcs[i], font->fid);
            }
         }
         x11_prev_font = TRUE;
      }
   } else {
      for (i=0; i < MAXWIND; i++) {
         if (x11_used[i]) {
#ifdef DEBUG
            printf("x11_set_font:  setting for window %d\n", i);
#endif
            XSetFont(display, x11_gcs[i], font->fid);
         }
      }
      x11_prev_font = TRUE;
   }
   x11.avail_fontscale = TRUE;
   if (x11_prev_font) {
      x11.char_height = (float)(font->ascent);
      x11.char_depth = (float)(font->descent);
      wmax = 0;
      for (i=(int)' '; i <= (int)'~'; i++) {
         c = (char)i;
         w = XTextWidth(font, &c, 1);
         if (w > wmax) {
            wmax = w;
         }
      }
      x11.char_width = (float)wmax;
      x11.char_space = 0.0;
#ifdef DEBUG
      printf("x11_set_font:  height=%f,width=%f,depth=%f\n",
        x11.char_height, x11.char_width, x11.char_depth);
#endif
   }
   return 0;
}

int x11_set_colours(offset, size, reds, greens, blues)
   int offset,size;
   float reds[],greens[],blues[];
{
   int i,k;

#ifdef TRUE24
   if (true_colour) {
   } else if (x11.avail_colour) {
#else
   if (x11.avail_colour) {
#endif
      k = offset;
      for (i=0; i < size; i++) {
         x11.reds[k] = reds[i];
         x11.greens[k] = greens[i];
         x11.blues[k] = blues[i];
         k++;
      }
      x11_dump_colour(offset, size);
   }
   return 0;
}

int x11_set_colour(index)
   int index;
{
#ifdef TRUE24
   if (true_colour) {
   } else if (x11.avail_colour) {
#else
   if (x11.avail_colour) {
#endif
      x11.colour_index = index;
#ifdef DEBUG
      printf("x11_set_colour:  index=%d\n", index);
#endif
      XSetForeground(display, gc, pixels[x11.colour_index]);
   }
   return 0;
}

int x11_set_rgb(red, green, blue)
   float red,green,blue;
{
   float rgbvec[3],hsv[3];
#ifdef TRUE24
   unsigned long pixel;
#endif

   x11.red = red;
   x11.green = green;
   x11.blue = blue;
#ifdef TRUE24
   if (true_colour) {
      pixel = x11_rgb_pixel(red, green, blue);
      XSetForeground(display, gc, pixel);
   } else if (x11_avail_rgbs) {
#else
   if (x11_avail_rgbs) {
#endif
      rgbvec[0] = red;
      rgbvec[1] = green;
      rgbvec[2] = blue;
      x11.colour_index = x11_rgb_index(rgbvec);
      XSetForeground(display, gc, pixels[x11.colour_index]);
   } else if (x11_avail_hsv) {
      rgbvec[0] = red;
      rgbvec[1] = green;
      rgbvec[2] = blue;
      x11_rgb_to_hsv(rgbvec, hsv);
      x11.colour_index = x11_hsv_index(hsv);
#ifdef DEBUG
      printf("x11_set_rgb:  rgb=[%5.2f,%5.2f,%5.2f]\n", rgbvec[0],
        rgbvec[1],rgbvec[2]);
      printf("x11_set_rgb:  hsv=[%5.2f,%5.2f,%5.2f]\n", hsv[0],
        hsv[1],hsv[2]);
      printf("..i=%d\n", x11.colour_index);
#endif
      XSetForeground(display, gc, pixels[x11.colour_index]);
   } else {
      x11_set_paint(x11.colour_index, red, green, blue);
   }
   return 0;
}

int x11_set_rgb_background(red, green, blue)
   float red,green,blue;
{
   x11.red_background = red;
   x11.green_background = green;
   x11.blue_background = blue;
#ifdef TRUE24
   if (true_colour) {
      x11_background_pixel = x11_rgb_pixel(red,green,blue);
      pixels[x11_background_index] = x11_background_pixel;
   } else if (x11.avail_colour) {
#else
   if (x11.avail_colour) {
#endif
      x11_set_paint(x11_background_index, red, green, blue);
   }
   return 0;
}

int x11_set_line_width(width)
   float width;
{
   x11.line_width = width;
   return 0;
}

int x11_colour_bitmap(ih0,iv0, nh,nv, cim)
   int ih0,iv0,nh,nv;
   int cim[];
{
   int k,ih,iv,i,kv;
   long nbyte;
   XImage *image;

   if (x11.avail_colour) {
#ifdef DEBUG
      printf("x11_colour_bitmap:  creating [%d,%d] image\n", nh,
        nv);
#endif
      if ((image=XCreateImage(display, DefaultVisual(display,screen),
       DefaultDepth(display,screen),
       ZPixmap, 0, (char*)NULL, nh, nv,
       XBitmapPad(display), 0)) == NULL) {
         fprintf(stderr, "error:  failed to get X Image\n");
      } else {
         nbyte = image->bytes_per_line * image->height;
         image->data = (char*)malloc(nbyte);
         if (image->data == NULL) {
            fprintf(stderr,"Error:  getting memory for x11 colour_bitmap\n");
         } else {
#ifdef DEBUG
            printf("x11_colour_bitmap:  allocated %d bytes\n", nbyte);
#endif
            k = 0;
            for (iv=0; iv < nv; iv++) {
              for (ih=0; ih < nh; ih++) {
                 kv = nv-iv-1;
                 i = cim[k];
                 XPutPixel(image, ih, kv, pixels[i]);
                 k++;
              }
           }
           iv = height - iv0 - nv;
           XPutImage(display, window, gc, image, 0, 0, ih0, iv, nh, nv);
#ifdef REPAIR
           XPutImage(display, pix, gc, image, 0, 0, ih0, iv, nh, nv);
#endif
         }
         XDestroyImage(image);
      }
   }
   return 0;
}

int x11_rgb_bitmap(ih0,iv0, nh,nv, rgb_image)
   int ih0,iv0,nh,nv;
   float rgb_image[];
{
   int k,ih,iv,i,kv;
   long nbyte;
   XImage *image;
   float rgb[3],hsv[3],dc;
#ifdef TRUE24
   unsigned int pixel;
#endif

#ifdef TRUE24
   if (true_colour || x11_avail_hsv || x11_avail_rgbs) {
#else
   if (x11_avail_hsv || x11_avail_rgbs) {
#endif
#ifdef DEBUG
      printf("x11_rgb_bitmap:  creating [%d,%d] image\n", nh,
        nv);
#endif
      if ((image=XCreateImage(display, DefaultVisual(display,screen),
       DefaultDepth(display,screen),
       ZPixmap, 0, (char*)NULL, nh, nv,
       XBitmapPad(display), 0)) == NULL) {
         fprintf(stderr, "error:  failed to get X Image\n");
      } else {
         nbyte = image->bytes_per_line * image->height;
         image->data = (char*)malloc(nbyte);
         if (image->data == NULL) {
            fprintf(stderr,"Error:  getting memory for x11 rgb_bitmap\n");
         } else {
#ifdef DEBUG
            printf("x11_rgb_bitmap:  allocated %d bytes\n", nbyte);
#endif
#ifdef TRUE24
            if (true_colour) {
               k = 0;
               for (iv=0; iv < nv; iv++) {
                  for (ih=0; ih < nh; ih++) {
                     rgb[0] = rgb_image[k++];
                     rgb[1] = rgb_image[k++];
                     rgb[2] = rgb_image[k++];
                     pixel = x11_rgb_pixel(rgb[0], rgb[1], rgb[2]);
                     kv = nv-iv-1;
                     XPutPixel(image, ih, kv, pixel);
                  }
               }
            } else if (x11_avail_rgbs) {
#else
            if (x11_avail_rgbs) {
#endif
               k = 0;
               for (iv=0; iv < nv; iv++) {
                  for (ih=0; ih < nh; ih++) {
                     rgb[0] = rgb_image[k++];
                     rgb[1] = rgb_image[k++];
                     rgb[2] = rgb_image[k++];
                     dc = randaf[rand_count++];
                     if (rand_count == MAXRAND) {
                        rand_count = 0;
                     }
                     rgb[0] = rgb[0] + dc*r_space;
                     dc = randaf[rand_count++];
                     if (rand_count == MAXRAND) {
                        rand_count = 0;
                     }
                     rgb[1] = rgb[1]+ dc*g_space;
                     dc = randaf[rand_count++];
                     if (rand_count == MAXRAND) {
                        rand_count = 0;
                     }
                     rgb[2] = rgb[2]+ dc*b_space;
                     kv = nv-iv-1;
                     i = x11_rgb_index(rgb);
                     XPutPixel(image, ih, kv, pixels[i]);
                  }
               }
            } else if (x11_avail_hsv) {
               k = 0;
               for (iv=0; iv < nv; iv++) {
                  for (ih=0; ih < nh; ih++) {
                     rgb[0] = rgb_image[k++];
                     rgb[1] = rgb_image[k++];
                     rgb[2] = rgb_image[k++];
                     x11_rgb_to_hsv(rgb, hsv);
                     dc = randaf[rand_count++];
                     if (rand_count == MAXRAND) {
                        rand_count = 0;
                     }
                     hsv[2] = hsv[2] + dc*v_space;
                     i = x11_hsv_index(hsv);
                     kv = nv-iv-1;
                     XPutPixel(image, ih, kv, pixels[i]);
#ifdef DDEBUG
                     printf("..hsv=[%5.2f,%5.2f,%5.2f],i=%d\n",
                       hsv[0],hsv[1],hsv[2],i);
#endif
                  }
               }
            }
            iv = height - iv0 - nv;
            XPutImage(display, window, gc, image, 0, 0, ih0, iv, nh, nv);
#ifdef REPAIR
            XPutImage(display, pix, gc, image, 0, 0, ih0, iv, nh, nv);
#endif
         }
         XDestroyImage(image);
      }
   }
   return 0;
}

int x11_rgbz_put_bitmap(ih0,iv0, nh,nv, rgb_image, d_image)
   int ih0,iv0,nh,nv;
   float rgb_image[],d_image[];
{
   return 0;
}

int x11_rgbz_get_bitmap(ih0,iv0, nh,nv, rgb_image, d_image)
   int ih0,iv0,nh,nv;
   float rgb_image[],d_image[];
{
   return 0;
}

/* RGBZ PRIMITIVES */

/* x11_rgbz_setup - setup the 3d projection etc. */
int x11_rgbz_setup(anglev, near, far, hmin, hmax, vmin, vmax)
   float anglev,near,far,hmin,hmax,vmin,vmax;
{
   return 0;
}

/* x11_rgbz_begin - begin 3d work */
int x11_rgbz_begin() {
   return 0;
}

/* x11_rgbz_end - end 3d work */
int x11_rgbz_end() {
   return 0;
}

/* x11_rgbz_clear - clear 3d picture */
int x11_rgbz_clear() {
   return 0;
}

int x11_rgbz_begin_polygon() {
   return 0;
}

int x11_rgbz_polygon_vertex(x, y, z)
   float x,y,z;
{
   return 0;
}

int x11_rgbz_end_polygon() {
   return 0;
}

int x11_rgbz_begin_polyline() {
   return 0;
}

int x11_rgbz_polyline_vertex(x, y, z)
   float x,y,z;
{
   return 0;
}

int x11_rgbz_end_polyline() {
   return 0;
}

int x11_flush() {
#ifdef REPAIR
#ifdef DEBUG
   printf("x11_flush:  invokes x11_pend_events\n");
#endif
   x11_pend_events();
#endif
   XFlush(display);
   return 0;
}

int x11_prepare()
{
   int i;

   x11_setup_rand();
/* window manager */
   for (i=0; i < MAXWIND; i++) {
      x11_used[i] = FALSE;
   }
/* window size-position */
   x11_pos_rel = TRUE;
   x11_pos_rh = 1.0;
   x11_pos_rv = 1.0;
   x11_size_rel = TRUE;
   x11_size_rh = 0.5;
   x11_size_rv = 0.5;
/* routine bindings (see "drp.h") */
   x11.prepare = x11_prepare;
   x11.open = x11_open;
   x11.close = x11_close;
   x11.select = x11_select;
   x11.erase = x11_erase;
   x11.rub = x11_rub;
   x11.rect = x11_rect;
   x11.polyline = x11_polyline;
   x11.polygon = x11_polygon;
   x11.text = x11_text;
   x11.text_box = x11_text_box;
   x11.cursor = x11_cursor;
   x11.sample = x11_sample;
   x11.set_font = x11_set_font;
   x11.set_colours = x11_set_colours;
   x11.set_colour = x11_set_colour;
   x11.set_rgb = x11_set_rgb;
   x11.set_rgb_background = x11_set_rgb_background;
   x11.set_line_width = x11_set_line_width;
   x11.colour_bitmap = x11_colour_bitmap;
   x11.rgb_bitmap = x11_rgb_bitmap;
/* rgbz bindings */
   x11.rgbz_setup = x11_rgbz_setup;
   x11.rgbz_begin = x11_rgbz_begin;
   x11.rgbz_end = x11_rgbz_end;
   x11.rgbz_clear = x11_rgbz_clear;
   x11.rgbz_begin_polygon = x11_rgbz_begin_polygon;
   x11.rgbz_polygon_vertex = x11_rgbz_polygon_vertex;
   x11.rgbz_end_polygon = x11_rgbz_end_polygon;
   x11.rgbz_begin_polyline = x11_rgbz_begin_polyline;
   x11.rgbz_polyline_vertex = x11_rgbz_polyline_vertex;
   x11.rgbz_end_polyline = x11_rgbz_end_polyline;
   x11.rgbz_put_bitmap = x11_rgbz_put_bitmap;
   x11.rgbz_get_bitmap = x11_rgbz_get_bitmap;
   x11.gflush = x11_flush;
   return 0;
}


/* PRIVATE FUNCTIONS */

int x11_setup_window(ixpmax,iypmax,ixmmax,iymmax, h0,v0, hsize,vsize)
   int ixpmax,iypmax,ixmmax,iymmax;
   float *h0,*v0,*hsize,*vsize;
{
   float xdpi,ydpi;
   float margin,hmin,vmin,hlen,vlen,hmax,vmax;

   xdpi = 25.4*ixpmax/ixmmax;
   ydpi = 25.4*iypmax/iymmax;

   margin = 0.025;
   hmin = margin*ixpmax;
   vmin = margin*ixpmax;
   hmax = (1.0-margin)*ixpmax;
   vmax = (1.0-margin*2)*iypmax;  /* increase margin to allow for WM frames */
   hlen = hmax-hmin;
   vlen = vmax-vmin;

   if (x11_size_rel) {
      *hsize = hlen*x11_size_rh;
      *vsize = vlen*x11_size_rv;
   } else {
      *hsize = x11_size_ah*xdpi;
      *vsize = x11_size_av*ydpi;
   }
   if (*hsize > 0.9*hlen) {
      *hsize = 0.9*hlen;
   } else if (*hsize < 0.1*hlen) {
      *hsize = 0.1*hlen;
   }
   if (*vsize > vlen*0.9) {
      *vsize = vlen*0.9;
   } else if (*vsize < 0.1*vlen) {
      *vsize = 0.1*vlen;
   }

   if (x11_pos_rel) {
      *h0 = hmin + x11_pos_rh * (hlen - *hsize);
      *v0 = vmin + x11_pos_rv * (vlen - *vsize);
   } else {
      *h0 = hmin + x11_pos_ah * xdpi;
      *v0 = hmin + x11_pos_av * ydpi;
   }
   if (*h0 < hmin) {
      *h0 = hmin;
   } else if (*h0 > hmax) {
      *h0 = 0.9*hlen;
   }
   if (*v0 < vmin) {
      *v0 = vmin;
   } else if (*v0 > vlen) {
      *v0 = 0.9*vlen;
   }

   if (*hsize+*h0 > hmax) {
      *hsize = hlen - *h0; }
   if (*vsize+*v0 > vmax) {
      *vsize = vlen - *v0;
   }
   return 0;
}


int x11_dump_colour(offset, size)
   int offset, size;
{
   float dc;
   int i,k;
   XColor c;

#ifdef DEBUG
   printf("x11_dump_colour:  offset=%d,size=%d\n", offset, size);
#endif
   dc = 65535.0;
   k = offset;
   for (i=0; i < size; i++) {
      c.flags = DoRed | DoGreen | DoBlue;
      c.pixel = pixels[k];
      c.red = (unsigned short)(x11.reds[k]*dc);
      c.green = (unsigned short)(x11.greens[k]*dc);
      c.blue = (unsigned short)(x11.blues[k]*dc);
      XStoreColor(display, cmap, &c);
#ifdef DEBUG
      printf("..i=%d,k=%d,pixel=%d,rgb=%5.2f,%5.2f,%5.2f\n",
         i, k, pixels[k], x11.reds[k], x11.greens[k], x11.blues[k]);
#endif
      k++;
   }
   return 0;
}

int x11_set_paint(offset, red, green, blue)
   int offset;
   float red, green, blue;
{
   x11.reds[offset] = red;
   x11.greens[offset] = green;
   x11.blues[offset] = blue;
   x11_dump_colour(offset, 1);
   return 0;
}

int x11_setup_rand() {
   int i;

   for (i=0; i < MAXRAND; i++) {
      randaf[i] = ranf*((float)rand()*divran-0.5);
   }
   rand_count = 0;
   norand = FALSE;
   return 0;
}

/* x11_mrs_videoc - MRS special to capture screen */
int x11_mrs_videoc(file, ios)
   char file[];
   int *ios;
{
   *ios = OK;
   return 0;
}

int x11_wm_hints(int argc, char *argv[], char name[], int h, int v) {
   XTextProperty w_name;
   XTextProperty i_name;
   XSizeHints sizehints;
   XWMHints wmhints;
   XClassHint classhints;
   //int status;
   (void)XStringListToTextProperty(&name, 1, &w_name);	// Cast to void BM
   (void)XStringListToTextProperty(&name, 1, &i_name);	// Cast to void BM
   sizehints.x = h;
   sizehints.y = v;
   sizehints.width = width;
   sizehints.height = height;
   sizehints.base_width = width;
   sizehints.base_height = height;
#ifdef OLD
   sizehints.flags = PPosition | PSize | PBaseSize;
#else
   sizehints.max_width = width;
   sizehints.max_height = height;
   sizehints.flags = PPosition | PSize | PBaseSize | PMaxSize;
#endif
   wmhints.initial_state = NormalState;
   wmhints.input = True;
   wmhints.flags = StateHint | InputHint;
   classhints.res_name = "MUVI";
   classhints.res_class = "MUVI";
   XSetWMProperties(display, window, &w_name, &i_name, argv, argc,
     &sizehints, &wmhints, &classhints);
   XFree(w_name.value);
   XFree(i_name.value);
   return 0;
}

/* Colour Model Support */

/* hsv_to_rgb:  convert HSV to RGB */
int x11_hsv_to_rgb(float hsv[3], float rgb[3]) {
   float h,s,v,r,g,b,f,p,q,t;
   int i;

   h = hsv[0]*6.0;
   s = hsv[1];
   v = hsv[2];
   if (s <= 0.0) {
      r = v;
      g = v;
      b = v;
   } else {
      i = (int)h;
      f = h - (float)i;
      p = v*(1.0 - s);
      q = v*(1.0 - (s*f));
      t = v*(1.0 - (s*(1.0-f)));
      switch (i) {
         case 0:
            r = v;
            g = t;
            b = p;
            break;
         case 1:
            r = q;
            g = v;
            b = p;
            break;
         case 2:
            r = p;
            g = v;
            b = t;
            break;
         case 3:
            r = p;
            g = q;
            b = v;
            break;
         case 4:
            r = t;
            g = p;
            b = v;
            break;
         default:
            r = v;
            g = p;
            b = q;
      }
   }
   rgb[0] = r;
   rgb[1] = g;
   rgb[2] = b;
   return 0;
}

/* x11_rgb_hsv:  convert RGB to HSV */
int x11_rgb_to_hsv(float rgb[3], float hsv[3]) {
   float min,max,h,s,v,r,g,b,rc,gc,bc;
   int i;

   min = 1.0;
   max = 0.0;
   for (i=0; i < 3; i++) {
      min = rgb[i] < min ? rgb[i] : min;
      max = rgb[i] > max ? rgb[i] : max;
   }
   v = max;
   if (max <= 0.0) {
      s = 0.0;
      h = 0.0;
   } else if (max <= min) {
      s = 0.0;
      h = 0.0;
   } else {
      r = rgb[0];
      g = rgb[1];
      b = rgb[2];
      s = (max - min)/max;
      rc = (max-r)/(max-min);
      gc = (max-g)/(max-min);
      bc = (max-b)/(max-min);
      if (r == max) {
         h = bc - gc;
      } else if (g == max) {
         h = 2.0 + rc - bc;
      } else {
         h = 4.0 + gc - rc;
      }
      if (h < 0.0) {
         h = h + 6.0;
      }
   }
   hsv[0] = h/6.0;
   hsv[1] = s;
   hsv[2] = v;
   return 0;
}

/* x11_setup_hsv:  setup HSV system cone */
int x11_setup_hsv(int nhue, int nsat) {
   int ih,is,iv,ig,k;
   XColor c;
   float hsv[3],rgb[3],g;

/* h[0.0,1.0] */
   h_count = nhue;
   h_space = 1.0/(float)h_count;
/* s[0.0,1.0] */
   s_count = nsat;
   s_space = 1.0/(float)s_count;
/* v[0.0,1.0] */
   v_count = (x11.colour_size-1)/(h_count*s_count+1);
   v_space = 1.0/(float)v_count;
/* index into main matrix */
   v_mult = 1;
   s_mult = v_count;
   h_mult = s_count*v_count;
/* Grey space */
   grey_space = 1.0/(float)v_count;
#ifdef DEBUG
   printf("x11_setup_hsv:  colour_size=%d\n", x11.colour_size);
   printf("x11_setup_hsv:  colour_offset=%d\n", x11.colour_offset);
   printf("x11_setup_hsv:  HSV[%d,%d,%d]\n", h_count, s_count, v_count);
#endif
/* setup hsv sub-cube */
   hsv_offset = x11.colour_offset;
#ifdef DEBUG
   printf("x11_setup_hsv:  main HSV sub-cube\n");
#endif
   for (ih=0; ih < h_count; ih++) {
      for (is=0; is < s_count; is++) {
         for (iv=0; iv < v_count; iv++) {
            hsv[0] = (float)ih*h_space;
            hsv[1] = (float)(is+1)*s_space;
            hsv[2] = (float)(iv+1)*v_space;
            x11_hsv_to_rgb(hsv, rgb);
            k = ih*h_mult + is*s_mult + iv*v_mult + hsv_offset;
            c.pixel = pixels[k];
            c.red = (unsigned short)(rgb[0]*65535.0);
            c.green = (unsigned short)(rgb[1]*65535.0);
            c.blue = (unsigned short)(rgb[2]*65535.0);
            c.flags = DoRed | DoGreen | DoBlue;
            XStoreColor(display, cmap, &c);
#ifdef DEBUG
            printf("..k=%d,i=[%d,%d,%d],hsv=[%5.2f,%5.2f,%5.2f]\n", k,
              ih,is,iv, hsv[0], hsv[1], hsv[2]);
#endif
         }
      }
   }
/* grey sub-vector */
   grey_offset = hsv_offset + (h_count*s_count*v_count);
#ifdef DEBUG
   printf("..grey_offset=%d\n", grey_offset);
#endif
   for (ig=0; ig <= v_count; ig++) {
      g = (float)ig*grey_space;
      k = ig + grey_offset;
      c.pixel = pixels[k];
      c.red = (unsigned short)(g*65535.0);
      c.green = (unsigned short)(g*65535.0);
      c.blue = (unsigned short)(g*65535.0);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(display, cmap, &c);
   }
   return 0;
}

/* x11_hsv_index:  get pixel index from HSV */
int x11_hsv_index(float hsv[3]) {
   int i;
   int ih,is,iv;
   float h,s,v;

   h = hsv[0];
   s = hsv[1];
   v = hsv[2];
   iv = (int)(v/v_space + 0.5);
   if (iv <= 0) {
      i = grey_offset;
   } else {
      if (iv > v_count) {
         iv = v_count;
      }
      is = (int)(s/s_space + 0.5);
      if (is <= 0) {
         i = iv + grey_offset;
      } else {
         if (is > s_count) {
            is = s_count;
         }
         ih = (int)(h/h_space + 0.5);
         ih = (ih+h_count)%h_count;
         i = h_mult*ih + s_mult*(is-1) + v_mult*(iv-1) + hsv_offset;
      }
   }
   return(i);
}

/* x11_setup_rgbs:  setup RGBS system cube */
int x11_setup_rgbs() {
   XColor c;
   int ir,ig,ib,k;
   float rgb[3];

   if (x11.colour_size >= 4096) {
      r_count = 16;
      g_count = 16;
      b_count = 16;
   } else if (x11.colour_size >= 216) {
      r_count = 6;
      g_count = 6;
      b_count = 6;
   } else if (x11.colour_size >=  125) {
      r_count = 5;
      g_count = 5;
      b_count = 5;
   } else if (x11.colour_size >=  64) {
      r_count = 4;
      g_count = 4;
      b_count = 4;
   } else {
      r_count = 3;
      g_count = 3;
      b_count = 3;
   }
#ifdef DEBUG
   printf("x11_setup_rgbs:  colour_size=%d\n", x11.colour_size);
   printf("..rgb[%d,%d,%d]\n", r_count, g_count, b_count);
#endif
/* spacings */
   r_space = 1.0/(float)(r_count-1);
   g_space = 1.0/(float)(g_count-1);
   b_space = 1.0/(float)(b_count-1);
/* index into main matrix */
   r_mult = 1;
   g_mult = b_count;
   b_mult = g_count*b_count;
/* setup rgb cube */
   rgb_offset = x11.colour_offset;
   for (ir=0; ir < r_count; ir++) {
      for (ig=0; ig < g_count; ig++) {
         for (ib=0; ib < b_count; ib++) {
            rgb[0] = (float)ir*r_space;
            rgb[1] = (float)ig*g_space;
            rgb[2] = (float)ib*b_space;
            k = ir*r_mult + ig*g_mult + ib*b_mult + rgb_offset;
            c.pixel = pixels[k];
            c.red = (unsigned short)(rgb[0]*65535.0);
            c.green = (unsigned short)(rgb[1]*65535.0);
            c.blue = (unsigned short)(rgb[2]*65535.0);
            c.flags = DoRed | DoGreen | DoBlue;
            XStoreColor(display, cmap, &c);
#ifdef DEBUG
            printf("..rgb=[%5.2f,%5.2f,%5.2f],k=%d\n", rgb[0],
              rgb[1], rgb[2], k);
#endif
         }
      }
   }
   return 0;
}

/* x11_rgb_index:  get pixel index from RGB */
int x11_rgb_index(float rgb[3]) {
   int i;
   int ir,ig,ib;
   float r,g,b;

   r = rgb[0];
   g = rgb[1];
   b = rgb[2];
   ir = (int)(r/r_space + 0.5);
   if (ir <= 0) {
      ir = 0;
   } else if (ir > r_count-1) {
      ir = r_count-1;
   }
   ig = (int)(g/g_space + 0.5);
   if (ig <= 0) {
      ig = 0;
   } else if (ig > g_count-1) {
      ig = g_count-1;
   }
   ib = (int)(b/b_space + 0.5);
   if (ib <= 0) {
      ib = 0;
   } else if (ib > b_count-1) {
      ib = b_count-1;
   }
   i = r_mult*ir + g_mult*ig + b_mult*ib + rgb_offset;
   return(i);
}
