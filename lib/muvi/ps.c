/* ps.c:  PostScript */
#include <stdio.h>
#include "drp.h"
static struct dev ps;
static bool used=FALSE;

#define NOT_MONO
#define NCMAP 34
#define MAXPT 4096
static FILE *fp;
static bool pending;

static int ps_foreground_index,ps_background_index;
static float ps_greys[NCMAP];
static bool ps_size_rel=TRUE;
static float ps_size_rh=0.5,ps_size_rv=0.5,ps_size_ah,ps_size_av;

/* EXPLICIT DECLARATIONS ADDED BY BM */
int ps_size_r(float rh, float rv);
int ps_size_a(float ah, float av);
int ps_open(char *name, char *type, struct dev *(*p), int *ios);
int ps_setup_window(float hlen,float vlen, float *hsize,float *vsize);
int ps_close(char name[], int *ios);
int ps_select(char name[], int *ios);
int ps_erase();
int ps_rub(float h1, float h2, float v1, float v2); 
int ps_rect(float h1, float h2, float v1, float v2);
int ps_polyline(int na, float ha[], float va[]);
int ps_polygon(int na, float ha[], float va[]);
int ps_text(float h, float v, char str[], char just[2]);
int ps_text_box(char str[], float *dh, float *dv);
int ps_cursor(float *h, float *v, int *button);
int ps_sample(float *h, float *v, int *button);
int ps_set_font();
int ps_set_colours(int offset, int size, float reds[], float greens[], float blues[]);
int ps_set_colour(int index);
int ps_set_rgb(float red, float green, float blue);
int ps_set_rgb_background(float red, float green, float blue);
int ps_set_line_width(float width);
int ps_colour_bitmap(int ih0,int iv0, int nh,int nv, int image[]);
int ps_rgb_bitmap(int ih0,int iv0, int nh,int nv, float rgbim[]);
int  ps_flush();
int ps_prepare();
int ps_setup_colour();
int ps_set_paint(int offset, float red, float green, float blue);
extern void ccopy(char a[], int offset, int n, char b[]);
extern bool csame(char a[], char b[]);
extern int slen(char s[]);
extern void slower(char a[]);






int ps_size_r(rh, rv)
   float rh,rv;
{
   ps_size_rh = rh;
   ps_size_rv = rv;
   ps_size_rel = TRUE;
   return 0;
}

int ps_size_a(ah, av)
   float ah,av;
{
   ps_size_ah = ah;
   ps_size_av = av;
   ps_size_rel = FALSE;
   return 0;
}

int ps_open(name, type, p, ios)
   char *name,*type;
   struct dev *(*p);
   int *ios;
{
   int offset;
   char file[SZLIN];
   float hlen,vlen;

/* return if repeat attempt */
   if (used) {
      printf("error:  postscript plot file already in use\n");
      *ios = ERR;
      return(0);
   }
/* return pointer to "tek" structure */
   ps.count = 1;
   *p = &ps;
/* record name and type */
   ccopy(name,0,SZNAM,ps.name);
   ccopy(type,0,SZNAM,ps.type);
/* fixed-width character dimensions */
   ps.avail_fontscale = TRUE;
   ccopy("Screen", 0, SZNAM, ps.font_style);
   ps.font_size = 10.0;
/* setup window */
   ps_setup_window(560.0, 792.0, &hlen, &vlen);
/* bitmap limits */
   ps.ihmin = 0;
   ps.ivmin = 0;
/*
   ps.ihmax = 72*(int)ps.rhu;
   ps.ivmax = 72*(int)ps.rvu;
*/
   ps.ihmax = (int)hlen;
   ps.ivmax = (int)vlen;
/* coordinate limits */
   ps.hmin = (float)ps.ihmin;
   ps.hmax = (float)ps.ihmax;
   ps.vmin = (float)ps.ivmin;
   ps.vmax = (float)ps.ivmax;
/* rectangles, polygons, double, bitmap */
   ps.avail_rect = TRUE;
   ps.avail_rub = FALSE;
   ps.avail_polygon = TRUE;
   ps.avail_bitmap = FALSE;
   ps.bitmap_dh = 1.0;
   ps.bitmap_dv = 1.0;
/* cursor */
   ps.avail_cursor = FALSE;
   ps.avail_sample = FALSE;
   ps.cursor_buttons = 1;
/* line width */
   ps.line_width = 0.16;
/* COLOUR */
   if (csame("ps_rgb",type)) {
      ps.avail_bitmap = TRUE;
      ps.avail_colour = FALSE;
      ps.change_colour = FALSE;
      ps.avail_rgb = TRUE;
   } else {
      ps.avail_colour = FALSE;
      ps.change_colour = FALSE;
      ps_setup_colour();
      ps.avail_rgb = FALSE;
   }
   ps.avail_rgbz = FALSE;
/* open spool file and plant intro */
   pending = FALSE;
/* construct helpful file name */
   ccopy(name,0,SZNAM,file);
   offset = slen(name);
   ccopy(".ps",offset,SZNAM,file);
   slower(file);
   if ((fp = fopen(file,"w")) == NULL) {
      printf("Error:  creating PostScript spool file.\n");
      *ios = ERR;
   } else {
      fprintf(fp, "/Xmin %d def\n",ps.ihmin);
      fprintf(fp, "/Xmax %d def\n",ps.ihmax);
      fprintf(fp, "/Ymin %d def\n",ps.ivmin);
      fprintf(fp, "/Ymax %d def\n",ps.ivmax);
/*
      fprintf(fp, "/Fh %d def\n",(int)ps.font_size);
*/
      ps_set_font();
      fprintf(fp, "BeginPlot\n");
      used = TRUE;
      *ios = OK;
   }
   return 0;
}

int ps_setup_window(hlen,vlen, hsize,vsize)
   float hlen,vlen, *hsize,*vsize;
{

   if (ps_size_rel) {
      *hsize = hlen*ps_size_rh;
      *vsize = vlen*ps_size_rv;
   } else {
      *hsize = ps_size_ah*72.0;
      *vsize = ps_size_av*72.0;
   }
   if (*hsize > 1.0*hlen) {
      *hsize = 1.0*hlen;
   } else if (*hsize < 0.1*hlen) {
      *hsize = 0.1*hlen;
   }
   if (*vsize > vlen*1.0) {
      *vsize = vlen*1.0;
   } else if (*vsize < 0.1*vlen) {
      *vsize = 0.1*vlen;
   }
   return 0;
}

int ps_close(name, ios)
   char name[];
   int *ios;
{
   if (!csame(name,ps.name)) {
      *ios = ERR;
   } else {
      if (pending) {
         fprintf(fp, "EndPlot\n");
      };
      fclose(fp);
      ps.used = FALSE;
      used = FALSE;
      *ios = OK;
   }
   return 0;
}

int ps_select(name, ios)
   char name[];
   int *ios;
{
   *ios = !csame(name,ps.name);
   return 0;
}

int ps_erase() {
   if (pending) {
      fprintf(fp, "EndPlot\n");
      pending = FALSE;
      fprintf(fp, "/Xmin %d def\n",ps.ihmin);
      fprintf(fp, "/Xmax %d def\n",ps.ihmax);
      fprintf(fp, "/Ymin %d def\n",ps.ivmin);
      fprintf(fp, "/Ymax %d def\n",ps.ivmax);
      fprintf(fp, "BeginPlot\n");
      fprintf(fp, "/Fh %d def\n", (int)ps.font_size);
   }
   return 0;
}

int ps_rub(h1, h2, v1, v2) 
   float h1,h2,v1,v2;
{
   return 0;
}

int ps_rect(h1, h2, v1, v2)
   float h1, h2, v1, v2;
{
   pending = TRUE;
   if (ps.avail_rgb) {
      fprintf(fp,"%.1f %.1f %.1f %.1f Block\n",h1,h2,v1,v2);
/*
   } else {
      fprintf(fp,"%.1f %.1f %.1f %.1f %.2f Rect\n",h1,h2,v1,v2,
        ps_greys[ps.colour_index]);
*/
   }
   return 0;
}

int ps_polyline(na, ha, va)
   float ha[], va[];
   int na;
{
   int i,count;

   pending = TRUE;
   fprintf(fp, "newpath\n");
   fprintf(fp, "  %.1f %.1f moveto\n", ha[0], va[0]);
   count = 0;
   for (i=1; i < na; i++) {
      fprintf(fp, "  %.1f %.1f lineto\n", ha[i], va[i]);
      count++;
      if ((count == 128) && ((i+1) < na)) {
         fprintf(fp, "stroke\n");
         fprintf(fp, "newpath\n");
         fprintf(fp, "  %.1f %.1f moveto\n", ha[i], va[i]);
         count = 0;
      }
   }
   fprintf(fp, "stroke\n");
   return 0;
}

int ps_polygon(na, ha, va)
   float ha[], va[];
   int na;
{
   int i;

   pending = TRUE;
   fprintf(fp, "newpath\n");
   fprintf(fp, "  %d %d moveto\n", (int)ha[0], (int)va[0]);
   for (i=1; i < na; i++) {
      fprintf(fp, "  %d %d lineto\n", (int)ha[i], (int)va[i]);
   }
   fprintf(fp, "closepath fill\n");
   return 0;
}

int ps_text(h, v, str, just)
   float h, v;
   char str[],just[2];
{
   int hj,vj;

   hj = 0;
   vj = 0;
   if (just[0] == 'r') {
      hj = 2;
   } else if (just[0] == 'm') {
      hj = 1;
   }
   if (just[1] == 't') {
      vj = 2;
   } else if (just[1] == 'm') {
      vj = 1;
   }
   fprintf(fp, "%d %d %d %d (%s) JustShow\n", (int)h, (int)v, hj, vj,
     str);
   return 0;
}

int ps_text_box(str, dh, dv)
   char str[];
   float *dh,*dv;
{
   int nchar;

   nchar = slen(str);
   *dh = (float)(nchar)*ps.char_width + (float)(nchar-1)*ps.char_space;
   *dv = ps.char_height;
   return 0;
}

int ps_cursor(h, v, button)
   float *h, *v;
   int *button;
{
   *h = ps.hmin;
   *v = ps.vmin;
   *button = 0;
   return 0;
}

int ps_sample(h, v, button)
   float *h,*v;
   int *button;
{
   *h = ps.hmin;
   *v = ps.vmin;
   *button = 0;
   return 0;
}

int ps_set_font() {
   fprintf(fp, "/Fh %d def\n",(int)ps.font_size);
   ps.char_height = ps.font_size;
   ps.char_width =  0.7*ps.char_height;
   ps.char_space = 0.0;
   if (csame("Screen", ps.font_style)) {
      ccopy("Courier", 0, SZNAM, ps.font_name);
   } else if (csame("Times", ps.font_style)) {
      ccopy("Times-Roman", 0, SZNAM, ps.font_name);
   } else {
      ccopy("Courier", 0, SZNAM, ps.font_name);
   }
   fprintf(fp, "/%s findfont Fh scalefont setfont\n", ps.font_name);
   return 0;
}

int ps_set_colours(offset, size, reds, greens, blues)
   int offset,size;
   float reds[],greens[],blues[];
{
   int i,k;

   if (ps.change_colour) {
      k = offset;
      for (i=0; i < size; i++) {
         ps.reds[k] = reds[i];
         ps.greens[k] = greens[i];
         ps.blues[k] = blues[i];
         k++;
      }
   }
   return 0;
}

int ps_set_colour(index)
   int index;
{
   if (ps.avail_colour) {
      ps.colour_index = index;
      ps.red = ps.reds[index];
      ps.green = ps.greens[index];
      ps.blue = ps.blues[index];
   }
   return 0;
}

int ps_set_rgb(red, green, blue)
   float red,green,blue;
{
   ps.red = red;
   ps.green = green;
   ps.blue = blue;
   if (ps.avail_rgb) {
      fprintf(fp, "%4.2f %4.2f %4.2f setrgbcolor\n", red, green, blue);
   } else if (ps.change_colour) {
      ps_set_paint(ps.colour_index, red, green, blue);
   }
   return 0;
}

int ps_set_rgb_background(red, green, blue)
   float red,green,blue;
{
   ps.red_background = red;
   ps.green_background = green;
   ps.blue_background = blue;
   if (ps.change_colour) {
      ps_set_paint(ps_background_index, red, green, blue);
   }
   return 0;
}

int ps_set_line_width(width)
   float width;
{
   ps.line_width = width;
   return 0;
}

int ps_colour_bitmap(ih0,iv0, nh,nv, image)
   int ih0,iv0,nh,nv;
   int image[];
{
   int kh,kv,k;
   float h0,v0,h1,h2,v1,v2,dh,dv;

   pending = TRUE;
   dh = ps.bitmap_dh;
   dv = ps.bitmap_dv;
   h0 = (float)(ih0-ps.ihmin)*dh + ps.hmin;
   v0 = (float)(iv0-ps.ivmin)*dv + ps.vmin;
   k = 0;
   for (kv=0; kv < nv; kv++) {
      for (kh=0; kh < nh; kh++) {
         h1 = (float)kh*dh + h0;
         v1 = (float)kv*dv + v0;
         h2 = h1 + dh;
         v2 = v1 + dv;
         fprintf(fp,"%.1f %.1f %.1f %.1f %.2f Rect\n",h1,h2,v1,v2,
           ps_greys[image[k]]);
         k++;
      }
   }
   return 0;
}

int ps_rgb_bitmap(ih0,iv0, nh,nv, rgbim)
   int ih0,iv0,nh,nv;
   float rgbim[];
{
   int kh,kv,ir,ig,ib,k;
   float h0,v0;
   if (ps.avail_rgb) {
      fprintf(fp, "gsave\n");
/*   A string to contain a line of rgb values */
#ifdef MONO
      fprintf(fp, "/picstr %d string def\n", nh);
#else
      fprintf(fp, "/picstr %d string def\n", nh*3);
#endif
/*   Translate to bottom-left and scale */
      h0 = (float)ih0;
      v0 = (float)iv0;
      fprintf(fp, "%f %f translate\n", h0, v0);
      fprintf(fp, "%d %d scale\n", nh, nv);
/*   The command:  nh nv nbit map proc multiproc ncolour colorimage */
      fprintf(fp, "%d %d %d\n", nh, nv, 8);
      fprintf(fp, "[%d 0 0 %d 0 0]\n", nh, nv);
      fprintf(fp, "{currentfile picstr readhexstring pop}\n");
#ifdef MONO
      fprintf(fp, "image\n");
#else
      fprintf(fp, "false 3\n");
      fprintf(fp, "colorimage\n");
#endif
/*   The image, one line at a time */
      k = 0;
      for (kv=0; kv < nv; kv++) {
         for (kh=0; kh < nh; kh++) {
            ir = (int)(255.0*rgbim[k++]);
            if (ir < 0) {
               ir = 0;
            } else if (ir > 255) {
               ir = 255;
            }
            ig = (int)(255.0*rgbim[k++]);
            if (ig < 0) {
               ig = 0;
            } else if (ig > 255) {
               ig = 255;
            }
            ib = (int)(255.0*rgbim[k++]);
            if (ib < 0) {
               ib = 0;
            } else if (ib > 255) {
               ib = 255;
            }
#ifdef MONO
            ir = (ir + ig  + ib)/3;
            fprintf(fp, "%02x", ir);
#else
            fprintf(fp, "%02x%02x%02x", ir, ig, ib);
#endif
         }
         fprintf(fp, "\n");
      }
      fprintf(fp, "grestore\n");
   }
   return 0;
}

int  ps_flush() {
   return 0;
}

int ps_prepare()
{
/* default window position-size */
   ps_size_rel = FALSE;
   ps_size_rh = 6.0;
   ps_size_rv = 4.0;
/* routine bindings (see "drp.h") */
   ps.prepare = ps_prepare;
   ps.open = ps_open;
   ps.close = ps_close;
   ps.select = ps_select;
   ps.erase = ps_erase;
   ps.rub = ps_rub;
   ps.rect = ps_rect;
   ps.polyline = ps_polyline;
   ps.polygon = ps_polygon;
   ps.text = ps_text;
   ps.text_box = ps_text_box;
   ps.cursor = ps_cursor;
   ps.sample = ps_sample;
   ps.set_font = ps_set_font;
   ps.set_colours = ps_set_colours;
   ps.set_colour = ps_set_colour;
   ps.set_rgb = ps_set_rgb;
   ps.set_rgb_background = ps_set_rgb_background;
   ps.set_line_width = ps_set_line_width;
   ps.colour_bitmap = ps_colour_bitmap;
   ps.rgb_bitmap = ps_rgb_bitmap;
   ps.gflush = ps_flush;
   return 0;
}


/* PRIVATE ROUTINES */

int ps_setup_colour() {
   int i,k;
   float di;

/* table position */
   ps.colour_offset = 2;
   ps.colour_size = NCMAP-2;
/* background (white) */
   ps_background_index = 0;
   ps_greys[0] = 1.0;
   ps.reds[ps_background_index] = 1.0;
   ps.greens[ps_background_index] = 1.0;
   ps.blues[ps_background_index] = 1.0;
/* foreground (black) */
   ps_greys[1] = 0.0;
   ps_foreground_index = 1;
   ps.reds[ps_foreground_index] = 0.0;
   ps.greens[ps_foreground_index] = 0.0;
   ps.blues[ps_foreground_index] = 0.0;
   ps.colour_index = ps_foreground_index;
/* remaining colours are a grey scale */
   k = ps.colour_offset;
   di = 1.0/(float)(ps.colour_size-1);
   for (i=0; i < ps.colour_size; i++) {
      ps_greys[k] = di*(float)i;
      ps.reds[k] = ps_greys[k];
      ps.greens[k] = ps_greys[k];
      ps.blues[k] = ps_greys[k];
      k++;
   }
   return 0;
}

int ps_set_paint(offset, red, green, blue)
   int offset;
   float red, green, blue;
{
   ps.reds[offset] = red;
   ps.greens[offset] = green;
   ps.blues[offset] = blue;
   return 0;
}

