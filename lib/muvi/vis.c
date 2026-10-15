/* vis.c:  visualisation package

     This module manages display surfaces and subwindows, and produces
     graphs in them from data supplied from above.  Scalar and vector
     fields are held in the grid/scatter_vis structures, whereas csg
     datasets are passed down from above?????
*/
#include <stdio.h>

#include <math.h>
#include <stdlib.h>
#include "const.h"
#include "mc.h"
#include "drp.h"
#include "csg.h"
#include "dr3.h"
#include "vis.h"
#define HSV_RAIN
#define CURV98CONT
#define MAXDEV 6
#define MAXSUB 7
#define MAXCAP 3
#define MAXCONTOUR 50
#define MAXTAG 256
#define MAXBAND 3
#define MAXKEY 32
#define MAXPAL 8

/* define MARIO for GRC, undef it for MANTIS */
#define MARIO

struct vis_sub {
   char name[SZNAM],kind[SZNAM];
   bool used,clear,drawn;
   float rhmin,rhmax,rvmin,rvmax,hmin,hmax,vmin,vmax;
   int ndc,nd;
   char caption[MAXCAP][SZLIN],gname[MAXDT][SZNAM],glabel[MAXDT][SZLIN],
     key_position[SZNAM],qlabel[SZLIN],qunits[SZLIN];
   float gcmin[MAXDT],gcmax[MAXDT],gasp[MAXDT],caption_h,caption_v,
     hgmin,hgmax,vgmin,vgmax,hdmin,hdmax,vdmin,vdmax,
     x0,dhdx,h0,y0,dvdy,v0,wcmin[MAXDT],wcmax[MAXDT],
     qmin,qmax,velmax,dldv,contour_levels[MAXCONTOUR],
     hsmin,hsmax,vsmin,vsmax,qscale,hkmin,hkmax,vkmin,vkmax,
     hxmin,hxmax,vxmin,vxmax, hzmin,hzmax,vzmin,vzmax,
     htmin,htmax,vtmin,vtmax;
   bool gscalp,gspace[MAXDT],log_scalar,nocmap,novel,nocont,noqr,nokey,
     noscales,scales_drawn,cmap_scales_drawn,contour_scales_drawn,
      key_drawn,mat_scales_drawn,vel_scales_drawn;
   int caption_count,caption_lines,ncontour,ncolour,key_size,
     key_label_field,key_number_field,key_count,qfield,qdp,qpower,
     scale_size,scale_count;
   bool noscene;
   scene_struct scene;
#ifdef CURV
   bool curv;
#endif
};

struct vis_surface {
   char name[SZNAM],type[SZNAM];
   bool used,avail_rub;
   struct vis_sub sub[MAXSUB];
   int sc;
   float hmin,hmax,vmin,vmax;
   float char_height,char_width,char_depth;
   bool colour_avail,cmap_avail,rgb_avail;
   int colour_black,colour_sky;
};

static bool vis_pos_rel,vis_size_rel=TRUE;
#ifndef NOGREY
static bool vis_grey_rain=FALSE;
#endif
static float vis_pos_rh,vis_pos_rv,vis_pos_ah,vis_pos_av,
  vis_size_rh,vis_size_rv,vis_size_ah,vis_size_av;
static struct vis_surface surface[MAXDEV], *wp;
static struct vis_sub *sp;
static scene_struct *scp;
static int wc=(-1);

static bool nobox,dolog,noqr,noqmin,noqmax,
  notr,notmin,notmax,noxr[MAXD],noxmin[MAXD],noxmax[MAXD],
  novmax,noscalp,noband,nocolour,nocontour,nocontour_levels,
  noarrow,nostream,nocaption,noxlabel[MAXD],noqlabel,notlabel,
  caption_below,label_centre,noline_marker,noline_connect,
  noarrow_marker,nokey,noscales,nocontour_label,
  noview_focus,noview_distance,nolight_distance,
  nobandmin[MAXBAND],nobandmax[MAXBAND],nobandr[MAXBAND],
  nokey_label,noview_near,noview_far,data_shade=TRUE,data_plane=FALSE,
  nodata_distance=TRUE,data_ortho=FALSE,
  nodata_coord=TRUE,wire_frame=FALSE,data_frame=TRUE;
static float qmin,qmax,tmin,tmax,xmin[MAXD],xmax[MAXD],vmax,
  aspect[MAXD],hborder,vborder,tick_scale,border,
  arrow_scale,contour_levels[MAXCONTOUR], arrow_ratio,
  stream_x[MAXD],stream_dx[MAXD],view_focus[MAXD],
  view_anglev,view_aspect,view_azim,view_altd,view_distance,
  light_fusion,light_azim,light_altd,light_distance,
  material_red[MAXMAT],material_green[MAXMAT],material_blue[MAXMAT],
  material_shine[MAXMAT],sky_red,sky_green,sky_blue,
  bandmin[MAXBAND],bandmax[MAXBAND],view_near,view_far,
  data_distance,data_coord,surface_height=0.5,
  palet_r[MAXPAL],palet_g[MAXPAL],palet_b[MAXPAL];
static int number_field,label_field,caption_lines,caption_inset,
  caption_count,line_marker,line_style,arrow_marker,arrow_rate,
  ncontour,ncolour,key_label_field,key_number_field,key_size,
  stream_size,nvmap,view_lines,material_count,
  tag_count,material_phong[MAXMAT],ngrey,data_ortho_id=2,nshade,
  npalet,line_palet;
static char caption[MAXCAP][SZLIN],xlabel[MAXD][SZLIN],tlabel[SZLIN],
  qlabel[SZLIN],key_position[SZNAM],material_name[MAXMAT][SZNAM],
  tag_name[MAXTAG][SZNAM], tag_material[MAXTAG][SZNAM],
  key_label[SZNAM], palet_name[MAXPAL][SZNAM];

static bool nocsg=TRUE,nogrid=TRUE,noscatter=TRUE,nopic=TRUE;
static csg_struct *vis_cp;
static grid_vis_struct *vis_gvp;
static scatter_vis_struct *vis_svp;
static int glevel=0;



/* EXPLICIT DECLARATIONS ADDED BY BM */
bool vis_material_find(char name[], int *material_index);
bool vis_tag_find(char name[], int *tag_index);
bool vis_palet_find(char name[], int *palet_index);
void vis_prepare();
void vis_reset();
void vis_flush();
void vis_print();
void vis_show();
void vis_open(char name[], char type[], int *ios);
void vis_select(char name[], int *ios);
void vis_close(char name[], int *ios);
void vis_clear();
void vis_active(bool *active);
void vis_open_sub(char name[], float  rhmin, float rhmax,
  float rvmin, float rvmax, int *ios);
void vis_select_sub(char name[], int *ios);
void vis_close_sub(char name[], int *ios);
void vis_erase();
void vis_free_sub(struct vis_sub *zp);
void vis_get_sub(char name[]);
void vis_get_hv(float *h1,float *h2,float *v1,float *v2);
void vis_begin();
void vis_end();
void vis_free();
void vis_csg(csg_struct *cp, int *ios);
void vis_grid(grid_vis_struct *gvp, int *ios);
void vis_scatter(scatter_vis_struct *svp, int *ios);
void vis_set0s(char name[], char value[], int *ios);
void vis_set0i(char name[], int value, int *ios);
void vis_set0r(char name[], float value, int *ios);
void vis_set0l(char name[], bool value, int *ios);
void vis_seta(char name[], int *ios);
void vis_set1r(char name[], int nvalue, float values[],
  int *ios);
void vis_get0s(char name[], int lenvalue, char value[], int *ios);
void vis_get0i(char name[], int *value, int *ios);
void vis_get0r(char name[], float *value, int *ios);
void vis_get0l(char name[], bool *value, int *ios);
void vis_view(int *ios);
#ifdef CURV
  float float_min(float a,float b);
  float float_max(float a,float b);
#endif
void vis_axes(int *ios);
void vis_design_2d_axes(int *ios);
void vis_design_3d_axes(int *ios);
void vis_divs(float vmin,float vmax,int max_div,float *voff,float *vdiv);
void vis_ndiv(float vmin,float vmax,float voff,float vdiv,int *ndiv);
void vis_design_scales();
void vis_scales();
void vis_design_key();
void vis_key();
void vis_draw(int *ios);
void vis_line_graph(int n, float xa[], int sa[], float qa[], int *ios);
void vis_contour_map(int nx, int ny, float xa[], float ya[], int sa[], float qa[],
  float qamin, float qamax, char qlabel[], char qunits[], int *ios);
void vis_colour_map(int nx, int ny, float xe[], float ye[], int sa[], float qa[],
  float qamin, float qamax, char qlabel[], char qunits[], int *ios);
void vis_band_map(int nx, int ny, float xe[], float ye[], float xa[], float ya[],
  int sa[], float qa[], float qamin, float qamax, char qlabel[], char qunits[], int *ios);
void vis_design_cmap(int nc, int *offset);
void vis_design_bmap(int *offset);
#ifdef HSV_RAIN
  void vis_h_to_rgb(float hue, float rgb[3]);
#endif
void vis_rain(int nc, int offset);
void vis_arrow_map(int nx, int ny, float xa[], float ya[], int sa[], int ndv, float va[], float velmax[], int *ios);
void vis_streamer(int nd, int nc[], float *xe[], int sa[], int ndv, float va[], float vamax[], int *ios);
void vis_track(int nd, int nc[], float *xe[], int sa[], int ndv, float va[], float vamin, float x[]);
void vis_colour_gram(int npt, int nd, float xpt[], int spt[], float qpt[],
  float qamin, float qamax, char qlabel[], char qunits[], int *ios);
void vis_marker_gram(int npt, int nd, float xpt[], int spt[], float qpt[], float qamin, float qamax, int *ios);
void vis_arrow_gram(int npt, int nd, float xpt[], int spt[], int ndv, float vpt[], float velmax[], int *ios);
void vis_arrow(float x, float y, float l, float ax, float ay);
void vis_bound_map(int nx, int ny, float xe[], float ye[], int sa[], int *ios);
void vis_grid_map(int nx, int ny, float xe[], float ye[], int sa[], int *ios);
void vis_design_3d_view(int *ios);
void vis_3d_grid_surface(int nd, int nc[], float *xe[], float *xa[],
  int sa[], float qa[], float qamin, float qamax, int *ios);
void vis_3d_draw_grid(int nc[], float *xe[], int sa[], int *ios);
void vis_grid_facet(int id, int side, int con[MAXD], float *xe[]);
void vis_left(int nc[], int k, int id, int con[], int *kl);
void vis_right(int nc[], int k, int id, int con[], int *kr);
void vis_3d_grid_vector(int nc[], float *xe[], int sa[], float va[], float vamax, int *ios);
void vis_3d_data_shade(int *ios);
void vis_3d_data_plane(int *ios);
void vis_3d_data_ortho(int nc[], float *xe[], int sa[], float qa[], float qamin, float qamax, int *ios);
void vis_3d_scatter_scalar(int npt,int nd,float xpt[],int spt[],float qpt[],float qamin,float qamax,int *ios);
void vis_3d_scatter_vector(int *ios);
void vis_3d_solid(int *ios);
void vis_mat_scales();
void vis_ld_scales();
void vis_zap_scene();
void vis_dump_pic(char file[], int *ios);
void vis_load_pic(char file[], int *ios);
void vis_zap_image();
void vis_bp();
void vis_ap(float x,float y);
void vis_ep();
void vis_marker(int mark, float x, float y);
void vis_text(float x, float y, char str[], char just[]);
#ifndef CURV
  void vis_cmap_rect(float x1, float x2, float y1, float y2, int ci);
#endif
#ifndef CURV
void vis_rect(float x1, float x2, float y1, float y2);
#endif
#ifndef CURV
void vis_black_rect(float x1, float x2, float y1, float y2);
#endif
#ifndef CURV
void vis_grid_rect(float x1, float x2, float y1, float y2);
#endif
#ifndef CURV
void vis_grid_cross(float x1, float x2, float y1, float y2);
#endif
void vis_box(float h1, float h2, float v1, float v2);
void vis_redesign_3d_axes();
void vis_rain_dance(float z, float *r, float *g, float *b);
void vis_default_palet();
void vis_select_palet();
#ifdef CURV
void vis_fill_rect(float x1, float x2, float y1, float y2);
#endif
#ifdef CURV
void vis_draw_rect(float x1, float x2, float y1, float y2);
#endif
#ifdef CURV
void vis_draw_cross(float x1, float x2, float y1, float y2);
#endif
void vis_2d_curv_colour_map(int nc[], float *xe[], int sa[], float qa[], float qamin, float qamax, int *ios);
void vis_2d_curv_grid_map(int nc[], float *xe[], int sa[], int *ios);
void vis_expand(int kad[MAXD],float z0, float z1, float z2, float xc[MAXD]);
#ifdef CURV98CONT
void vis_2d_curv_contour_map(int nx, int ny, float xa[], float ya[], int sa[], float qa[],
  float qamin, float qamax, char qlabel[], char qunits[], int *ios);
#endif

/* --- EXTERNAL FUNCTIONS --- */
extern void ccopy(char a[], int offset, int n, char b[]);
extern int slen(char s[]);
extern int rem(int a,int b);
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
extern int drp_gflush();
extern int dr3_setup_rand();
extern void dr3_create_scene(float anglev, float azim, float altd, float twist, float distance, float f[MAXD], float near, float far, float hmin, float hmax, float vmin, float vmax);
extern void dr3_fresh_scene();
extern void dr3_zap_scene();
extern void dr3_view(float azim, float altd, float twist, float distance, float f[MAXD]);
extern void dr3_project(float anglev, float near, float far, float hmin, float hmax, float vmin, float vmax);
extern void dr3_put_scene(scene_struct *scp);
extern void dr3_get_scene(scene_struct *scp);
extern void dr3_pop_image(scene_struct *scp);
extern void dr3_push_image(scene_struct *scp);
extern void dr3_create_image(int *ios);
extern void dr3_alloc_image(int *ios);
extern void dr3_free_image();   
extern void dr3_clear_image();
extern void dr3_shade(char name[], int *ios);
extern void dr3_get_nl(int *n);
extern void dr3_get_nd(int *n);
extern void dr3_sky(float red, float green, float blue);
extern void dr3_lights(float azim, float altd, float distance, float fusion);
extern void dr3_zap_materials();   
extern void dr3_define_material(char material[], float red, float green, float blue, float rough, float glow);
extern void dr3_select_material(char material[]);
extern void dr3_begin_solid();
extern int dr3_end_solid();
#ifdef STL
  extern void dr3_cut_solid(char shape[], float par[], int npar);
#else
  extern void dr3_cut_solid(char shape[], float par[]);
#endif
#ifdef STL
  extern void dr3_add_solid(char shape[], float par[], int npar);
#else
  extern void dr3_add_solid(char shape[], float par[]);
#endif
//extern void dr3_alloc_nibs();
//extern void dr3_add_nibs();
//extern void dr3_cut_nibs();
//extern void dr3_fix_nib(int n, float d[], int op[], int c[], int *nnn);
//extern void dr3_render_nibs();
//extern void dr3_sort_nib(int n, float d[], int op[], int c[]);
//extern void dr3_pop_nib(nib_struct *nib, int *n, float d[], int op[], int c[]); 
//extern void dr3_push_nib(int n, float d[], int op[], int c[], nib_struct *(*nib));
//extern void dr3_free_nibs();   
extern int dr3_begin_polygon();
extern void dr3_polygon_vertex(float x,float y,float z);
extern void dr3_polygon_normal(float nx,float ny,float nz);
extern void dr3_polygon_scalar(float s);
extern void dr3_end_polygon();
extern void dr3_zap_polylines();
extern void dr3_display_polylines();
extern void dr3_hide_polylines();
extern void dr3_begin_polyline();
extern void dr3_polyline_vertex(float x,float y,float z);
extern void dr3_end_polyline();
extern void dr3_push_matrix();
extern void dr3_pop_matrix();
extern void dr3_get_matrix(float a[MAXD4][MAXD4]);
extern void dr3_put_matrix(float a[MAXD4][MAXD4]);
extern void dr3_zap_matrix();
extern void dr3_rotran_matrix(float x, float y, float z, float azim, float altd, float twist);
extern void dr3_rot_matrix(int id, float angled);
extern void dr3_tran_matrix(float tx, float ty, float tz);
extern void dr3_study_polygon();
extern void dr3_begin_surface();
extern void dr3_end_surface();
extern void dr3_render_surface();
extern void dr3_render_solid();
extern void dr3_shade_cm(float h, float v, float d, float c[], float *dimk, int *cimk);
extern void dr3_shade_rgb(float h, float v, float d, float c[], float *dimk, float rgbimk[3]);
extern void dr3_shade_rgbz(float xv[MAXD4], float rgb[3]);
extern void dr3_light(float xv[MAXD4], float normal[MAXD], float *light);
extern void dr3_wtod(float x,float y,float z, float *h, float *v, float *d);
extern void dr3_get_class(int *class);
extern void dr3_rain(float phi, float *red, float *green, float *blue);
extern void dr3_scalar(float xw[MAXD4], float *s, bool *found);
extern void dr3_vtod(float xv[MAXD4], float *h, float *v, float *d);
extern void dr3_wtov(float xo[MAXD4], float xv[MAXD4]);
extern void dr3_dtov(float h,float v,float d, float xv[MAXD4]);
extern void dr3_vtow(float xv[MAXD4], float xw[MAXD4]);
extern void dr3_clock(float *cp);
extern void dr3_loads(int nc[], float *xe[MAXD4], int sa[], float qa[], float qmin, float qmax, int *ios);
extern void dr3_zaps();
extern void dr3_data_plane(float distance, int *ios);
extern void dr3_data_shade(int *ios);
extern void dr3_coerce_wire(int *ios);
extern void dr3_coerce_ld(int *ios);
extern void dr3_coerce_lm(int *ios);
extern void dr3_coerce_d(int *ios);
extern void dr3_seta(char name[]);
extern void dr3_set0i(char name[], int i);
extern void dr3_get_shade(scene_struct *scp);
extern void dr3_get_image_type(scene_struct *scp);
extern void dr3_default_limits(int nh,int nv, int *ivmin,int *ivmax,int ihmin[], int ihmax[]);
extern void dr3_extend_limits(int mh, int mv, int jvmin, int jvmax, int jhmin[], int jhmax[],
  int *ivmin, int *ivmax, int ihmin[], int ihmax[]);
extern void dr3_dump_scene(int unit, int *ios);
extern void dr3_load_scene(int unit, int *ios);
extern void dr3_wire_icon(float size, float x, float y, float z);
extern void dr3_solid_icon(float s, float size, float x, float y, float z);
extern void dr3_begin();
extern void dr3_end();
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
int ds_get1c(int unit, char name[], int nvalue, char value[], int *ios);
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
extern void user_setup_model(char geom[], int nd, char co[], char un[], 
  char phys[], int nq, char vl[], char vu[], int vmlt[][3], int *ios);
extern void user_check_scalar(char name[], int *ios);
extern void user_get_scalar(float q[], char name[], float *value);
extern bool geo_met(char geom[], int nd);
extern void geo_x(int lumd, float *z, float *x);
extern void contour_map(int nx, int ny, float xa[], float ya[], int sa[], float qa[], float qlevel, char qlab[], int *ios);
extern void dr3_wire_block(float xmin, float xmax, float ymin, float ymax,
  float zmin, float zmax);




bool vis_material_find(name, material_index)
   char name[];
   int *material_index;
{
   int i;

   *material_index = -1;
   for (i=0; i < material_count; i++) {
      if (csame(name,material_name[i])) {
         *material_index = i;
         break;
      }
   }
   return(*material_index >= 0);
}

bool vis_tag_find(name, tag_index)
   char name[];
   int *tag_index;
{
   int i;

   *tag_index = -1;
   for (i=0; i < tag_count; i++) {
      if (csame(name,tag_name[i])) {
         *tag_index = i;
         break;
      }
   }
   return(*tag_index >= 0);
}

bool vis_palet_find(name, palet_index)
   char name[];
   int *palet_index;
{
   int i;

   *palet_index = -1;
   for (i=0; i < npalet; i++) {
      if (csame(name,palet_name[i])) {
         *palet_index = i;
         break;
      }
   }
   return(*palet_index >= 0);
}


/* SETUP */

void vis_prepare() {
   int i;

   for (i=0; i < MAXDEV; i++) {
      surface[i].used = FALSE;
   }
   vis_pos_rel = TRUE;
   vis_pos_rh = 1.0;
   vis_pos_rv = 1.0;
   vis_size_rel = TRUE;
   vis_size_rh = 0.5;
   vis_size_rv = 0.5;
   vis_reset();
   drp_prepare();
}

void vis_reset() {
   int id,i;

/* overlay control */
   nobox = TRUE;
/* coordinate,time and value range control */
   dolog = FALSE;
   noqr = TRUE;
   noqmin = TRUE;
   noqmax = TRUE;
   noqlabel = TRUE;
   notr = TRUE;
   notmin = TRUE;
   notmax = TRUE;
   notlabel = TRUE;
   nocaption = TRUE;
   for (id=0; id < MAXD; id++) {
      noxr[id] = TRUE;
      noxmin[id] = TRUE;
      noxmax[id] = TRUE;
      aspect[id] = 1.0;
      noxlabel[id] = TRUE;
   }
   novmax = TRUE;
/* control over the kind of display methods to be used */
   noscalp = TRUE;
/* colour bands */
   noband = TRUE;
   for (i=0; i < MAXBAND; i++) {
      nobandmin[i] = TRUE;
      nobandmax[i] = TRUE;
      nobandr[i] = TRUE;
   }
/* colour map */
   nocolour = FALSE;
   ncolour = 10;
   nshade = 10;
/* greyscales (lights) */
   ngrey = 10;
/* contour map */
   nocontour = FALSE;
   ncontour = 12;
   nocontour_levels = TRUE;
   nocontour_label = FALSE;
/* arrow map */
   noarrow = FALSE;
   nvmap = 7;
/* streamline */
   nostream = TRUE;
/* annotation control */
   hborder = 0.0;
   vborder = 0.0;
   tick_scale = 0.5;
   border = 0.02;
   caption_below = TRUE;
   label_centre = TRUE;
   caption_lines = 1;
   caption_inset = 8;
   number_field = 8;
   label_field = 12;
   nokey = FALSE;
   key_size = 4;
   nokey_label = TRUE;
   ccopy(" ", 0, SZNAM, key_label);
   key_label_field = 8;
   key_number_field = 6;
   ccopy("default", 0, SZNAM, key_position);
   noscales = FALSE;
/* line-graphs */
   noline_marker = TRUE;
   line_marker = 1;
   noline_connect = FALSE;
   line_style = 1;
   line_palet=0;
/* arrow-map\/gram */
   noarrow_marker = TRUE;
   arrow_scale = 1.0;
   arrow_rate = 1;
   arrow_ratio = 0.025;
   arrow_marker = 4;
/* streamline-map */
   stream_size = 1;
   for (id=0; id < MAXD; id++) {
      stream_x[id] = 0.0;
      stream_dx[id] = 0.0;
   }
/* viewing in 3-dimensions */
   noview_focus = TRUE;
   noview_distance = TRUE;
   noview_near = TRUE;
   noview_far = TRUE;
   view_anglev = 10.0;
   view_aspect = 1.01;
   view_lines = 21;
   view_azim = 0.0;
   view_altd = 0.0;
   light_fusion = 0.5;
   light_altd = 90.0;
   light_azim = 0.0;
   nolight_distance = TRUE;
/* sky colour */
   sky_red = 0.8;
   sky_green = 0.8;
   sky_blue = 1.0;
/* materials */
   material_count = 1;
   ccopy("DEFAULT", 0, SZNAM, material_name[0]);
   material_red[0] = 1.0;
   material_green[0] = 1.0;
   material_blue[0] = 1.0;
   material_shine[0] = 0.0;
   material_phong[0] = 2;
/* tag-material maps */
   tag_count = 0;
/* palet (generic) */
   vis_default_palet();
}

void vis_flush() {
   int i,ios;

   for (i=0; i < MAXDEV; i++) {
      if (surface[i].used) {
         vis_close(surface[i].name, &ios);
      }
   }
}

void vis_print() {
   struct vis_sub *zp;
   char font_style[SZNAM];
   float font_size;
   int i,line_styles;

   printf("Request:-\n");
   if (vis_pos_rel) {
      printf("  Relative position in screen [%f,%f]\n", vis_pos_rh,
        vis_pos_rv);
   } else {
      printf("  Absolute position in screen [%f,%f]\n", vis_pos_ah,
        vis_pos_av);
      printf("  ([left,bottom] in inches)\n");
   }
   if (vis_size_rel) {
      printf("  Relative size in screen [%f,%f]\n", vis_size_rh,
        vis_size_rv);
   } else {
      printf("  Absolute size in screen [%f,%f]\n", vis_size_ah,
        vis_size_av);
      printf("  (in inches)\n");
   }
   printf("Worksurfaces:-\n");
   for (i=0; i < MAXDEV; i++) {
      if (surface[i].used) {
         printf("  Name='%s',Type='%s'\n", surface[i].name,
           surface[i].type);
      }
   }
   if (wc < 0) {
      printf("  No worksurface selected\n");
   } else {
      printf("  Selected worksurface=%s,type=%s\n", wp->name, wp->type);
      printf("    Subwindows:-\n");
      for (i=0; i < MAXSUB; i++) {
         zp = &(wp->sub[i]);
         if (zp->used) {
            printf("      Name='%s',Port=[%4.2f,%4.2f,%4.2f,%4.2f]\n",
               zp->name, zp->rhmin, zp->rhmax, zp->rvmin, zp->rvmax);
         }
      }
/*   worksurface capabilities */
      printf("    Screen hmin=%f,hmax=%f,vmin=%f,vmax=%f\n",
        wp->hmin, wp->hmax, wp->vmin, wp->vmax);
      drp_get_font_style(font_style);
      drp_get_font_size(&font_size);
      printf("    Font style=%s,size=%f\n", font_style, font_size);
      printf("    Character height=%f,width=%f,depth=%f\n",
        wp->char_height, wp->char_width, wp->char_depth);
      drp_get_line_styles(&line_styles);
      printf("    Line styles=%d\n", line_styles);
      if (wp->avail_rub) {
         printf("    Block erase available.\n");
      } else {
         printf("    No block erase available\n");
      }
      if (wp->colour_avail) {
         printf("    Colour available");
         printf("    Black colour index=%d\n", wp->colour_black);
         printf("    Sky colour index=%d\n", wp->colour_sky);
         if (wp->cmap_avail) {
            printf("    Colour-map available\n");
         } else {
            printf("    No colour-map available.\n");
         }
      } else {
         printf("    No colour available\n");
      }
   }
}

void vis_show() {
   char xnames[MAXD];
   int i,id;

   printf("Visualisation details:-\n");
   printf("  Options:-\n");
   if (nobox) {
      printf("    nobox");
   } else {
      printf("    box");
   }
   if (noscalp) {
      printf(",noscalp");
   } else {
      printf(",scalp");
   }
   if (noband) {
      printf(",noband");
   } else {
      printf(",band");
   }
   if (nocontour) {
      printf(",nocontour");
   } else {
      printf(",contour");
   }
   if (nocolour) {
      printf(",nocolour");
   } else {
      printf(",colour");
   }
   if (noarrow) {
      printf(",noarrow");
   } else {
      printf(",arrow");
   }
   if (nostream) {
      printf(",nostream");
   } else {
      printf(",stream");
   }
   if (noscales) {
      printf(",noscales");
   } else {
      printf(",scales");
   }
   if (nokey) {
      printf(",nokey");
   } else {
      printf(",key");
   }
   printf("\n");
   if (nocaption) {
      printf("    nocaption\n");
   } else {
      printf("    Caption:-\n");
      for (i=0; i < caption_count; i++) {
         printf("      %d: %s\n", i, caption[i]);
      }
      printf("    Max. caption lines %d\n", caption_lines);
   }
   if (noqr) {
      printf("    noqr");
   } else {
      printf("    qmin=%E,qmax=%E",qmin,qmax);
   }
   if (noqlabel) {
      printf(",noqlabel\n");
   } else {
      printf(",qlabel='%s'\n",qlabel);
   }
   if (novmax) {
      printf("    novmax\n");
   } else {
      printf("    vmax=%E\n", vmax);
   }
   if (notr) {
      printf("    notr");
   } else {
      printf("    tmin=%E,tmax=%E",tmin,tmax);
   }
   if (notlabel) {
      printf(",notlabel\n");
   } else {
      printf(",tlabel='%s'\n",tlabel);
   }
   xnames[0] = 'x';
   xnames[1] = 'y';
   xnames[2] = 'z';
   for (id=0; id < MAXD; id++) {
      printf("    %c",xnames[id]);
      if (noxr[id]) {
         printf(",norange");
      } else {
         printf(",min=%E,max=%E", xmin[id], xmax[id]);
      }
      printf(",aspect=%f", aspect[id]);
      if (noxlabel[id]) {
         printf(",nolabel\n");
      } else {
         printf(",label='%s'\n",xlabel[id]);
      }
   }
   printf("    label_field=%d", label_field);
   printf(",number_field=%d\n", label_field);
   printf("    line_graph");
   if (noline_connect) {
      printf(",noline_connect");
   } else {
      printf(",line_connect");
   }
   printf(",line_style=%d", line_style);
   printf(",line_palet=%d", line_palet+1);
   if (noline_marker) {
      printf(",noline_marker\n");
   } else {
      printf(",line_marker=%d\n", line_marker);
   }
   if (!nokey) {
      printf("    key,size=%d,label_field=%d,number_field=%d", key_size,
        key_label_field, key_number_field);
      if (!nokey_label) {
         printf(",label=%s", key_label);
      }
      printf("\n");
   }
   if (!nocolour) {
      printf("    colour_map");
      printf(",colours=%d\n",ncolour);
   }
   if (!nocontour) {
      printf("    contour_map");
      printf(",contours=%d",ncontour);
      if (nocontour_label) {
         printf(",with labels\n");
      } else {
         printf(",without labels\n");
      }
      if (!nocontour_levels) {
         printf("      contour levels:-\n");
         for (i=0; i < ncontour; i++) {
            printf("      %3d %E\n", i, contour_levels[i]);
         }
      }
   }
   if (!noband) {
      printf("    bands:\n");
      for (i=0; i < MAXBAND; i++) {
         if (!nobandr[i]) {
            printf("      %d(min=%E,max=%E)\n",i,bandmin[i],bandmax[i]);
         }
      }
   }
   if (!noarrow) {
      printf("    arrow");
      printf(",scale=%f", arrow_scale);
      printf(",rate=%d", arrow_rate);
      printf(",marker=%d", arrow_marker);
      printf(",ratio=%f\n", arrow_ratio);
   }
   if (!nostream) {
      printf("    streamlines");
      printf(",size=%d,x0=[", stream_size);
      for (id=0; id < MAXD; id++) {
         printf(" %E", stream_x[id]);
      }
      printf(" ]\n");
      printf("      dx=[");
      for (id=0; id < MAXD; id++) {
         printf(" %E", stream_dx[id]);
      }
      printf(" ]\n");
   }
/* viewing in 3-dimensions stuff */
   printf("    views,anglev=%f,aspect=%f,lines=%d\n",
     view_anglev, view_aspect, view_lines);
   printf("      azim=%f,altd=%f", view_azim, view_altd);
   if (noview_distance) {
      printf(",nodistance");
   } else {
      printf(",distance=%E", view_distance);
   }
   if (noview_near) {
      printf(",nonear");
   } else {
      printf(",near=%E", view_near);
   }
   if (noview_far) {
      printf(",nofar");
   } else {
      printf(",far=%E", view_far);
   }
   printf("\n");
   if (noview_focus) {
      printf("      nofocus\n");
   } else {
      printf("      focus=[");
      for (id=0; id < MAXD; id++) {
         printf(" %E", view_focus[id]);
      }
      printf(" ]\n");
   }
/* light model */
   printf("    light,fusion=%4.2f,azim=%f,altd=%f",
     light_fusion, light_azim, light_altd);
   if (nolight_distance) {
      printf(",nolight_distance\n");
   } else {
      printf(",distance=%E\n", light_distance);
   }
/* materials */
   printf("    materials,count=%d:-\n",material_count);
   for (i=0; i < material_count; i++) {
      printf("      name=%s,red=%4.2f,green=%4.2f,blue=%4.2f",
        material_name[i], material_red[i], material_green[i],
        material_blue[i]);
      printf(",shine=%4.2f", material_shine[i]);
      if (material_shine[i] > 0.0) {
         printf(",phong=%d", material_phong[i]);
      }
      printf("\n");
   }
/* tags */
   printf("    tags,count=%d:-\n", tag_count);
   for (i=0; i < tag_count; i++) {
      printf("      tag=%s,material=%s\n",tag_name[i],tag_material[i]);
   }
/* 3d scalars */
   printf("    3d_scalar");
   if (data_ortho) {
      printf(",ortho,id=%d", data_ortho_id);
      if (nodata_coord) {
         printf(",nocoord\n");
      } else {
         printf(",coord=%E\n", data_coord);
      }
   } else if (data_shade) {
      printf(",shade\n");
   } else if (data_plane) {
      printf(",plane");
      if (nodata_distance) {
         printf(",nodistance\n");
      } else {
         printf(",distance=%E\n", data_distance);
      }
   } else {
      printf("\n");
   }
/*subwindow stuff */
   if (wc < 0) {
      printf("  No worksurface selected\n");
   } else if (wp->sc < 0) {
      printf("  No subwindow selected\n");
   } else {
      printf("  Worksurface=%s", wp->name);
      printf(",subwindow=%s", sp->name);
      if (sp->clear) {
         printf(",clear");
      } else {
         printf(",unclear");
      }
      if (sp->drawn) {
         printf(",drawn\n");
      } else {
         printf(",undrawn\n");
      }
   }
/* palet */
   printf("  Palet(generic):-\n");
   for (i=0; i < npalet; i++) {
      printf("    %s,red=%5.2f,green=%5.2f,blue=%5.2f\n",palet_name[i],
        palet_r[i], palet_g[i], palet_b[i]);
   }
}

/* WORKSURFACE MANAGEMENT */

void vis_open(char name[], char type[], int *ios)
{
   int i,wd,colour_size,colour_offset;
   bool colour_change;

   if (!*ios) {
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
      } else {
         if (vis_pos_rel) {
            drp_pos_r(vis_pos_rh, vis_pos_rv);
         } else {
            drp_pos_a(vis_pos_ah, vis_pos_av);
         }
         if (vis_size_rel) {
            drp_size_r(vis_size_rh, vis_size_rv);
         } else {
            drp_size_a(vis_size_ah, vis_size_av);
         }
         drp_open(name, type, ios);
         if (*ios) {
            printf("error:  failed to open worksurface (vis_open)\n");
         } else {
/*         admin */
            wc = wd;
            wp = &surface[wc];
            wp->used = TRUE;
            ccopy(name, 0, SZNAM, wp->name);
            ccopy(type, 0, SZNAM, wp->type);
/*         hardware coordinate limits */
            drp_get_hv(&(wp->hmin),&(wp->hmax),&(wp->vmin),&(wp->vmax));
            drp_avail_rub(&(wp->avail_rub));
/*         default characted dimensions */
            drp_get_char(&(wp->char_height), &(wp->char_width), 
               &(wp->char_depth));
/*         colour availability and default configuration */
            drp_avail_colour(&(wp->colour_avail));
            if (wp->colour_avail) {
               drp_change_colour(&colour_change);
               drp_get_colour_shape(&colour_offset, &colour_size);
               wp->cmap_avail = colour_change && (colour_size > 7);
               if (wp->cmap_avail) {
                  drp_put_paint("BLACK", &(wp->colour_black));
                  if (wp->colour_black < 0) {
                     printf("error:  defining BLACK paint\n");
                  } else {
                     drp_set_paint(wp->colour_black, 0.0, 0.0, 0.0);
#ifndef BUG
                     drp_set_colour(wp->colour_black);
#endif
                  }
                  drp_put_paint("SKY", &(wp->colour_sky));
                  if (wp->colour_sky < 0) {
                     printf("error:  defining SKY paint\n");
                     wp->colour_sky = wp->colour_black;
                  } else {
                     drp_set_paint(wp->colour_sky, 0.5, 0.5, 1.0);
                  }
               }
            }
            drp_avail_rgb(&wp->rgb_avail);
            if (wp->rgb_avail) {
               wp->cmap_avail = TRUE;
               drp_set_rgb(0.0, 0.0, 0.0);
            }
/*         setup subwindows */
            for (i=0; i < MAXSUB; i++) {
               sp = &(wp->sub[i]);
               sp->used = FALSE;
            }
            sp = &(wp->sub[0]);
            sp->used = TRUE;
            ccopy("whole",0,SZNAM,sp->name);
            sp->rhmin = 0.0;
            sp->rhmax = 1.0;
            sp->rvmin = 0.0;
            sp->rvmax = 1.0;
            sp->hmin = wp->hmin;
            sp->hmax = wp->hmax;
            sp->vmin = wp->vmin;
            sp->vmax = wp->vmax;
            sp->clear = FALSE;
            sp->drawn = FALSE;
/*
            vis_clear();
*/
            drp_erase();
            drp_get_count(&i);
            if (i == 1) {
               drp_zap_ramps();
            }
            sp->noscene = TRUE;
            scp = &(sp->scene);
         }
      }
   }
}

void vis_select(char name[], int *ios)
{
   int i,wd,sd;

   if (!*ios) {
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
         drp_select(name, ios);
         if (*ios) {
            printf("error:  selecting worksurface\n");
         } else {
            wc = wd;
            wp = &surface[wc];
            sd = wp->sc;
            sp = &(wp->sub[sd]);
            scp = &(sp->scene);
         }
      }
   }
}

void vis_close(name, ios)
   char name[];
   int *ios;
{
   int i,wd;
   struct vis_sub *zp;

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
      for (i=0; i < MAXSUB; i++) {
         if (surface[wd].sub[i].used) {
/*          free any resources for this subwindow */
            zp = &(surface[wd].sub[i]);
            vis_free_sub(zp);
            surface[wd].sub[i].used = FALSE;
         }
      }
      drp_close(name, ios);
      if (*ios) {
         printf("error:  closing worksurface\n");
      }
      surface[wd].used = FALSE;
      if (wd == wc) {
         wc = -1;
      }
   }
}

void vis_clear() {
   struct vis_sub *zp;
   int i;

   if (wc < 0) {
   } else {
      drp_erase();
      for (i=0; i<MAXSUB; i++) {
         if (wp->sub[i].used) {
/*          free any resources for this subwindow */
            zp = &(wp->sub[i]);
            vis_free_sub(zp);
            wp->sub[i].clear = TRUE;
            wp->sub[i].drawn = FALSE;
         }
      }
      drp_get_count(&i);
      if (i > 0) {
/*      clear all associated surfaces */
      }
      drp_zap_ramps();
   }
}

void vis_active(active)
   bool *active;
{
   *active = (wc >= 0);
}

/* SUBWINDOWS */

void vis_open_sub(char name[], float  rhmin, float rhmax,
  float rvmin, float rvmax, int *ios)
{
   float hlen,vlen;
   int i,sd;

   if (!*ios) {
      if (wc < 0) {
         printf("error: no current worksurface (vis_open_sub)\n");
         *ios = ERR;
      } else {
         sd = -1;
         for (i=0; i < MAXSUB; i++) {
            if (*ios) {
            } else if (!surface[wc].sub[i].used) {
               if (sd < 0) {
                  sd = i;
               };
            } else if (csame(name, surface[wc].sub[i].name)) {
               printf("Error:  duplicate subwindow named %s\n", name);
               *ios = ERR;
            }
         }
         if (*ios) {
         } else if (sd < 0) {
            printf("error:  too many subwindows (vis_open_sub)\n");
            *ios = ERR;
         } else {
            wp->sc = sd;
            sp = &(wp->sub[sd]);
            scp = &(sp->scene);
            sp->used = TRUE;
            ccopy(name, 0, SZNAM, sp->name);
            sp->rhmin = rhmin;
            sp->rhmax = rhmax;
            sp->rvmin = rvmin;
            sp->rvmax = rvmax;
            hlen = wp->hmax - wp->hmin;
            sp->hmin = (sp->rhmin)*hlen + wp->hmin;
            sp->hmax = (sp->rhmax)*hlen + wp->hmin;
            vlen = wp->vmax - wp->vmin;
            sp->vmin = (sp->rvmin)*vlen + wp->vmin;
            sp->vmax = (sp->rvmax)*vlen + wp->vmin;
            sp->clear = FALSE;
            sp->drawn = FALSE;
/*          maybe initialise any resources */
         }
      }
   }
}

void vis_select_sub(char name[], int *ios)
{
   int i,sd;

   if (!*ios) {
      if (wc < 0) {
         printf("error:  no current worksurface (vis_select_sub)\n");
         *ios = ERR;
     } else {
         sd = -1;
         for (i=0; i < MAXSUB; i++) {
            if (surface[wc].sub[i].used) {
               if (csame(name, surface[wc].sub[i].name)) {
                  sd = i;
               }
            }
         }
         if (sd < 0) {
            printf("error:  unknown subwindow named %s\n", name);
            *ios = ERR;
         } else {
            wp->sc = sd;
            sp = &(wp->sub[sd]);
            scp = &(sp->scene);
         }
      }
   }
}

void vis_close_sub(char name[], int *ios)
{
   int i,sd;

   if (wc < 0) {
      printf("error:  no current worksurface (drp_select_sub)\n");
      *ios = ERR;
  } else {
      sd = -1;
      for (i=0; i < MAXSUB; i++) {
         if (surface[wc].sub[i].used) {
            if (csame(name, surface[wc].sub[i].name)) {
               sd = i;
            }
         }
      }
      if (sd < 0) {
         printf("error:  unknown subwindow named %s\n", name);
         *ios = ERR;
      } else {
         if (sd > 0) {
            vis_free_sub(&(wp->sub[sd]));
            wp->sub[sd].used = FALSE;
            if (sd == wp->sc) {
               wp->sc = 0;
               sp = &(wp->sub[0]);
               scp = &(sp->scene);
            }
         }
      }
   }
}

void vis_erase() {
   struct vis_sub *zp;
   int i;

   if (wc < 0) {
   } else if (wp->avail_rub && (wp->sc > 0)) {
      vis_free_sub(sp);
      drp_rub(sp->hmin, sp->hmax, sp->vmin, sp->vmax);
      sp->clear = TRUE;
      sp->drawn = FALSE;
/*   set clear any subwindows entirely hidden by this one */
   } else {
/*
      vis_clear();
*/
      drp_erase();
      for (i=0; i<MAXSUB; i++) {
         if (wp->sub[i].used) {
/*          free any resources for this subwindow */
            zp = &(wp->sub[i]);
            vis_free_sub(zp);
            wp->sub[i].clear = TRUE;
            wp->sub[i].drawn = FALSE;
         }
      }
      drp_get_count(&i);
      if (i == 1) {
         drp_zap_ramps();
      }
   }
}

void vis_free_sub(zp)
   struct vis_sub *zp;
{
   scene_struct *zcp;

   if (!zp->noscene) {
      zcp = scp;
      scp = &(zp->scene);
      vis_zap_scene();
      scp = zcp;
   }
}

void vis_get_sub(name)
   char name[];
{
   if (wc < 0) {
      ccopy("?", 0, SZNAM, name);
   } else {
      ccopy((sp->name),0,SZNAM,name);
   }
}

void vis_get_hv(h1,h2,v1,v2)
   float *h1,*h2,*v1,*v2;
{
   float hlen,vlen,h0,v0;

   drp_get_hv(h1, h2, v1, v2);
   if (wp->sc > 0) {
      hlen = *h2 - *h1;
      h0 = *h1;
      *h1 = (sp->rhmin)*hlen + h0;
      *h2 = (sp->rhmax)*hlen + h0;
      vlen = *v2 - *v1;
      v0 = *v1;
      *v1 = (sp->rvmin)*vlen + v0;
      *v2 = (sp->rvmax)*vlen + v0;
   }
}

/* DATA ACCEPTANCES */

void vis_begin()
{
   if (glevel == 0) {
      nocsg = TRUE;
      nogrid = TRUE;
      noscatter = TRUE;
      nopic = TRUE;
   }
   glevel = glevel + 1;
}

void vis_end()
{
   int ios=OK;
   glevel = glevel - 1;
   if (glevel < 0) {
      printf("error:  vis_end called too often\n");
   } else {
      if (glevel == 0) {
         vis_view(&ios);
         vis_free();
         nogrid = TRUE;
         noscatter = TRUE;
         nocsg = TRUE;
         nopic = TRUE;
      }
   }
}


void vis_free() {
   //int n,id,ns,nv;
   int id,ns;	// Removed 'n' and 'nv' declarations BM

/* free any grid model */
   if (!nogrid) {
      ns = 1;
      for (id=0; id < vis_gvp->nd; id++) {
         ns = ns*vis_gvp->nc[id];
         if (vis_gvp->xe[id] != NULL) {
            //cfree(vis_gvp->xe[id], (vis_gvp->nc[id]+1), sizeof(float));
            free(vis_gvp->xe[id]);	// Swapped cfree for standard free BM
            vis_gvp->xe[id] = NULL;
         }
         if (vis_gvp->xa[id] != NULL) {
            //cfree(vis_gvp->xa[id], vis_gvp->nc[id], sizeof(float));
            free(vis_gvp->xa[id]);
            vis_gvp->xa[id] = NULL;
         }
      }
      //nv = ns*vis_gvp->nd;
      if (vis_gvp->sa != NULL) {
         //cfree(vis_gvp->sa, ns, sizeof(int));
         free(vis_gvp->sa);
         vis_gvp->sa = NULL;
      }
      if (vis_gvp->qa != NULL) {
         //cfree(vis_gvp->qa, ns, sizeof(float));
         free(vis_gvp->qa);
         vis_gvp->qa = NULL;
      }
      if (vis_gvp->va != NULL) {
         //cfree(vis_gvp->va, nv, sizeof(float));
         free(vis_gvp->va);
         vis_gvp->va = NULL;
      }
   }
/* free any scatter model */
   if (!noscatter) {
      if (vis_svp->xpt != NULL) {
         //n = vis_svp->npt*vis_svp->nd;
         //cfree(vis_svp->xpt, n, sizeof(float));
         free(vis_svp->xpt);
         vis_svp->xpt = NULL;
      }
      if (vis_svp->spt != NULL) {
         //cfree(vis_svp->spt, vis_svp->npt, sizeof(int));
         free(vis_svp->spt);
         vis_svp->spt = NULL;
      }
      if (vis_svp->qpt != NULL) {
         //cfree(vis_svp->qpt, vis_svp->npt, sizeof(int));
         free(vis_svp->qpt);
         vis_svp->qpt = NULL;
      }
      if (vis_svp->vpt != NULL) {
         //n = vis_svp->npt*vis_svp->nd;
         //cfree(vis_svp->vpt, n, sizeof(float));
         free(vis_svp->vpt);
         vis_svp->vpt = NULL;
      }
   }
}

void vis_csg(cp, ios)
   csg_struct *cp;
   int *ios;
{
   if (*ios) {
      nocsg = TRUE;
   } else {
      vis_cp = cp;
      nocsg = FALSE;
   }
}

void vis_grid(gvp, ios)
   grid_vis_struct *gvp;
   int *ios;
{
   if (*ios) {
      nogrid = TRUE;
   } else {
      vis_gvp = gvp;
      nogrid = FALSE;
   }
}

void vis_scatter(svp, ios)
   scatter_vis_struct *svp;
   int *ios;
{
   if (*ios) {
      noscatter = TRUE;
   } else {
      vis_svp = svp;
      noscatter = FALSE;
   }
}

void vis_set0s(char name[], char value[], int *ios)
{
   int offset,i;

   if (*ios) {
   } else if (csame("font_style",name)) {
      if (wc >= 0) {
         drp_set_font_style(value);
      }
   } else if (csame("caption",name)) {
      if (nocaption) {
         caption_count = 0;
         nocaption = FALSE;
      }
      if (caption_count >= MAXCAP) {
         printf("error:  too many caption lines\n");
         *ios = ERR;
      } else {
         ccopy(value, 0, SZLIN, caption[caption_count]);
         caption_count = caption_count + 1;
      }
   } else if (csame("xlabel",name)) {
      ccopy(value, 0, SZLIN, xlabel[0]);
      noxlabel[0] = FALSE;
   } else if (csame("ylabel",name)) {
      ccopy(value, 0, SZLIN, xlabel[1]);
      noxlabel[1] = FALSE;
   } else if (csame("zlabel",name)) {
      ccopy(value, 0, SZLIN, xlabel[2]);
      noxlabel[2] = FALSE;
   } else if (csame("tlabel",name)) {
      ccopy(value, 0, SZLIN, tlabel);
      notlabel = FALSE;
   } else if (csame("qlabel",name)) {
      ccopy(value, 0, SZLIN, qlabel);
      noqlabel = FALSE;
   } else if (cabrev("tag",name,&offset)) {
      if (vis_tag_find(&(name[offset]), &i)) {
         ccopy(value, 0, SZNAM, tag_material[i]);
      } else if (tag_count >= MAXTAG) {
         printf("error:  too many tags\n");
         *ios = ERR;
      } else {
         ccopy(&(name[offset]), 0, SZNAM, tag_name[tag_count]);
         ccopy(value, 0, SZNAM, tag_material[tag_count]);
         tag_count = tag_count + 1;
      }
   } else if (csame("key_label",name)) {
      ccopy(value, 0, SZNAM, key_label);
      nokey_label = FALSE;
   } else {
      printf("error:  unknown string attribute '%s' (vis_set0s)\n",
        name);
      *ios = ERR;
   }
}

void vis_set0i(char name[], int value, int *ios)
{
   int offset,i;

   if (*ios) {
   } else if (csame("caption_lines",name)) {
      if ((value < 0) || (value > MAXCAP)) {
         printf("error:  caption_lines invalid (vis_set0i)\n");
         *ios = ERR;
      } else {
         caption_lines = value;
      }
   } else if (csame("caption_inset",name)) {
      if (value < 0) {
         printf("error:  caption_inset invalid (vis_set0i)\n");
         *ios = ERR;
      } else {
         caption_inset = value;
      }
   } else if (csame("label_field",name)) {
      if (value < 0) {
         printf("error:  label_field invalid (vis_set0i)\n");
         *ios = ERR;
      } else {
         label_field = value;
      }
   } else if (csame("number_field",name)) {
      if (value < 0) {
         printf("error:  number_field invalid (vis_set0i)\n");
         *ios = ERR;
      } else {
         number_field = value;
      }
   } else if (csame("arrow_marker",name)) {
      arrow_marker = value;
      noarrow_marker = FALSE;
   } else if (csame("line_marker",name)) {
      line_marker = value;
      noline_marker = FALSE;
   } else if (csame("line_style",name)) {
      line_style = value;
      noline_connect = FALSE;
   } else if (csame("line_palet",name)) {
      if ((value < 0) || (value >= npalet)) {
         printf("error:  line palet out of range\n");
         *ios = ERR;
      } else {
         line_palet = value;
      }
   } else if (csame("ncontour",name)) {
      if (value < 1) {
         printf("error:  too few contours\n");
         *ios = ERR;
      } else if (value > MAXCONTOUR) {
         printf("error:  too many contours\n");
         *ios = ERR;
      } else {
         ncontour = value;
         nocontour_levels = TRUE;
      }
   } else if (csame("ncolour",name)) {
      if (value < 1) {
         printf("error:  too few colours\n");
      } else if (value > MAXCM) {
         printf("error:  too many colours\n");
         *ios = ERR;
      } else {
         ncolour = value;
      }
   } else if (csame("nshade",name)) {
      if (value < 1) {
         printf("error:  too few shades\n");
      } else {
         nshade = value;
      }
   } else if (csame("arrow_rate",name)) {
      if (value < 1) {
         printf("error:  invalid arrow_rate\n");
         *ios = ERR;
      } else {
         arrow_rate = value;
      }
   } else if (csame("stream_size",name)) {
      stream_size = value;
   } else if (csame("view_lines",name)) {
      if (value < 2) {
         printf("error:  view_lines too small\n");
         *ios = ERR;
      } else {
         view_lines = value;
      }
   } else if (cabrev("material_phong",name,&offset)) {
      if ((value < 1) || (value > 32)) {
         printf("error:  material phong out of range\n");
         *ios = ERR;
      } else if (!vis_material_find(&(name[offset]), &i)) {
         if (material_count >= MAXMAT) {
            printf("error:  too many materials\n");
            *ios = ERR;
         } else {
            ccopy(&(name[offset]),0,SZNAM,material_name[material_count]);
            material_red[material_count] = 0.5;
            material_green[material_count] = 0.5;
            material_blue[material_count] = 0.5;
            material_shine[material_count] = 0.0;
            material_phong[material_count] = value;
            material_count = material_count + 1;
         }
      } else {
         material_phong[i] = value;
      }
   } else if (csame("key_size",name)) {
      if ((value < 0) || (value > MAXKEY)) {
         printf("error:  illegal key size\n");
         *ios = ERR;
      } else {
         key_size = value;
      }
   } else if (csame("key_label_field",name)) {
      if ((value < 0) || (value > SZNAM)) {
         printf("error:  illegal key label field\n");
         *ios = ERR;
      } else {
         key_label_field = value;
      }
   } else if (csame("key_number_field",name)) {
      if ((value < 0) || (value > SZNAM)) {
         printf("error:  illegal key number field\n");
         *ios = ERR;
      } else {
         key_number_field = value;
      }
   } else if (csame("data_ortho_id",name)) {
      if ((value < 0) || (value >= MAXD)) {
         printf("error:  data_ortho_id ignored\n");
      } else {
         data_ortho_id = value;
      }
   } else {
      printf("error:  unknown integer attribute '%s' (vis_set0i)\n",
        name);
   }
}

void vis_set0r(char name[], float value, int *ios)
{
   int offset,i,iband;

   if (*ios) {
   } else if (csame("font_size",name)) {
      if (wc >= 0) {
         drp_set_font_size(value);
      }
   } else if (csame("surface_height",name)) {
      surface_height = value;
   } else if (csame("pos_rh",name)) {
      vis_pos_rh = value;
      vis_pos_rel = TRUE;
   } else if (csame("pos_rv",name)) {
      vis_pos_rv = value;
      vis_pos_rel = TRUE;
   } else if (csame("pos_ah",name)) {
      vis_pos_ah = value;
      vis_pos_rel = FALSE;
   } else if (csame("pos_av",name)) {
      vis_pos_av = value;
      vis_pos_rel = FALSE;
   } else if (csame("size_rh",name)) {
      vis_size_rh = value;
      vis_size_rel = TRUE;
   } else if (csame("size_rv",name)) {
      vis_size_rv = value;
      vis_size_rel = TRUE;
   } else if (csame("size_ah",name)) {
      vis_size_ah = value;
      vis_size_rel = FALSE;
   } else if (csame("size_av",name)) {
      vis_size_av = value;
      vis_size_rel = FALSE;
   } else if (csame("qmin",name)) {
      noqmin = FALSE;
      qmin = value;
      noqr = noqmax;
   } else if (csame("qmax",name)) {
      noqmax = FALSE;
      qmax = value;
      noqr = noqmin;
   } else if (csame("vmax",name)) {
      novmax = FALSE;
      vmax = value;
   } else if (csame("tmin",name)) {
      notmin = FALSE;
      tmin = value;
      notr = notmax;
   } else if (csame("tmax",name)) {
      notmax = FALSE;
      tmax = value;
      notr = notmin;
   } else if (csame("xaspect",name)) {
      aspect[0] = value;
   } else if (csame("yaspect",name)) {
      aspect[1] = value;
   } else if (csame("zaspect",name)) {
      aspect[2] = value;
   } else if (csame("xmin",name)) {
      noxmin[0] = FALSE;
      xmin[0] = value;
      noxr[0] = noxmax[0];
   } else if (csame("xmax",name)) {
      noxmax[0] = FALSE;
      xmax[0] = value;
      noxr[0] = noxmin[0];
   } else if (csame("ymin",name)) {
      noxmin[1] = FALSE;
      xmin[1] = value;
      noxr[1] = noxmax[1];
   } else if (csame("ymax",name)) {
      noxmax[1] = FALSE;
      xmax[1] = value;
      noxr[1] = noxmin[1];
   } else if (csame("zmin",name)) {
      noxmin[2] = FALSE;
      xmin[2] = value;
      noxr[2] = noxmax[2];
   } else if (csame("zmax",name)) {
      noxmax[2] = FALSE;
      xmax[2] = value;
      noxr[2] = noxmin[2];
   } else if (csame("vmax",name)) {
      vmax = value;
   } else if (csame("hborder",name)) {
      hborder = value;
   } else if (csame("vborder",name)) {
      vborder = value;
   } else if (csame("tick_scale",name)) {
      tick_scale = value;
   } else if (csame("marker_scale",name)) {
      drp_set_marker_scale(value);
   } else if (csame("arrow_scale",name)) {
      arrow_scale = value;
   } else if (csame("arrow_ratio",name)) {
      arrow_ratio = value;
   } else if (csame("contour_level",name)) {
      nocontour_levels = FALSE;
      ncontour = 1;
      contour_levels[0] = value;
   } else if (csame("xstream",name)) {
      stream_x[0] = value;
   } else if (csame("ystream",name)) {
      stream_x[1] = value;
   } else if (csame("zstream",name)) {
      stream_x[2] = value;
   } else if (csame("dxstream",name)) {
      stream_dx[0] = value;
   } else if (csame("dystream",name)) {
      stream_dx[1] = value;
   } else if (csame("dzstream",name)) {
      stream_dx[2] = value;
   } else if (csame("view_distance",name)) {
      if (value <= 0.0) {
         printf("error:  invalid view_distance (rejected)\n");
         *ios = ERR;
      } else {
         view_distance = value;
         noview_distance = FALSE;
      }
   } else if (csame("view_near",name)) {
      if (value <= 0.0) {
         printf("error:  invalid view_near\n");
         *ios = ERR;
      } else {
         view_near = value;
         noview_near = FALSE;
      }
   } else if (csame("view_far",name)) {
      if (value <= 0.0) {
         printf("error:  invalid view_far\n");
         *ios = ERR;
      } else {
         view_far = value;
         noview_far = FALSE;
      }
   } else if (csame("view_azim",name)) {
      if ((value < 0.0) || (value >= 360.0)) {
         printf("error:  invalid view_azim\n");
         *ios = ERR;
      } else {
         view_azim = value;
      }
   } else if (csame("view_altd",name)) {
      if ((value < -90.0) || (value > 90.0)) {
         printf("error:  invalid view_altd\n");
         *ios = ERR;
      } else {
         view_altd = value;
      }
   } else if (csame("view_aspect",name)) {
      if ((value < 0.3) || (value > 3.0)) {
         printf("error:  invalid view_aspect\n");
         *ios = ERR;
      } else {
         view_aspect = value;
      }
   } else if (csame("view_anglev",name)) {
      if ((value < 5.0) || (value > 25.0)) {
         printf("error:  invalid view_anglev\n");
         *ios = ERR;
      } else {
         view_anglev = value;
      }
   } else if (csame("light_azim",name)) {
      if ((value < 0.0) || (value >= 360.0)) {
         printf("error:  invalid light_azim\n");
         *ios = ERR;
      } else {
         light_azim = value;
      }
   } else if (csame("light_altd",name)) {
      if ((value < -90.0) || (value > 90.0)) {
         printf("error:  invalid light_altd\n");
         *ios = ERR;
      } else {
         light_altd = value;
      }
   } else if (csame("light_distance",name)) {
      if (value <= 0.0) {
         printf("error:  invalid light_distance\n");
         *ios = ERR;
      } else {
         light_distance = value;
         nolight_distance = FALSE;
      }
   } else if (csame("light_fusion",name)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  invalid light_fusion\n");
         *ios = ERR;
      } else {
         light_fusion = value;
      }
   } else if (cabrev("material_red",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  material colour out of range\n");
         *ios = ERR;
      } else if (!vis_material_find(&(name[offset]), &i)) {
         if (material_count >= MAXMAT) {
            printf("error:  too many materials\n");
            *ios = ERR;
         } else {
            ccopy(&(name[offset]),0,SZNAM,material_name[material_count]);
            material_red[material_count] = value;
            material_green[material_count] = 0.5;
            material_blue[material_count] = 0.5;
            material_shine[material_count] = 0.0;
            material_phong[material_count] = 2;
            material_count = material_count + 1;
         }
      } else {
         material_red[i] = value;
      }
   } else if (cabrev("material_green",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  material colour out of range\n");
         *ios = ERR;
      } else if (!vis_material_find(&(name[offset]), &i)) {
         if (material_count >= MAXMAT) {
            printf("error:  too many materials\n");
            *ios = ERR;
         } else {
            ccopy(&(name[offset]),0,SZNAM,material_name[material_count]);
            material_red[material_count] = 0.5;
            material_green[material_count] = value;
            material_blue[material_count] = 0.5;
            material_shine[material_count] = 0.0;
            material_phong[material_count] = 2;
            material_count = material_count + 1;
         }
      } else {
         material_green[i] = value;
      }
   } else if (cabrev("material_blue",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  material colour out of range\n");
         *ios = ERR;
      } else if (!vis_material_find(&(name[offset]), &i)) {
         if (material_count >= MAXMAT) {
            printf("error:  too many materials\n");
            *ios = ERR;
         } else {
            ccopy(&(name[offset]),0,SZNAM,material_name[material_count]);
            material_red[material_count] = 0.5;
            material_green[material_count] = 0.5;
            material_blue[material_count] = value;
            material_shine[material_count] = 0.0;
            material_phong[material_count] = 2;
            material_count = material_count + 1;
         }
      } else {
         material_blue[i] = value;
      }
   } else if (cabrev("material_shine",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  material shine out of range\n");
         *ios = ERR;
      } else if (!vis_material_find(&(name[offset]), &i)) {
         if (material_count >= MAXMAT) {
            printf("error:  too many materials\n");
            *ios = ERR;
         } else {
            ccopy(&(name[offset]),0,SZNAM,material_name[material_count]);
            material_red[material_count] = 0.5;
            material_green[material_count] = 0.5;
            material_blue[material_count] = 0.5;
            material_shine[material_count] = value;
            material_phong[material_count] = 2;
            material_count = material_count + 1;
         }
      } else {
         material_shine[i] = value;
      }
   } else if (csame("sky_red",name)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  invalid sky red\n");
         *ios = ERR;
      } else {
         sky_red = value;
      }
   } else if (csame("sky_green",name)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  invalid sky green\n");
         *ios = ERR;
      } else {
         sky_green = value;
      }
   } else if (csame("sky_blue",name)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  invalid sky blue\n");
         *ios = ERR;
      } else {
         sky_blue = value;
      }
   } else if (cabrev("band_min",name,&offset)) {
      if (csame("0",&(name[offset]))) {
         iband = 0;
      } else if (csame("1",&(name[offset]))) {
         iband = 1;
      } else if (csame("2",&(name[offset]))) {
         iband = 2;
      } else {
         printf("error:  illegal band index\n");
         *ios = ERR;
      }
      if (!*ios) {
         bandmin[iband] = value;
         nobandmin[iband] = FALSE;
         nobandr[iband] = nobandmin[iband] || nobandmax[iband];
      }
   } else if (cabrev("band_max",name,&offset)) {
      if (csame("0",&(name[offset]))) {
         iband = 0;
      } else if (csame("1",&(name[offset]))) {
         iband = 1;
      } else if (csame("2",&(name[offset]))) {
         iband = 2;
      } else {
         printf("error:  illegal band index\n");
         *ios = ERR;
      }
      if (!*ios) {
         bandmax[iband] = value;
         nobandmax[iband] = FALSE;
         nobandr[iband] = nobandmin[iband] || nobandmax[iband];
      }
   } else if (csame("data_distance",name)) {
      if (value <= 0.0) {
         printf("error:  data_distance must be positive\n");
      } else {
         data_distance = value;
         nodata_distance = FALSE;
      }
   } else if (csame("data_coord",name)) {
      data_coord = value;
      nodata_coord = FALSE;
   } else if (cabrev("palet_red",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  palet_red colour out of range\n");
         *ios = ERR;
      } else if (!vis_palet_find(&(name[offset]), &i)) {
         printf("error:  palet=%s not known\n", &(name[offset]));
      } else {
         palet_r[i] = value;
      }
   } else if (cabrev("palet_green",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  palet_green colour out of range\n");
         *ios = ERR;
      } else if (!vis_palet_find(&(name[offset]), &i)) {
         printf("error:  palet=%s not known\n", &(name[offset]));
      } else {
         palet_g[i] = value;
      }
   } else if (cabrev("palet_blue",name,&offset)) {
      if ((value < 0.0) || (value > 1.0)) {
         printf("error:  palet_blue colour out of range\n");
         *ios = ERR;
      } else if (!vis_palet_find(&(name[offset]), &i)) {
         printf("error:  palet=%s not known\n", &(name[offset]));
      } else {
         palet_b[i] = value;
      }
   } else {
      printf("error:  unknown real attribute '%s' (vis_set0r)\n",
        name);
   }
}

void vis_set0l(name, value, ios)
   char name[];
   int *ios;
   bool value;
{

   if (*ios) {
   } else {
      printf("error:  unknown logical attribute '%s' (vis_set0l)\n",
        name);
   }
}

void vis_seta(char name[], int *ios)
{
   int offset,iband;

   if (*ios) {
   } else if (csame("font",name)) {
      if (wc >= 0) {
         drp_set_font();
         drp_get_char(&(wp->char_height), &(wp->char_width), 
            &(wp->char_depth));
      }
#ifndef NOGREY
   } else if (csame("nogrey_rain",name)) {
      vis_grey_rain = FALSE;
   } else if (csame("grey_rain",name)) {
      vis_grey_rain = TRUE;
#endif
   } else if (csame("nobox",name)) {
      nobox = TRUE;
   } else if (csame("box",name)) {
      nobox = FALSE;
   } else if (csame("noqr",name)) {
      noqr = TRUE;
   } else if (csame("novmax",name)) {
      novmax = TRUE;
      noqmin = TRUE;
      noqmax = TRUE;
   } else if (csame("notr",name)) {
      notr = TRUE;
      notmin = TRUE;
      notmax = TRUE;
   } else if (csame("noxr",name)) {
      noxr[0] = TRUE;
      noxmin[0] = TRUE;
      noxmax[0] = TRUE;
   } else if (csame("noyr",name)) {
      noxr[1] = TRUE;
      noxmin[1] = TRUE;
      noxmax[1] = TRUE;
   } else if (csame("nozr",name)) {
      noxr[2] = TRUE;
      noxmin[2] = TRUE;
      noxmax[2] = TRUE;
   } else if (csame("noscalp",name)) {
      noscalp = TRUE;
   } else if (csame("scalp",name)) {
      noscalp = FALSE;
   } else if (csame("noband",name)) {
      noband = TRUE;
   } else if (csame("band",name)) {
      noband = FALSE;
      nocolour = TRUE;
      nocontour = TRUE;
   } else if (csame("nocolour",name)) {
      nocolour = TRUE;
   } else if (csame("colour",name)) {
      nocolour = FALSE;
      noband = TRUE;
   } else if (csame("nocontour",name)) {
      nocontour = TRUE;
   } else if (csame("contour",name)) {
      nocontour = FALSE;
      noband = TRUE;
   } else if (csame("nocontour_label",name)) {
      nocontour_label = TRUE;
   } else if (csame("contour_label",name)) {
      nocontour_label = FALSE;
   } else if (csame("noarrow",name)) {
      noarrow = TRUE;
   } else if (csame("arrow",name)) {
      noarrow = FALSE;
      nostream = TRUE;
   } else if (csame("nostream",name)) {
      nostream = TRUE;
   } else if (csame("stream",name)) {
      nostream = FALSE;
      noarrow = TRUE;
   } else if (csame("noxlabel",name)) {
     noxlabel[0] = TRUE;
   } else if (csame("noylabel",name)) {
     noxlabel[1] = TRUE;
   } else if (csame("nozlabel",name)) {
     noxlabel[2] = TRUE;
   } else if (csame("notlabel",name)) {
     notlabel = TRUE;
   } else if (csame("noqlabel",name)) {
     noqlabel = TRUE;
   } else if (csame("caption_below",name)) {
      caption_below = TRUE;
   } else if (csame("caption_above",name)) {
      caption_below = FALSE;
   } else if (csame("nocaption",name)) {
      nocaption = TRUE;
   } else if (csame("noarrow_marker",name)) {
      noarrow_marker = TRUE;
   } else if (csame("noline_marker",name)) {
      noline_marker = TRUE;
   } else if (csame("noline_connect",name)) {
      noline_connect = TRUE;
   } else if (csame("scales",name)) {
      noscales = FALSE;
   } else if (csame("noscales",name)) {
      noscales = TRUE;
   } else if (csame("key",name)) {
      nokey = FALSE;
   } else if (csame("nokey",name)) {
      nokey = TRUE;
   } else if (csame("nokey_label",name)) {
      nokey_label = TRUE;
   } else if (csame("noview_focus",name)) {
      noview_focus = TRUE;
   } else if (csame("noview_distance",name)) {
      noview_distance = TRUE;
   } else if (csame("noview_near",name)) {
      noview_near = TRUE;
   } else if (csame("noview_far",name)) {
      noview_far = TRUE;
   } else if (csame("nolight_distance",name)) {
      nolight_distance = TRUE;
   } else if (csame("material_reset",name)) {
      material_count = 1;
      material_red[0] = 1.0;
      material_green[0] = 1.0;
      material_blue[0] = 1.0;
      material_shine[0] = 0.0;
      material_phong[0] = 2;
   } else if (csame("tag_reset",name)) {
      tag_count = 0;
   } else if (cabrev("nobandr",name,&offset)) {
      if (csame("0",&(name[offset]))) {
         iband = 0;
      } else if (csame("1",&(name[offset]))) {
         iband = 1;
      } else if (csame("2",&(name[offset]))) {
         iband = 2;
      } else {
         printf("error:  illegal band index\n");
         *ios = ERR;
      }
      if (!*ios) {
         nobandmin[iband] = TRUE;
         nobandmax[iband] = TRUE;
         nobandr[iband] = TRUE;
      }
   } else if (csame("log",name)) {
      dolog = TRUE;
   } else if (csame("nolog",name)) {
      dolog = FALSE;
   } else if (csame("data_shade",name)) {
      data_shade = TRUE;
      data_plane = FALSE;
      data_ortho = FALSE;
   } else if (csame("data_plane",name)) {
      data_plane = TRUE;
      data_shade = FALSE;
      data_ortho = FALSE;
   } else if (csame("data_ortho",name)) {
      data_ortho = TRUE;
      data_shade = FALSE;
      data_plane = FALSE;
   } else if (csame("nodata_distance",name)) {
      nodata_distance = TRUE;
   } else if (csame("nodata_coord",name)) {
      nodata_coord = TRUE;
   } else if (csame("noview_near",name)) {
      noview_near = TRUE;
   } else if (csame("noview_far",name)) {
      noview_far = TRUE;
   } else if (csame("wire_frame",name)) {
      wire_frame = TRUE;
   } else if (csame("nowire_frame",name)) {
      wire_frame = FALSE;
   } else if (csame("data_frame",name)) {
      data_frame = TRUE;
   } else if (csame("nodata_frame",name)) {
      data_frame = FALSE;
   } else {
      printf("error:  unknown assertion '%s' (vis_seta)\n",
        name);
   }
}

void vis_set1r(char name[], int nvalue, float values[],
  int *ios)
{
   int i,id;

   if (*ios) {
   } else if (csame("contour_levels",name)) {
      if (nvalue < 1) {
         printf("error:  too few contour levels\n");
         *ios = ERR;
      } else if (nvalue > MAXCONTOUR) {
         printf("error:  too many contour levels\n");
         *ios = ERR;
      } else {
         nocontour_levels = FALSE;
         ncontour = nvalue;
         for (i=0; i < nvalue; i++) {
            contour_levels[i] = values[i];
         }
      }
   } else if (csame("view_focus",name)) {
      if (nvalue != 3) {
         printf("error:  too few view_focus coordinates\n");
         *ios = ERR;
      } else {
         noview_focus = FALSE;
         for (id=0; id < MAXD; id++) {
            view_focus[id] = values[id];
         }
      }
   } else {
      printf("error:  unknown real attribute '%s' (vis_set1r)\n",
        name);
   }
}

void vis_get0s(name, lenvalue, value, ios)
   char name[];
   int lenvalue,*ios;
   char value[];
{

   if (*ios) {
   } else if (wc < 0) {
      printf("error:  attribute '%s' unknown or undefined\n", name);
   } else if (csame("subwindow/name",name)) {
      ccopy((sp->name), 0, lenvalue, value);
   } else {
      printf("error:  unknown string attribute '%s' (vis_get0s)\n",
        name);
   }
}

void vis_get0i(name, value, ios)
   char name[];
   int *ios;
   int *value;
{

   if (*ios) {
   } else {
      printf("error:  unknown integer attribute '%s' (vis_get0i)\n",
        name);
   }
}

void vis_get0r(name, value, ios)
   char name[];
   int *ios;
   float *value;
{

   if (*ios) {
   } else if (wc < 0) {
      printf("error:  attribute '%s' unknown or undefined\n", name);
   } else {
      printf("error:  unknown real attribute '%s' (vis_get0r)\n",
        name);
   }
}

void vis_get0l(name, value, ios)
   char name[];
   int *ios;
   bool *value;
{

   if (*ios) {
   } else {
      printf("error:  unknown string attribute '%s' (vis_get0l)\n",
        name);
   }
}

/* VISUALISATION ROUTINES (INTERNAL) */

void vis_view(ios) 
   int *ios;
{
   bool draw_axes;

   if (wc < 0) {
      printf("error:  no active worksurface (vis_view)\n");
   } else {
/*   decide on whether we NEED to design and draw axes */
      if (!sp->drawn) {
         draw_axes = TRUE;
      } else if (!nobox) {
         draw_axes = TRUE;
      } else {
         draw_axes = FALSE;
      }
/*   if so then clear surface and design and draw axes */
      if (draw_axes) {
/*       zap any scene associated with the subwindow */
         if (!(sp->noscene)) {
            vis_zap_scene();
         }
/*       clear the subwindow surface */
         if (!(sp->clear)) {
            vis_erase();
         }
/*       design and draw axes (tell dr3) */
         vis_axes(ios);
/*    otherwise capture scene (tell dr3) */
      } else if (sp->nd == 3) {
         printf("vis_view:  re-use axes\n");
         dr3_put_scene(scp);
         dr3_push_image(scp);
      }
/*   attempt to draw each available dataset on current axes */
      if (!*ios) {
         vis_draw(ios);
         if (sp->nd == 3) {
            dr3_fresh_scene();
         }
      }
/*    Get what dr3 has got */
      if (sp->nd == 3) {
         dr3_get_scene(scp);
         dr3_pop_image(scp);
         if (*ios) {
            vis_zap_scene();
         } else {
             sp->noscene = FALSE;
          }
      }
#ifndef GFLUSH_BUG
      drp_gflush();
#endif
   }
#ifdef GFLUSH_BUG
   drp_gflush();
#endif
}

#ifdef CURV
#define MWC 40
float float_min(a,b)
   float a,b;
{
   if (a > b) {
      return(b);
   } else {
      return(a);
   }
}

float float_max(a,b)
   float a,b;
{
   if (a > b) {
      return(a);
   } else {
      return(b);
   }
}
#endif

void vis_axes(ios) 
   int *ios;
{
   float gclen,wclen;
   int id,i;
#ifdef CURV
   int j,k,jd,kad[MAXD],ndc,ndom;
   float dc0,dc1,dc2,u[MAXD],w[MAXD],xo2,zm[MAXD];
   xo2 = 0.0;
   for (i=0; i < MAXD; i++)
     kad[i] = 0;
#endif

/* check error */
   if (*ios) {
/* design axes based on grid dataset? */
   } else if (!nogrid) {
/*   caption */
      if (nocaption) {
         sp->caption_count = 1;
         ccopy(vis_gvp->title, 0, SZLIN, sp->caption[0]);
      } else {
         sp->caption_count = caption_count;
         for (i=0; i < caption_count; i++) {
            ccopy(caption[i], 0, SZLIN, sp->caption[i]);
         }
      }
#ifdef CURV
/*   First the remaining space coordinates */
      sp->nd = vis_gvp->nd;
      sp->ndc = vis_gvp->ndc;
#ifdef DEBUG
      printf("vis_axes:  gvp->nd=%d,gvp->ndc=%d\n", vis_gvp->nd,
        vis_gvp->ndc);
#endif
      for (id=0; id < vis_gvp->ndc; id++) {
         ccopy(vis_gvp->coords[id], 0, SZNAM, sp->gname[id]);
         sp->gspace[id] = TRUE;
         sp->gasp[id] = aspect[id];
         if (noxlabel[id]) {
            ccopy(vis_gvp->clabel[id], 0, SZLIN, sp->glabel[id]);
         } else {
            ccopy(xlabel[id], 0, SZLIN, sp->glabel[id]);
         }
         if (noxr[id]) {
            sp->gcmin[id] = vis_gvp->xmin[id];
            sp->gcmax[id] = vis_gvp->xmax[id];
         } else {
            sp->gcmin[id] = xmin[id];
            sp->gcmax[id] = xmax[id];
         }
#ifdef DEBUG
         printf("vis_axes:  id=%d,gcmin=%f,gcmax=%f\n", id,
           sp->gcmin[id], sp->gcmax[id]);
#endif
      }
/*   Add a dummy space coordinate to escalate curv 2D to 3D */
      sp->curv = vis_gvp->curv && (vis_gvp->ndc >= 2);
#ifdef DEBUG
      printf("vis_axes:  gvp->curv=%d,sp->curv=%d\n", vis_gvp->curv,
        sp->curv);
#endif
      if (vis_gvp->curv && (vis_gvp->ndc==2) &&
        ((vis_gvp->nd + vis_gvp->ndo))==4) {
#ifdef DEBUG
         printf("vis_axes:  add dummy space coordinate\n");
#endif
/*       capture kad and the masked space coordinate value */
         xo2 = vis_gvp->xo[0];
/*       original model space dimensions */
         ndom = (vis_gvp->nd+vis_gvp->ndo-1);
/*       masked model space dimensions */
         ndc = vis_gvp->ndc;
         for (jd=0; jd < MAXD; jd++) {
            kad[jd] = vis_gvp->kad[jd];
#ifdef DEBUG
            printf("vis_axes:  kad[%d]=%d\n", jd, kad[jd]);
#endif
         }
         id = sp->ndc;
         ccopy("?", 0, SZNAM, sp->gname[id]);
         sp->gspace[id] = TRUE;
         sp->gasp[id] = aspect[id];
         ccopy("?", 0, SZNAM, sp->glabel[id]);
         sp->gcmin[id] = 0.0;
         sp->gcmax[id] = 1.0;
         sp->ndc++;
         sp->nd++;
      }
/*   Now add the optional time coordinate */
      jd = sp->ndc;
      for (id=vis_gvp->ndc; id < vis_gvp->nd; id++) {
#ifdef DEBUG
         printf("vis_axes:  adding time coordinate,id=%d\n", id);
#endif
         ccopy(vis_gvp->coords[id], 0, SZNAM, sp->gname[jd]);
         sp->gspace[jd] = FALSE;
         if (notlabel) {
            ccopy(sp->gname[jd], 0, SZLIN, sp->glabel[jd]);
         } else {
            ccopy(tlabel, 0, SZLIN, sp->glabel[jd]);
         }
         if (notr) {
            sp->gcmin[jd] = vis_gvp->xmin[id];
            sp->gcmax[jd] = vis_gvp->xmax[id];
         } else {
            sp->gcmin[jd] = tmin;
            sp->gcmax[jd] = tmax;
         }
         jd++;
      }
#else
/*   basic coordinates (space,time) */
      sp->nd = vis_gvp->nd;
      sp->ndc = vis_gvp->ndc;
      for (id=0; id < sp->nd; id++) {
         ccopy(vis_gvp->coords[id], 0, SZNAM, sp->gname[id]);
         if (id < vis_gvp->ndc) {
            sp->gspace[id] = TRUE;
            sp->gasp[id] = aspect[id];
            if (noxlabel[id]) {
               ccopy(vis_gvp->clabel[id], 0, SZLIN, sp->glabel[id]);
            } else {
               ccopy(xlabel[id], 0, SZLIN, sp->glabel[id]);
            }
            if (noxr[id]) {
               sp->gcmin[id] = vis_gvp->xmin[id];
               sp->gcmax[id] = vis_gvp->xmax[id];
            } else {
               sp->gcmin[id] = xmin[id];
               sp->gcmax[id] = xmax[id];
            }
         } else {
            sp->gspace[id] = FALSE;
            if (notlabel) {
               ccopy(sp->gname[id], 0, SZLIN, sp->glabel[id]);
            } else {
               ccopy(tlabel, 0, SZLIN, sp->glabel[id]);
            }
            if (notr) {
               sp->gcmin[id] = vis_gvp->xmin[id];
               sp->gcmax[id] = vis_gvp->xmax[id];
            } else {
               sp->gcmin[id] = tmin;
               sp->gcmax[id] = tmax;
            }
         }
      }
#endif
/*   promote scalar to coordinate? */
      sp->gscalp = FALSE;
      if (vis_gvp->noscalar) {
      } else if (sp->nd > 2) {
      } else if ((sp->nd == 2) && noscalp) {
      } else {
         sp->gscalp = TRUE;
         sp->nd = sp->nd + 1;
         id = sp->nd-1;
         ccopy(vis_gvp->qlabel, 0, SZNAM, sp->gname[id]);
         sp->gspace[id] = FALSE;
         if (noqlabel) {
            ccopy(sp->gname[id], 0, SZLIN, sp->glabel[id]);
         } else {
            ccopy(qlabel, 0, SZLIN, sp->glabel[id]);
         }
         if (noqr) {
            sp->gcmin[id] = vis_gvp->qmin;
            sp->gcmax[id] = vis_gvp->qmax;
         } else {
            sp->gcmin[id] = qmin;
            sp->gcmax[id] = qmax;
         }
      }
/* otherwise, based on scatter dataset? */
   } else if (!noscatter) {
/*   caption */
      if (nocaption) {
         sp->caption_count = 1;
         ccopy(vis_svp->title, 0, SZLIN, sp->caption[0]);
      } else {
         sp->caption_count = caption_count;
         for (i=0; i < caption_count; i++) {
            ccopy(caption[i], 0, SZLIN, sp->caption[i]);
         }
      }
#ifdef CURV
/*   First the remaining space coordinates */
      sp->nd = vis_svp->nd;
      sp->ndc = vis_svp->ndc;
      for (id=0; id < vis_svp->ndc; id++) {
         ccopy(vis_svp->coords[id], 0, SZNAM, sp->gname[id]);
         sp->gspace[id] = TRUE;
         sp->gasp[id] = aspect[id];
         if (noxlabel[id]) {
            ccopy(vis_svp->clabel[id], 0, SZLIN, sp->glabel[id]);
         } else {
            ccopy(xlabel[id], 0, SZLIN, sp->glabel[id]);
         }
         if (noxr[id]) {
            sp->gcmin[id] = vis_svp->xmin[id];
            sp->gcmax[id] = vis_svp->xmax[id];
         } else {
            sp->gcmin[id] = xmin[id];
            sp->gcmax[id] = xmax[id];
         }
      }
/*   Add a dummy space coordinate to escalate curv 2D to 3D */
      sp->curv = vis_svp->curv && (vis_svp->ndc >= 2);
      if (vis_svp->curv && (vis_svp->ndc==2) &&
        ((vis_svp->nd + vis_svp->ndo))==4) {
/*       capture kad and the masked space coordinate value */
         xo2 = vis_svp->xo[0];
/*       original model space dimensions */
         ndom = (vis_svp->nd+vis_svp->ndo-1);
/*       masked model space dimensions */
         ndc = vis_svp->ndc;
         for (jd=0; jd < MAXD; jd++) {
            kad[jd] = vis_svp->kad[jd];
         }
         id = sp->ndc;
         ccopy("?", 0, SZNAM, sp->gname[id]);
         sp->gspace[id] = TRUE;
         sp->gasp[id] = aspect[id];
         ccopy("?", 0, SZNAM, sp->glabel[id]);
         sp->gcmin[id] = 0.0;
         sp->gcmax[id] = 1.0;
         sp->ndc++;
         sp->nd++;
      }
/*   No time coordinate for scatter grams */
#else
/*   optional space coordinates */
      sp->ndc = vis_svp->ndc;
      sp->nd = vis_svp->nd;
      if (vis_svp->nd > 0) {
         for (id=0; id < sp->ndc; id++) {
            ccopy(vis_svp->coords[id], 0, SZNAM, sp->gname[id]);
            sp->gspace[id] = TRUE;
            sp->gasp[id] = aspect[id];
            if (noxlabel[id]) {
               ccopy(vis_svp->clabel[id], 0, SZLIN, sp->glabel[id]);
            } else {
               ccopy(xlabel[id], 0, SZLIN, sp->glabel[id]);
            }
            if (noxr[id]) {
               sp->gcmin[id] = vis_svp->xmin[id];
               sp->gcmax[id] = vis_svp->xmax[id];
            } else {
               sp->gcmin[id] = xmin[id];
               sp->gcmax[id] = xmax[id];
            }
         }
      }
#endif
/*   promote scalar to coordinate? */
      sp->gscalp = FALSE;
      if (vis_svp->noscalar) {
      } else if (sp->nd > 2) {
      } else if ((sp->nd == 2) && noscalp) {
      } else {
         sp->gscalp = TRUE;
         sp->nd = sp->nd + 1;
         id = sp->nd - 1;
         ccopy(vis_svp->qlabel, 0, SZNAM, sp->gname[id]);
         sp->gspace[id] = FALSE;
         if (noqlabel) {
            ccopy(sp->gname[id], 0, SZNAM, sp->glabel[id]);
         } else {
            ccopy(qlabel, 0, SZLIN, sp->glabel[id]);
         }
         if (noqr) {
            sp->gcmin[id] = vis_svp->qmin;
            sp->gcmax[id] = vis_svp->qmax;
         } else {
            sp->gcmin[id] = qmin;
            sp->gcmax[id] = qmax;
         }
      }
/* otherwise, based on CSG? */
   } else if (!nocsg) {
/*   caption */
      if (nocaption) {
         sp->caption_count = 1;
         ccopy(vis_cp->title, 0, SZLIN, sp->caption[0]);
      } else {
         sp->caption_count = caption_count;
         for (i=0; i < caption_count; i++) {
            ccopy(caption[i], 0, SZLIN, sp->caption[i]);
         }
      }
#ifdef CURV
      sp->curv = FALSE;
#endif
/*   fixed 3-dimensional cartesians */
      sp->ndc = 3;
      sp->nd = sp->ndc;
      ccopy("x", 0, SZNAM, sp->gname[0]);
      ccopy("y", 0, SZNAM, sp->gname[1]);
      ccopy("z", 0, SZNAM, sp->gname[2]);
      for (id=0; id < MAXD; id++) {
         if (noxlabel[id]) {
            ccopy(sp->gname[id], 0, SZLIN, sp->glabel[id]);
         } else {
            ccopy(xlabel[id], 0, SZLIN, sp->glabel[id]);
         }
         if (noxr[id]) {
            if (vis_cp->noxr[id]) {
               printf("error:  no x[%d] coord range (vis_axes)\n",id);
               *ios = ERR;
               break;
            } else {
               sp->gcmin[id] = vis_cp->xmin[id];
               sp->gcmax[id] = vis_cp->xmax[id];
            }
         } else {
            sp->gcmin[id] = xmin[id];
            sp->gcmax[id] = xmax[id];
         }
      }
      sp->gscalp = FALSE;
   } else {
      printf("error:  nothing there (vis_axes)\n");
      *ios = ERR;
   }
/* check for valid graph coordinate ranges and then draw axes */
#ifdef DEBUG
   printf("vis_axes:  now construct world coords\n");
#endif
   if (!*ios) {
#ifdef CURV
/*   Modify graph coordinates is curvilinear */
#ifdef DEBUG
      printf("vis_axes:  modify world coordinates?\n");
#endif
      if (sp->curv) {
#ifdef DEBUG
         printf("vis_axes:  sp->curv=TRUE\n");
         printf("vis_axes:  sp->ndc=%d\n", sp->ndc);
#endif
         if (sp->ndc == 2) {
            ccopy("x", 0, SZNAM, sp->gname[0]);
            ccopy("x", 0, SZNAM, sp->glabel[0]);
            ccopy("y", 0, SZNAM, sp->gname[1]);
            ccopy("y", 0, SZNAM, sp->glabel[1]);
         } else if (sp->ndc==3) {
            ccopy("x", 0, SZNAM, sp->gname[0]);
            ccopy("x", 0, SZNAM, sp->glabel[0]);
            ccopy("y", 0, SZNAM, sp->gname[1]);
            ccopy("y", 0, SZNAM, sp->glabel[1]);
            ccopy("z", 0, SZNAM, sp->gname[2]);
            ccopy("z", 0, SZNAM, sp->glabel[2]);
         }
      }
#endif
/*   construct world coordinates */
#ifdef CURV
/*    Curvilinear case, requires cartesian conversion */
      if (sp->curv) {
         wclen = 0.0;
/*      First add a border to graph coordinates */
#ifdef DEBUG
         printf("vis_axes:  curvilinear,ndc=%d\n", sp->ndc);
#endif
         for (id=0; id < sp->nd; id++) {
            gclen = sp->gcmax[id] - sp->gcmin[id];
            sp->gcmin[id] = sp->gcmin[id] - border*gclen;
            sp->gcmax[id] = sp->gcmax[id] + border*gclen;
#ifdef DEBUG
            printf("..id=%d,gcmin=%f,gcmax=%f\n", id, sp->gcmin[id],
              sp->gcmax[id]);
#endif
         }
/*      Now convert 2D curv to 2D cartesian */
         if (sp->ndc == 2) {
#ifdef DEBUG
            printf("vis_axes:  convert 2D curv to 2D cartesian\n");
#endif
            dc0 = (sp->gcmax[0] - sp->gcmin[0])/(float)MWC;
            dc1 = (sp->gcmax[1] - sp->gcmin[1])/(float)MWC;
            sp->wcmin[0] = BIG;
            sp->wcmax[0] = (-BIG);
            sp->wcmin[1] = BIG;
            sp->wcmax[1] = (-BIG);
            for (i=0; i < MWC; i++) {
               u[0] = (float)(i)*dc0 + sp->gcmin[0];
               for (j=0; j < MWC; j++) {
                  u[1] = (float)(j)*dc1 + sp->gcmin[1];
                  geo_x(2, u, w);
                  sp->wcmin[0] = float_min(sp->wcmin[0],w[0]);
                  sp->wcmax[0] = float_max(sp->wcmax[0],w[0]);
                  sp->wcmin[1] = float_min(sp->wcmin[1],w[1]);
                  sp->wcmax[1] = float_max(sp->wcmax[1],w[1]);
               }
            }
/*      Or else convert 2D curv to 3D cartesians */
         } else if ((ndom == 3) && (ndc==2)) {
#ifdef DEBUG
            printf("vis_axes:  convert 2D curv to 3D cartesians\n");
#endif
            dc0 = (sp->gcmax[0] - sp->gcmin[0])/(float)MWC;
            dc1 = (sp->gcmax[1] - sp->gcmin[1])/(float)MWC;
            sp->wcmin[0] = BIG;
            sp->wcmax[0] = (-BIG);
            sp->wcmin[1] = BIG;
            sp->wcmax[1] = (-BIG);
            sp->wcmin[2] = BIG;
            sp->wcmax[2] = (-BIG);
            zm[2] = xo2;
#ifdef DEBUG
            printf("vis_axes:  xo2=%f\n", xo2);
#endif
            for (i=0; i < MWC; i++) {
               zm[0] = (float)(i)*dc0 + sp->gcmin[0];
               for (j=0; j < MWC; j++) {
                  zm[1] = (float)(j)*dc1 + sp->gcmin[1];
                  for (id=0; id < 3; id++) {
                     jd = kad[id];
                     u[jd] = zm[id];
                  }
                  geo_x(3, u, w);
                  sp->wcmin[0] = float_min(sp->wcmin[0],w[0]);
                  sp->wcmax[0] = float_max(sp->wcmax[0],w[0]);
                  sp->wcmin[1] = float_min(sp->wcmin[1],w[1]);
                  sp->wcmax[1] = float_max(sp->wcmax[1],w[1]);
                  sp->wcmin[2] = float_min(sp->wcmin[2],w[2]);
                  sp->wcmax[2] = float_max(sp->wcmax[2],w[2]);
               }
            }
/*      Or else deal with 3D curvilinears */
         } else {
#ifdef DEBUG
            printf("vis_axes:  deal with 3D curv to 3D cartesians\n");
#endif
            dc0 = (sp->gcmax[0] - sp->gcmin[0])/(float)MWC;
            dc1 = (sp->gcmax[1] - sp->gcmin[1])/(float)MWC;
            dc2 = (sp->gcmax[2] - sp->gcmin[2])/(float)MWC;
            sp->wcmin[0] = BIG;
            sp->wcmax[0] = (-BIG);
            sp->wcmin[1] = BIG;
            sp->wcmax[1] = (-BIG);
            sp->wcmin[2] = BIG;
            sp->wcmax[2] = (-BIG);
            for (i=0; i < MWC; i++) {
               u[0] = (float)(i)*dc0 + sp->gcmin[0];
               for (j=0; j < MWC; j++) {
                  u[1] = (float)(j)*dc1 + sp->gcmin[1];
                  for (k=0; k < MWC; k++) {
                     u[2] = (float)(j)*dc2 + sp->gcmin[2];
                     geo_x(3, u, w);
                     sp->wcmin[0] = float_min(sp->wcmin[0],w[0]);
                     sp->wcmax[0] = float_max(sp->wcmax[0],w[0]);
                     sp->wcmin[1] = float_min(sp->wcmin[1],w[1]);
                     sp->wcmax[1] = float_max(sp->wcmax[1],w[1]);
                     sp->wcmin[2] = float_min(sp->wcmin[2],w[2]);
                     sp->wcmax[2] = float_max(sp->wcmax[2],w[2]);
                  }
               }
            }
         }
/*      Add the scalar? */
         if (sp->gscalp && (sp->ndc==2)) {
#ifdef DEBUG
            printf("vis_axes:  scalar ranges\n");
#endif
            id = sp->ndc;
            wclen = surface_height*sqrt(wclen);
            sp->wcmin[id] = 0.0;
            sp->wcmax[id] = wclen;
         }
/*      Now a quick check on degeneracy */
         for (id=0; id < sp->nd; id++) {
#ifdef DEBUG
            printf("vis_axes:  id=%d,wcmin=%f,wcmax=%f\n", id,
              sp->wcmin[id],sp->wcmax[id]);
#endif
            if (sp->wcmax[id] <= sp->wcmin[id]) {
               printf("error:  world coordinate %s range degenerate\n",
                 sp->gname[id]);
               *ios = ERR;
            }
         }
/*    Previous cartesian-like case */
      } else {
#ifdef DEBUG
         printf("vis_axes:  thinks it is cartesian\n");
#endif
         wclen = 0.0;
         for (id=0; id < sp->nd; id++) {
            gclen = sp->gcmax[id] - sp->gcmin[id];
            sp->gcmin[id] = sp->gcmin[id] - border*gclen;
            sp->gcmax[id] = sp->gcmax[id] + border*gclen;
            if (sp->gscalp && (sp->ndc == 2) && (id == sp->ndc)) {
               wclen = surface_height*sqrt(wclen);
               sp->wcmin[id] = 0.0;
               sp->wcmax[id] = wclen;
            } else {
               sp->wcmin[id] = sp->gcmin[id];
               sp->wcmax[id] = sp->gcmax[id];
#ifdef WCLEN_BUG
               wclen = wclen + gclen*gclen;
#endif
            }
            if (sp->wcmax[id] <= sp->wcmin[id]) {
               printf("error:  world coordinate %s range degenerate\n",
                 sp->gname[id]);
               *ios = ERR;
            }
         }
      }
#else
      wclen = 0.0;
      for (id=0; id < sp->nd; id++) {
         gclen = sp->gcmax[id] - sp->gcmin[id];
         sp->gcmin[id] = sp->gcmin[id] - border*gclen;
         sp->gcmax[id] = sp->gcmax[id] + border*gclen;
         if (sp->gscalp && (sp->ndc == 2) && (id == sp->ndc)) {
            wclen = surface_height*sqrt(wclen);
            sp->wcmin[id] = 0.0;
            sp->wcmax[id] = wclen;
         } else {
            sp->wcmin[id] = sp->gcmin[id];
            sp->wcmax[id] = sp->gcmax[id];
#ifdef WCLEN_BUG
            wclen = wclen + gclen*gclen;
#endif
         }
         if (sp->wcmax[id] <= sp->wcmin[id]) {
            printf("error:  world coordinate %s range degenerate\n",
              sp->gname[id]);
            *ios = ERR;
         }
      }
#endif
      if (*ios) {
         printf("error:  not drawing axes!\n");
      } else if (sp->nd == 2) {
         vis_design_2d_axes(ios);
      } else if (sp->nd == 3) {
         vis_design_3d_axes(ios);
      } else {
         printf("error:  cannot draw %d-dimensional graph axes\n",
           sp->nd);
         *ios = ERR;
      }
      if (!*ios) {
         sp->clear = FALSE;
         sp->drawn = TRUE;
         sp->noqr = TRUE;
         sp->nocont = TRUE;
         sp->nocmap = TRUE;
         sp->novel = TRUE;
      }
   }
   drp_set_hv(sp->hmin, sp->hmax, sp->vmin, sp->vmax);
}

void vis_design_2d_axes(ios)
   int *ios;
{
   float unit,hclen[MAXD],min_unit,hl[MAXD],hr[MAXD],wclen[MAXD],db,
     h1,h2,v1,v2,h0,dhdx,v0,dvdy,x0,y0,xoff,xdiv,yoff,ydiv,
     h,v,dvt,dvn,dvl,x,y,dvc,scale,dve;
   char str[SZLIN],just[SZNAM];
   int i,id,n,ndiv,max_div=4,ipower,field,ndp,nspace;
   float r,g,b;
   v = 0.0;
   for (i=0; i < MAXD; i++)
     hclen[i] = 0.0;

   drp_get_rgb(&r,&g,&b);
/* determine number of "space" dimensions */
   nspace = 0;
   for (id=0; id < 2; id++) {
      if (sp->gspace[id]) {
         nspace++;
      }
   }
/* graph extent (subwindow minus borders) */
   hclen[0] = sp->hmax -  sp->hmin;
   hclen[1] = sp->vmax -  sp->vmin;
   sp->hgmin = sp->hmin + hborder*hclen[0];
   sp->hgmax = sp->hmax - hborder*hclen[0];
   sp->vgmin = sp->vmin + vborder*hclen[1];
   sp->vgmax = sp->vmax - vborder*hclen[1];
/* some attributes */
   sp->key_number_field = key_number_field;
   sp->key_label_field = key_label_field;
/* design scalesbox and let it restrict graph box */
   if (nspace < 2) {
      sp->noscales = TRUE;
      sp->nokey = nokey;
   } else if (!noscales || !nokey)  {
      vis_design_scales();
      sp->hgmax = sp->hxmin;
   } else {
      sp->noscales = TRUE;
      sp->nokey = TRUE;
   }
/* set clipping viewport */
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
/* initialise borders within graph box */
   hl[0] = 0.0;
   hr[0] = 0.0;
/* left h-border (ylabel|ynumbers) */
   if (label_field > number_field) {
      n = label_field;
   } else {
      n = number_field;
   }
   hl[0] = hl[0] + ((float)n+0.5)*wp->char_width;
/* safety margins */
   hl[0] = hl[0] + 0.5*wp->char_width;
   hr[0] = hr[0] + 0.5*wp->char_width;
/* tick marks outside? */
   if (tick_scale < 0.0) {
      hl[0] = hl[0] - tick_scale*wp->char_height;
      hr[0] = hr[0] - tick_scale*wp->char_height;
   }
/* construct h-extent */
   hclen[0] = (sp->hgmax - sp->hgmin) - (hl[0]+hr[0]);
/* initialise vertical borders within graph box */
   hl[1] = 0.0;
   hr[1] = 0.0;
/* clip to graph box */
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
/* half character internal border */
   hl[1] = hl[1] + 0.5*wp->char_height;
   hr[1] = hr[1] + 0.5*wp->char_height;
/* restrict number of caption lines (really) */
   sp->caption_lines = caption_lines;
/* caption below  */
   if (caption_below) {
      hl[1] = hl[1] + (1.5*sp->caption_lines+0.5)*wp->char_height;
/* caption above */
   } else {
      hr[1] = hr[1] + (1.5*sp->caption_lines+0.5)*wp->char_height;
   }
/* ticks (if outside) */
   if (tick_scale < 0.0) {
      hl[1] = hl[1] - tick_scale*wp->char_height;
      hr[1] = hr[1] - tick_scale*wp->char_height;
   }
/* axis numbers and label below */
   hl[1] = hl[1] + 3.0*wp->char_height;
/* key above */
   hr[1] = hr[1] + 2.0*wp->char_height;
/* construct v-extent for this */
   hclen[1] = (sp->vgmax - sp->vgmin) - (hl[1]+hr[1]);
/* determine any space-units consistent with aspect ratios */
   min_unit = BIG;
   nspace = 0;
   for (id=0; id < sp->nd; id++) {
      wclen[id] = sp->wcmax[id] - sp->wcmin[id];
      if (sp->gspace[id]) {
         nspace++;
         unit = hclen[id]/(sp->gasp[id]*wclen[id]);
         if (unit < min_unit) {
            min_unit = unit;
         }
      }
   }
/* apply any aspect ratio effects to space coords */
   for (id=0; id < sp->nd; id++) {
      if (sp->gspace[id]) {
         db = hclen[id] - min_unit*sp->gasp[id]*wclen[id];
         hl[id] = hl[id] + 0.5*db;
         hr[id] = hr[id] + 0.5*db;
         hclen[id] = hclen[id] - db;
      }
   }
/* define data window */
   sp->hdmin = sp->hgmin + hl[0];
   sp->hdmax = sp->hgmax - hr[0];
   sp->vdmin = sp->vgmin + hl[1];
   sp->vdmax = sp->vgmax - hr[1];
/* design inner key box */
   if (!nokey && (nspace == 1)) {
      vis_design_key();
   }
   h1 = sp->hdmin;
   h2 = sp->hdmax;
   v1 = sp->vdmin;
   v2 = sp->vdmax;
/* capture the (x,y) <> (h,v) relations */
   sp->x0 = sp->wcmin[0];
   sp->dhdx = hclen[0]/wclen[0];
   sp->h0 = h1;
   sp->y0 = sp->wcmin[1];
   sp->dvdy = hclen[1]/wclen[1];
   sp->v0 = v1;
/* design x-axis divisions */
   vis_divs(sp->wcmin[0], sp->wcmax[0], max_div, &xoff, &xdiv);
   vis_ndiv(sp->wcmin[0], sp->wcmax[0], xoff, xdiv, &ndiv);
/* x-axes, ticks and numbers */
/* tick extent */
   dvt = -tick_scale*wp->char_height;
/* gap for number */
   dvn = 0.5*wp->char_height;
/* increase number gap if tick on outside */
   if (dvt > 0.0) {
      dvn = dvn + dvt;
   }
   dhdx = hclen[0]/(sp->wcmax[0] - sp->wcmin[0]);
   x0 = sp->wcmin[0];
   cvector(ndiv, (x0+xoff), xdiv, &ipower, &scale, &field, &ndp);
   h0 = h1;
   drp_line(h1, v1, h2, v1);
   if (nspace == 2) {
      drp_line(h1, v2, h2, v2);
   }
   for (i=0; i < (ndiv+1); i++) {
      x = (int)i*xdiv + (x0+xoff);
      h = (x-x0)*dhdx + h0;
      drp_line(h, v1, h, v1-dvt);
      x = x/scale;
      cfloat(x, field, ndp, SZLIN, str);
      srmblk(str);
      drp_text(h, v1-dvn, str, "mt");
      if (nspace == 2) {
         drp_line(h, v2, h, v2+dvt);
      }
   }
/* x-label */
   dvl = dvn + 1.5*wp->char_height;
   if (label_centre) {
      h = 0.5*(h1+h2);
      drp_text(h, v1-dvl, sp->glabel[0], "mt");
   } else {
      drp_text(h2, v1-dvl, sp->glabel[0], "rt");
   }
/* xpower? */
   if (ipower != 0) {
      h = h2 - 6.0*wp->char_width;
      drp_text(h, v1-dvl, "x10", "lt");
      dvc = dvl - 0.5*wp->char_height;
      sprintf(str, "%d", ipower);
      drp_text(h+3.0*wp->char_width, v1-dvc, str, "lt");
   }
/* caption below */
   if (caption_below) {
      dvc = dvl + 2.0*wp->char_height;
      sp->caption_h = sp->hgmin + (float)caption_inset*wp->char_width;
      sp->caption_v = sp->vdmin-dvc;
      ccopy("lt", 0, SZNAM, just);
   } else {
      dvc = (1.5*(float)sp->caption_lines+0.5)*wp->char_height;
      if (tick_scale < 0.0) {
         dvc = dvc - tick_scale*wp->char_height;
      }
      sp->caption_h = sp->hgmin + (float)caption_inset*wp->char_width;
      sp->caption_v = sp->vdmax+dvc;
      ccopy("lb", 0, SZNAM, just);
   }
   for (i=0; i < sp->caption_count; i++) {
      if (i < sp->caption_lines) {
         drp_text(sp->caption_h, sp->caption_v, sp->caption[i], just);
         sp->caption_v = sp->caption_v - 1.5*wp->char_height;
      } else {
         printf("note:  lost caption line\n");
      }
   }
/* design y-axis divisions */
   vis_divs(sp->wcmin[1], sp->wcmax[1], max_div, &yoff, &ydiv);
   vis_ndiv(sp->wcmin[1], sp->wcmax[1], yoff, ydiv, &ndiv);
/* y-axes, ticks and numbers */
/* tick extent */
   dvt = -tick_scale*wp->char_height;
/* gap for number */
   dvn = 0.5*wp->char_width;
/* increase number gap if tick on outside */
   if (dvt > 0.0) {
      dvn = dvn + dvt;
   }
   dvdy = hclen[1]/(sp->wcmax[1] - sp->wcmin[1]);
   y0 = sp->wcmin[1];
   cvector(ndiv, (y0+yoff), ydiv, &ipower, &scale, &field, &ndp);
   drp_line(h1,v1,h1,v2);
   if (nspace == 2) {
      drp_line(h2, v1, h2, v2);
   }
   v0 = v1;
   for (i=0; i < (ndiv+1); i++) {
      y = (int)i*ydiv + (y0+yoff);
      v = (y-y0)*dvdy + v0;
      drp_line(h1, v, h1-dvt, v);
      y = y/scale;
      cfloat(y, field, ndp, SZLIN, str);
      drp_text(h1-dvn, v, str, "rm");
      if (nspace == 2) {
         drp_line(h2, v, h2+dvt, v);
      }
   }
/* y-axis label */
   dvl = 0.5*wp->char_width;
   if (dvt > 0.0) {
      dvl = dvl + dvt;
   }
   v = v - 0.5*ydiv*dvdy;
   drp_text(h1-dvl, v, sp->glabel[1], "rm");
/* y-number power */
   if (ipower != 0) {
      dvc = 6.0*wp->char_width;
      dve = 0.5*wp->char_height;
      drp_text(h1-dvc, v2+dve, "x10", "lb");
      dvc = dvc - 3.0*wp->char_width;
      sprintf(str, "%d", ipower);
      dve = 2.0*dve;
      drp_text(h1-dvc, v2+dve, str, "lb");
   }
   drp_get_rgb(&r,&g,&b);
}

void vis_design_3d_axes(ios)
   int *ios;
{
   float unit,hclen[MAXD],min_unit,hl[MAXD],hr[MAXD],wclen[MAXD],db,
     h1,h2,v1,v2;
   char just[SZNAM];
   int i,id;

/* graph box less borders */
   hclen[0] = (sp->hmax - sp->hmin);
   hclen[1] = (sp->vmax - sp->vmin);
   sp->hgmin = sp->hmin + hborder*hclen[0];
   sp->hgmax = sp->hmax - hborder*hclen[0];
   sp->vgmin = sp->vmin + vborder*hclen[1];
   sp->vgmax = sp->vmax - vborder*hclen[0];
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
/* some attributes */
   sp->key_number_field = key_number_field;
   sp->key_label_field = key_label_field;
/* account for any scales-keys */
   if (!noscales || !nokey) {
      vis_design_scales();
      sp->hgmax = sp->hxmin;
   } else {
      sp->noscales = TRUE;
      sp->nokey = TRUE;
   }
/* proceed to map this (X,Y) view-surface space onto (h,v) */
   wclen[0] = view_aspect;
   wclen[1] = 1.0;
/* safety margins */
   hl[0] = 0.5*wp->char_width;
   hr[0] = 0.5*wp->char_width;
   hl[1] = 0.5*wp->char_height;
   hr[1] = 0.5*wp->char_height;
/* caption below  */
   sp->caption_lines = caption_lines;
   if (caption_below) {
      hl[1] = hl[1] + (1.5*sp->caption_lines+0.5)*wp->char_height;
/* caption above */
   } else {
      hr[1] = hr[1] + (1.5*sp->caption_lines+0.5)*wp->char_height;
   }
/* construct extents subject to these */
   hclen[0] = (sp->hgmax - sp->hgmin) - (hl[0]+hr[0]);
   hclen[1] = (sp->vgmax - sp->vgmin) - (hl[1]+hr[1]);
/* aspect rations in [h,v] */
   min_unit = BIG;
   for (id=0; id < 2; id++) {
      unit = hclen[id]/wclen[id];
      if (unit < min_unit) {
         min_unit = unit;
      }
   }
/* apply aspect rations */
   for (id=0; id < 2; id++) {
      db = hclen[id] - min_unit*wclen[id];
      hl[id] = hl[id] + 0.5*db;
      hr[id] = hr[id] + 0.5*db;
      hclen[id] = hclen[id] - db;
   }
/* draw data window */
   h1 = sp->hgmin + hl[0];
   h2 = h1 + hclen[0];
   v1 = sp->vgmin + hl[1];
   v2 = v1 + hclen[1];
   sp->hdmin = h1;
   sp->hdmax = h2;
   sp->vdmin = v1;
   sp->vdmax = v2;
/* caption below */
   sp->caption_h = sp->hgmin + (float)caption_inset*wp->char_width;
   ccopy("lb", 0, SZNAM, just);
   if (caption_below) {
/*
      dvc = 0.5*wp->char_height;
      sp->caption_v = sp->vdmin-dvc;
      ccopy("lt", 0, SZNAM, just);
*/
      sp->caption_v = sp->vgmin + ((float)(caption_lines-1)+0.5)*
         wp->char_height;
   } else {
/*
      dvc = (1.5*(float)sp->caption_lines+0.5)*wp->char_height;
      sp->caption_v = sp->vdmax+dvc;
      ccopy("lb", 0, SZNAM, just);
*/
      sp->caption_v = sp->vgmax - 2.0*wp->char_height;
   }
   for (i=0; i < sp->caption_count; i++) {
      if (i < sp->caption_lines) {
         drp_text(sp->caption_h, sp->caption_v, sp->caption[i], just);
         sp->caption_v = sp->caption_v - 1.5*wp->char_height;
      } else {
         printf("note:  caption line lost\n");
      }
   }
/* store reference data for key and colour scale */
   vis_design_3d_view(ios);
}

void vis_divs(vmin,vmax,max_div,voff,vdiv)
   float vmin,vmax,*voff,*vdiv;
   int max_div;
{
   int idiv,ipower,i1;
   float vd,power,d1;

/* stab at division length */
   if (vmax > vmin) {
      vd = (vmax-vmin)/(float)max_div;
   } else {
      vd = (vmin-vmax)/(float)max_div;
   }

/* get suitable base-10 exponent */
/*
   glog10(vd,&power);
*/
   power = log10(vd);
   ipower = (int)power;
   power = (float)ipower;
/*
   gpow(10.0, power, &power);
*/
   power = pow(10.0, power);
   if (power >= vd) {
      ipower = ipower - 1;
      power = 0.1*power;
   }

/* choose a division size (0.2,0.5,1.0)*10**ipower */
   idiv = (int)(vd/power);
   if (idiv > 5) {
      idiv = 10;
   } else if (idiv > 2) {
      idiv = 5;
   } else {
      idiv = 2;
   }
   *vdiv = (float)idiv*power;
   if (vmax < vmin) {
      *vdiv = -(*vdiv);
   }

/* suitable offset to first tick */
   d1 = vmin/(*vdiv);
   i1 = (int)d1;
   if (d1 > 0.0) {
      i1 = i1 + 1;
   }
   *voff = (int)i1*(*vdiv) - vmin;
}

void vis_ndiv(vmin,vmax,voff,vdiv,ndiv)
   float vmin,vmax,voff,vdiv;
   int *ndiv;
{
   float vfree,vgap;

   vfree = vmax - (vmin+voff);
   *ndiv = (int)(vfree/vdiv);
   vgap = (float)(*ndiv)*vdiv - vfree;
   if (vgap*vgap < 1.0E-6*vdiv*vdiv) {
      *ndiv = *ndiv -1;
   }
}

void vis_design_scales() {
   float dhc,dhkl,dhs,dh,dv, dvk,dvz;
   int n;

/* copy some bits */
   sp->key_number_field = key_number_field;
   sp->key_label_field = key_label_field;
/* width required for scalar contour/colourmap labels */
   dhc = ((float)key_number_field + 6.0)*wp->char_width;
/* width required for key label field */
   dhkl = ((float)key_label_field + 3.5)*wp->char_width;
/* decide which width dominates */
   dhs = dhc;
   if (dhkl > dhs) {
      dhs = dhkl;
   }
/* horizontal and vertical extent */
   dh = (sp->hmax - sp->hmin);
   dv = (sp->vmax - sp->vmin);
   sp->hxmax = sp->hmax - hborder*dh;
   sp->hxmin = sp->hxmax - dhs;
   sp->vxmin = sp->vmin + vborder*dv;
   sp->vxmax = sp->vmax - vborder*dv;
/* horizontal extents are boring */
   sp->hkmin = sp->hxmin;
   sp->hkmax = sp->hxmax;
   sp->hzmin = sp->hxmin;
   sp->hzmax = sp->hxmax;
/* total vertical extent available */
   dv = sp->vxmax - sp->vxmin;
   sp->nokey = nokey;
   sp->noscales = noscales;
/* no key means must have scales only */
   if (nokey) {
      if (dv < 12.0*wp->char_height) {
         printf("error:  no room for scales\n");
         sp->noscales = TRUE;
      } else {
         sp->vzmin = sp->vxmin;
         sp->vzmax = sp->vxmax;
      }
/* noscales means we must have key alone */
   } else if (noscales) {
      sp->vkmin = sp->vxmin;
      sp->vkmax = sp->vxmax;
      sp->key_size = (int)(dv/(1.5*wp->char_height));
/* otherwise got key AND scales */
   } else {
      n = key_size;
      dvk = 1.5*wp->char_height*(float)n;
      dvz = dv - dvk;
      if (dv < 12.0*wp->char_height) {
         printf("error:  no room for scales\n");
         sp->noscales = TRUE;
         dvk = dv;
      } else if (dvz < 12.0*wp->char_height) {
         printf("error:  losing key entries for scales\n");
         dvz = 12.0*wp->char_height;
         dvk = dv - dvz;
      }
      sp->vzmax = sp->vxmax;
      sp->vzmin = sp->vxmax - dvz;
      sp->key_size = (int)(dvk/(1.5*wp->char_height));
      if (sp->key_size < 0) {
         sp->nokey = TRUE;
      }
      sp->vkmin = sp->vxmin;
      sp->vkmax = sp->vkmin + dvk;
   }
   sp->key_count = 0;
/* the scale area itself */
   if (!sp->noscales) {
      sp->hsmin = sp->hzmin + 3.0*wp->char_width;
      sp->hsmax = sp->hsmin + 2.0*wp->char_width;
      sp->vsmax = sp->vzmax - 4.5*wp->char_height;
      sp->vsmin = sp->vzmin + 4.5*wp->char_height;
   }
   sp->key_drawn = FALSE;
   sp->scales_drawn = FALSE;
   sp->cmap_scales_drawn = FALSE;
   sp->contour_scales_drawn = FALSE;
   sp->vel_scales_drawn = FALSE;
   sp->mat_scales_drawn = FALSE;
}

/* vis_scales - actually draw the scales */
void vis_scales() {
   float v1,v2,h,v,dh,dv;
   char str[SZLIN];
   int n;

/* set clipping port */
   drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
/* fill it with "sky" */
/*
   if (wp->colour_avail) {
      drp_get_colour(&colour_index);
      drp_set_colour(wp->colour_sky);
      drp_rect(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      drp_set_colour(colour_index);
   }
*/
/* draw the information outer box box */
/*
   vis_box(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
*/
/* draw the scales outer box box */
   vis_box(sp->hzmin, sp->hzmax, sp->vzmin, sp->vzmax);
/* initial outer box for scale */
/*
   vis_box(sp->hsmin, sp->hsmax, sp->vsmin, sp->vsmax);
*/
/* select number format for scalars */
   cvector(1, sp->qmin, (sp->qmax-sp->qmin), &sp->qpower, &sp->qscale,
     &sp->qfield, &sp->qdp);
   n = (sp->qfield - sp->key_number_field);
   if (n > 0) {
      if (n > sp->qdp) {
         n = sp->qdp;
      }
      sp->qfield = sp->qfield - n;
      sp->qdp = sp->qdp - n;
   }
/* min scale limit */
   drp_bp();
   h = sp->hsmax;
   v = sp->vsmin;
   drp_ap(h, v);
   v = v - 1.5*wp->char_height;
   drp_ap(h, v);
   h = h + 0.5*wp->char_width;
   drp_ap(h, v);
   drp_ep();
   cfloat(sp->qmin/sp->qscale, sp->qfield, sp->qdp, SZLIN, str);
   drp_text(h, v, str, "lm");
/* max scale limit */
   drp_bp();
   h = sp->hsmax;
   v = sp->vsmax;
   drp_ap(h, v);
   v = v + 1.5*wp->char_height;
   drp_ap(h, v);
   h = h + 0.5*wp->char_width;
   drp_ap(h, v);
   drp_ep();
   cfloat(sp->qmax/sp->qscale, sp->qfield, sp->qdp, SZLIN, str);
   drp_text(h, v, str, "lm");
/* scale factor and units */
   v2 = sp->vzmin + 2.0*wp->char_height;
/*
   vis_box(sp->hzmin, sp->hzmax, sp->vzmin, v2);
*/
   h = sp->hzmin + 0.5*wp->char_width;
   v = sp->vzmin + 0.25*wp->char_height;
   if (sp->qpower != 0) {
      drp_text(h, v, "x10", "lb");
      drp_text_box("x10", &dh, &dv);
      h = h + dh;
      v = v + 0.5*wp->char_height;
      sprintf(str, "%d", sp->qpower);
      drp_text(h, v, str, "lb");
      drp_text_box(str, &dh, &dv);
      h = h + dh + wp->char_width;
      v = v - 0.5*wp->char_height;
   }
   drp_text(h, v, sp->qunits, "lb");
/* label at top */
   v2 = sp->vzmax;
   v1 = v2 - 2.0*wp->char_height;
/*
   vis_box(sp->hzmin, sp->hzmax, v1, v2);
*/
   h = sp->hzmin + 0.5*wp->char_width;
   v = v1 + 0.25*wp->char_height;
   drp_text(h, v, sp->qlabel, "lb");
/* we have drawn scales structure */
   sp->scales_drawn = TRUE;
   sp->cmap_scales_drawn = FALSE;
   sp->contour_scales_drawn = FALSE;
}

void vis_design_key() {
   float dhk;

/* copy some bits */
   sp->key_label_field = key_label_field;
   sp->key_number_field = key_number_field;
   sp->key_size = key_size;
/* width required for key label field */
   dhk = ((float)key_label_field + 3.5)*wp->char_width;
/* tuck it into corner of data box */
   sp->hkmax = sp->hdmax;
   sp->hkmin = sp->hkmax - dhk;
   sp->vkmax = sp->vdmax;
   sp->vkmin = sp->vkmax - 1.5*wp->char_height*sp->key_size;
   sp->key_count = 0;
   ccopy(key_position, 0, SZNAM, sp->key_position);
   sp->key_drawn = FALSE;
}

void vis_key() {
   drp_set_hv(sp->hmin, sp->hmax, sp->vmin, sp->vmax);
   vis_box(sp->hkmin, sp->hkmax, sp->vkmin, sp->vkmax);
   sp->key_drawn = TRUE;
}

void vis_draw(ios)
   int *ios;
{
   bool temp_noline_connect,temp_noline_marker;
   char temp_qlabel[SZLIN];
   int id;

/* check errors */
   if (*ios) {
/* 2-dimensional world coordinates */
   } else if (sp->nd == 2) {
/*   scalar-promoted */
      if (sp->gscalp) {
         if (!nogrid) {
#ifdef CURV_OVER
            sp->curv = vis_gvp->curv && (vis_gvp->ndc >= 2);
#endif
            if (!vis_gvp->noscalar && (vis_gvp->nd == 1)) {
/*
               if (vis_gvp->crunched) {
                  temp_noline_connect = noline_connect;
                  temp_noline_marker = noline_marker;
                  noline_connect = TRUE;
                  noline_marker = FALSE;
               }
*/
               vis_line_graph(vis_gvp->nc[0],vis_gvp->xa[0],vis_gvp->sa,
                 vis_gvp->qa, ios);
/*
               if (vis_gvp->crunched) {
                  noline_connect = temp_noline_connect;
                  noline_marker = temp_noline_marker;
               }
*/
            }
         }
         if (!noscatter) {
#ifdef CURV_OVER
            sp->curv = vis_svp->curv && (vis_svp->ndc >= 2);
#endif
            if (!vis_svp->noscalar && (vis_svp->nd == 1)) {
               temp_noline_connect = noline_connect;
               temp_noline_marker = noline_marker;
               noline_connect = TRUE;
               noline_marker = FALSE;
               vis_line_graph(vis_svp->npt,vis_svp->xpt,vis_svp->spt,
                 vis_svp->qpt, ios);
               noline_connect = temp_noline_connect;
               noline_marker = temp_noline_marker;
            }
         }
/*   unscalar promoted */
      } else {
         if (!nogrid) {
#ifdef CURV_OVER
            sp->curv = vis_gvp->curv && (vis_gvp->ndc >= 2);
#endif
            if (!vis_gvp->noscalar && (vis_gvp->nd == 2)) {
               if (!nocolour && wp->cmap_avail) {
                  if (noqlabel) {
                     ccopy(vis_gvp->qlabel, 0, SZLIN, temp_qlabel);
                  } else {
                     ccopy(qlabel, 0, SZLIN, temp_qlabel);
                  }
#ifdef CURV
                  if (((vis_gvp->nd+vis_gvp->ndo)==4) &&
                    vis_gvp->curv) {
#ifdef DEBUG
                     printf("vis_draw:  not calling colour_map\n");
#endif
                  } else {
#endif
                  vis_colour_map(vis_gvp->nc[0], vis_gvp->nc[1],
                    vis_gvp->xe[0], vis_gvp->xe[1], vis_gvp->sa,
                    vis_gvp->qa, vis_gvp->qmin, vis_gvp->qmax, 
                    temp_qlabel, vis_gvp->qunits, ios);
#ifdef CURV
                  }
#endif
               }
               if (!noband && wp->cmap_avail) {
                  if (noqlabel) {
                     ccopy(vis_gvp->qlabel, 0, SZLIN, temp_qlabel);
                  } else {
                     ccopy(qlabel, 0, SZLIN, temp_qlabel);
                  }
                  vis_band_map(vis_gvp->nc[0], vis_gvp->nc[1],
                    vis_gvp->xe[0], vis_gvp->xe[1],
                    vis_gvp->xa[0], vis_gvp->xa[1], vis_gvp->sa,
                    vis_gvp->qa, vis_gvp->qmin, vis_gvp->qmax, 
                    temp_qlabel, vis_gvp->qunits, ios);
               }
               if (!nocontour) {
                  if (noqlabel) {
                     ccopy(vis_gvp->qlabel, 0, SZLIN, temp_qlabel);
                  } else {
                     ccopy(qlabel, 0, SZLIN, temp_qlabel);
                  }
                  vis_contour_map(vis_gvp->nc[0], vis_gvp->nc[1],
                    vis_gvp->xa[0], vis_gvp->xa[1], vis_gvp->sa,
                    vis_gvp->qa, vis_gvp->qmin, vis_gvp->qmax, 
                    temp_qlabel, vis_gvp->qunits, ios);
                     vis_bound_map(vis_gvp->nc[0], vis_gvp->nc[1],
                       vis_gvp->xe[0], vis_gvp->xe[1], vis_gvp->sa,
                       ios);
               }
            }
            if (!vis_gvp->nogrid && (vis_gvp->nd == 2)) {
               vis_grid_map(vis_gvp->nc[0], vis_gvp->nc[1],
                 vis_gvp->xe[0], vis_gvp->xe[1], vis_gvp->sa,
                 ios);
            }
            if (!vis_gvp->novector && (vis_gvp->nd == 2)) {
               if (!noarrow) {
                  vis_arrow_map(vis_gvp->nc[0], vis_gvp->nc[1],
                    vis_gvp->xa[0], vis_gvp->xa[1], vis_gvp->sa,
                    vis_gvp->ndv, vis_gvp->va, vis_gvp->vmax, ios);
               }
               if (!nostream) {
#ifdef CURV
                  if (sp->curv) {
                     printf("vis_draw:  no 2D curv streamlines\n");
                  } else {
#endif
                  vis_streamer(vis_gvp->nd, vis_gvp->nc, vis_gvp->xe,
                    vis_gvp->sa, vis_gvp->ndv, vis_gvp->va,
                    vis_gvp->vmax, ios);
#ifdef CURV
                  }
#endif
               }
            }
         }
         if (!noscatter) {
#ifdef CURV_OVER
            sp->curv = vis_svp->curv && (vis_svp->ndc >= 2);
#endif
            if (!vis_svp->noscalar && (vis_svp->nd == 2)) {
               if (!nocolour && wp->cmap_avail) {
                  if (noqlabel) {
                     ccopy(vis_svp->qlabel, 0, SZLIN, temp_qlabel);
                  } else {
                     ccopy(qlabel, 0, SZLIN, temp_qlabel);
                  }
                  vis_colour_gram(vis_svp->npt,vis_svp->nd,vis_svp->xpt,
                    vis_svp->spt, vis_svp->qpt, vis_svp->qmin,
                    vis_svp->qmax, temp_qlabel, vis_svp->qunits, ios);
               } else {
                  vis_marker_gram(vis_svp->npt,vis_svp->nd,vis_svp->xpt,
                    vis_svp->spt, vis_svp->qpt, vis_svp->qmin,
                    vis_svp->qmax, ios);
               }
            }
            if (!vis_svp->novector && (vis_svp->nd == 2)) {
               vis_arrow_gram(vis_svp->npt, vis_svp->nd, vis_svp->xpt,
                 vis_svp->spt, vis_svp->ndv, vis_svp->vpt,
                 vis_svp->vmax, ios);
            }
         }
      }
/* 3-dimensional world coordinates */
   } else if (sp->nd == 3) {
/*    Dump some defaults */
      dr3_set0i("ncolour", ncolour);
      dr3_set0i("nshade", nshade);
/*    Scalar promoted */
      if (sp->gscalp) {
/*       Grid model (gvp) */
         if (!nogrid) {
#ifdef CURV_OVER
            sp->curv = vis_gvp->curv && (vis_gvp->ndc >= 2);
#endif
            if ((!vis_gvp->noscalar) && (vis_gvp->nd == 2)) {
               if (sp->noqr) {
                  if (noqr) {
                     sp->qmin = vis_gvp->qmin;
                     sp->qmax = vis_gvp->qmax;
                  } else {
                     sp->qmin = qmin;
                     sp->qmax = qmax;
                  }
                  ccopy(vis_gvp->qlabel, 0, SZNAM, sp->qlabel);
                  ccopy(vis_gvp->qunits, 0, SZNAM, sp->qunits);
                  sp->noqr = sp->qmin >= sp->qmax;
               }
               if (sp->noqr) {
                  printf("warning:  no scalar range\n");
               } else {
                  printf("Scalar range %E to %E\n", sp->qmin, sp->qmax);
                  vis_3d_grid_surface(vis_gvp->nd, vis_gvp->nc,
                    vis_gvp->xe, vis_gvp->xa, vis_gvp->sa, vis_gvp->qa,
                    sp->qmin, sp->qmax, ios);
               }
            }
         }
/*    Unscalar promoted */
      } else {
/*       Solid model (cp) */
         if (!nocsg) {
#ifdef CURV_OVER
            sp->curv = FALSE;
#endif
            if (wire_frame) {
               dr3_seta("wire_frame");
            }
            vis_3d_solid(ios);
            dr3_seta("nowire_frame");
         }
/*       Grid model (gvp) */
         if (!nogrid) {
#ifdef CURV_OVER
            sp->curv = vis_gvp->curv && (vis_gvp->ndc >= 2);
#endif
/*          Grid field in 3D */
            if ((!vis_gvp->nogrid) && (vis_gvp->nd == 3)) {
               vis_3d_draw_grid(vis_gvp->nc,vis_gvp->xe,vis_gvp->sa,
                 ios);
            }
/*          Scalar field in 3D */
            if ((!vis_gvp->noscalar) && (vis_gvp->nd == 3)) {
               if (sp->noqr) {
                  if (noqr) {
                     sp->qmin = vis_gvp->qmin;
                     sp->qmax = vis_gvp->qmax;
                  } else {
                     sp->qmin = qmin;
                     sp->qmax = qmax;
                  }
                  ccopy(vis_gvp->qlabel, 0, SZNAM, sp->qlabel);
                  ccopy(vis_gvp->qunits, 0, SZNAM, sp->qunits);
                  sp->noqr = sp->qmin >= sp->qmax;
               }
               if (sp->noqr) {
                  printf("warning:  no scalar range\n");
               } else {
                  printf("Scalar range %E to %E\n", sp->qmin, sp->qmax);
                  if (data_shade) {
                     dr3_loads(vis_gvp->nc,vis_gvp->xe,vis_gvp->sa,
                       vis_gvp->qa,sp->qmin,sp->qmax,ios);
                     vis_3d_data_shade(ios);
                     dr3_zaps();
                  } else if (data_plane) {
                     dr3_loads(vis_gvp->nc,vis_gvp->xe,vis_gvp->sa,
                       vis_gvp->qa,sp->qmin,sp->qmax,ios);
                     vis_3d_data_plane(ios);
                     dr3_zaps();
                  } else if (data_ortho) {
                     vis_3d_data_ortho(vis_gvp->nc, vis_gvp->xe,
                       vis_gvp->sa, vis_gvp->qa, sp->qmin, sp->qmax,
                       ios);
                  }
               }
#ifdef CURV
/*          Grid field in 2D (3D perspective) */
            } else if ((!vis_gvp->nogrid) && (vis_gvp->nd==2)) {
               if (((vis_gvp->nd+vis_gvp->ndo)==4) &&
                 vis_gvp->curv) {
#ifdef DEBUG
                  printf("vis_draw:  calling 2d_curv_grid_map\n");
#endif
                  vis_2d_curv_grid_map(vis_gvp->nc, vis_gvp->xe,
                    vis_gvp->sa, ios);
               }
/*          Scalar field in 2D (3D perspective) */
            } else if ((!vis_gvp->noscalar) && (vis_gvp->nd==2)) {
               if (sp->noqr) {
                  if (noqr) {
                     sp->qmin = vis_gvp->qmin;
                     sp->qmax = vis_gvp->qmax;
                  } else {
                     sp->qmin = qmin;
                     sp->qmax = qmax;
                  }
                  ccopy(vis_gvp->qlabel, 0, SZNAM, sp->qlabel);
                  ccopy(vis_gvp->qunits, 0, SZNAM, sp->qunits);
                  sp->noqr = sp->qmin >= sp->qmax;
               }
               if (sp->noqr) {
                  printf("warning:  no scalar range\n");
               } else {
                  printf("Scalar range %E to %E\n", sp->qmin, sp->qmax);
                  if (((vis_gvp->nd+vis_gvp->ndo)==4) &&
                    vis_gvp->curv) {
#ifdef DEBUG
                     printf("vis_draw:  calling 2d_curv_colour_map\n");
#endif
#ifdef CURV98CONT
                     if (!nocolour) {
#endif
                        vis_2d_curv_colour_map(vis_gvp->nc, vis_gvp->xe,
                          vis_gvp->sa, vis_gvp->qa, sp->qmin, sp->qmax,
                          ios);
#ifdef CURV98CONT
                     }
                     if (!nocontour) {
                        if (noqlabel) {
                           ccopy(vis_gvp->qlabel, 0, SZLIN, temp_qlabel);
                        } else {
                           ccopy(qlabel, 0, SZLIN, temp_qlabel);
                        }
                        vis_2d_curv_contour_map(vis_gvp->nc[0],vis_gvp->nc[1],
                          vis_gvp->xa[0], vis_gvp->xa[1], vis_gvp->sa,
                          vis_gvp->qa, vis_gvp->qmin, vis_gvp->qmax, 
                          temp_qlabel, vis_gvp->qunits, ios);
                     }
#endif
                  }
               }
#endif
            }
/*          Vector field in 3D */
            if ((!vis_gvp->novector) && (vis_gvp->nd == 3)) {
               if (sp->novel) {
                  if (novmax) {
                     sp->velmax = 0.0;
                     for (id=0; id < vis_gvp->ndv; id++) {
                        if (vis_gvp->vmax[id] > sp->velmax) {
                           sp->velmax = vis_gvp->vmax[id];
                        }
                     }
                  } else {
                     sp->velmax = vmax;
                  }
                  sp->novel = (sp->velmax <= 0.0);
               }
               if (sp->novel) {
                  printf("warning:  vector range zero\n");
               } else {
                  if (!noarrow) {
                     printf("VectorMax %E\n", sp->velmax);
                     vis_3d_grid_vector(vis_gvp->nc, vis_gvp->xe,
                       vis_gvp->sa, vis_gvp->va, sp->velmax, ios);
                  }
                  if (!nostream) {
                     vis_streamer(vis_gvp->nd, vis_gvp->nc, vis_gvp->xe,
                       vis_gvp->sa, vis_gvp->ndv, vis_gvp->va,
                       vis_gvp->vmax, ios);
                  }
               }
            }
         }
/*       Scatter model (svp) */
         if (!noscatter) {
#ifdef CURV_OVER
            sp->curv = vis_svp->curv && (vis_svp->ndc >= 2);
#endif
/*          Scalar points in 3D */
            if ((!vis_svp->noscalar) && (vis_svp->nd == 3)) {
               if (sp->noqr) {
                  if (noqr) {
                     sp->qmin = vis_svp->qmin;
                     sp->qmax = vis_svp->qmax;
                  } else {
                     sp->qmin = qmin;
                     sp->qmax = qmax;
                  }
                  sp->noqr = sp->qmin >= sp->qmax;
               }
               if (sp->noqr) {
                  printf("warning:  no scalar range\n");
               } else {
                  printf("Scalar range %E to %E\n", sp->qmin, sp->qmax);
                  vis_3d_scatter_scalar(vis_svp->npt,vis_svp->nd,
                    vis_svp->xpt,vis_svp->spt,vis_svp->qpt,
                    sp->qmin,sp->qmax,ios);
               }
            }
/*          Vector points in 3D */
            if ((!vis_svp->novector) && (vis_svp->nd == 3)) {
               if (sp->novel) {
                  if (novmax) {
                     sp->velmax = 0.0;
                     for (id=0; id < vis_svp->ndv; id++) {
                        if (vis_svp->vmax[id] > sp->velmax) {
                           sp->velmax = vis_svp->vmax[id];
                        }
                     }
                  } else {
                     sp->velmax = vmax;
                  }
                  sp->novel = (sp->velmax <= 0.0);
               }
               if (sp->novel) {
                  printf("warning:  vector range zero\n");
               } else {
                  printf("VectorMax %E\n", sp->velmax);
                  vis_3d_scatter_vector(ios);
               }
            }
         }
      }
/* unconventional number of dimensions (weird) */
   } else {
      printf("error:  not 2 or 3-dimensional (vis_draw)\n");
      *ios = ERR;
   }
}

void vis_line_graph(n, xa, sa, qa, ios)
   int n,*ios,sa[];
   float xa[],qa[];
{
   bool on,bad;
   int i,istyle,colour_index;
   float r,g,b;
   float x,q,h,h1,h2,v;

/* grab current colour settings */
   if (wp->colour_avail) {
      drp_get_colour(&colour_index);
   } else if (wp->rgb_avail) {
      drp_get_rgb(&r, &g, &b);
   }
/* produce key box? */
   if (!sp->nokey && !sp->key_drawn) {
      vis_key();
   }
/* line style select */
   if (line_style > 1) {
      drp_get_line_style(&istyle);
      drp_set_line_style(line_style);
   }
/* key entry? */
   if (!sp->nokey && !nokey_label && (sp->key_drawn)) {
      if (sp->key_count >= sp->key_size) {
         printf("warning:  no room for further key entries\n");
      } else {
         drp_set_hv(sp->hkmin, sp->hkmax, sp->vkmin, sp->vkmax);
         h1 = sp->hkmin + 0.5*wp->char_width;
         h2 = h1 + 2.0*wp->char_width;
         v = sp->vkmax - wp->char_height*(1.5*sp->key_count + 0.75);
/*      Select palet for key */
         vis_select_palet();
         if (!noline_connect) {
            drp_bp();
            drp_ap(h1, v);
            drp_ap(h2, v);
            drp_ep();
         }
         if (!noline_marker) {
            drp_marker(line_marker, 0.5*(h1+h2), v);
         }
/*       Restore default line colour */
         if (wp->colour_avail)  {
            drp_set_colour(colour_index);
         } else if (wp->rgb_avail) {
            drp_set_rgb(r, g, b);
         }
         h = h2 + 0.5*wp->char_width;
         drp_text(h, v, key_label, "lm");
         sp->key_count = sp->key_count + 1;
      }
   }
/* set data clipping viewport */
   drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
/* select palet for data */
   vis_select_palet();
/*draw as polyline with optionbal markers */
   on = FALSE;
   for (i=0; i < n; i++) {
      if (sa[i] == FLUID) {
         x = xa[i];
         q = qa[i];
         bad = FALSE;
         if (bad) {
            if (on) {
               vis_ep();
               on = FALSE;
            }
         } else if ((x < sp->gcmin[0]) || (x > sp->gcmax[0])) {
            if (on) {
               vis_ep();
               on = FALSE;
            }
         } else if ((q < sp->gcmin[1]) || (q > sp->gcmax[1])) {
            if (on) {
               vis_ep();
               on = FALSE;
            }
         } else {
            if (on) {
               vis_ap(x, q);
            } else if (!noline_connect) {
               vis_bp();
               vis_ap(x, q);
               on = TRUE;
            }
         }
      } else if (on) {
         vis_ep();
         on = FALSE;
      }
   }
   if (on) {
      vis_ep();
   }
   if (!noline_marker) {
      for (i=0; i < n; i++) {
         if (sa[i] == FLUID) {
            x = xa[i];
            q = qa[i];
            vis_marker(line_marker, x, q);
         }
      }
   }
   if (line_style > 1) {
      drp_set_line_style(istyle);
   }
/* revert to graph clipping viewport */
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
/* restore current colour settings */
   if (wp->colour_avail)  {
      drp_set_colour(colour_index);
   } else if (wp->rgb_avail) {
      drp_set_rgb(r, g, b);
   }
}

void vis_contour_map(nx, ny, xa, ya, sa, qa, qamin, qamax, qlabel, qunits,
  ios)
   int nx,ny,*ios,sa[];
   float xa[],ya[],qa[],qamin,qamax;
   char qlabel[],qunits[];
{
   int i;
   float dqdi,dvdq,h,v,vtop;
   char str[SZLIN];

/* design contours */
   if (sp->nocont) {
/*   establish scalar range if none */
      if (sp->noqr) {
         if (noqr) {
            sp->qmin = qamin;
            sp->qmax = qamax;
         } else {
            sp->qmin = qmin;
            sp->qmax = qmax;
         }
         sp->noqr = sp->qmin >= sp->qmax;
      }
/*   number of contours */
      sp->ncontour = ncontour;
/*   user-supplied contour limits */
      if (!nocontour_levels) {
         for (i=0; i < ncontour; i++) {
            sp->contour_levels[i] = contour_levels[i];
         }
         sp->nocont = FALSE;
      } else if (sp->noqr) {
         printf("error:  no scalar range to determine contours\n");
         *ios = ERR;
      } else {
         if (ncontour == 1) {
            sp->contour_levels[0] = 0.5*(sp->qmin + sp->qmax);
         } else {
            dqdi = (sp->qmax - sp->qmin)/(float)(ncontour+1);
            for (i=0; i < ncontour; i++) {
               sp->contour_levels[i] = (float)(i+1)*dqdi + sp->qmin;
            }
         }
         sp->nocont = FALSE;
      }
   }
/* continue if we really have some contour levels */
   if (sp->nocont) {
      printf("error:  no contour levels\n");
      *ios = ERR;
   } else {
/*    Draw scales if not done already */
      if (!sp->noscales) {
         if (!sp->scales_drawn) {
            ccopy(qlabel, 0, SZNAM, sp->qlabel);
            ccopy(qunits, 0, SZNAM, sp->qunits);
            vis_scales();
         }
/*      draw contour levels? */
         if (!sp->contour_scales_drawn) {
            drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
            dvdq = (sp->vsmax-sp->vsmin)/(sp->qmax-sp->qmin);
            vtop = sp->vmin;
            for (i=0; i < sp->ncontour; i++) {
               if (sp->contour_levels[i] < sp->qmin) {
               } else if (sp->contour_levels[i] > sp->qmax) {
               } else {
                  v = (sp->contour_levels[i]-sp->qmin)*dvdq + sp->vsmin;
/*               check for crowding */
                  if (v-wp->char_height > vtop) {
/*                  value on right */
                     h = sp->hsmax + 0.5*wp->char_width;
                     cfloat(sp->contour_levels[i]/sp->qscale,sp->qfield,
                       sp->qdp, SZLIN, str);
                     drp_text(h, v, str, "lm");
/*                  index on left */
                     sprintf(str, "%d", i);
                     h = sp->hsmin - 0.5*wp->char_width;
                     drp_text(h, v, str, "rm");
                     vtop = v + 0.5*wp->char_height;
                  }
/*               contour level as a line */
                  drp_bp();
                  drp_ap(sp->hsmin, v);
                  drp_ap(sp->hsmax, v);
                  drp_ep();
               }
            }
            sp->contour_scales_drawn = TRUE;
         }
      }
/*   draw contours */
      for (i=0; i < sp->ncontour; i++) {
         if (nocontour_label) {
            sprintf(str, "%c", '\0');
         } else {
            sprintf(str, "%d", i);
         }
         drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
         contour_map(nx, ny, xa, ya, sa, qa, 
           sp->contour_levels[i], str, ios);
         drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
      }
   }
}

void vis_colour_map(nx, ny, xe, ye, sa, qa, qamin, qamax, qlabel, qunits, 
  ios)
   int nx,ny,*ios,sa[];
   float xe[],ye[],qa[],qamin,qamax;
   char qlabel[],qunits[];
{
   int colour_index,i,k,ix,iy,offset,ndi,nci[MAXD];
   float didq,q0,v1,v2,dvdi,phi,dphi,r,g,b,red,green,blue,dpdq;
   dpdq = 0.0;

/* design colour ranges */
   if (sp->nocmap) {
/*   establish scalar range if none */
      if (sp->noqr) {
         if (noqr) {
            sp->qmin = qamin;
            sp->qmax = qamax;
         } else {
            sp->qmin = qmin;
            sp->qmax = qmax;
         }
         sp->noqr = sp->qmin >= sp->qmax;
      }
/*   colour ranges */
      if (sp->noqr) {
         printf("error:  no scalar range to determine contours\n");
         *ios = ERR;
      } else {
         if (wp->colour_avail) {
            vis_design_cmap(ncolour, &offset);
            if (offset < 0) {
               printf("error:  no scalar colour ramp possible\n");
               *ios = ERR;
             }
         } else if (wp->rgb_avail) {
            *ios = OK;
         } else {
            printf("error:  no colour??? \n");
            *ios = ERR;
         }
         sp->nocmap = *ios;
      }
/* pick up existing colour ramp */
   } else if (wp->colour_avail) {
      vis_design_cmap(ncolour, &offset);
      if (offset < 0) {
         printf("error:  scalar colour ramp lost, not recovered\n");
         *ios = ERR;
          sp->nocmap = TRUE;
      }
   } else if (wp->rgb_avail) {
      *ios = OK;
   }
/* continue only is we have a colour map */
   if (sp->nocmap) {
      printf("error:  no colour map defined\n");
      *ios = ERR;
   } else {
/*    Draw scales if not done already */
      if (wp->colour_avail) {
         drp_get_ramp("S", &ndi, nci, &offset);
         sp->ncolour = nci[0];
      } else if (wp->rgb_avail) {
         sp->ncolour = 64;
      }
      if (!sp->noscales) {
         if (!sp->scales_drawn) {
            ccopy(qlabel, 0, SZNAM, sp->qlabel);
            ccopy(qunits, 0, SZNAM, sp->qunits);
            vis_scales();
         }
/*      draw colours bits? */
         if (!sp->cmap_scales_drawn) {
            if (wp->colour_avail) {
               drp_get_colour(&colour_index);
            } else if (wp->rgb_avail) {
               drp_get_rgb(&r, &g, &b);
            }
            dvdi = (sp->vsmax - sp->vsmin)/(float)(sp->ncolour);
            drp_set_hv(sp->hsmin, sp->hsmax, sp->vsmin, sp->vsmax);
            phi = 0.0;
            dphi = 1.0/(float)(sp->ncolour - 1);
            for (i=0; i < sp->ncolour; i++) {
               v1 = (float)i*dvdi + sp->vsmin;
               v2 = v1 + dvdi;
               if (wp->colour_avail) {
                  drp_set_colour(offset + i);
               } else if (wp->rgb_avail) {
                  vis_rain_dance(phi, &red, &green, &blue);
                  phi = phi + dphi;
                  drp_set_rgb(red, green, blue);
               }
               drp_rect(sp->hsmin, sp->hsmax, v1, v2);
            }
            if (wp->colour_avail)  {
               drp_set_colour(colour_index);
            } else if (wp->rgb_avail) {
               drp_set_rgb(r, g, b);
            }
/*         outer box for scale */
            drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
/*
            vis_box(sp->hsmin, sp->hsmax, sp->vsmin, sp->vsmax);
*/
            sp->cmap_scales_drawn = TRUE;
         }
      }
/*   draw image */
      drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
      didq = (float)(sp->ncolour)/(sp->qmax - sp->qmin);
      q0 = sp->qmin;
      if (wp->colour_avail) {
         drp_get_colour(&colour_index);
      } else if (wp->rgb_avail) {
         drp_get_rgb(&r, &g, &b);
         dpdq = 1.0/(sp->qmax - sp->qmin);
      }
      k = 0;
      for (iy=0; iy < ny; iy++) {
         for (ix=0; ix < nx; ix++) {
/*
            if (xe[ix] < sp->wcmin[0]) {
            } else if (xe[ix] > sp->wcmax[0]) {
            } else if (xe[ix+1] < sp->wcmin[0]) {
            } else if (xe[ix+1] > sp->wcmax[0]) {
            } else if (ye[iy] < sp->wcmin[1]) {
            } else if (ye[iy] > sp->wcmax[1]) {
            } else if (ye[iy+1] < sp->wcmin[1]) {
            } else if (ye[iy+1] > sp->wcmax[1]) {
            } else {
*/
               if (sa[k] == FLUID) {
                  if (wp->colour_avail) {
                     i = (qa[k]-q0)*didq;
                     if (i < 0) {
                        i = 0;
                     } else if (i >= sp->ncolour) {
                        i = sp->ncolour-1;
                     }
#ifdef CURV
                     drp_set_colour(offset+i);
                     vis_fill_rect(xe[ix],xe[ix+1],ye[iy],ye[iy+1]);
#else
                     vis_cmap_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1],
                       offset+i);
#endif
                  } else if (wp->rgb_avail) {
                     phi = (qa[k]-q0)*dpdq;
                     if (phi < 0.0) {
                        phi = 0.0;
                     } else if (phi > 1.0) {
                        phi = 1.0;
                     }
                     vis_rain_dance(phi, &red, &green, &blue);
                     drp_set_rgb(red, green, blue);
#ifdef CURV
                     vis_fill_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#else
                     vis_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#endif
                  }
               } else if (sa[k] == SOLID) {
#ifdef CURV
                  if (wp->colour_avail) {
                     drp_set_colour(wp->colour_black);
                  } else if (wp->rgb_avail) {
                     drp_set_rgb(0.0, 0.0, 0.0);
                  }
                  vis_fill_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#else
                  vis_black_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#endif
               }
/*
            }
*/
            k = k + 1;
         }
      }
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
      if (wp->colour_avail) {
         drp_set_colour(colour_index);
      } else if (wp->rgb_avail) {
         drp_set_rgb(r, g, b);
      }
   }
}


void vis_band_map(nx, ny, xe, ye, xa, ya, sa, qa, qamin, qamax, qlabel, 
  qunits, ios)
   int nx,ny,*ios,sa[];
   float xe[],ye[],xa[],ya[],qa[],qamin,qamax;
   char qlabel[],qunits[];
{
   int colour_index,i,k,ix,iy,offset,ndi,nci[MAXD];
   float v1,v2,dvdq;
   char str[SZLIN];

/* design colour ranges */
   if (sp->nocmap) {
/*   establish scalar range if none */
      if (sp->noqr) {
         if (noqr) {
            sp->qmin = qamin;
            sp->qmax = qamax;
         } else {
            sp->qmin = qmin;
            sp->qmax = qmax;
         }
         sp->noqr = sp->qmin >= sp->qmax;
      }
/*    setup scalar colour ramp */
      vis_design_bmap(&offset);
/*   colour ranges */
      if (sp->noqr) {
         printf("error:  no scalar range to determine contours\n");
         *ios = ERR;
      } else if (offset < 0) {
         printf("error:  unable to create scalr range ramp\n");
         *ios = ERR;
      } else {
         sp->nocmap = FALSE;
      }
   } else {
      vis_design_bmap(&offset);
      if (offset < 0) {
         printf("error:  unable to create scalar range ramp\n");
         *ios = ERR;
      }
   }
/* continue only is we have a colour map */
   if (sp->nocmap) {
      printf("error:  no colour map defined\n");
      *ios = ERR;
   } else {
/*    Draw scales if not done already */
      drp_get_ramp("S", &ndi, nci, &offset);
      sp->ncolour = nci[0];
      if (!sp->noscales) {
         if (!sp->scales_drawn) {
            ccopy(qlabel, 0, SZNAM, sp->qlabel);
            ccopy(qunits, 0, SZNAM, sp->qunits);
            vis_scales();
         }
/*      draw colours bits? */
         if (!sp->cmap_scales_drawn) {
            drp_get_colour(&colour_index);
            dvdq = (sp->vsmax - sp->vsmin)/(sp->qmax - sp->qmin);
            for (i=0; i < MAXBAND; i++) {
               if (!nobandr[i]) {
                  v1 = (bandmin[i]-sp->qmin)*dvdq + sp->vsmin;
                  v2 = (bandmax[i]-sp->qmin)*dvdq + sp->vsmin;
                  drp_set_colour(offset + i);
                  drp_rect(sp->hsmin, sp->hsmax, v1, v2);
               }
            }
            drp_set_colour(colour_index);
/*         outer box for scale */
            vis_box(sp->hsmin, sp->hsmax, sp->vsmin, sp->vsmax);
            sp->cmap_scales_drawn = TRUE;
         }
      }
/*   draw image */
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
      drp_get_colour(&colour_index);
      k = 0;
      for (iy=0; iy < ny; iy++) {
         for (ix=0; ix < nx; ix++) {
/*
            if (xe[ix] < sp->wcmin[0]) {
            } else if (xe[ix] > sp->wcmax[0]) {
            } else if (xe[ix+1] < sp->wcmin[0]) {
            } else if (xe[ix+1] > sp->wcmax[0]) {
            } else if (ye[iy] < sp->wcmin[1]) {
            } else if (ye[iy] > sp->wcmax[1]) {
            } else if (ye[iy+1] < sp->wcmin[1]) {
            } else if (ye[iy+1] > sp->wcmax[1]) {
            } else {
*/
               if (sa[k] == FLUID) {
                  for (i=0; i < MAXBAND; i++) {
                     if (nobandr[i]) {
                     } else if (qa[k] < bandmin[i]) {
                     } else if (qa[k] > bandmax[i]) {
                     } else {
#ifdef CURV
                        drp_set_colour(offset+i);
                        vis_fill_rect(xe[ix],xe[ix+1],ye[iy],ye[iy+1]);
#else
                        vis_cmap_rect(xe[ix],xe[ix+1],ye[iy],ye[iy+1],
                          offset+i);
#endif
                     }
                  }
               } else if (sa[k] == SOLID) {
#ifdef CURV
                  drp_set_colour(wp->colour_black);
                  vis_fill_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#else
                  vis_black_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#endif
               }
/*
            }
*/
            k = k + 1;
         }
      }
      for (i=0; i < MAXBAND; i++) {
         if (!nobandr[i]) {
            str[0] = '\0';
            drp_set_colour(offset+i);
            contour_map(nx, ny, xa, ya, sa, qa, 
              bandmin[i], str, ios);
            contour_map(nx, ny, xa, ya, sa, qa, 
              bandmax[i], str, ios);
         }
      }
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
      drp_set_colour(colour_index);
   }
}

/* vis_design_cmap - design (scalar) colourmap */
void vis_design_cmap(nc, offset)
   int nc,*offset;
{
   int ndi,nci[MAXD],n;

/* see if there is one already there */
   drp_get_ramp("S", &ndi, nci, offset);
/* if not try to create one */
   if (*offset < 0) {
      drp_get_colour_space(&n);
/*    only viable if some solour space available */
      if (n > 4) {
         if (n > nc) {
            n = nc;
         }
         ndi = 1;
         nci[0] = n;
         drp_put_ramp("S", ndi, nci, offset);
         if (*offset < 0) {
            printf("error:  trouble creating scalar colour ramp\n");
         } else {
            vis_rain(n, *offset);
         }
      }
/* check that it obeys the convention */
   } else if (ndi != 1) {
      printf("error:  colour ramp S breaks convention (report)\n");
   }
}

/* vis_design_bmap - design (scalar-band) colourmap */
void vis_design_bmap(offset)
   int *offset;
{
   int ndi,nci[MAXD],n;

/* see if there is one already there */
   drp_get_ramp("B", &ndi, nci, offset);
/* if not try to create one */
   if (*offset < 0) {
      drp_get_colour_space(&n);
/*    only viable if some solour space available */
      if (n >= MAXBAND) {
         ndi = 1;
         nci[0] = MAXBAND;
         drp_put_ramp("B", ndi, nci, offset);
         if (*offset < 0) {
            printf("error:  trouble creating scalar colour band\n");
         } else {
            vis_rain(MAXBAND, *offset);
         }
      }
/* check that it obeys the convention */
   } else if ((ndi != 1) || (nci[0] != MAXBAND)) {
      printf("error:  colour ramp B breaks convention (report)\n");
      *offset = -1;
   }
}

#ifdef HSV_RAIN
void vis_h_to_rgb(float hue, float rgb[3]) {
   float h,f;
   int i;

   h = hue*6.0;
   i = (int)h;
   f = h - (float)i;
   switch (i) {
      case 0:
         rgb[0] = 1.0;
         rgb[1] = f;
         rgb[2] = 0.0;
         break;
      case 1:
         rgb[0] = 1.0-f;
         rgb[1] = 1.0;
         rgb[2] = 0.0;
         break;
      case 2:
         rgb[0] = 0.0;
         rgb[1] = 1.0;
         rgb[2] = f;
         break;
      case 3:
         rgb[0] = 0.0;
         rgb[1] = 1.0-f;
         rgb[2] = 1.0;
         break;
      case 4:
         rgb[0] = f;
         rgb[1] = 0.0;
         rgb[2] = 1.0;
         break;
      default:
         rgb[0] = 1.0;
         rgb[1] = 0.0;
         rgb[2] = 1.0 - f;
   }
}
#endif

void vis_rain(nc, offset)
   int offset,nc;
{
   int i;
   float dzdi,z,red[MAXCM],green[MAXCM],blue[MAXCM];
#ifdef HSV_RAIN
   float h,rgb[3];
#endif

   dzdi = 1.0/(float)(nc-1);
   for (i=0; i < nc; i++) {
      z = (float)i*dzdi;
#ifndef NOGREY
      if (vis_grey_rain) {
         red[i] = 1.0-z;
         green[i] = 1.0-z;
         blue[i] = 1.0-z;
      } else {
#endif
#ifdef HSV_RAIN
      h = 0.66666*(1.0-z);
      vis_h_to_rgb(h, rgb);
      red[i] = rgb[0];
      green[i] = rgb[1];
      blue[i] = rgb[2];
#else
      red[i] = z;
      if (z < 0.5) {
         green[i] = 2.0*z;
      } else {
         green[i] = 2.0 - 2.0*z;
      }
      blue[i] = 1.0-z;
#endif
#ifndef NOGREY
      }
#endif
   }
   drp_set_colours(offset, nc, red, green, blue);
}

void vis_arrow_map(nx, ny, xa, ya, sa, ndv, va, velmax, ios)
   int nx,ny,ndv,*ios,sa[];
   float xa[],ya[],va[],velmax[];
{
   int ix,iy,k,id,kv;
   float x,y,vx,vy,dldv,lmin,l,v,ax,ay;

/* determine velocity range and arrow scaling */
   if (sp->novel) {
      if (novmax) {
         sp->velmax = 0.0;
         for (id=0; id < ndv; id++) {
            if (velmax[id] > sp->velmax) {
               sp->velmax = velmax[id];
            }
         }
      } else {
         sp->velmax = vmax;
      }
      sp->novel = (sp->velmax == 0.0);
      if (sp->novel) {
         printf("warning:  vector range zero\n");
      } else {
         printf("VectorMax %E\n", sp->velmax);
         lmin = BIG;
         for (id=0; id < 2; id++) {
            l = sp->wcmax[id] - sp->wcmin[id];
            if (l < lmin) {
               lmin = l;
            }
         }
         if (nx > ny) {
            lmin = lmin/(float)nx;
         } else {
            lmin = lmin/(float)nx;
         }
         sp->dldv = 0.5*arrow_rate*lmin/sp->velmax;
      }
   }
   if (!sp->novel) {
      drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
      dldv = sp->dldv;
      k = 0;
      kv = 0;
      for (iy=0; iy < ny; iy++) {
         y = ya[iy];
         for (ix=0; ix < nx; ix++) {
            if (sa[k] != FLUID) {
/*
            } else if (xa[ix] < sp->wcmin[0]) {
            } else if (xa[ix] > sp->wcmax[0]) {
            } else if (ya[iy] < sp->wcmin[1]) {
            } else if (ya[iy] > sp->wcmax[1]) {
*/
            } else if (ix != ((ix/arrow_rate)*arrow_rate)) {
            } else if (iy != ((iy/arrow_rate)*arrow_rate)) {
            } else {
               x = xa[ix];
               vx = va[kv];
               y = ya[iy];
               vy = va[kv+1];
               v = sqrt(vx*vx + vy*vy);
               if (v > 0.02*sp->velmax) {
                  ax = vx/v;
                  ay = vy/v;
                  if (v > sp->velmax) {
                     v = sp->velmax;
                  }
                  l = dldv*v*arrow_scale;
                  vis_arrow(x, y, l, ax, ay);
               }
            }
            kv = kv + ndv;
            k = k + 1;
         }
      }
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
   }
}

void vis_streamer(nd, nc, xe, sa, ndv, va, vamax, ios)
   int nd,nc[],sa[],ndv,*ios;
   float *xe[],va[],vamax[];
{
   int id,i;
   float x[MAXD],vamin;
   bool zap;

   vamin = 0.0;
   for (id=0; id < nd; id++) {
      vamin = vamin + vamax[id]*vamax[id];
   }
   vamin = 0.01*sqrt(vamin);
   if (sp->nd == 2) {
      drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
      *ios = OK;
   } else {
      dr3_coerce_wire(ios);
      if (*ios) {
         printf("error:  unable to obtain image\n");
      }
   }
   if (!*ios) {
      for (i=0; i < stream_size; i++) {
         zap = FALSE;
         for (id=0; id < nd; id++) {
            x[id] = (float)i*stream_dx[id] + stream_x[id];
            zap = zap || (x[id] < sp->wcmin[id]) ||
              (x[id] > sp->wcmax[id]);
         }
         if (!zap) {
            vis_track(nd, nc, xe, sa, ndv, va, vamin, x);
         }
      }
   }
   if (sp->nd == 2) {
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
   }
}

void vis_track(nd, nc, xe, sa, ndv, va, vamin, x)
   int nd,nc[],sa[],ndv;
   float *xe[],va[],vamin,x[];
{
   float dt,dxbar,vbar,divdx[MAXD],x0[MAXD],cn,dx[MAXD];
   int id,mstep,istep,j,n,nm[MAXD],con[MAXD],k;
   bool out;

   cn = 0.3;
   mstep = 0;
   dxbar = 0.0;
   for (id=0; id < nd; id++) {
      nm[id] = nc[id] - 1;
      n = nc[id];
      x0[id] = xe[id][0];
      divdx[id] = (float)(nc[id])/(xe[id][n] - x0[id]);
      dx[id] = 1.0/divdx[id];
      dxbar = dxbar + dx[id]*dx[id];
      mstep = mstep + nc[id];
   }
   mstep = mstep*50;
   dxbar = sqrt(dxbar);
   if (nd == 3) {
      dr3_wire_icon(0.1*dxbar, x[0], x[1], x[2]);
   } else if (nd == 2) {
#ifdef CURV
      vis_draw_rect(x[0]-0.1*dxbar, x[0]+0.1*dxbar, x[1]-0.1*dxbar,
        x[1]+0.1*dxbar);
#else
      vis_grid_rect(x[0]-0.1*dxbar, x[0]+0.1*dxbar, x[1]-0.1*dxbar,
        x[1]+0.1*dxbar);
#endif
   }
   istep = 0;
   while (istep < mstep) {
      out = FALSE;
      for (id=0; id < nd; id++) {
         j = (x[id] - x0[id])*divdx[id];
         if ((j < 0) || (j > nm[id])) {
            out = TRUE;
            break;
         } else if ((x[id]<sp->wcmin[id]) || (x[id]>sp->wcmax[id])) {
            out = TRUE;
            break;
         } else {
            con[id] = j;
         }
      }
      if (out) {
         break;
      } else {
         j = 0;
         for (id=nd-1; id >=0 ; id--) {
            j = j*nc[id] + con[id];
         }
         if (sa[j] == FLUID) {
            vbar = 0.0;
            k = j*ndv;
            for (id=0; id < nd; id++) {
               vbar = vbar + fabs(va[k++]);
            }
            if (vbar > vamin) {
               if (istep == 0) {
                  if (nd == 2) {
                     vis_marker(line_marker, x[0], x[1]);
                     vis_bp();
                     vis_ap(x[0], x[1]);
                  } else {
                     dr3_begin_polyline();
                     dr3_polyline_vertex(x[0], x[1], x[2]);
                  }
               }
               dt = cn*dxbar/vbar;
               k = j*ndv;
               for (id=0; id < nd; id++) {
                  x[id] = x[id] + dt*va[k++];
               }
               if (nd == 2) {
                  vis_ap(x[0], x[1]);
               } else {
                  dr3_polyline_vertex(x[0], x[1], x[2]);
               }
            } else {
               break;
            }
         } else {
            break;
         }
      }
      istep++;
   }
   if (istep > 0) {
      if (nd == 2) {
         vis_ep();
      } else {
         dr3_end_polyline();
      }
   }
}

void vis_colour_gram(npt, nd, xpt, spt, qpt, qamin, qamax, qlabel, qunits,
  ios)
   int npt,nd,*ios,spt[];
   float xpt[],qpt[],qamin,qamax;
   char qlabel[],qunits[];
{
   int colour_index,i,k,offset,ndi,nci[MAXD],id,j;
   float didq,q0,v1,v2,dvdi,phi,dphi,r,g,b,red,green,blue,dpdq,l,lmin,
     x,y;
   dpdq = 0.0;
/* design colour ranges */
   if (sp->nocmap) {
/*   establish scalar range if none */
      if (sp->noqr) {
         if (noqr) {
            sp->qmin = qamin;
            sp->qmax = qamax;
         } else {
            sp->qmin = qmin;
            sp->qmax = qmax;
         }
         sp->noqr = sp->qmin >= sp->qmax;
      }
/*   colour ranges */
      if (sp->noqr) {
         printf("error:  no scalar range to determine contours\n");
         *ios = ERR;
      } else {
         if (wp->colour_avail) {
            vis_design_cmap(ncolour, &offset);
            if (offset < 0) {
               printf("error:  no scalar colour ramp possible\n");
               *ios = ERR;
             }
         } else if (wp->rgb_avail) {
            *ios = OK;
         } else {
            printf("error:  no colour??? \n");
            *ios = ERR;
         }
         sp->nocmap = *ios;
      }
/* pick up existing colour ramp */
   } else if (wp->colour_avail) {
      vis_design_cmap(ncolour, &offset);
      if (offset < 0) {
         printf("error:  scalar colour ramp lost, not recovered\n");
         *ios = ERR;
          sp->nocmap = TRUE;
      }
   } else if (wp->rgb_avail) {
      *ios = OK;
   }
/* continue only is we have a colour map */
   if (sp->nocmap) {
      printf("error:  no colour map defined\n");
      *ios = ERR;
   } else {
/*    Draw scales if not done already */
      if (wp->colour_avail) {
         drp_get_ramp("S", &ndi, nci, &offset);
         sp->ncolour = nci[0];
      } else if (wp->rgb_avail) {
         sp->ncolour = 64;
      }
      if (!sp->noscales) {
         if (!sp->scales_drawn) {
            ccopy(qlabel, 0, SZNAM, sp->qlabel);
            ccopy(qunits, 0, SZNAM, sp->qunits);
            vis_scales();
         }
/*      draw colours bits? */
         if (!sp->cmap_scales_drawn) {
            if (wp->colour_avail) {
               drp_get_colour(&colour_index);
            } else if (wp->rgb_avail) {
               drp_get_rgb(&r, &g, &b);
            }
            dvdi = (sp->vsmax - sp->vsmin)/(float)(sp->ncolour);
            drp_set_hv(sp->hsmin, sp->hsmax, sp->vsmin, sp->vsmax);
            phi = 0.0;
            dphi = 1.0/(float)(sp->ncolour - 1);
            for (i=0; i < sp->ncolour; i++) {
               v1 = (float)i*dvdi + sp->vsmin;
               v2 = v1 + dvdi;
               if (wp->colour_avail) {
                  drp_set_colour(offset + i);
               } else if (wp->rgb_avail) {
                  vis_rain_dance(phi, &red, &green, &blue);
                  phi = phi + dphi;
                  drp_set_rgb(red, green, blue);
               }
               drp_rect(sp->hsmin, sp->hsmax, v1, v2);
            }
            if (wp->colour_avail)  {
               drp_set_colour(colour_index);
            } else if (wp->rgb_avail) {
               drp_set_rgb(r, g, b);
            }
/*         outer box for scale */
            drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
/*
            vis_box(sp->hsmin, sp->hsmax, sp->vsmin, sp->vsmax);
*/
            sp->cmap_scales_drawn = TRUE;
         }
      }
/*   draw image */
      drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
      didq = (float)(sp->ncolour)/(sp->qmax - sp->qmin);
      q0 = sp->qmin;
      if (wp->colour_avail) {
         drp_get_colour(&colour_index);
      } else if (wp->rgb_avail) {
         drp_get_rgb(&r, &g, &b);
         dpdq = 1.0/(sp->qmax - sp->qmin);
      }
      lmin = BIG;
      for (id=0; id < 2; id++) {
         l = sp->wcmax[id] - sp->wcmin[id];
         if (l < lmin) {
            lmin = l;
         }
      }
      lmin = arrow_ratio*lmin;
      k = 0;
      for (j=0; j < npt; j++) {
         x = xpt[k];
         y = xpt[k+1];
         k = k + nd;
         if (x > sp->wcmax[0]) {
         } else if (x > sp->wcmax[0]) {
         } else if (y < sp->wcmin[1]) {
         } else if (y > sp->wcmax[1]) {
         } else if (qpt[j] < qamin) {
         } else if (qpt[j] > qamax) {
         } else {
            if (spt[j] == FLUID) {
               if (wp->colour_avail) {
                  i = (qpt[j]-q0)*didq;
                  if (i < 0) {
                     i = 0;
                  } else if (i >= sp->ncolour) {
                     i = sp->ncolour-1;
                  }
#ifdef CURV
                  drp_set_colour(offset+i);
                  vis_fill_rect(x-lmin,x+lmin,y-lmin,y+lmin);
                  drp_set_colour(colour_index);
                  vis_draw_rect(x-lmin,x+lmin,y-lmin,y+lmin);
#else
                  vis_cmap_rect(x-lmin,x+lmin,y-lmin,y+lmin,offset+i);
                  drp_set_colour(colour_index);
                  vis_grid_rect(x-lmin, x+lmin, y-lmin, y+lmin);
#endif
               } else if (wp->rgb_avail) {
                  phi = (qpt[j]-q0)*dpdq;
                  vis_rain_dance(phi, &red, &green, &blue);
                  drp_set_rgb(red, green, blue);
#ifdef CURV
                  vis_fill_rect(x-lmin, x+lmin, y-lmin, y+lmin);
#else
                  vis_rect(x-lmin, x+lmin, y-lmin, y+lmin);
#endif
                  drp_set_rgb(r, g, b);
#ifdef CURV
                  vis_draw_rect(x-lmin, x+lmin, y-lmin,y+lmin);
#else
                  vis_grid_rect(x-lmin, x+lmin, y-lmin,y+lmin);
#endif
               }
/*
            } else {
               vis_black_rect(x-lmin, x+lmin, y-lmin, y+lmin);
*/
            }
         }
      }
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
      if (wp->colour_avail) {
         drp_set_colour(colour_index);
      } else if (wp->rgb_avail) {
         drp_set_rgb(r, g, b);
      }
   }
}

void vis_marker_gram(npt, nd, xpt, spt, qpt, qamin, qamax, ios)
   int npt,nd,*ios,spt[];
   float xpt[],qpt[],qamin,qamax;
{
   char str[SZLIN];
   int j,k,id;
   float x,y,l,lmin;

   lmin = BIG;
   for (id=0; id < 2; id++) {
      l = sp->wcmax[id] - sp->wcmin[id];
      if (l < lmin) {
         lmin = l;
      }
   }
   lmin = arrow_ratio*lmin;
   drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
   k = 0;
   for (j=0; j < npt; j++) {
      x = xpt[k];
      y = xpt[k+1];
      k = k + nd;
      if (x > sp->wcmax[0]) {
      } else if (x > sp->wcmax[0]) {
      } else if (y < sp->wcmin[1]) {
      } else if (y > sp->wcmax[1]) {
      } else if (qpt[j] < qamin) {
      } else if (qpt[j] > qamax) {
      } else {
         sprintf(str, " %9.2E", qpt[j]);
         vis_marker(line_marker, x, y);
/*
         vis_text(x, y, str, "lm");
*/
      }
   }
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
}


void vis_arrow_gram(npt, nd, xpt, spt, ndv, vpt, velmax, ios)
   int npt,nd,ndv,*ios,spt[];
   float xpt[],vpt[],velmax[];
{
   int j,k,kv,id;
   float x,y,vx,vy,dldv,lmin,l,v,ax,ay;

/* determine velocity range and arrow scaling */
   if (sp->novel) {
      if (novmax) {
         sp->velmax = 0.0;
         for (id=0; id < ndv; id++) {
            if (velmax[id] > sp->velmax) {
               sp->velmax = velmax[id];
            }
         }
      } else {
         sp->velmax = vmax;
      }
      sp->novel = (sp->velmax == 0.0);
      if (sp->novel) {
         printf("warning:  velocity range zero\n");
      } else {
         lmin = BIG;
         for (id=0; id < 2; id++) {
            l = sp->wcmax[id] - sp->wcmin[id];
            if (l < lmin) {
               lmin = l;
            }
         }
         lmin = arrow_ratio*lmin;
         dldv = lmin/sp->velmax;
         sp->dldv = dldv;
      }
   }
   if (!sp->novel) {
      drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
      dldv = sp->dldv;
      k = 0;
      kv = 0;
      for (j=0; j < npt; j++) {
         x = xpt[k];
         y = xpt[k+1];
         k = k + nd;
/*
         if (x > sp->wcmax[0]) {
         } else if (x > sp->wcmax[0]) {
         } else if (y < sp->wcmin[1]) {
         } else if (y > sp->wcmax[1]) {
         } else {
*/
            vx = vpt[kv];
            vy = vpt[kv+1];
            v = sqrt(vx*vx + vy*vy);
            if (v > 0.02*sp->velmax) {
               ax = vx/v;
               ay = vy/v;
               if (v > sp->velmax) {
                  v = sp->velmax;
               }
               l = dldv*v*arrow_scale;
               vis_arrow(x, y, l, ax, ay);
            }
/*
         }
*/
         kv = kv + ndv;
      }
      drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
   }
}

void vis_arrow(x, y, l, ax, ay)
   float x,y,l,ax,ay;
{
   float bx,by,xp,yp,gamma=0.3;
#ifdef CURV
   float u[MAXD],w[MAXD],hh,vv;
#endif

   if (!noarrow_marker) {
      vis_marker(arrow_marker, x, y);
   }
#ifdef CURV
   if (sp->curv) {
      u[0] = x;
      u[1] = y;
      geo_x(2, u, w);
      drp_bp();
      hh = (w[0]-sp->x0)*sp->dhdx + sp->h0;
      vv = (w[1]-sp->y0)*sp->dvdy + sp->v0;
      drp_ap(hh,vv);
      hh = hh + (sp->dhdx)*l*ax;
      vv = vv + (sp->dvdy)*l*ay;
      drp_ap(hh,vv);
      drp_ep();
   } else {
/*    arrow stem */
      vis_bp();
      vis_ap(x, y);
      xp = x + l*ax;
      yp = y + l*ay;
      vis_ap(xp, yp);
/*    anticlockwise arrow fin */
      bx = -ay;
      by = ax;
      xp = x + l*((1.0-gamma)*ax + gamma*bx);
      yp = y + l*((1.0-gamma)*ay + gamma*by);
      vis_ap(xp, yp);
/*    clockwise to right fin */
      xp = x + l*((1.0-gamma)*ax - gamma*bx);
      yp = y + l*((1.0-gamma)*ay - gamma*by);
      vis_ap(xp, yp);
/*    return to tip */
      xp = x + ax*l;
      yp = y + ay*l;
      vis_ap(xp, yp);
      vis_ep();
   }
#else
/* arrow stem */
   vis_bp();
   vis_ap(x, y);
   xp = x + l*ax;
   yp = y + l*ay;
   vis_ap(xp, yp);
/* anticlockwise arrow fin */
   bx = -ay;
   by = ax;
   xp = x + l*((1.0-gamma)*ax + gamma*bx);
   yp = y + l*((1.0-gamma)*ay + gamma*by);
   vis_ap(xp, yp);
/* clockwise to right fin */
   xp = x + l*((1.0-gamma)*ax - gamma*bx);
   yp = y + l*((1.0-gamma)*ay - gamma*by);
   vis_ap(xp, yp);
/* return to tip */
   xp = x + ax*l;
   yp = y + ay*l;
   vis_ap(xp, yp);
   vis_ep();
#endif
}

void vis_bound_map(nx, ny, xe, ye, sa, ios)
   int nx,ny,*ios,sa[];
   float xe[],ye[];
{
   int k,ix,iy,kk;
   bool dowt;

   drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
   k = 0;
   for (iy=0; iy < ny; iy++) {
      for (ix=0; ix < nx; ix++) {
/*
         if (xe[ix] < sp->wcmin[0]) {
         } else if (xe[ix] > sp->wcmax[0]) {
         } else if (xe[ix+1] < sp->wcmin[0]) {
         } else if (xe[ix+1] > sp->wcmax[0]) {
         } else if (ye[iy] < sp->wcmin[1]) {
         } else if (ye[iy] > sp->wcmax[1]) {
         } else if (ye[iy+1] < sp->wcmin[1]) {
         } else if (ye[iy+1] > sp->wcmax[1]) {
         } else {
*/
            if (sa[k] == SOLID) {
               dowt = FALSE;
               if (ix > 0) {
                  kk = k - 1;
                  dowt = (sa[kk] == FLUID);
               }
               if (ix < (nx-1)) {
                  kk = k + 1;
                  dowt = dowt || (sa[kk] == FLUID);
               }
               if (iy > 0) {
                  kk = (iy-1)*nx + ix;
                  dowt = dowt || (sa[kk] == FLUID);
               }
               if (iy < (ny-1)) {
                  kk = (iy+1)*nx + ix;
                  dowt = dowt || (sa[kk] == FLUID);
               }
               if (dowt) {
#ifdef CURV
                  vis_draw_rect(xe[ix],xe[ix+1],ye[iy],ye[iy+1]);
                  vis_draw_cross(xe[ix],xe[ix+1],ye[iy],ye[iy+1]);
#else
                  vis_grid_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
                  vis_grid_cross(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#endif
               }
            }
/*
         }
*/
         k = k + 1;
      }
   }
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
}

void vis_grid_map(nx, ny, xe, ye, sa, ios)
   int nx,ny,*ios,sa[];
   float xe[],ye[];
{
   int k,ix,iy;

   drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
   k = 0;
   for (iy=0; iy < ny; iy++) {
      for (ix=0; ix < nx; ix++) {
/*
         if (xe[ix] < sp->wcmin[0]) {
         } else if (xe[ix] > sp->wcmax[0]) {
         } else if (xe[ix+1] < sp->wcmin[0]) {
         } else if (xe[ix+1] > sp->wcmax[0]) {
         } else if (ye[iy] < sp->wcmin[1]) {
         } else if (ye[iy] > sp->wcmax[1]) {
         } else if (ye[iy+1] < sp->wcmin[1]) {
         } else if (ye[iy+1] > sp->wcmax[1]) {
         } else {
*/
            if (sa[k] == SOLID) {
#ifdef CURV
               vis_draw_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
               vis_draw_cross(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#else
               vis_grid_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
               vis_grid_cross(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#endif
            } else if (sa[k] != BLANK) {
#ifdef CURV
               vis_draw_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#else
               vis_grid_rect(xe[ix], xe[ix+1], ye[iy], ye[iy+1]);
#endif
            }
/*
         }
*/
         k = k + 1;
      }
   }
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
}


void vis_design_3d_view(ios)
   int *ios;
{
   float dl,scene_size,betav[MAXD],view_point[MAXD],dc;
   int id;

/* construct view-focus is absent */
   if (noview_focus) {
      printf("default:  view_focus=[");
      for (id=0; id < MAXD; id++) {
         view_focus[id] = 0.5*(sp->wcmin[id]+sp->wcmax[id]);
         printf(" %E", view_focus[id]);
      }
      printf(" ]\n");
   }
/* construct view_distance if absent */
   if (noview_distance) {
      scene_size = 0.0;
      for (id=0; id < MAXD; id++) {
         dc = 0.5*(sp->wcmax[id] - sp->wcmin[id]);
         if (fabs(dc) > scene_size) {
            scene_size = fabs(dc);
         }
      }
      if (scene_size <= 0.0) {
         printf("warning:  cannot use scene size for view_distance\n");
         view_distance = 1.0;
      } else {
         view_distance = 2.0*scene_size/tan(view_anglev*DEGRAD);
      }
      printf("default:  view_distance=%E\n", view_distance);
   }
/* construct view_near if missing */
   if (noview_near) {
      view_near = 0.5*view_distance;
      printf("default:  view_near=%E\n", view_near);
   } else if ((view_near <= 0.0) || (view_near >= view_distance)) {
      view_near = 0.5*view_distance;
      printf("warning:  adjust view_near=%E\n", view_near);
   }
/* construct view_far if missing */
   if (noview_far) {
      view_far = 2.0*view_distance;
      printf("default:  view_far=%E\n", view_far);
   } else if (view_far <= view_distance) {
      view_far = 2.0*view_distance;
      printf("warning:  adjust view_far=%E\n", view_far);
   }
/* construct view direction vector and hence view-point */
   betav[0] = cos(view_azim*DEGRAD)*cos(view_altd*DEGRAD);
   betav[1] = sin(view_azim*DEGRAD)*cos(view_altd*DEGRAD);
   betav[2] = sin(view_altd*DEGRAD);
   for (id=0; id < MAXD; id++) {
      view_point[id] = view_focus[id] + betav[id]*view_distance;
   }
/* calculate missing light distance */
   if (nolight_distance) {
      light_distance = 0.0;
      for (id=0; id < MAXD; id++) {
         dl = view_point[id] - view_focus[id];
         light_distance = light_distance + dl*dl;
      }
      light_distance = sqrt(light_distance);
      printf("default:  light_distance=%E\n",light_distance);
   }
/* introduce this to the DR3 implementation */
   dr3_create_scene(view_anglev, view_azim, view_altd, 0.0,
     view_distance, view_focus, view_near, view_far,
      sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
   dr3_lights(light_azim, light_altd, light_distance, light_fusion);
/* draw domain extents */
   if (data_frame) {
      dr3_begin();
      dr3_wire_block(sp->wcmin[0],sp->wcmax[0],sp->wcmin[1],
        sp->wcmax[1],sp->wcmin[2],sp->wcmax[2]);
      dr3_end();
   }
/* get a copy of the scene so far */
   dr3_get_scene(scp);
}

void vis_3d_grid_surface(nd, nc, xe, xa, sa, qa, qamin, qamax, ios)
   int nd,nc[],sa[],*ios;
   float *xe[],*xa[],qa[],qamin,qamax;
{
   int ix,iy,k0,k1,k2,k3,k;
   float x0,x1,y0,y1,dsdq,s,dzds;
   float par[3];

   dr3_coerce_ld(ios);
   if (*ios) {
      dr3_coerce_d(ios);
      if (*ios) {
         dr3_coerce_wire(ios);
      }
   }
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
      if (wire_frame) {
         dr3_seta("wire_frame");
      }
      dr3_begin();
      dr3_begin_surface();
      dsdq = 1.0/(qamax-qamin);
      dzds = (sp->wcmax[2] - sp->wcmin[2]);
      for (iy=1; iy < nc[1]; iy++) {
         y0 = xa[1][iy-1];
         y1 = xa[1][iy];
         for (ix=1; ix < nc[0]; ix++) {
            x0 = xa[0][ix-1];
            x1 = xa[0][ix];
            k0 = (iy-1)*nc[0] + (ix-1);
            k1 = (iy-1)*nc[0] + ix;
            k2 = iy*nc[0] + ix;
            k3 = iy*nc[0] + (ix-1);
            dr3_begin_polygon();
            dr3_polygon_normal(0.0, 0.0, 1.0);
/*
            s = 0.0;
            count = 0;
            if (sa[k0] == FLUID) {
               s = s + (qa[k0] - qamin)*dsdq;
               count++;
            }
            if (sa[k1] == FLUID) {
               s = s + (qa[k1] - qamin)*dsdq;
               count++;
            }
            if (sa[k2] == FLUID) {
               s = s + (qa[k2] - qamin)*dsdq;
               count++;
            }
            if (sa[k3] == FLUID) {
               s = s + (qa[k3] - qamin)*dsdq;
               count++;
            }
            if (count > 0) {
               s = s/(float)count;
            }
*/
            if (sa[k0] == FLUID) {
               s = (qa[k0] - qamin)*dsdq;
               dr3_polygon_scalar(s);
               dr3_polygon_vertex(x0, y0, (dzds*s));
            }
            if (sa[k1] == FLUID) {
               s = (qa[k1] - qamin)*dsdq;
               dr3_polygon_scalar(s);
               dr3_polygon_vertex(x1, y0, (dzds*s));
            }
            if (sa[k2] == FLUID) {
               s = (qa[k2] - qamin)*dsdq;
               dr3_polygon_scalar(s);
               dr3_polygon_vertex(x1, y1, (dzds*s));
            }
            if (sa[k3] == FLUID) {
               s = (qa[k3] - qamin)*dsdq;
               dr3_polygon_scalar(s);
               dr3_polygon_vertex(x0, y1, (dzds*s));
            }
            dr3_end_polygon();
         }
      }
      dr3_end_surface();
      dr3_begin_surface();
      k = 0;
      for (iy=0; iy < nc[1]; iy++) {
         y0 = xe[1][iy];
         y1 = xe[1][iy+1];
         for (ix=0; ix < nc[0]; ix++) {
            x0 = xe[0][ix];
            x1 = xe[0][ix+1];
            if (sa[k] == SOLID) {
               dr3_push_matrix();
               dr3_tran_matrix(x0, y0, sp->wcmin[2]);
               par[0] = (x1-x0);
               par[1] = (y1-y0);
               par[2] = (sp->wcmax[2] - sp->wcmin[2]);
#ifdef STL
               dr3_add_solid("block",par,3);
#else
               dr3_add_solid("block",par);
#endif
               dr3_pop_matrix();
            }
            k++;
         }
      }
      dr3_end_surface();
      dr3_end();
      dr3_seta("nowire_frame");
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         if (!scp->vector_bitmap) {
            vis_ld_scales();
         }
      }
   }
}

void vis_3d_draw_grid(nc, xe, sa, ios)
   int nc[],sa[],*ios;
   float *xe[];
{
   int id,ix,iy,iz,k,kl,kr,con[MAXD];

   if (wire_frame) {
      dr3_seta("wire_frame");
   }
/* zap materials to save colour space */
   dr3_zap_materials();
/* coerce an LM or WIRE image */
   dr3_coerce_lm(ios);
   if (*ios) {
      dr3_coerce_wire(ios);
   }
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
/*    do surfaces? */
      dr3_select_material("default");
/*    work through grid looking at solid cells with fluid neighbours */
      k = 0;
      for (iz=0; iz < nc[2]; iz++) {
         con[2] = iz;
         for (iy=0; iy < nc[1]; iy++) {
            con[1] = iy;
            for (ix=0; ix < nc[0]; ix++) {
               con[0] = ix;
               if (sa[k] == SOLID) {
                  for (id=0; id < MAXD; id++) {
                     vis_left(nc, k, id, con, &kl);
                     if (kl != k) {
                        if (sa[kl] == FLUID) {
/*                         draw the kl:k facet */
                           vis_grid_facet(id,LEFT,con,xe);
                        }
                     }
                     vis_right(nc, k, id, con, &kr);
                     if (kr != k) {
                        if (sa[kr] == FLUID) {
/*                         draw the k:kr facet */
                           vis_grid_facet(id,RIGHT,con,xe);
                        }
                     }
                  }
               }
               k++;
            }
         }
      }
   }
   dr3_seta("nowire_frame");
}

void vis_grid_facet(id, side, con, xe)
   int id,side,con[MAXD];
   float *xe[];
{
   float xl[MAXD],xr[MAXD],vn,x;
   int jd,j;

   for (jd=0; jd < MAXD; jd++) {
      j = con[jd];
/*
      xl[jd] = (float)con[jd]*dxdj_s[jd] + x0_s[jd];
      xr[jd] = xl[jd] + dxdj_s[jd];
*/
      xl[jd] = xe[jd][j];
      xr[jd] = xe[jd][j+1];
   }
   if (side == LEFT) {
      x = xl[id];
   } else {
      x = xr[id];
   }
   vn = (float)(2*side-1);
   dr3_begin_polygon();
   if (id == 0) {
      dr3_polygon_normal(vn, 0.0, 0.0);
      dr3_polygon_vertex(x, xl[1], xl[2]);
      dr3_polygon_vertex(x, xr[1], xl[2]);
      dr3_polygon_vertex(x, xr[1], xr[2]);
      dr3_polygon_vertex(x, xl[1], xr[2]);
   } else if (id == 1) {
      dr3_polygon_normal(0.0, vn, 0.0);
      dr3_polygon_vertex(xl[0], x, xl[2]);
      dr3_polygon_vertex(xr[0], x, xl[2]);
      dr3_polygon_vertex(xr[0], x, xr[2]);
      dr3_polygon_vertex(xl[0], x, xr[2]);
   } else {
      dr3_polygon_normal(0.0, 0.0, vn);
      dr3_polygon_vertex(xl[0], xl[1], x);
      dr3_polygon_vertex(xr[0], xl[1], x);
      dr3_polygon_vertex(xr[0], xr[1], x);
      dr3_polygon_vertex(xl[0], xr[1], x);
   }
   dr3_end_polygon();
}

void vis_left(nc, k, id, con, kl)
   int nc[],k,id,con[],*kl;
{
   int com;

   if (con[id] == 0) {
      *kl = k;
   } else {
      com = con[id];
      con[id] = com-1;
      *kl = con[0] + nc[0]*(con[1] + nc[1]*con[2]);
      con[id] = com;
   }
}

void vis_right(nc, k, id, con, kr)
   int nc[],k,id,con[],*kr;
{
   int com;

   if (con[id] < nc[id]) {
      com = con[id];
      con[id] = com+1;
      *kr = con[0] + nc[0]*(con[1] + nc[1]*con[2]);
      con[id] = com;
   } else {
      *kr = k;
   }
}

void vis_3d_grid_vector(nc, xe, sa, va, vamax, ios)
   int nc[],sa[],*ios;
   float *xe[],va[],vamax;
{
   int ix,nx,iy,ny,iz,nz,k,l,n,jd,id,i;
   float x,y,z,xl,xr,yl,yr,zl,zr, dx,dy,dz,dldv,dlmin,dl,dk2,dk;

/* default data coordinate */
   if (nodata_coord) {
      id = data_ortho_id;
      data_coord = 0.5*(sp->gcmin[id]+sp->gcmax[id]);
      printf("default:  data_coord=%f\n", data_coord);
   }
/* short-hands */
   nx = nc[0];
   ny = nc[1];
   nz = nc[2];
/* Coerce WIRE only */
   dr3_coerce_wire(ios);
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
/*    arrow length */
      dlmin = BIG;
      for (jd=0; jd < MAXD; jd++) {
         n = nc[jd];
         dl = (xe[jd][n] - xe[jd][0])/(float)n;
         if (dl < dlmin) {
            dlmin = dl;
         }
      }
      dlmin = dlmin*(float)arrow_rate;
      dldv = 0.4*(dlmin/vamax);
      sp->dldv = dldv;
/*    x-planes */
      if (data_ortho_id == 0) {
/*       find closest x-plane */
         for (ix=0; ix < nx; ix++) {
            if (data_coord < xe[0][ix+1]) {
               break;
            }
         }
         if (ix > (nx-1)) {
            ix = nx-1;
         }
/*       A sequence of y-z planes */
         x = 0.5*(xe[0][ix] + xe[0][ix+1]);
         for (iz=0; iz < nz; iz++) {
            zl = xe[2][iz];
            zr = xe[2][iz+1];
            for (iy=0; iy < ny; iy++) {
               yl = xe[1][iy];
               yr = xe[1][iy+1];
               k = ix + nx*(iy + ny*iz);
               if (sa[k] == FLUID) {
                  l = k*MAXD;
                  dx = va[l]*dldv*arrow_scale;
                  dy = va[l+1]*dldv*arrow_scale;
                  dz = va[l+2]*dldv*arrow_scale;
                  dk2 = dx*dx + dy*dy + dz*dz;
                  if (dk2 > (0.0001*dlmin*dlmin)) {
                     y = 0.5*(yl+yr);
                     z = 0.5*(zl+zr);
                     dr3_begin_polyline();
                     dr3_polyline_vertex(x, y, z);
                     dr3_polyline_vertex(x+dx, y+dy, z+dz);
                     dr3_end_polyline();
                     dk = sqrt(dk2);
                     dr3_wire_icon(0.1*dk, x, y, z);
                  }
               }
               for(i=1; i < arrow_rate; i++) {
                  iy++;
               }
            }
            for(i=1; i < arrow_rate; i++) {
               iz++;
            }
         }
/*    y-planes */
      } else if (data_ortho_id == 1) {
/*       find closest y-plane */
         for (iy=0; iy < ny; iy++) {
            if (data_coord < xe[1][iy+1]) {
               break;
            }
         }
         if (iy > (ny-1)) {
            iy = ny-1;
         }
/*       A sequence of x-z planes */
         y = 0.5*(xe[1][iy] + xe[1][iy+1]);
         for (iz=0; iz < nz; iz++) {
            zl = xe[2][iz];
            zr = xe[2][iz+1];
            for (ix=0; ix < nx; ix++) {
               xl = xe[0][ix];
               xr = xe[0][ix+1];
               k = ix + nx*(iy + ny*iz);
               if (sa[k] == FLUID) {
                  l = k*MAXD;
                  dx = va[l]*dldv*arrow_scale;
                  dy = va[l+1]*dldv*arrow_scale;
                  dz = va[l+2]*dldv*arrow_scale;
                  dk2 = dx*dx + dy*dy + dz*dz;
                  if (dk2 > (0.0001*dlmin*dlmin)) {
                     x = 0.5*(xl+xr);
                     z = 0.5*(zl+zr);
                     dr3_begin_polyline();
                     dr3_polyline_vertex(x, y, z);
                     dr3_polyline_vertex(x+dx, y+dy, z+dz);
                     dr3_end_polyline();
                     dk = sqrt(dk2);
                     dr3_wire_icon(0.1*dk, x, y, z);
                  }
               }
               for(i=1; i < arrow_rate; i++) {
                  ix++;
               }
            }
            for(i=1; i < arrow_rate; i++) {
               iz++;
            }
         }
/*    z-planes */
      } else if (data_ortho_id == 2) {
/*       find closest z-plane */
         for (iz=0; iz < nz; iz++) {
            if (data_coord < xe[2][iz+1]) {
               break;
            }
         }
         if (iz > (nz-1)) {
            iz = nz-1;
         }
/*       A sequence of x-y planes */
         z = 0.5*(xe[2][iz] + xe[2][iz+1]);
         for (iy=0; iy < ny; iy++) {
            yl = xe[1][iy];
            yr = xe[1][iy+1];
            for (ix=0; ix < nx; ix++) {
               xl = xe[0][ix];
               xr = xe[0][ix+1];
               k = ix + nx*(iy + ny*iz);
               if (sa[k] == FLUID) {
                  l = k*MAXD;
                  dx = va[l]*dldv*arrow_scale;
                  dy = va[l+1]*dldv*arrow_scale;
                  dz = va[l+2]*dldv*arrow_scale;
                  dk2 = dx*dx + dy*dy + dz*dz;
                  if (dk2 > (0.0001*dlmin*dlmin)) {
                     x = 0.5*(xl+xr);
                     y = 0.5*(yl+yr);
                     dr3_begin_polyline();
                     dr3_polyline_vertex(x, y, z);
                     dr3_polyline_vertex(x+dx, y+dy, z+dz);
                     dr3_end_polyline();
                     dk = sqrt(dk2);
                     dr3_wire_icon(0.1*dk, x, y, z);
                  }
               }
               for(i=1; i < arrow_rate; i++) {
                  ix++;
               }
            }
            for(i=1; i < arrow_rate; i++) {
               iy++;
            }
         }
      }
   }
}

void vis_3d_data_shade(ios)
   int *ios;
{

   dr3_data_shade(ios);
   if (!*ios) {
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         vis_ld_scales();
      }
   }
}

void vis_3d_data_plane(ios)
   int *ios;
{

   if (nodata_distance) {
      data_distance = scp->view_distance;
      printf("default:  data_distance=%f\n",
        data_distance);
   }
   dr3_data_plane(data_distance, ios);
   if (!*ios) {
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         vis_ld_scales();
      }
   }
}

void vis_3d_data_ortho(nc, xe, sa, qa, qamin, qamax, ios)
   int nc[],sa[],*ios;
   float *xe[],qa[],qamin,qamax;
{
   int ix,nx,iy,ny,iz,nz,k,id;
   float x,y,z,xl,xr,yl,yr,zl,zr,s,dpdq;

/* default the coordinate */
   if (nodata_coord) {
      id = data_ortho_id;
      data_coord = 0.5*(sp->gcmin[id]+sp->gcmax[id]);
      printf("default:  data_coord=%f\n", data_coord);
   }
/* short-hand */
   id = data_ortho_id;
   nx = nc[0];
   ny = nc[1];
   nz = nc[2];
   dpdq = 1.0/(qamax - qamin);
/* arrange an LD image */
#ifdef CURV_PLAY
   dr3_coerce_d(ios);
#else
   dr3_coerce_ld(ios);
#endif
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
      if (wire_frame) {
         dr3_seta("wire_frame");
      }
      dr3_seta("nohard_surface");
      dr3_seta("lucid_surface");
/*    x-planes */
      if (data_ortho_id == 0) {
/*       find closest x-plane */
         for (ix=0; ix < nx; ix++) {
            if (data_coord < xe[0][ix+1]) {
               break;
            }
         }
         if (ix > (nx-1)) {
            ix = nx-1;
         }
/*       A sequence of y-z planes */
         x = 0.5*(xe[0][ix] + xe[0][ix+1]);
         for (iz=0; iz < nz; iz++) {
            zl = xe[2][iz];
            zr = xe[2][iz+1];
            for (iy=0; iy < ny; iy++) {
               yl = xe[1][iy];
               yr = xe[1][iy+1];
               k = ix + nx*(iy + ny*iz);
#ifdef MARIO
               if (zl > sp->gcmax[2]) {
               } else if (zr < sp->gcmin[2]) {
               } else if (yl > sp->gcmax[1]) {
               } else if (yr < sp->gcmin[1]) {
               } else if (sa[k] == FLUID) {
#else
               if (sa[k] == FLUID) {
#endif
                  dr3_begin_polygon();
                  dr3_polygon_normal(0.0, 1.0, 1.0);
                  s = (qa[k]-qamin)*dpdq;
                  dr3_polygon_scalar(s);
                  dr3_polygon_vertex(x, yl, zl);
                  dr3_polygon_vertex(x, yr, zl);
                  dr3_polygon_vertex(x, yr, zr);
                  dr3_polygon_vertex(x, yl, zr);
                  dr3_polygon_vertex(x, yl, zl);
                  dr3_end_polygon();
               }
            }
         }
/*    y-planes */
      } else if (data_ortho_id == 1) {
/*       find closest y-plane */
         for (iy=0; iy < ny; iy++) {
            if (data_coord < xe[1][iy+1]) {
               break;
            }
         }
         if (iy > (ny-1)) {
            iy = ny-1;
         }
/*       A sequence of x-z planes */
         y = 0.5*(xe[1][iy] + xe[1][iy+1]);
         for (iz=0; iz < nz; iz++) {
            zl = xe[2][iz];
            zr = xe[2][iz+1];
            for (ix=0; ix < nx; ix++) {
               xl = xe[0][ix];
               xr = xe[0][ix+1];
               k = ix + nx*(iy + ny*iz);
#ifdef MARIO
               if (zl > sp->gcmax[2]) {
               } else if (zr < sp->gcmin[2]) {
               } else if (xl > sp->gcmax[0]) {
               } else if (xr < sp->gcmin[0]) {
               } else if (sa[k] == FLUID) {
#else
               if (sa[k] == FLUID) {
#endif
                  dr3_begin_polygon();
                  s = (qa[k]-qamin)*dpdq;
                  dr3_polygon_scalar(s);
                  dr3_polygon_vertex(xl, y, zl);
                  dr3_polygon_vertex(xr, y, zl);
                  dr3_polygon_vertex(xr, y, zr);
                  dr3_polygon_vertex(xl, y, zr);
                  dr3_polygon_vertex(xl, y, zl);
                  dr3_end_polygon();
               }
            }
         }
/*    z-planes */
      } else if (data_ortho_id == 2) {
/*       find closest z-plane */
         for (iz=0; iz < nz; iz++) {
            if (data_coord < xe[2][iz+1]) {
               break;
            }
         }
         if (iz > (nz-1)) {
            iz = nz-1;
         }
/*       A sequence of x-y planes */
         z = 0.5*(xe[2][iz] + xe[2][iz+1]);
         for (iy=0; iy < ny; iy++) {
            yl = xe[1][iy];
            yr = xe[1][iy+1];
            for (ix=0; ix < nx; ix++) {
               xl = xe[0][ix];
               xr = xe[0][ix+1];
               k = ix + nx*(iy + ny*iz);
#ifdef MARIO
               if (yl > sp->gcmax[1]) {
               } else if (yr < sp->gcmin[1]) {
               } else if (xl > sp->gcmax[0]) {
               } else if (xr < sp->gcmin[0]) {
               } else if (sa[k] == FLUID) {
#else
               if (sa[k] == FLUID) {
#endif
                  dr3_begin_polygon();
                  s = (qa[k]-qamin)*dpdq;
                  dr3_polygon_scalar(s);
                  dr3_polygon_vertex(xl, yl, z);
                  dr3_polygon_vertex(xr, yl, z);
                  dr3_polygon_vertex(xr, yr, z);
                  dr3_polygon_vertex(xl, yr, z);
                  dr3_end_polygon();
               }
            }
         }
      }
      dr3_seta("hard_surface");
      dr3_seta("nolucid_surface");
      dr3_seta("nowire_frame");
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         vis_ld_scales();
      }
   }
}

void vis_3d_scatter_scalar(npt,nd,xpt,spt,qpt,qamin,qamax,ios)
   int npt,nd,*ios,spt[];
   float xpt[],qpt[],qamin,qamax;
{
   int j,id,k;
   float dpdq,x,y,l,lmin,s,z;

   dpdq = 1.0/(qamax - qamin);
   dr3_coerce_ld(ios);
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
      dr3_seta("wire_frame");
      dr3_seta("hard_surface");
      dr3_seta("nolucid_surface");
      dr3_begin_surface();
      lmin = BIG;
      for (id=0; id < 2; id++) {
         l = sp->wcmax[id] - sp->wcmin[id];
         if (l < lmin) {
            lmin = l;
         }
      }
      lmin = arrow_ratio*lmin;
      k = 0;
      for (j=0; j < npt; j++) {
         x = xpt[k];
         y = xpt[k+1];
         z = xpt[k+2];
         k = k + nd;
         if (x > sp->wcmax[0]) {
         } else if (x > sp->wcmax[0]) {
         } else if (y < sp->wcmin[1]) {
         } else if (y > sp->wcmax[1]) {
         } else if (z < sp->wcmin[2]) {
         } else if (z > sp->wcmax[2]) {
         } else {
            if (spt[j] == FLUID) {
               s = dpdq*(qpt[j]-qamin);
               dr3_solid_icon(s, lmin, x, y, z);
            } else {
               dr3_wire_icon(lmin, x, y, z);
            }
         }
      }
      dr3_end_surface();
      dr3_seta("hard_surface");
      dr3_seta("nolucid_surface");
      dr3_seta("nowire_frame");
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         vis_ld_scales();
      }
   }
}

void vis_3d_scatter_vector(ios)
   int *ios;
{
   printf("vis_3d_scatter_vector:  not implemented yet\n");
   *ios = OK;
}

void vis_3d_solid(ios)
   int *ios;
{
   int iob,iob0,nobi,icl,i,colour_index;
   float x,y,z,azim,altd,twistd,r,g,b;
   bool cut,visible;

/* save current colour */
   if (wp->colour_avail) {
      drp_get_colour(&colour_index);
   } else if (wp->rgb_avail) {
      drp_get_rgb(&r, &g, &b);
   }
/* load materials */
   dr3_zap_materials();
   for (i=0; i < material_count; i++) {
      dr3_define_material(material_name[i], material_red[i],
        material_green[i], material_blue[i], material_shine[i],
        0.0);
   }
/* indicate sky colour */
   dr3_sky(sky_red, sky_green, sky_blue);
/* Coerce LM or WIRE */
   dr3_coerce_lm(ios);
   if (*ios) {
      dr3_coerce_wire(ios);
   }
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
      dr3_get_shade(scp);
/*    have we already displayed the material key? */
      if (!sp->nokey) {
         if (!sp->key_drawn) {
            vis_key();
         }
         if (!sp->mat_scales_drawn) {
            if (!scp->vector_bitmap) {
               vis_mat_scales();
            }
         }
      }
/*    now go through the csg data base */
      dr3_begin();
      for (icl=0; icl < vis_cp->ncl; icl++) {
/*      find material index for cluster */
         visible = TRUE;
         if (vis_tag_find(vis_cp->cltag[icl], &i)) {
            if (csame("invisible",tag_material[i])) {
               visible = FALSE;
            } else {
               dr3_select_material(tag_material[i]);
            }
         } else {
            printf("warning:  tag=%s not found\n", vis_cp->cltag[icl]);
            dr3_select_material(material_name[0]);
         }
/*      bypass invisible material? */
         if (visible) {
/*         go through primitives to find cuts */
            iob0 = vis_cp->clob0[icl];
            nobi = vis_cp->clnob[icl];
            cut = FALSE;
            for (iob=0; iob < nobi; iob++) {
               i = iob+iob0;
               cut = cut || (vis_cp->obop[i] == SUB);
            }
            if (cut) {
               dr3_begin_solid();
            } else {
               dr3_begin_surface();
            }
/*         go through primitives drawing */
            for (iob=0; iob < nobi; iob++) {
               i = iob+iob0;
               dr3_push_matrix();
               x = vis_cp->at[i][0];
               y = vis_cp->at[i][1];
               z = vis_cp->at[i][2];
               azim = vis_cp->azim[i];
               altd = vis_cp->altd[i];
               twistd = vis_cp->twistd[i];
               dr3_rotran_matrix(x, y, z, azim, altd, twistd);
               if (vis_cp->obop[i] == SUB) {
#ifdef STL
                  dr3_cut_solid(vis_cp->obnam[i], vis_cp->obdat[i],
                    vis_cp->obndat[i]);
#else
                  dr3_cut_solid(vis_cp->obnam[i], vis_cp->obdat[i]);
#endif
               } else if (vis_cp->obop[i] == ADD) {
#ifdef STL
                  dr3_add_solid(vis_cp->obnam[i], vis_cp->obdat[i],
                     vis_cp->obndat[i]);
#else
                  dr3_add_solid(vis_cp->obnam[i], vis_cp->obdat[i]);
#endif
               }
               dr3_pop_matrix();
            }
            if (cut) {
               dr3_end_solid();
            } else {
               dr3_end_surface();
            }
         }
      }
      dr3_end();
/*   Tidy up */
      if (wp->colour_avail) {
         drp_set_colour(colour_index);
      } else if (wp->rgb_avail) {
         drp_set_rgb(r, g, b);
      }
   }
}

void vis_mat_scales() {
   int i,offset,nj,j,k,colour_index;
   float h1,h2,dhdj,dv,v,v1,v2,ha,hb,h,s,dsdj,r,g,b,
     keep_r,keep_g,keep_b;
   offset = 0;

   if (wp->colour_avail) {
      drp_get_colour(&colour_index);
   } else if (wp->rgb_avail) {
      drp_get_rgb(&keep_r, &keep_g, &keep_b);
   }
   if (scp->colour_bitmap) {
      offset = scp->shade_offset;
      nj = scp->shade_nc[0];
   } else {
      nj = 16;
      dsdj = 1.0/(float)(nj-1);
   }
   h1 = sp->hkmin + 0.5*wp->char_width;
   h2 = h1 + 2.0*wp->char_width;
   dhdj = (h2-h1)/(float)(nj-1);
   drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
   for (i=0; i <  scp->shade_nc[1]; i++) {
      if (sp->key_count < sp->key_size) {
         r = scp->shade_mat_r[i];
         g = scp->shade_mat_g[i];
         b = scp->shade_mat_b[i];
/*
         r = material_red[i];
         g = material_green[i];
         b = material_blue[i];
*/
         dv = wp->char_height;
         v = sp->vkmax - 1.5*dv*((float)sp->key_count + 0.5);
         v1 = v - 0.5*dv;
         v2 = v + 0.5*dv;
         for (j=1; j < nj; j++) {
            ha = dhdj*(float)(j-1) + h1;
            hb = ha + dhdj;
            if (scp->colour_bitmap) {
               k = (j + nj*i) + offset;
               drp_set_colour(k);
            } else {
               s = dsdj*(float)j;
               drp_set_rgb(s*r, s*g, s*b);
            }
            drp_rect(ha, hb, v1, v2);
         }
         if (wp->colour_avail) {
            drp_set_colour(wp->colour_black);
         } else if (wp->rgb_avail) {
            drp_set_rgb(0.0, 0.0, 0.0);
         }
         vis_box(h1, h2, v1, v2);
         h = h2 + 0.5*wp->char_width;
         drp_text(h, v, scp->shade_mat[i], "lm");
         sp->key_count = sp->key_count + 1;
      } else {
         printf("warning:  some material key lost\n");
         break;
      }
   }
   if (wp->colour_avail) {
      drp_set_colour(colour_index);
   } else if (wp->rgb_avail) {
      drp_set_rgb(keep_r, keep_g, keep_b);
   }
   sp->mat_scales_drawn = TRUE;
}

void vis_ld_scales() {
   int i,j,k,ni,nj,colour_index,offset;
   float dhdj,dvdi,dsdj,dsdi,s,r,g,b,v1,v2,h1,h2,
     keep_r,keep_g,keep_b;
   offset = 0;

   if (!sp->scales_drawn) {
      vis_scales();
   }
   if (wp->colour_avail) {
      drp_get_colour(&colour_index);
   } else if (wp->rgb_avail) {
      drp_get_rgb(&keep_r, &keep_g, &keep_b);
   }
   if (scp->colour_bitmap) {
      offset = scp->shade_offset;
      nj = scp->shade_nc[0];
      ni = scp->shade_nc[1];
   } else {
      nj = 64;
      ni = 16;
   }
   dvdi = (sp->vsmax - sp->vsmin)/(float)(ni-1);
   dhdj = (sp->hsmax - sp->hsmin)/(float)(nj-1);
   dsdi = 1.0/(float)(ni-1);
   dsdj = 1.0/(float)(nj-1);
   for (i=1; i < ni; i++) {
      if (!scp->colour_bitmap) {
         s = dsdi*(float)i;
         dr3_rain(s, &r, &g, &b);
      }
      v1 = dvdi*(float)(i-1) + sp->vsmin;
      v2 = v1 + dvdi;
      for (j=1; j < nj; j++) {
         h1 = dhdj*(float)(j-1) + sp->hsmin;
         h2 = h1 + dhdj;
         if (scp->colour_bitmap) {
            k = (j + nj*i) + offset;
            drp_set_colour(k);
         } else {
            s = dsdj*(float)j;
            drp_set_rgb(s*r, s*g, s*b);
         }
         drp_rect(h1, h2, v1, v2);
      }
   }
   if (scp->colour_bitmap) {
      drp_set_colour(colour_index);
   } else if (wp->rgb_avail) {
      drp_set_rgb(keep_r, keep_g, keep_b);
   }
}

/* SCENE PRIMITIVES */

/* vis_zap_scene - zap scene associated with current subwindow */
void vis_zap_scene() {
/* zap image if present */
   if (!(scp->noimage)) {
      vis_zap_image();
   }
   sp->noscene = TRUE;
}


/* vis_dump_pic - dump picture(scene) to file */
void vis_dump_pic(file, ios)
   char file[];
   int *ios;
{
   int unit,i,id;

   if (*ios) {
   } else if (wc < 0) {
      printf("error:  no worksurface selected\n");
      *ios = ERR;
   } else if (sp->noscene) {
      printf("error:  subwindow does not contain scene\n");
      *ios = ERR;
   } else {
      ds_create(file, &unit, ios);
      if (*ios) {
         printf("error:  creating scene file '%s'\n", file);
      } else {
         ds_put_begin(unit, "scene", ios);
         if (*ios) {
            printf("error:  dumping type for scene '%s' \n",file);
         } else {
/*         graph bits */
            ds_put0i(unit, "sp_caption_count", sp->caption_count, ios);
            for (i=0; i < sp->caption_count; i++) {
               ds_put0s(unit, "sp_caption", sp->caption[i], ios);
            }
            ds_put0i(unit, "sp_nd", sp->nd, ios);
            ds_put0i(unit, "sp_ndc", sp->ndc, ios);
            for (id=0; id < sp->nd; id++) {
               ds_put0s(unit, "sp_gname", sp->gname[id], ios);
               ds_put0s(unit, "sp_glabel", sp->glabel[id], ios);
               ds_put0r(unit, "sp_gcmin", sp->gcmin[id], ios);
               ds_put0r(unit, "sp_gcmax", sp->gcmax[id], ios);
               ds_put0r(unit, "sp_wcmin", sp->wcmin[id], ios);
               ds_put0r(unit, "sp_wcmax", sp->wcmax[id], ios);
            }
/*         let dr3 do the scene structure */
            dr3_put_scene(scp);
            dr3_push_image(scp);
            dr3_dump_scene(unit, ios);
            dr3_pop_image(scp);
            ds_put_begin(unit, "other", ios);
            ds_put_end(unit, "other",ios);
            ds_put_end(unit, "scene", ios);
         }
         ds_close(unit, ios);
      }
   }
}

void vis_load_pic(file, ios)
   char file[];
   int *ios;
{
   int unit,id,i;
   char caption[SZLIN],name[SZNAM],type[SZNAM];

   if (*ios) {
   } else if (wc < 0) {
      printf("error:  no worksurface selected\n");
      *ios = ERR;
   } else {
      vis_erase();
      ds_open(file, &unit, ios);
      if (*ios) {
         printf("error:  opening scene file '%s'\n", file);
      } else {
         ds_get_begin(unit, "scene", ios);
         if (*ios) {
            printf("error:  loading type for scene '%s' \n",file);
         } else {
/*         graph bits */
            ds_get0i(unit, "sp_caption_count", &(sp->caption_count),
              ios);
            for (i=0; i < sp->caption_count; i++) {
               ds_get0s(unit, "sp_caption", caption, ios);
               if (*ios) {
               } else if (i < MAXCAP) {
                  ccopy(caption, 0, SZLIN, sp->caption[i]);
               }
            }
            if (sp->caption_count > MAXCAP) {
               sp->caption_count = MAXCAP;
            }
            ds_get0i(unit, "sp_nd", &(sp->nd), ios);
            ds_get0i(unit, "sp_ndc", &(sp->ndc), ios);
            if (*ios) {
            } else if (sp->nd != MAXD) {
               printf("error:  picture has wrong number of dims\n");
               *ios = ERR;
            } else {
               for (id=0; id < sp->nd; id++) {
                  ds_get0s(unit, "sp_gname", sp->gname[id], ios);
                  ds_get0s(unit, "sp_glabel", sp->glabel[id], ios);
                  ds_get0r(unit, "sp_gcmin", &(sp->gcmin[id]), ios);
                  ds_get0r(unit, "sp_gcmax", &(sp->gcmax[id]), ios);
                  ds_get0r(unit, "sp_wcmin", &(sp->wcmin[id]), ios);
                  ds_get0r(unit, "sp_wcmax", &(sp->wcmax[id]), ios);
               }
            }
/*         Create a graph space in which dr3 can load its stuff */
            vis_redesign_3d_axes();
/*         get dr3 to load scene structure */
            dr3_load_scene(unit, ios);
/*          Get any additional stuff */
            while (ds_get(unit, name, type, ios)) {
/*             other */
               if (!csame("begin",name)) {
                  printf("error:  %s unexpected in scene-record\n", name);
                  *ios = ERR;
               } else if (csame("other",type)) {
                  ds_get_end(unit, "other", ios);
               } else {
                  printf("error:  %s unexpected in scene-record\n", type);
                  *ios = ERR;
               }
            }
/*          Get the scene? */
            if (*ios) {
               printf("error:  loading (3d) scene\n");
            } else {
               dr3_get_scene(scp);
               sp->noscene = FALSE;
               dr3_pop_image(scp);
            }
         }
         ds_close(unit, ios);
      }
      if (!*ios) {
         sp->clear = FALSE;
         sp->drawn = TRUE;
/*       Captions */
         drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin, sp->vgmax);
         for (i=0; i < sp->caption_count; i++) {
            if (i < sp->caption_lines) {
               drp_text(sp->caption_h, sp->caption_v, sp->caption[i], "lb");
               sp->caption_v = sp->caption_v - 1.5*wp->char_height;
            } else {
               printf("note:  lost caption line\n");
            }
         }
/*       key-scales */
         if (!sp->nokey) {
            if (!sp->nokey) {
               vis_key();
            }
            if (!sp->mat_scales_drawn) {
               if (scp->shade_index == SHADE_LM) {
                  if (!scp->vector_bitmap) {
                     vis_mat_scales();
                  }
               }
            }
         }
         sp->noqr = TRUE;
         sp->nocont = TRUE;
         sp->nocmap = TRUE;
         sp->novel = TRUE;
      }
#ifndef GFLUSH_BUG
      drp_gflush();
#endif
   }
#ifdef GFLUSH_BUG
   drp_gflush();
#endif
}

/* IMAGE PRIMITIVES */

/* vis_zap_image - zap image if any */
void vis_zap_image() {
   if (scp->dim != NULL) {
      free((scp->dim));
      scp->dim = NULL;
   }
/* Colourmap images */
   if (scp->cim != NULL) {
      free((scp->cim));
      scp->cim = NULL;
/* RGB images */
   } else if (scp->rgbim != NULL) {
      free((scp->rgbim));
      scp->rgbim = NULL;
   }
   scp->noimage = TRUE;
}

/* UTILITY PRIMITIVES */

void vis_bp() {
#ifdef CURV98CONT
   if (sp->curv && (vis_gvp->nd+vis_gvp->ndo==4)) {
      dr3_begin_polyline();
   } else {
      drp_bp();
   }
#else
   drp_bp();
#endif
}

void vis_ap(x,y)
   float x,y;
{
   float h,v;
#ifdef CURV
   float u[MAXD],w[MAXD];
#endif

#ifdef CURV98CONT
   if (sp->curv && (vis_gvp->nd+vis_gvp->ndo==4)) {
      vis_expand(vis_gvp->kad, x, y, vis_gvp->xo[0], u);
      dr3_polyline_vertex(u[0], u[1], u[2]);
   } else {
#endif
#ifdef CURV
   if (sp->curv) {
      u[0] = x;
      u[1] = y;
      geo_x(2, u, w);
      h = (w[0]-sp->x0)*sp->dhdx + sp->h0;
      v = (w[1]-sp->y0)*sp->dvdy + sp->v0;
   } else {
      h = (x-sp->x0)*sp->dhdx + sp->h0;
      v = (y-sp->y0)*sp->dvdy + sp->v0;
   }
#else
   h = (x-sp->x0)*sp->dhdx + sp->h0;
   v = (y-sp->y0)*sp->dvdy + sp->v0;
#endif
   drp_ap(h, v);
#ifdef CURV98CONT
   }
#endif
}

void vis_ep() {
#ifdef CURV98CONT
   if (sp->curv && (vis_gvp->nd+vis_gvp->ndo==4)) {
      dr3_end_polyline();
   } else {
      drp_ep();
   }
#else
   drp_ep();
#endif
}

void vis_marker(mark, x, y)
   int mark;
   float x,y;
{
   float h,v;
#ifdef CURV
   float u[MAXD],w[MAXD];
#endif

#ifdef CURV
   if (sp->curv) {
      u[0] = x;
      u[1] = y;
      geo_x(2, u, w);
      h = (w[0]-sp->x0)*sp->dhdx + sp->h0;
      v = (w[1]-sp->y0)*sp->dvdy + sp->v0;
   } else {
      h = (x-sp->x0)*sp->dhdx + sp->h0;
      v = (y-sp->y0)*sp->dvdy + sp->v0;
   }
#else
   h = (x-sp->x0)*sp->dhdx + sp->h0;
   v = (y-sp->y0)*sp->dvdy + sp->v0;
#endif
/*
   if ((h < sp->hdmin) || (h > sp->hdmax)) {
   } else if ((v < sp->vdmin) || (v > sp->vdmax)) {
   } else {
      drp_marker(mark, h, v);
   }
*/
   drp_marker(mark, h, v);
}

void vis_text(x, y, str, just)
   float x,y;
   char str[],just[];
{
   float h,v;
#ifdef CURV
   float u[MAXD],w[MAXD];
#endif

#ifdef CURV
   if (sp->curv) {
      u[0] = x;
      u[1] = y;
      geo_x(2, u, w);
      h = (w[0]-sp->x0)*sp->dhdx + sp->h0;
      v = (w[1]-sp->y0)*sp->dvdy + sp->v0;
   } else {
      h = (x-sp->x0)*sp->dhdx + sp->h0;
      v = (y-sp->y0)*sp->dvdy + sp->v0;
   }
#else
   h = (x-sp->x0)*sp->dhdx + sp->h0;
   v = (y-sp->y0)*sp->dvdy + sp->v0;
#endif
/*
   if ((h < sp->hdmin) || (h > sp->hdmax)) {
   } else if ((v < sp->vdmin) || (v > sp->vdmax)) {
   } else {
      drp_text(h, v, str, just);
   }
*/
   drp_text(h, v, str, just);
}

#ifndef CURV
void vis_cmap_rect(x1, x2, y1, y2, ci)
   float x1,x2,y1,y2;
   int ci;
{
   float h1,h2,v1,v2;

   h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
   v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
   h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
   v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
   drp_set_colour(ci);
   drp_rect(h1, h2, v1, v2);
/*
   printf("cmap_rect:  ci=%d,h1=%f,h2=%f,v1=%f,v2=%f\n", ci,
     h1,h2,v1,v2);
*/
}
#endif

#ifndef CURV
void vis_rect(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float h1,h2,v1,v2;

   h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
   v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
   h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
   v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
   drp_rect(h1, h2, v1, v2);
}
#endif

#ifndef CURV
void vis_black_rect(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float h1,h2,v1,v2;

   h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
   v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
   h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
   v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
   if (wp->colour_avail) {
      drp_set_colour(wp->colour_black);
   } else if (wp->rgb_avail) {
      drp_set_rgb(0.0, 0.0, 0.0);
   }
   drp_rect(h1, h2, v1, v2);
}
#endif

#ifndef CURV
void vis_grid_rect(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float h1,h2,v1,v2;

   h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
   v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
   h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
   v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
   drp_bp();
   drp_ap(h1,v1);
   drp_ap(h2,v1);
   drp_ap(h2,v2);
   drp_ap(h1,v2);
   drp_ap(h1,v1);
   drp_ep();
}
#endif

#ifndef CURV
void vis_grid_cross(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float h1,h2,v1,v2;

   h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
   v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
   h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
   v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
   drp_bp();
   drp_ap(h1,v1);
   drp_ap(h2,v2);
   drp_ep();
   drp_bp();
   drp_ap(h2,v1);
   drp_ap(h1,v2);
   drp_ep();
}
#endif

void vis_box(h1, h2, v1, v2)
   float h1,h2,v1,v2;
{
   drp_bp();
   drp_ap(h1,v1);
   drp_ap(h2,v1);
   drp_ap(h2,v2);
   drp_ap(h1,v2);
   drp_ap(h1,v1);
   drp_ep();
}

void vis_redesign_3d_axes()
{
   float hlen,vlen,dv,dvtop,dvbot;

/* graph box less borders */
   hlen = (sp->hmax - sp->hmin);
   vlen = (sp->vmax - sp->vmin);
   sp->hgmin = sp->hmin + hborder*hlen;
   sp->hgmax = sp->hmax - hborder*hlen;
   sp->vgmin = sp->vmin + vborder*vlen;
   sp->vgmax = sp->vmax - vborder*hlen;
/* some attributes */
   sp->key_number_field = key_number_field;
   sp->key_label_field = key_label_field;
/* account for any scales-keys */
   if (!noscales || !nokey) {
      vis_design_scales();
      sp->hgmax = sp->hxmin;
   } else {
      sp->noscales = TRUE;
      sp->nokey = TRUE;
   }
/* Account for caption box */
   sp->caption_lines = caption_lines;
   dvtop = 0.5*wp->char_height;
   dvbot = 0.5*wp->char_height;
   dv = 1.5*(float)caption_lines*wp->char_height;
   sp->caption_h = sp->hgmin + (float)caption_inset*wp->char_width;
   if (caption_below) {
      dvbot = dvbot + dv;
      sp->caption_v = sp->vgmin + ((float)(caption_lines-1)+0.5)*
         wp->char_height;
   } else {
      dvtop = dvtop + dv;
      sp->caption_v = sp->vgmax - 2.0*wp->char_height;
   }
   drp_set_hv(sp->hgmin, sp->hgmax, sp->vgmin+dvbot, sp->vgmax-dvtop);
/* set some flags */
}

void vis_rain_dance(z, r, g, b)
   float z,*r,*g,*b;
{
#ifdef HSV_RAIN
   float h,rgb[3];
#endif

#ifndef NOGREY
   if (vis_grey_rain) {
      *r = 1.0-z;
      *g = 1.0-z;
      *b = 1.0-z;
   } else {
#endif
#ifdef HSV_RAIN
   h = 0.66666*(1.0-z);
   vis_h_to_rgb(h, rgb);
   *r = rgb[0];
   *g = rgb[1];
   *b = rgb[2];
#else
   *r = z;
   if (z < 0.5) {
      *g = 2.0*z;
   } else {
      *g = 2.0 - 2.0*z;
   }
   *b = 1.0-z;
#endif
#ifndef NOGREY
   }
#endif
}

/* vis_default_palet - define target colours for line-graph palet */
void vis_default_palet() {
   npalet = 0;
/* black */
   ccopy("1", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 0.0;
   palet_g[npalet] = 0.0;
   palet_b[npalet] = 0.0;
   npalet++;
/* red */
   ccopy("2", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 1.0;
   palet_g[npalet] = 0.0;
   palet_b[npalet] = 0.0;
   npalet++;
/* green */
   ccopy("3", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 0.0;
   palet_g[npalet] = 1.0;
   palet_b[npalet] = 0.0;
   npalet++;
/* yellow */
   ccopy("4", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 1.0;
   palet_g[npalet] = 1.0;
   palet_b[npalet] = 0.0;
   npalet++;
/* blue */
   ccopy("5", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 0.0;
   palet_g[npalet] = 0.0;
   palet_b[npalet] = 1.0;
   npalet++;
/* magenta */
   ccopy("6", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 1.0;
   palet_g[npalet] = 0.0;
   palet_b[npalet] = 1.0;
   npalet++;
/* cyan */
   ccopy("7", 0, SZNAM, palet_name[npalet]);
   palet_r[npalet] = 0.0;
   palet_g[npalet] = 1.0;
   palet_b[npalet] = 1.0;
   npalet++;
}

void vis_select_palet() {
   int i,offset,ndi,nci[4];

   i = line_palet;
   if (wp->colour_avail) {
      drp_get_ramp("PALET", &ndi, nci, &offset);
      if (offset < 0) {
         nci[0] = npalet;
         drp_put_ramp("PALET", 1, nci, &offset);
         if (offset < 0) {
            printf("warning:  unable to implement palet\n");
         }
      }
      if (offset >= 0) {
         drp_set_colour(offset+line_palet);
         drp_set_rgb(palet_r[i], palet_g[i], palet_b[i]);
      }
   } else if (wp->rgb_avail) {
      drp_set_rgb(palet_r[i], palet_g[i], palet_b[i]);
   }
}

#ifdef CURV
void vis_fill_rect(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float v[4][MAXD],w[MAXD],hh,vv;
   int i;
   float h1,h2,v1,v2;

   if (sp->curv) {
/*    bottom-left vertex */
      v[0][0] = x1;
      v[0][1] = y1;
/*    bottom-right vertex */
      v[1][0] = x2;
      v[1][1] = y1;
/*    top-right vertex */
      v[2][0] = x2;
      v[2][1] = y2;
/*    top-left vertex */
      v[3][0] = x1;
      v[3][1] = y2;
      drp_bg();
      for (i=0; i < 4; i++) {
         geo_x(2, v[i], w);
         hh = (w[0]-sp->x0)*sp->dhdx + sp->h0;
         vv = (w[1]-sp->y0)*sp->dvdy + sp->v0;
         drp_ag(hh, vv);
      }
      drp_eg();
   } else {
      h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
      v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
      h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
      v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
      drp_rect(h1, h2, v1, v2);
   }
}
#endif

#ifdef CURV
void vis_draw_rect(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float v[5][MAXD],w[MAXD],hh,vv;
   int i;
   float h1,h2,v1,v2;

   if (sp->curv) {
/*    bottom-left vertex */
      v[0][0] = x1;
      v[0][1] = y1;
/*    bottom-right vertex */
      v[1][0] = x2;
      v[1][1] = y1;
/*    top-right vertex */
      v[2][0] = x2;
      v[2][1] = y2;
/*    top-left vertex */
      v[3][0] = x1;
      v[3][1] = y2;
/*    and back */
      v[4][0] = v[0][0];
      v[4][1] = v[0][1];
      drp_bp();
      for (i=0; i < 5; i++) {
         geo_x(2, v[i], w);
         hh = (w[0]-sp->x0)*sp->dhdx + sp->h0;
         vv = (w[1]-sp->y0)*sp->dvdy + sp->v0;
         drp_ap(hh, vv);
      }
      drp_ep();
   } else {
      h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
      v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
      h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
      v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
      drp_bp();
      drp_ap(h1,v1);
      drp_ap(h2,v1);
      drp_ap(h2,v2);
      drp_ap(h1,v2);
      drp_ap(h1,v1);
      drp_ep();
   }
}
#endif

#ifdef CURV
void vis_draw_cross(x1, x2, y1, y2)
   float x1,x2,y1,y2;
{
   float v[4][MAXD],w[4][MAXD],hh,vv;
   int i;
   float h1,h2,v1,v2;

   if (sp->curv) {
/*    bottom-left vertex */
      v[0][0] = x1;
      v[0][1] = y1;
/*    bottom-right vertex */
      v[1][0] = x2;
      v[1][1] = y1;
/*    top-right vertex */
      v[2][0] = x2;
      v[2][1] = y2;
/*    top-left vertex */
      v[3][0] = x1;
      v[3][1] = y2;
      for (i=0; i < 4; i++) {
         geo_x(2, v[i], w[i]);
      }
/*    [0,0] to [1,1] */
      drp_bp();
      hh = (w[0][0]-sp->x0)*sp->dhdx + sp->h0;
      vv = (w[0][1]-sp->y0)*sp->dvdy + sp->v0;
      drp_ap(hh,vv);
      hh = (w[2][0]-sp->x0)*sp->dhdx + sp->h0;
      vv = (w[2][1]-sp->y0)*sp->dvdy + sp->v0;
      drp_ap(hh,vv);
      drp_ep();
/*    [0,1] to [1,0]) */
      drp_bp();
      hh = (w[1][0]-sp->x0)*sp->dhdx + sp->h0;
      vv = (w[1][1]-sp->y0)*sp->dvdy + sp->v0;
      drp_ap(hh,vv);
      hh = (w[3][0]-sp->x0)*sp->dhdx + sp->h0;
      vv = (w[3][1]-sp->y0)*sp->dvdy + sp->v0;
      drp_ap(hh,vv);
      drp_ep();
   } else {
      h1 = (x1-sp->x0)*sp->dhdx + sp->h0;
      v1 = (y1-sp->y0)*sp->dvdy + sp->v0;
      h2 = (x2-sp->x0)*sp->dhdx + sp->h0;
      v2 = (y2-sp->y0)*sp->dvdy + sp->v0;
      drp_bp();
      drp_ap(h1,v1);
      drp_ap(h2,v2);
      drp_ep();
      drp_bp();
      drp_ap(h2,v1);
      drp_ap(h1,v2);
      drp_ep();
   }
}
void vis_2d_curv_colour_map(nc, xe, sa, qa, qamin, qamax, ios)
   int nc[],sa[],*ios;
   float *xe[],qa[],qamin,qamax;
{
   int ix,nx,iy,ny,k;
   float xl,xr,yl,yr,s,dpdq,z,xc[MAXD];

/* short-hand */
   nx = nc[0];
   ny = nc[1];
   dpdq = 1.0/(qamax - qamin);
/* The frozen coordinate */
   z = vis_gvp->xo[0];
/* arrange an LD image */
   dr3_coerce_ld(ios);
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
      if (wire_frame) {
         dr3_seta("wire_frame");
      }
      dr3_seta("nohard_surface");
      dr3_seta("lucid_surface");
      for (iy=0; iy < ny; iy++) {
         yl = xe[1][iy];
         yr = xe[1][iy+1];
         for (ix=0; ix < nx; ix++) {
            xl = xe[0][ix];
            xr = xe[0][ix+1];
            k = ix + nx*iy;
            if (sa[k] == FLUID) {
               dr3_begin_polygon();
               s = (qa[k]-qamin)*dpdq;
               dr3_polygon_scalar(s);
/*            Bottom left vertex */
               vis_expand(vis_gvp->kad, xl, yl, z, xc);
               dr3_polygon_vertex(xc[0], xc[1], xc[2]);
/*            Bottom right vertex */
               vis_expand(vis_gvp->kad, xr, yl, z, xc);
               dr3_polygon_vertex(xc[0], xc[1], xc[2]);
/*            Top right vertex */
               vis_expand(vis_gvp->kad, xr, yr, z, xc);
               dr3_polygon_vertex(xc[0], xc[1], xc[2]);
/*            Top left vertex */
               vis_expand(vis_gvp->kad, xl, yr, z, xc);
               dr3_polygon_vertex(xc[0], xc[1], xc[2]);
               dr3_end_polygon();
            }
         }
      }
      dr3_seta("hard_surface");
      dr3_seta("nolucid_surface");
      dr3_seta("nowire_frame");
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         vis_ld_scales();
      }
   }
}

void vis_2d_curv_grid_map(nc, xe, sa, ios)
   int nc[],sa[],*ios;
   float *xe[];
{
   int ix,nx,iy,ny,k;
   float xl,xr,yl,yr,z,xc[MAXD];

/* short-hand */
   nx = nc[0];
   ny = nc[1];
/* The frozen coordinate */
   z = vis_gvp->xo[0];
/* arrange an LD image */
   dr3_coerce_ld(ios);
   if (*ios) {
      printf("error:  unable to obtain image\n");
   } else {
      if (wire_frame) {
         dr3_seta("wire_frame");
      }
      dr3_seta("nohard_surface");
      dr3_seta("lucid_surface");
      for (iy=0; iy < ny; iy++) {
         yl = xe[1][iy];
         yr = xe[1][iy+1];
         for (ix=0; ix < nx; ix++) {
            xl = xe[0][ix];
            xr = xe[0][ix+1];
            k = ix + nx*iy;
            if (sa[k] == FLUID) {
               dr3_begin_polyline();
/*            Bottom left vertex */
               vis_expand(vis_gvp->kad, xl, yl, z, xc);
               dr3_polyline_vertex(xc[0], xc[1], xc[2]);
/*            Bottom right vertex */
               vis_expand(vis_gvp->kad, xr, yl, z, xc);
               dr3_polyline_vertex(xc[0], xc[1], xc[2]);
/*            Top right vertex */
               vis_expand(vis_gvp->kad, xr, yr, z, xc);
               dr3_polyline_vertex(xc[0], xc[1], xc[2]);
/*            Top left vertex */
               vis_expand(vis_gvp->kad, xl, yr, z, xc);
               dr3_polyline_vertex(xc[0], xc[1], xc[2]);
/*            Bottom left vertex */
               vis_expand(vis_gvp->kad, xl, yl, z, xc);
               dr3_polyline_vertex(xc[0], xc[1], xc[2]);
               dr3_end_polyline();
            }
         }
      }
      dr3_seta("hard_surface");
      dr3_seta("nolucid_surface");
      dr3_seta("nowire_frame");
      dr3_get_shade(scp);
      drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
      if (!sp->noscales) {
         vis_ld_scales();
      }
   }
}
void vis_expand(kad,z0, z1, z2, xc)
   int kad[MAXD];
   float z0,z1,z2,xc[MAXD];
{
   float z[MAXD],zm[MAXD];
   int id,jd;

/* masked curv coordinates (0:2) */
   zm[0] = z0;
   zm[1] = z1;
   zm[2] = z2;
/* put curv coordinates into proper order */
   for (id=0; id < MAXD; id++) {
      jd = kad[id];
      z[jd] = zm[id];
   }
   geo_x(3, z, xc);
#ifdef VERBOSE
   printf("vis_expand:  geo_x z=");
   for (id=0; id < MAXD; id++) {
      printf(" %f ", z[id]);
   }
   printf(" -> xc=");
   for (id=0; id < MAXD; id++) {
      printf(" %f ", xc[id]);
   }
   printf("\n");
#endif
}
#endif

#ifdef CURV98CONT
void vis_2d_curv_contour_map(nx, ny, xa, ya, sa, qa, qamin, qamax, qlabel,
  qunits, ios)
   int nx,ny,*ios,sa[];
   float xa[],ya[],qa[],qamin,qamax;
   char qlabel[],qunits[];
{
   int i;
   float dqdi,dvdq,h,v,vtop;
   char str[SZLIN];

/* design contours */
   if (sp->nocont) {
/*   establish scalar range if none */
      if (sp->noqr) {
         if (noqr) {
            sp->qmin = qamin;
            sp->qmax = qamax;
         } else {
            sp->qmin = qmin;
            sp->qmax = qmax;
         }
         sp->noqr = sp->qmin >= sp->qmax;
      }
/*   number of contours */
      sp->ncontour = ncontour;
/*   user-supplied contour limits */
      if (!nocontour_levels) {
         for (i=0; i < ncontour; i++) {
            sp->contour_levels[i] = contour_levels[i];
         }
         sp->nocont = FALSE;
      } else if (sp->noqr) {
         printf("error:  no scalar range to determine contours\n");
         *ios = ERR;
      } else {
         if (ncontour == 1) {
            sp->contour_levels[0] = 0.5*(sp->qmin + sp->qmax);
         } else {
            dqdi = (sp->qmax - sp->qmin)/(float)(ncontour+1);
            for (i=0; i < ncontour; i++) {
               sp->contour_levels[i] = (float)(i+1)*dqdi + sp->qmin;
            }
         }
         sp->nocont = FALSE;
      }
   }
/* continue if we really have some contour levels */
   if (sp->nocont) {
      printf("error:  no contour levels\n");
      *ios = ERR;
   } else {
/*    Draw scales if not done already */
      if (!sp->noscales) {
         if (!sp->scales_drawn) {
            ccopy(qlabel, 0, SZNAM, sp->qlabel);
            ccopy(qunits, 0, SZNAM, sp->qunits);
            vis_scales();
         }
/*      draw contour levels? */
         if (!sp->contour_scales_drawn) {
            drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
            dvdq = (sp->vsmax-sp->vsmin)/(sp->qmax-sp->qmin);
            vtop = sp->vmin;
            for (i=0; i < sp->ncontour; i++) {
               if (sp->contour_levels[i] < sp->qmin) {
               } else if (sp->contour_levels[i] > sp->qmax) {
               } else {
                  v = (sp->contour_levels[i]-sp->qmin)*dvdq + sp->vsmin;
/*               check for crowding */
                  if (v-wp->char_height > vtop) {
/*                  value on right */
                     h = sp->hsmax + 0.5*wp->char_width;
                     cfloat(sp->contour_levels[i]/sp->qscale,sp->qfield,
                       sp->qdp, SZLIN, str);
                     drp_text(h, v, str, "lm");
/*                  index on left */
                     sprintf(str, "%d", i);
                     h = sp->hsmin - 0.5*wp->char_width;
                     drp_text(h, v, str, "rm");
                     vtop = v + 0.5*wp->char_height;
                  }
/*               contour level as a line */
                  drp_bp();
                  drp_ap(sp->hsmin, v);
                  drp_ap(sp->hsmax, v);
                  drp_ep();
               }
            }
            sp->contour_scales_drawn = TRUE;
         }
      }
/*   draw contours */
      for (i=0; i < sp->ncontour; i++) {
         if (nocontour_label) {
            sprintf(str, "%c", '\0');
         } else {
            sprintf(str, "%d", i);
         }
         drp_set_hv(sp->hdmin, sp->hdmax, sp->vdmin, sp->vdmax);
/*       arrange an LD image */
         dr3_coerce_ld(ios);
         if (*ios) {
            printf("error:  unable to obtain image\n");
         } else {
            if (wire_frame) {
               dr3_seta("wire_frame");
            }
            dr3_seta("nohard_surface");
            dr3_seta("lucid_surface");
            contour_map(nx, ny, xa, ya, sa, qa, 
              sp->contour_levels[i], str, ios);
            dr3_seta("hard_surface");
            dr3_seta("nolucid_surface");
            dr3_seta("nowire_frame");
            dr3_get_shade(scp);
            drp_set_hv(sp->hxmin, sp->hxmax, sp->vxmin, sp->vxmax);
            if (!sp->noscales) {
               vis_ld_scales();
            }
         }
      }
   }
}
#endif
