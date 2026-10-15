/* dr3.h:  dr3 module definitions */
#define MAXCD 9 /* colour model dimensions */
#define OLD
#ifdef OLD
#define MAXSEG 8192 /* polyline segments (was 4096) */
#define MAXPOL 65536 /* polyline points */
#else
#define MAXSEG 4096 /* polyline segments (was 4096) */
#define MAXPOL 32768 /* polyline points */
#endif
#define MAXMAT 32 /* materials */

#define SHADE_WIRE 1 /* wire-frame only */
#define SHADE_LM 2 /* light-material shading */
#define SHADE_LD 3 /* light-data shading */
#define SHADE_L 4 /* light shading */
#define SHADE_D 5 /* data shading only */

#define DEV_RGBZ 4 /* RGB+zbuffer class */
#define DEV_RGB 3 /* RGB class */
#define DEV_CM 2 /* ColourMap class */
#define DEV_VEC 1 /* VECtor class */

/* SceneStructure */
typedef struct scene_type {
/* Viewing parameters */
   float view_azim,view_altd,view_twist,view_distance,
     view_focus[MAXD4],view_point[MAXD4],view_matrix[MAXD4][MAXD4],
     view_world[MAXD4][MAXD4],view_nx[MAXD],view_ny[MAXD],
     view_nz[MAXD],view_w0[MAXD],view_near,view_far;
   bool noview;
/* Projection parameters */
   float project_anglev,project_near,project_far,project_hmin,
     project_hmax,project_vmin,project_vmax,project_dhdxz,
     project_h0,project_dvdyz,project_v0,project_dddz,project_d0;
   bool noproject;
/* Lighting parameters */
   float light_azim,light_altd,light_distance,light_point[MAXD4],
     light_ambient,light_candle,light_view[MAXD4];
   bool nolight;
/* Shading model parameters */
   float shade_red[MAXCM],shade_green[MAXCM],shade_blue[MAXCM],
     shade_didc[MAXCD],sky_red,sky_green,sky_blue,
     shade_mat_r[MAXMAT],shade_mat_g[MAXMAT],shade_mat_b[MAXMAT];
   int shade_index,shade_nd,shade_nc[MAXCD],shade_ncm,shade_offset;
   char shade_name[SZNAM],shade_mat[MAXMAT][SZNAM];
   bool noshade;
/* Image data */
   float dhidh,hi0,dvidv,vi0,*dim,*rgbim;
   int ih0,iv0,nh,nv,*cim,dev_class;
   bool vector_bitmap,colour_bitmap,rgb_bitmap,rgbz_bitmap,noimage,
     image_owned;
/* Polyline data */
   int np,ip0[MAXSEG],nps[MAXSEG],ns;
   float hp[MAXPOL],vp[MAXPOL],dp[MAXPOL];
/* Scene control */
   int scene_index;
} scene_struct;
