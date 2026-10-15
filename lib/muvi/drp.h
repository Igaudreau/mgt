/* drp.h:  definitions/declarations common to all backends */
#include "const.h"
#define MAXRAMP 256
#define MAXPAINT 16
#define MAXCM 4096
struct dev {
   int (*prepare)();
   int (*open)();
   int (*close)();
   int (*select)();
   int (*erase)();
   int (*rub)();
   int (*rect)();
   int (*polyline)();
   int (*polygon)();
   int (*text)();
   int (*text_box)();
   int (*cursor)();
   int (*sample)();
   int (*set_font)();
   int (*set_colours)();
   int (*set_colour)();
   int (*set_rgb)();
   int (*set_rgb_background)();
   int (*set_line_width)();
   int (*colour_bitmap)();
   int (*rgb_bitmap)();
   int (*rgbz_setup)();
   int (*rgbz_begin)();
   int (*rgbz_end)();
   int (*rgbz_clear)();
   int (*rgbz_begin_polygon)();
   int (*rgbz_polygon_vertex)();
   int (*rgbz_end_polygon)();
   int (*rgbz_begin_polyline)();
   int (*rgbz_polyline_vertex)();
   int (*rgbz_end_polyline)();
   int (*rgbz_put_bitmap)();
   int (*rgbz_get_bitmap)();
   int (*gflush)();
   char name[SZNAM],type[SZNAM],font_style[SZNAM],font_name[SZNAM],
     ramp_names[MAXRAMP][SZNAM],paint_names[MAXPAINT][SZNAM];
   float hmin,hmax,vmin,vmax,
     red,green,blue,red_background,green_background,blue_background,
     reds[MAXCOLOUR],greens[MAXCOLOUR],blues[MAXCOLOUR],
     font_size,char_height,char_width,char_depth,char_space,
     line_width,bitmap_dh,bitmap_dv;
   bool used,avail_rect,avail_polygon,avail_colour,avail_rgb,
     change_colour,avail_cursor,avail_bitmap,avail_sample,avail_rub,
     paint_used[MAXPAINT],ramp_used[MAXRAMP],avail_fontscale,
     avail_rgbz;
   int colour_index,colour_offset,colour_size,cursor_buttons,
     ihmin,ihmax,ivmin,ivmax, ramp_count,ramp_top,ramp_offset[MAXRAMP],
     ramp_nd[MAXRAMP],ramp_nc[MAXRAMP][MAXD],ramp_size[MAXRAMP],
     paint_count,paint_top,paint_offset[MAXPAINT],count;
};
