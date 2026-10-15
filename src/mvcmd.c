/*
 *Copyright
 *
 * (c) 1989 mantis numerics limited
 *     46 the calls, leeds ls2 7ey, uk.
 */ 

#include <string.h>
#include <stdio.h>
#define MAXCON 50
#define SZNAM 1000
#ifdef DOUBLE
#define float double
#endif

extern bool get_r(float &value);
extern bool get_i(int &value);
extern bool get_n(char key[]);
extern void report();
extern "C"{
  void model_seta(char name[], int *ios);
  void model_set0r(char name[], float value, int *ios);
  void model_set0i(char name[], int value, int *ios);
  void model_set0s(char name[], char value[], int *ios);
  void vis_print();
  void vis_show();
  void vis_seta(char name[], int *ios);
  void vis_set0r(char name[], float value, int *ios);
  void vis_set0i(char name[], int value, int *ios);
  void vis_open(char name[], char type[], int *ios);
  void vis_select(char name[], int *ios);
  void vis_close(char name[], int *ios);
  void vis_open_sub(char name[], float  rhmin, float rhmax,
		    float rvmin, float rvmax, int *ios);
  void vis_close_sub(char name[], int *ios); 
  void vis_erase();
  void vis_clear();
  void vis_select_sub(char name[], int *ios);
  void vis_set0s(char name[], char value[], int *ios);
  void vis_set1r(char name[], int nvalue, float values[],
		 int *ios);
  void  parc_rank(int *mype, int *limpe, int *ios);
}

//file mvcmd.c
/*
  + mvcmd - Muvi command interpreter.
  Description
  This  function checks  to see  if  it can  implement the  command
  supplied in  the key  string, and if  so it reads  any additional
  parameters and invokes the appropriate muvi routine.
*/
bool mvcmd(char key[])
{
  const int maxcon = MAXCON;
  int ios, ix, i, nline, nchar, mark, istyle, n, mype, limpe;
  float x, h, v, rhmin, rhmax, rvmin, rvmax, bord, red, green, blue, r,
    qmin, qmax, asp, xmin, xmax, ymin, ymax, zmin, zmax, cont[MAXCON],
    y, z, dx, dy, dz, xc[3], phong, shine;
  char tname[SZNAM], tlabel[SZNAM], name[SZNAM], type[SZNAM],
    label[SZNAM];
  bool cmd;
  /* Get process rank */
  ios = 0;
  parc_rank(&mype, &limpe, &ios);
  cmd = true;
  ios = 0;
  /* unmask coordinate */
  if (strcmp(key, "UNMASK") == 0){
    if (!get_n(tname)){
      printf("Error:  reading coord name. \n");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "unmask/");
	strcat(name, tname);
	model_seta(name, &ios);
      }
    }
  }
  /* mask coordinate by nearest value */
  else if (strcmp(key, "MASK") == 0){
    if (!get_n(tname)){
      printf("Error:  reading coord name.");
      report();
    }
    else if (!get_r(x)){
      printf("Error:  reading mask coordinate value.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "mask/");
	strcat(name, tname);
	model_set0r(name, x, &ios);
      }
    }
  }
  /* mask coordinate by index */
  else if (strcmp(key, "MASKI") == 0){
    if (!get_n(tname)){
      printf("Error:  reading coord name.");
      report();
    }
    else if (!get_i(ix)){
      printf("Error:  reading mask coordinate index.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "maski/");
	strcat(name, tname);
	model_set0i(name, ix, &ios);
      }
    }
  }
  /* vprint - print general status of vis */
  else if (strcmp(key, "VPRINT") == 0){
    if (mype == 0)
      vis_print();
  }
  /* vshow - show details of current vis worksurface/window */
  else if (strcmp(key, "VSHOW") == 0){
    if (mype == 0)
      vis_show();
  }
  /* paper hinch vinch - specify paper size (inches) */
  else if (strcmp(key, "PAPER") == 0){
    if (!get_r(h)){
      printf("Error:  reading worksurface width.");
      report();
    }
    else if (!get_r(v)){
      printf("Error:  reading worksurface height.");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"size_ah", h, &ios);
	vis_set0r((char*)"size_av", v, &ios);
      }
    }
  }
  /* rpaper hrel vrel - specify paper size (normalised) */
  else if (strcmp(key, "RPAPER") == 0){
    if (!get_r(h)){
      printf("Error:  reading worksurface width.");
      report();
    }
    else if (!get_r(v)){
      printf("Error:  reading worksurface height.");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"size_rh", h, &ios);
	vis_set0r((char*)"size_rv", v, &ios);
      }
    }
  }
  /* position hinch vinch - specify paper position (inches) */
  else if (strcmp(key, "POSITION") == 0){
    if (!get_r(h)){
      printf("Error:  reading worksurface h-position.");
      report();
    }
    else if (!get_r(v)){
      printf("Error:  reading worksurface v-position.");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"pos_ah", h, &ios);
	vis_set0r((char*)"pos_av", v, &ios);
      }
    }
  }
  /* rposition hrel vrel - specify paper position (normalised) */
  else if (strcmp(key, "RPOSITION") == 0){
    if (!get_r(h)){
      printf("Error:  reading worksurface h-position.");
      report();
    }
    else if (!get_r(v)){
      printf("Error:  reading worksurface v-position.");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"pos_rh", h, &ios);
	vis_set0r((char*)"pos_rv", v, &ios);
      }
    }
  }
  /* open name type - open worksurface */
  else if (strcmp(key, "OPEN") == 0){
    if (!get_n(tname)){
      printf("Error:  reading worksurface name.");
      report();
    }
    else{
      strcpy(name, tname);
      if (!get_n(tname)){
	printf("Error:  reading worksurface type.");
	report();
      }
      else if (mype ==0){
	strcpy(type, tname);
	vis_open(name, type, &ios);
      }
      if (ios != 0){
	printf("Error:  opening worksurface.");
	report();
      }
    }
  }
  /* select name - select already open graphics device as active */
  else if (strcmp(key, "SELECT") == 0){
    if (!get_n(tname)){
      printf("Error:  reading worksurface name.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, tname);
	vis_select(name,&ios);
	if (ios != 0){
	  printf("Error:  unable to select worksurface.");
	  report();
	}
      }
    }
  }
  /* close name - close graphics device */
  else if (strcmp(key, "CLOSE") == 0){
    if (!get_n(tname)){
      printf("Error:  reading worksurface name.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, tname);
	vis_close(name, &ios);
	if (ios != 0){
	  printf("Error:  closing worksurface.");
	  report();
	}
      }
    }
  }
  /* openw name type - open subwindow */
  else if (strcmp(key, "OPENW") == 0){
    if (!get_n(tname)){
      printf("Error:  reading subwindow name.");
      report();
    }
    else if (!get_r(rhmin)){
      printf("Error:  reading subwindow left position.");
      report();
    }
    else if (!get_r(rhmax)){
      printf("Error:  reading subwindow right position.");
      report();
    }
    else if (!get_r(rvmin)){
      printf("Error:  reading subwindow bottom position.");
      report();
    }
    else if (!get_r(rvmax)){
      printf("Error:  reading subwindow top position.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, tname);
	vis_open_sub(name, rhmin, rhmax, rvmin, rvmax, &ios);
	if (ios != 0){
	  printf("Error:  opening subwindow.");
	  report();
	}
      }
    }
  }
  /* selectw name - select subwindow */
  else if (strcmp(key, "SELECTW") == 0){
    if (!get_n(tname)){
      printf("Error:  reading subwindow name.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, tname);
	vis_select_sub(name,&ios);
	if (ios != 0){
	  printf("Error:  unable to select that subwindow.");
	  report();
	}
      }
    }
  }
  /* closew name - close subwindow */
  else if (strcmp(key, "CLOSEW") == 0){
    if (!get_n(tname)){
      printf("Error:  reading subwindow name.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, tname);
	vis_close_sub(name, &ios);
	if (ios != 0){
	  printf("Error:  closing subwindow.");
	  report();
	}
      }
    }
  }
  /* er - clear graphics subwindow */
  else if (strcmp(key,"ER") == 0){
    if (mype == 0)
      vis_erase();
  }
  /* clear - clear graphics worksurface */
  else if (strcmp(key, "CLEAR") == 0){
    if (mype == 0)
      vis_clear();
  }
  /* [no]scales - whether contour/colour scales are drawn */
  else if (strcmp(key, "NOSCALES") == 0){
    if (mype == 0)
      vis_seta((char*)"noscales", &ios);
  }
  else if (strcmp(key, "SCALES") == 0){
    if (mype == 0)
      vis_seta((char*)"scales", &ios);
  }
  /* [no]key - whether key is drawn */
  else if (strcmp(key, "NOKEY") == 0){
    if (mype == 0)
      vis_seta((char*)"nokey", &ios);
  }
  else if (strcmp(key, "KEY") == 0){
    if (mype == 0)
      vis_seta((char*)"key", &ios);
  }
  /* key_size size - specify number of key elements */
  else if (strcmp(key, "KEY_SIZE") == 0){
    if (!get_i(i)){
      printf("Error:  reading key size.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"key_size", i, &ios);
    }
  }
  /* key_label_field size - specify key label field size */
  else if (strcmp(key, "KEY_LABEL_FIELD") == 0){
    if (!get_i(i)){
      printf("Error:  reading key label field size.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"key_label_field", i, &ios);
    }
  }
  /* key_number_field size - specify key number field size */
  else if (strcmp(key, "KEY_NUMBER_FIELD") == 0){
    if (!get_i(i)){
      printf("Error:  reading key number field size.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"key_number_field", i, &ios);
    }
  }

  /* key_label label - specify current lkey label */
  else if (strcmp(key, "KEY_LABEL") == 0){
    if (!get_n(tname)){
      printf("Error:  reading key label.");
      report();
    }
    else{
      if (mype == 0)
	strcpy(label, tname);
      vis_set0s((char*)"key_label", label, &ios);
    }
  }
  /* nokey_label */
  else if (strcmp(key, "NOKEY_LABEL") == 0){
    if (mype == 0)
      vis_seta((char*)"nokey_label", &ios);
  }
  /* [no]box - whether axes are drawn afresh each time */
  else if (strcmp(key, "NOBOX") == 0){
    if (mype == 0)
      vis_seta((char*)"nobox", &ios);
  }
  else if (strcmp(key, "BOX") == 0){
    if (mype == 0)
      vis_seta((char*)"box", &ios);
  }
  /* [no]qlabel [label] - set q-label string */
  else if (strcmp(key, "NOQLABEL") == 0){
    if (mype == 0)
      vis_seta((char*)"noqlabel", &ios);
  }
  else if (strcmp(key, "QLABEL") == 0){
    if (!get_n(label)){
      printf("Error:  reading qlabel");
      report();
    }
    else{
      if (mype == 0)
	vis_set0s((char*)"qlabel", label, &ios);
    }
  }
  /*  [no]xlabel [label] - set x-label string */
  else if (strcmp(key, "NOXLABEL") ==0){
    if (mype == 0)
      vis_seta((char*)"noxlabel", &ios);
  }
  else if (strcmp(key, "XLABEL") == 0){
    if (!get_n(label)){
      printf("Error:  reading xlabel");
      report();
    }
    else{
      if (mype == 0)
	vis_set0s((char*)"xlabel", label, &ios);
    }
  }
  /* [no]ylabel [label] - set y-label string */
  else if (strcmp(key, "NOYLABEL") == 0){
    if (mype == 0)
      vis_seta((char*)"noylabel", &ios);
  }
  else if (strcmp(key, "YLABEL") == 0){
    if (!get_n(label)){
      printf("Error:  reading ylabel");
      report();
    }
    else{
      if (mype == 0)
	vis_set0s((char*)"ylabel", label, &ios);
    }
  }
  /* [no]zlabel [label] - set z-label string */
  else if (strcmp(key, "NOZLABEL") == 0){
    if (mype == 0)
      vis_seta((char*)"nozlabel", &ios);
  }
  else if (strcmp(key, "ZLABEL") == 0){
    if (!get_n(label)){
      printf("Error:  reading zlabel");
      report();
    }
    else{
      if (mype == 0)
	vis_set0s((char*)"zlabel", label, &ios);
    }
  }
  /* [no]caption [caption] - set caption string */
  else if (strcmp(key, "NOCAPTION") == 0){
    if (mype == 0)
      vis_seta((char*)"nocaption", &ios);
  }
  else if (strcmp(key, "CAPTION") == 0){
    if (!get_n(label)){
      printf("Error:  reading caption");
      report();
    }
    else{
      if (mype == 0)
	vis_set0s((char*)"caption", label, &ios);
    }
  }
  /* caption_lines nlines - number of caption lines reserved */
  else if (strcmp(key, "CAPTION_LINES") == 0){
    if (!get_i(nline)){
      printf("Error:  reading caption_lines.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"caption_lines", nline, &ios);
    }
  }
  /* caption_inset nchar - number of indented caption chars */
  else if (strcmp(key, "CAPTION_INSET") == 0){
    if (!get_i(nchar)){
      printf("Error:  reading caption_inset.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"caption_inset", nchar, &ios);
    }
  }
  /* label_field nchar - max. number of characters in label fields */
  else if (strcmp(key, "LABEL_FIELD") == 0){
    if (!get_i(nchar)){
      printf("Error:  reading label_field.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"label_field", nchar, &ios);
    }
  }
  /* number_field nchar - max. number of characters in number fields */
  else if (strcmp(key, "NUMBER_FIELD") == 0){
    if (!get_i(nchar)){
      printf("Error:  reading number_field.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"number_field", nchar, &ios);
    }
  }
  /* caption_below/above - where is caption box? */
  else if (strcmp(key, "CAPTION_BELOW") == 0){
    if (mype == 0)
      vis_seta((char*)"caption_below", &ios);
  }
  else if (strcmp(key, "CAPTION_ABOVE") ==0){
    if (mype == 0)
      vis_seta((char*)"caption_above", &ios);
  }
  /* tick_scale scale - set tick size in charheight units */
  else if (strcmp(key, "TICK_SCALE") == 0){
    if (!get_r(bord)){
      printf("Error;:  reading tick_scale.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"tick_scale", bord, &ios);
    }
  }
  /* vscale scale - set arrow scale factor */
  else if (strcmp(key, "VSCALE") == 0){
    if (!get_r(bord)){
      printf("Error:  reading arrow_scale.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"arrow_scale", bord, &ios);
    }
  }
  /* vratio ratio - set arrow ratio for scatter grams */
  else if (strcmp(key, "VRATIO") == 0){
    if (!get_r(bord)){
      printf("Error:  reading arrow_ratio.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"arrow_ratio", bord, &ios);
    }
  }
  /* vrate rate - set arrow rate */
  else if (strcmp(key, "VRATE") == 0){
    if (!get_i(i)){
      printf("Error:  reading arrow_rate.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"arrow_rate", i, &ios);
    }
  }
  /* novmark - revoke markers on arrow plots */
  else if (strcmp(key, "NOVMARK") == 0){
    if (mype == 0)
      vis_seta((char*)"noarrow_marker", &ios);
  }
  /*  vmark marker - set arrow marker */
  else if (strcmp(key, "VMARK") == 0){
    if (!get_i(mark)){
      printf("Error:  reading arrow_marker.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"arrow_marker", mark, &ios);
    }
  }
  /* font name size - set font name and size */
  else if (strcmp(key, "FONT") == 0){
    if (!get_n(name)){
      printf("Error:  reading font.");
      report();
    }
    else{
      if (!get_r(bord)){
	printf("Error:  reading font_size.");
	report();
      }
      else{
	if (mype == 0){
	  vis_set0s((char*)"font_style", name, &ios);
	  vis_set0r((char*)"font_size", bord, &ios);
	  vis_seta((char*)"font", &ios);
	}
      }
    }
  }
  /* mark_scale scale - set marker size in charheight units */
  else if (strcmp(key, "MARK_SCALE") == 0){
    if (!get_r(bord)){
      printf("Error:  reading marker_scale.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"marker_scale", bord, &ios);
    }
  }
  /* marker marker - set marker */
  else if (strcmp(key, "MARKER") == 0){
    if (!get_i(mark)){
      printf("Error:  reading marker.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"line_marker", mark, &ios);
    }
  }
  /* nomarker - revoke marker */
  else if (strcmp(key, "NOMARKER") == 0){
    if (mype == 0)
      vis_seta((char*)"noline_marker", &ios);
  }
  /* line style - set line style */
  else if (strcmp(key, "LINE") == 0){
    if (!get_i(istyle)){
      printf("Error:  reading style.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"line_style", istyle, &ios);
    }
  }
  /* line_palet palet - set line palet */
  else if (strcmp(key, "LINE_PALET") == 0){
    if (!get_i(istyle)){
      printf("Error:  reading palet index.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"line_palet", istyle-1, &ios);
    }
  }
  /* noline - revoke line connects */
  else if (strcmp(key, "NOLINE") == 0){
    if (mype == 0)
      vis_seta((char*)"noline_connect", &ios);
  }
  /* palet name red green blue - define palet */
  else if (strcmp(key, "PALET") ==0){
    if (!get_n(tname)){
      printf("Error:  reading palet name.");
      report();
    }
    else if (!get_r(red)){
      printf("Error:  reading palet red");
      report();
    }
    else if (!get_r(green)){
      printf("Error:  reading palet green");
      report();
    }
    else if (!get_r(blue)){
      printf("Error:  reading palet blue");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name,"palet_red/");
	strcat(name, tname);
	vis_set0r(name, red, &ios);
	strcpy(name,"palet_green/");
	strcat(name, tname);
	vis_set0r(name, green, &ios);
	strcpy(name,"palet_blue/");
	strcat(name, tname);
	vis_set0r(name, red, &ios);
      }
    }
  }
  /* hborder fraction - set fractional border in subwindow */
  else if (strcmp(key, "HBORDER") == 0){
    if (!get_r(bord)){
      printf("Error:  reading hborder.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"hborder", bord, &ios);
    }
  }
  /* vborder fraction - set fractional border in subwindow */
  else if (strcmp(key, "VBORDER") ==0){
    if (!get_r(bord)){
      printf("Error:  reading vborder."); 
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"vborder", bord, &ios);
    }
  }
  /* novmax - quantity range */
  else if (strcmp(key, "NOVMAX") == 0){
    if (mype == 0)
      vis_seta((char*)"novmax", &ios);
  }
  /* vmax vmax - set velocity range */
  else if (strcmp(key, "VMAX") == 0){
    if (!get_r(r)){
      printf("Error:  reading vmax.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"vmax", r, &ios);
    }
  }
  /* noqr - quantity range */
  else if (strcmp(key, "NOQR") == 0){
    if (mype == 0)
      vis_seta((char*)"noqr", &ios);
  }
  /* qr qmin qmax - set quantity range */
  else if (strcmp(key, "QR") == 0){
    if (!get_r(qmin)){
      printf("Error:  reading qmin.");
      report();
    }
    else if (!get_r(qmax)){
      printf("Error:  reading qmax.");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"qmin", qmin, &ios);
	vis_set0r((char*)"qmax", qmax, &ios);
      }
    }
  }
  /* xaspect asp - make x-axis longer by xaspect times */
  else if (strcmp(key, "XASPECT") == 0){
    if (!get_r(asp)){
      printf("Error:  reading xaspect");
      report();
    }
    else if (asp <= 0.0){
      printf("Error:  xaspect invalid.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"xaspect", asp, &ios);
    }
  }
  /* yaspect asp - make y-axis longer by xaspect times */
  else if (strcmp(key, "YASPECT") == 0){
    if (!get_r(asp)){
      printf("Error:  reading yaspect");
      report();
    }
    else if (asp <= 0.0){
      printf("Error:  yaspect invalid.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"yaspect", asp, &ios);
    }
  }
  /* zaspect asp - make z-axis longer by xaspect times */
  else if (strcmp(key, "ZASPECT") == 0){
    if (!get_r(asp)){
      printf("Error:  reading zaspect");
      report();
    }
    else if (asp <= 0.0){
      printf("Error:  zaspect invalid.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"zaspect", asp, &ios);
    }
  }
  /* noxr - no x-range */
  else if (strcmp(key, "NOXR") == 0){
    if (mype == 0)
      vis_seta((char*)"noxr", &ios);
  }
  /* xr xmin xmax - set x-range */
  else if (strcmp(key, "XR") == 0){
    if (!get_r(xmin)){
      printf("Error:  reading xmin");
      report();
    }
    else if (!get_r(xmax)){
      printf("Error:  reading xmax");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"xmin", xmin, &ios);
	vis_set0r((char*)"xmax", xmax, &ios);
      }
    }
  }
  /* noyr - y-range */
  else if (strcmp(key, "NOYR") == 0){
    if (mype == 0)
      vis_seta((char*)"noyr", &ios);
  }
  /* yr ymin ymax - set y-range */
  else if (strcmp(key, "YR") == 0){
    if (!get_r(ymin)){
      printf("Error:  reading ymin");
      report();
    }
    else if (!get_r(ymax)){
      printf("Error:  reading ymax");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"ymin", ymin, &ios);
	vis_set0r((char*)"ymax", ymax, &ios);
      }
    }
  }
  /* nozr - z-range */
  else if (strcmp(key, "NOZR") == 0){
    if (mype == 0)
      vis_seta((char*)"nozr", &ios);
  }
  /* zr zmin zmax - set z-range */
  else if (strcmp(key, "zr") == 0){
    if (!get_r(zmin)){
      printf("Error:  reading zmin");
      report();
    }
    else if (!get_r(zmax)){
      printf("Error:  reading zmax");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"zmin", zmin, &ios);
	vis_set0r((char*)"zmax", zmax, &ios);
      }
    }
  }
  /* [no]scalp - whether scalar is promoted to coordinate */
  else if (strcmp(key, "NOSCALP") == 0){
    if (mype == 0)
      vis_seta((char*)"noscalp", &ios);
  }
  else if (strcmp(key, "SCALP") == 0){
    if (mype == 0)
      vis_seta((char*)"scalp", &ios);
  }
  /* [no]log - whether scalars are converted to log10 */
  else if (strcmp(key, "LOG") == 0){
    if (mype == 0)
      model_seta((char*)"log", &ios);
  }
  else if (strcmp(key, "NOLOG") == 0){
    if (mype == 0)
      model_seta((char*)"nolog", &ios);
  }
  /* [no]band - whether colour band (range) */
  else if (strcmp(key, "NOBAND") == 0){
    if (mype == 0)
      vis_seta((char*)"noband", &ios);
  }
  else if (strcmp(key, "BAND") == 0){
    if (mype == 0)
      vis_seta((char*)"band", &ios);
  }
  /* bandr index min max - colour range */
  else if (strcmp(key, "BANDR") == 0){
    if (!get_n(tname)){
      printf("Error:  reading band index");
      report();
    }
    else if (!get_r(qmin)){
      printf("Error:  reading band min");
      report();
    }
    else if (!get_r(qmax)){
      printf("Error:  reading band max");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "band_min/");
	strcat(name, tname);
	vis_set0r(name, qmin, &ios);
	strcpy(name, "band_max/");
	strcat(name, tname);
	vis_set0r(name, qmax, &ios);
      }
    }
  }
  /* nobandr index - colour range */
  else if (strcmp(key, "NOBANDR") == 0){
    if (!get_n(tname)){
      printf("Error:  reading band index");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "nobandr/");
	strcat(name, tname);
	vis_seta(name, &ios);
      }
    }
  }
  /* [no]colour - whether colour map */
  else if (strcmp(key, "NOCOLOUR") == 0){
    if (mype == 0)
      vis_seta((char*)"nocolour", &ios);
  }
  else if (strcmp(key, "COLOUR") == 0){
    if (mype == 0)
      vis_seta((char*)"colour", &ios);
  }
  /* ncolour ncolour - number of colour bands */
  else if (strcmp(key, "NCOLOUR") == 0){
    if (!get_i(i)){
      printf("Error:  reading number of colour bands.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"ncolour", i, &ios);
    }
  }
  /* nshade nshade - number of shade bands */
  else if (strcmp(key, "NSHADE") == 0){
    if (!get_i(i)){
      printf("Error:  reading number of shades.");
      report();
    }
    else
      vis_set0i((char*)"nshade", i, &ios);
  }
  /* data_shade - shade surfaces with data */
  else if (strcmp(key, "DATA_SHADE") == 0){
    if (mype == 0)
      vis_seta((char*)"data_shade", &ios);
  }
  /* data_ortho id - shade ortho surfaces with data */
  else if (strcmp(key, "DATA_ORTHO") == 0){
    if (!get_i(i)){
      printf("Error:  reading ortho direction");
      report();
    }
    else{
      if (mype == 0){
	vis_seta((char*)"data_ortho", &ios);
	vis_set0i((char*)"data_ortho_id", i, &ios);
      }
    }
  }
  /* [no]grey_rain - whether grey rainbow */
  else if (strcmp(key, "NOGREY_RAIN") == 0){
    if (mype == 0)
      vis_seta((char*)"nogrey_rain", &ios);
  }
  else if (strcmp(key, "GREY_RAIN") == 0){
    if (mype == 0)
      vis_seta((char*)"grey_rain", &ios);
  }
  /* [no]contour - whether contour map */
  else if (strcmp(key, "NOCONTOUR") == 0){
    if (mype == 0)
      vis_seta((char*)"nocontour", &ios);
  }
  else if (strcmp(key, "CONTOUR") == 0){
    if (mype == 0)
      vis_seta((char*)"contour", &ios);
  }
  /* ncont ncont - number of contour levels */
  else if (strcmp(key, "NCONT") == 0){
    if (!get_i(i)){
      printf("Error:  reading number of contour levels.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0i((char*)"ncontour", i, &ios);
    }
  }
  /* cont contour_level - single contour level */
  else if (strcmp(key, "CONT") == 0){
    if (!get_r(r)){
      printf("Error:  reading contour level.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"contour_level", r, &ios);
    }
  }
  /* conts ncontour cont1 cont2 ... - set of contour levels */
  else if (strcmp(key, "CONTS") == 0){
    if (!get_i(n)){
      printf("Error:  reading number of contour levels.");
      report();
    }
    else if ((n < 1) || (n > maxcon)){
      printf("Error:  invalid number of contours.");
      report();
    }
    else{
      for (i = 0; i < n; i++){
	if (ios != 0);
	else if (!get_r(cont[i])){
	  printf("Error:  reading contour level %d ", i + 1);
	  report();
	  ios = 1;
	}
      }
      if (mype == 0)
	vis_set1r((char*)"contour_levels", n, cont, &ios);
    }
  }
  /* nocontour_label */
  else if (strcmp(key, "NOCONTOUR_LABEL") == 0){
    if (mype == 0)
      vis_seta((char*)"nocontour_label", &ios);
  }
  /* contour_label */
  else if (strcmp(key, "CONTOUR_LABEL") == 0){
    if (mype == 0)
      vis_seta((char*)"contour_label", &ios);
  }
  /* [no]arrow - whether velocity arrows drawn */
  else if (strcmp(key, "NOARROW") == 0){
    if (mype == 0)
      vis_seta((char*)"noarrow", &ios);
  }
  else if (strcmp(key, "ARROW") == 0){
    if (mype == 0)
      vis_seta((char*)"arrow", &ios);
  }
  /* stream n x y z dx dy dz - velocity streamlines drawn */
  else if (strcmp(key, "STREAM") == 0){
    if (!get_i(i)){
      printf("Error:  reading stream_size.");
      report();
    }
    else if (!get_r(x)){
      printf("Error:  reading xstream.");
      report();
    }
    else if (!get_r(y)){
      printf("Error:  reading ystream.");
      report();
    }
    else if (!get_r(z)){
      printf("Error:  reading zstream.");
      report();
    }
    else if (!get_r(dx)){
      printf("Error:  reading dxstream.");
      report();
    }
    else if (!get_r(dy)){
      printf("Error:  reading dystream.");
      report();
    }
    else if (!get_r(dz)){
      printf("Error:  reading dzstream.");
      report();
    }
    else{
      if (mype == 0){
	vis_seta((char*)"stream", &ios);
	vis_set0i((char*)"stream_size", i, &ios);
	vis_set0r((char*)"xstream", x, &ios);
	vis_set0r((char*)"ystream", y, &ios);
	vis_set0r((char*)"zstream", z, &ios);
	vis_set0r((char*)"dxstream", dx, &ios);
	vis_set0r((char*)"dystream", dy, &ios);
	vis_set0r((char*)"dzstream", dz, &ios);
      }
    }
  }
  /* view_azim view_azim - azimuthal direction of view-point */
  else if (strcmp(key, "VIEW_AZIM") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_azim");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_azim", r, &ios);
    }
  }
  /* view_altd view_altd - altitudinal direction of view-point */
  else if (strcmp(key, "VIEW_ALTD") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_altd");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_altd", r, &ios);
    }
  }
  /* view_anglev view_anglev - vertical viewing angle */
  else if (strcmp(key, "VIEW_ANGLEV") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_anglev");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_anglev", r, &ios);
    }
  }
  /* view_aspect view_aspect - view-surface aspect ratio */
  else if (strcmp(key, "VIEW_ASPECT") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_aspect");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_aspect", r, &ios);
    }
  }
  /* noview_distance - revoke view distance */
  else if (strcmp(key, "NOVIEW_DISTANCE") == 0){
    if (mype == 0)
      vis_seta((char*)"noview_distance", &ios);
  }
  /* view_distance view_distance - view-surface distance */
  else if (strcmp(key, "VIEW_DISTANCE") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_distance");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_distance", r, &ios);
    }
  }
  /* noview_near - revoke near view */
  else if (strcmp(key, "NOVIEW_NEAR") == 0){
    if (mype == 0)
      vis_seta((char*)"noview_near", &ios);
  }
  /* view_near front clipping distance */
  else if (strcmp(key, "VIEW_NEAR") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_near");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_near", r, &ios);
    }
  }
  /* noview_far - revoke far view */
  else if (strcmp(key, "NOVIEW_FAR") == 0){
    if (mype == 0)
      vis_seta((char*)"noview_far", &ios);
  }
  /* view_far back clipping distance */
  else if (strcmp(key, "VIEW_FAR") == 0){
    if (!get_r(r)){
      printf("Error:  reading view_far");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"view_far", r, &ios);
    }
  }
  /* nodata_coord - revoke data coord */
  else if (strcmp(key, "NODATA_COORD") == 0){
    if (mype == 0)
      vis_seta((char*)"nodata_coord", &ios);
  }
  /* data_coord back clipping distance */
  else if (strcmp(key, "DATA_COORD") == 0){
    if (!get_r(r)){
      printf("Error:  reading data_coord");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"data_coord", r, &ios);
    }
  }
  /* noview_focus - revoke it */
  else if (strcmp(key, "noview_focus") == 0){
    if (mype == 0)
      vis_seta((char*)"noview_focus", &ios);
  }
  /* view_focus view_[x,y,z] - view-focus position */
  else if (strcmp(key, "VIEW_FOCUS") == 0){
    if (!get_r(xc[0])){
      printf("Error:  reading view_focus_x");
      report();
    }
    else if (!get_r(xc[1])){
      printf("Error:  reading view_focus_y");
      report();
    }
    else if (!get_r(xc[2])){
      printf("Error:  reading view_focus_z");
      report();
    }
    else{
      if (mype == 0)
	vis_set1r((char*)"view_focus", 3, xc, &ios);
    }
  }
  /* [no]wire_frame */
  else if (strcmp(key, "WIRE_FRAME") == 0){
    if (mype == 0)
      vis_seta((char*)"wire_frame", &ios);
  }
  else if (strcmp(key, "NOWIRE_FRAME") == 0){
    if (mype == 0)
      vis_seta((char*)"nowire_frame", &ios);
  }
  /* [no]data_frame */
  else if (strcmp(key, "DATA_FRAME") == 0){
    if (mype == 0)
      vis_seta((char*)"data_frame", &ios);
  }
  else if (strcmp(key, "NODATA_FRAME") == 0){
    if (mype == 0)
      vis_seta((char*)"nodata_frame", &ios);
  }
  /* hsurf */
  else if (strcmp(key, "hsurf") == 0){
    if (!get_r(r)){
      printf("Error:  reading hsurf.");
      report();
    }
    else if (r <= 0.0){
      printf("Error:  hsurf must be positive.");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"surface_height", r, &ios);
    }
  } 
  /* light_fusion light_fusion - light-fusion intensity */
  else if (strcmp(key, "LIGHT_FUSION") == 0){
    if (!get_r(r)){
      printf("Error:  reading light_fusion");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"light_fusion", r, &ios);
    }
  }
  /* light_azim light_azim - azimuthal direction of light-point */
  else if (strcmp(key, "LIGHT_AZIM") == 0){
    if (!get_r(r)){
      printf("Error:  reading light_azim");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"light_azim", r, &ios);
    }
  }
  /* light_altd light_altd - altitudinal direction of light-point */
  else if (strcmp(key, "LIGHT_ALTD") == 0){
    if (!get_r(r)){
      printf("Error:  reading light_altd");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"light_altd", r, &ios);
    }
  }
  /* light_distance light_distance of light-point from focus */
  else if (strcmp(key, "LIGHT_DISTANCE") == 0){
    if (!get_r(r)){
      printf("Error:  reading light_distance");
      report();
    }
    else{
      if (mype == 0)
	vis_set0r((char*)"light_distance", r, &ios);
    }
  }
  /* nolight_distance */
  else if (strcmp(key, "NOLIGHT_DISTANCE") == 0){
    vis_seta((char*)"nolight_distance", &ios);
  }
  /* sky red green blue - set sky colour */
  else if (strcmp(key, "SKY") == 0){
    if (!get_r(red)){
      printf("Error:  reading sky red");
      report();
    }
    else if (!get_r(green)){
      printf("Error:  reading sky green");
      report();
    }
    else if (!get_r(blue)){
      printf("Error:  reading sky blue");
      report();
    }
    else{
      if (mype == 0){
	vis_set0r((char*)"sky_red", red, &ios);
	vis_set0r((char*)"sky_green", green, &ios);
	vis_set0r((char*)"sky_blue", blue, &ios);
      }
    }
  }
  /* material name rgb red green blue end - define material */
  else if (strcmp(key, "MATERIAL") == 0){
    if (!get_n(tname)){
      printf("Error:  reading material name.");
      report();
    }
    else{
      if (!get_n(key)){
	printf("Error:  reading material key.");
	report();
      }
      else if (strcmp(key, "rgb") != 0){
	printf("Error:  expecting material rgb key.");
	report();
      }
      else if (!get_r(red)){
	printf("Error:  reading material red");
	report();
      }
      else if (!get_r(green)){
	printf("Error:  reading material green");
	report();
      }
      else if (!get_r(blue)){
	printf("Error:  reading material blue");
	report();
      }
      else if (!get_n(key)){
	printf("Error:  reading material key.");
	report();
      }
      else if (strcmp(key, "end") != 0){
	printf("Error:  expecting material end key.");
	report();
      }
      else{
	if (mype == 0){
	  strcpy(name, "material_red");
	  strcat(name, tname);
	  vis_set0r(name, red, &ios);
	  strcpy(name, "material_green");
	  strcat(name, tname);
	  vis_set0r(name, green, &ios);
	  strcpy(name, "material_blue");
	  strcat(name, tname);
	  vis_set0r(name, blue, &ios);
	}
      }
    }
  }
  /* material_red name red - define material red reflectance */
  else if (strcmp(key, "MATERIAL_RED") == 0){
    if (!get_n(tname)){
      printf("Error:  reading material name.");
      report();
    }
    else if (!get_r(red)){
      printf("Error:  reading material red");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "material_red/");
	strcat(name, tname);
	vis_set0r(name, red, &ios);
      }
    }
  }
  /* material_green name green - define material green reflectance */
  else if (strcmp(key, "MATERIAL_GREEN") == 0){
    if (!get_n(tname)){
      printf("Error:  reading material name.");
      report();
    }
    else if (!get_r(green)){
      printf("Error:  reading material green");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "material_green/");
	strcat(name, tname);
	vis_set0r(name, green, &ios);
      }
    }
  }
  /* material_blue name blue - define material blue reflectance */
  else if (strcmp(key, "MATERIAL_BLUE") == 0){
    if (!get_n(tname)){
      printf("Error:  reading material name.");
      report();
    }
    else if (!get_r(blue)){
      printf("Error:  reading material blue");
      report();
    }
    else{
      if (mype == 0){
	strcpy(name, "material_blue/");
	strcat(name, tname);
	vis_set0r(name, blue, &ios);
      }
    }
  }
  /* material_shine name shine - define material shine reflectance */
  else if (strcmp(key, "MATERIAL_SHINE") == 0){
    if (!get_n(tname)){
      printf("Error:  reading material name.");
      report();
    }
    else if (!get_r(shine)){
      printf("Error:  reading material shine");
      report();
    }
    else{
      if (mype == 0){
	printf("Warning; material shine not implemented.");
	/*
	  strcpy(name, "material_shine/");
	  strcat(name, tname);
	  vis_set0r(name, shine, &ios);
	*/
      }
    }
  }
  /* material_phong name phong - define material phong reflectance */
  else if (strcmp(key, "MATERIAL_PHONG") == 0){
    if (!get_n(tname)){
      printf("Error:  reading material name.");
      report();
    }
    else if (!get_r(phong)){
      printf("Error:  reading material phong");
      report();
    }
    else{
      if (mype == 0){
	printf("Warning; material shine not implemented.");
	/*
	  strcpy(name, "material_phong/");
	  strcat(name, tname);
	  vis_set0r(name, phong, &ios);
	*/
      }
    }
  }
  /* nomaterials - remove them */
  else if (strcmp(key, "NOMATERIALS") == 0){
    vis_seta((char*)"material_reset", &ios);
  }
  /* tag name material - assign a material to a csg tag */
  else if (strcmp(key, "TAG") == 0){
    if (!get_n(tname)){
      printf("Error:  reading tag name.");
      report();
    }
    else if (!get_n(name)){
      printf("Error:  reading material name.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(tlabel, "tag/");
	strcat(tlabel, tname);
	vis_set0s(tlabel, name, &ios);
      }
    }
  }
  /* notags - remove them */
  else if (strcmp(key, "NOTAGS") == 0){
    vis_seta((char*)"tag_reset", &ios);
  }
  /* vplot coord scalar - select coord-scalar pair for vector */
  else if (strcmp(key, "VPLOT") == 0){
    if (!get_n(tname)){
      printf("Error:  reading vector coord.");
      report();
    }
    else if (!get_n(name)){
      printf("Error:  reading vector scalar.");
      report();
    }
    else{
      if (mype == 0){
	strcpy(type, "vplot/");
	strcat(type, tname);
	model_set0s(type, name, &ios);
      }
    }
  }
  /* unknown command */
  else
    cmd = false;
  return (cmd);
}
//end mvcmd.c
