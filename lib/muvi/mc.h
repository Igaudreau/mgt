/* mc.h:  define which machine we are running on

      SUN      Sun Microsystems (SunOS+Sunview)
      SGI      Silicon Graphics Iris (unix+GraphicsLibrary)
      DEC      Digital Alpha
      VAX      DEC VAX (VMS)
      FX       Alliant FX-series

and which back-ends are supported:

      GL       Old SG Graphics Library
      GLX      OpenGL + X11
      X11      X11
*/

#include "machine.h"

#ifdef SUN
#undef sun
#define BY_VALUE
#endif

#ifdef SGI
#undef sgi
#define BY_VALUE
#endif

#ifdef IBM
#undef ibm
#define BY_VALUE
#endif

#ifdef DEC
#undef dec
#define BY_VALUE
#endif

#ifdef HP
#undef hp
#define BY_VALUE
#endif


#ifdef FX
#define BY_ADDRESS
#endif

#ifdef VAX
#define BY_ADDRESS
#endif

#ifdef BY_ADDRESS
#define LOCV(char_len) *char_len
#define ADRV(char_len) &char_len
#endif

#ifdef BY_VALUE
#define LOCV(char_len) char_len
#define ADRV(char_len) char_len
#endif

