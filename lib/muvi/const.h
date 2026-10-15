/* const.h:  general constants */
#define SZNAM 32
#define SZLIN 80
#define MAXGREY 4096
#define MAXCOLOUR 4096
#define bool int
#define TRUE 1
#define FALSE 0
#define OK FALSE
#define ERR TRUE
#define MAXD 3
#define NIL 0
#define BIG 1.0E36
/*
 define rem(N,M) (N - M*(N/M))
*/
#define BLANK 0
#define FLUID 1
#define SOLID 2
#define SZDATE 64
#define PI 3.1415927
#define LEFT 0
#define RIGHT 1
#define CURV
#define CURV_MUVI
#define CURV_OVER
#undef CURV_PLAY
#define PLOTG_SLICE
#define STL

#ifdef DOUBLE
#define float double
#endif
