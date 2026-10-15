/* Physical parameters */
#define MAXPAR 10
/* Steady solutions structure */
#ifdef USER
int iqd, iqu0, iqe, iqk, iqeps, iqal0, iqie, delpart, iqpot, iqte, iqsf,
  niso, nrat, jaccall, reactcall, nsub, niter, solev, curlev;
float g, pfloor, dfloor, avisp, avise, tcond, sdiff, vis, ce1,
  ce2, cs, ct, cd, ce, ceps, kmin, epsmin, mtl, gconst, rat, T_e, T_0,
  rho_0, T_s, rho_s, cns, dtp, tp_old, curv, depth, cns1, cns2, D, D_l,
  D_r, par[MAXPAR], amplitude, period, gslrtol, gslatol, gslhmin, gslhmax;
/* Species indices */
int ih1, ih2, ih3, ihe3, ihe4, ili6, ili7, ili8, ibe7, ibe8, ibe9,
ibe10, ib7, ib8, ib10, ib11, ic11, ic12, ic13, ic14, in13, in14, in15,
io14, io15, io16, io17, io18, io19, if17, if18, if19, ine18, ine19,
ine20, ine21, ine22,  ine23, ina18, ina19, ina20, ina21, ina22,ina23,
ina24, img22, img23, img24, img25, img26, img27, ial25, ial26, ial27,
ial28, isi27, isi28, isi29, isi30, isi31, isi32, ip29, ip30, ip31, ip32,
ip33, ip34, is30, is31, is32, is33, is34, is35,  is36, is37, icl33,
icl34, icl35, icl36, icl37, icl38, iar35, iar36, iar37, iar38, iar39,
iar40, iar41, ik37, ik38, ik39, ik40, ik41, ik42, ica39, ica40, ica41,
ica42, ica43, ica44, ica45, ica46, ica47, ica48, ica49,  isc41, isc42,
isc43, isc44, isc45, isc46, isc47, isc48, isc49, iti43, iti44, iti45,
iti46, iti47, iti48, iti49, iti50, iti51, iti52, iti53, iti54, iti55,
iti56, iti57, iti58, iv45, iv46, iv47, iv48, iv49,iv50, iv51, iv52,
iv53, iv54, iv55, iv56, iv57, iv58, icr47, icr48, icr49, icr50, icr51,
icr52, icr53, icr54, icr55,  icr56, icr57, icr58, icr59, icr60, imn51,
imn52, imn53, imn54, imn55, imn56, imn57, imn58, imn59, imn60, imn61,
ife52, ife53, ife54, ife55, ife56, ife57, ife58, ife59, ife60, ife61,
ico55, ico56, ico57,ico58, ico59, ico60, ico61, ico62, ico63, ico64,
ico65, ini56, ini57, ini58, ini59, ini60, ini61, ini62, ini63, ini64,
ini65, icu57, icu58, icu59, icu60, icu61, icu62, icu63, icu64, icu65,
icu66, izn60, izn61, izn62, izn63, izn64, izn65, izn66, izn67, izn68,
izn69, ineut,iprot;
int ir3a, irg3a, ircag, iroga, iroag, irnega, irneag, irmgga, irmgag,
  irsiga, ircaag, irtiga, ir1212, ir1216, ir1616, irsi2ni, irni2si;
/* torch switches */
int doweak, doneut;
bool mon, fxp, inten, hll, god, kt, tim7, tim13, tim21, timtorch,
  rrscreen, revscreen, react, hydrostatic, isothermal, cons, init, sootfoil;
char ode[SZNAM];
char monpref[SZNAM], monmod[SZNAM];
/* Helmholtz free energy and derivatives and local copy */ 
float f[JMAX+1][IMAX+1], fd[JMAX+1][IMAX+1], ft[JMAX+1][IMAX+1],
  fdd[JMAX+1][IMAX+1], ftt[JMAX+1][IMAX+1], fdt[JMAX+1][IMAX+1],
  fddt[JMAX+1][IMAX+1], fdtt[JMAX+1][IMAX+1], fddtt[JMAX+1][IMAX+1],
  fl[37];
/* Range and spacings for EOS */
const float   tlo = 4.0, dlo = -10.0, thi = 11.0, dhi = 11.0;
const float  tstp  = (thi - tlo)/float(JMAX-1),
  dstp  = (dhi - dlo)/float(IMAX-1);
/* Physical constants */
const float k_b = 1.380649e-16, N_A = 6.0221409e23,
  phygconst = 6.6743e-8, stefan = 7.565914e-15, mn = 1.67492721184e-24,
  mp = 1.67262163783e-24, me = 9.1093821545e-28,
  mev2gr = 1.7826617585133994e-27,
  enuc_conv2 = -5.4124311235461752e+44;
/* Density and temperature points */
float d[IMAX+1], te[JMAX+1];
/* Nuclear properties */
float Zbar, Abar, Z2bar;
/* Arrays  for atomic  mass, atomic  number, binding  energy (ergs/g),
   total binding energy (erg/gm) of each isotope
*/   
float Ai[MAXSPECIES], Zi[MAXSPECIES], Qi[MAXSPECIES], aion[MAXSPECIES],
  zion[MAXSPECIES], bion[MAXSPECIES];
char ionam [MAXSPECIES][5];
/* Array for screen functions */
float zs13[25], zhat[25], zhat2[25], lzav[25], aznut[25], zs13inv[25];
#else
extern int iqd, iqu0, iqe, iqk, iqeps, iqal0, iqie, delpart, iqpot, iqsf,
  iqte, niso, jaccall, reactcall, nsub, niter, solev, curlev;
extern float g, pfloor, dfloor, avisp, avise, tcond, sdiff, vis, ce1,
  ce2, cs, ct, cd, ce, ceps, kmin, epsmin, mtl, gconst, rat, T_e, T_0,
  rho_0, cns, dtp, tp_old, curv, depth, cns1, cns2, rho_0, T_s, rho_s,
  cns, dtp, tp_old, curv, depth, cns1, cns2, D, D_l, D_r, par[MAXPAR],
  amplitude, period, gslrtol, gslatol, gslhmin, gslhmax;
/* Species indices */
extern int ih1, ih2, ih3, ihe3, ihe4, ili6, ili7, ili8, ibe7, ibe8, ibe9,
ibe10, ib7, ib8, ib10, ib11, ic11, ic12, ic13, ic14, in13, in14, in15,
io14, io15, io16, io17, io18, io19, if17, if18, if19, ine18, ine19,
ine20, ine21, ine22,  ine23, ina18, ina19, ina20, ina21, ina22,ina23,
ina24, img22, img23, img24, img25, img26, img27, ial25, ial26, ial27,
ial28, isi27, isi28, isi29, isi30, isi31, isi32, ip29, ip30, ip31, ip32,
ip33, ip34, is30, is31, is32, is33, is34, is35,  is36, is37, icl33,
icl34, icl35, icl36, icl37, icl38, iar35, iar36, iar37, iar38, iar39,
iar40, iar41, ik37, ik38, ik39, ik40, ik41, ik42, ica39, ica40, ica41,
ica42, ica43, ica44, ica45, ica46, ica47, ica48, ica49,  isc41, isc42,
isc43, isc44, isc45, isc46, isc47, isc48, isc49, iti43, iti44, iti45,
iti46, iti47, iti48, iti49, iti50, iti51, iti52, iti53, iti54, iti55,
iti56, iti57, iti58, iv45, iv46, iv47, iv48, iv49,iv50, iv51, iv52,
iv53, iv54, iv55, iv56, iv57, iv58, icr47, icr48, icr49, icr50, icr51,
icr52, icr53, icr54, icr55,  icr56, icr57, icr58, icr59, icr60, imn51,
imn52, imn53, imn54, imn55, imn56, imn57, imn58, imn59, imn60, imn61,
ife52, ife53, ife54, ife55, ife56, ife57, ife58, ife59, ife60, ife61,
ico55, ico56, ico57,ico58, ico59, ico60, ico61, ico62, ico63, ico64,
ico65, ini56, ini57, ini58, ini59, ini60, ini61, ini62, ini63, ini64,
ini65, icu57, icu58, icu59, icu60, icu61, icu62, icu63, icu64, icu65,
icu66, izn60, izn61, izn62, izn63, izn64, izn65, izn66, izn67, izn68,
izn69, ineut,iprot;
extern bool mon, fxp, inten, hll, god, kt, iso7, iso13, react,
  hydrostatic, isothermal, sootfoil;
extern char ode[SZNAM];
extern char monpref[SZNAM], monmod[SZNAM];
extern float Ai[MAXSPECIES], Zi[MAXSPECIES], Qi[MAXSPECIES],
  Mi[MAXSPECIES];
#endif
