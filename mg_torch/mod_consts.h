/* Model constants */
#ifdef USER
int iqdn,iqx,iqt;
bool helium,cons_he,cons_max_rate,logstep,log_on,prim_only,carbonoxygen,isomix,purecarbon;

bool log_shk_rels = false;
double mass_flux = 1.0; // Default for hydrostatic runs to behave
double mom_flux = 1.0; // Momentum flux
const int FILE_NAME_PREC_D = 15; // Precision for det speed in file names and output
const int VAR_PREC_OUT = 15;

int idrho = 14;
int idte = 13;

#else
extern int iqdn,iqx,iqt;
extern bool helium,cons_he,cons_max_rate,logstep,log_on,prim_only,carbonoxygen,isomix,purecarbon;
extern bool log_shk_rels;
extern double mass_flux; // Default for hydrostatic runs to behave
extern double mom_flux; // Momentum flux

extern const int FILE_NAME_PREC_D; // Precision for det speed in file names and output
extern const int VAR_PREC_OUT;

extern int idrho;
extern int idte;

/* Variabls user.h does not export */
extern int ih1, ihe3, ihe4, ic12, in14, io16, ine20, img24, isi28, is32, iar36,
  ica40, iti44, icr48, icr56, ife52, ife54, ife56, ini56, ineut, iprot;
extern float Zi[MAXSPECIES], Qi[MAXSPECIES],
  Mi[MAXSPECIES];
extern int nsub;
/* Variabls user.h does not export */

#endif
