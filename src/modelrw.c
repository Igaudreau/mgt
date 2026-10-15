#undef MAIN
#include <iostream>
#include<iomanip>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <new>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "funcs.h"
#include "para.h"
#include "user.h"
extern void mk_level();
extern void mk_plevel();
extern void mk_memory();
extern struct part_list *mov_part(bool all, struct part_list *start);
extern void mk_partlist(struct part_list *&start, 
  struct part_list *&end);
extern void clear();
extern float cell_vol(int ilev, float xc[]);
extern struct cell *find_fcell(int &ilev, float xp[]);
extern struct part_list *assign_part(struct part_list *start);
extern void usrwr();
extern bool usrrd(char vname[]);
extern void usrpath(char path[]);
extern void broadcast(int *value);
extern void broadcast(int buf[], int nvalues);
extern void broadcast(float *value);
extern void broadcast(float buf[], int nvalues);
extern void broadcast(bool *value);
extern void broadcast(char value[]);
extern void merge_array(float buf[], int nvalues[]);
extern float global_max(float value);
extern struct cell *new_cell();
extern void connect(int ilev);
extern void connect_par(int ilev, struct cell *cptr);
extern void de_frag();
extern void create_halo(int ilev);
extern void mk_locarr();
extern void set_bound_flag(int ilev);
extern int re_part();
extern bool init_part();
bool writem(char name[]);
void write_global();
void write_var(const char vname[], float value);
void write_var(const char vname[], int size, float value[]); 
void write_var(const char vname[], int value);
void write_var(const char vname[], int size, int value[]);
void write_var(const char vname[], bool value);
void write_var(const char vname[], char value[]);
bool readm(char name[]);
bool distrib();
bool re_distrib(int tnlevs, int ilev, float *&inbufr, 
  float *&outbufr);
bool open_input(int npe);
bool find_data_files(int npe, int ilev, int &lnpe, int &hnpe);
void rec_cell_data(int ilev, float *&inbufr);
bool read_global(); 
void read_var(float *value);
void read_var(int *value);
void read_var(bool *value);
void read_var(char *value);
void safe_read(int fd, void *buf, int count);		// BM
void safe_write(int fd, const void *buf, int count);	// BM

int unit, maxcom, maxpcom, olimpe, bufsiz;
char head[SZNAM], input[SZNAM];

/* NOTE:  This version of modelrw.c is updated for the normal particle system (not LEM)
      Saves particle information on every CPU */


//file writem.c
/*
+ writem - write model to file.
     Description :
        This writes the information  that is necessary to recreate the
        model to files.  The global data is put in filename.h.dsb, and
        the data  for each  processor in files  filename.npe.dsb where
        npe is  the processor number.  It returns  true if successful,
        false if not.
*/   
bool writem(char name[])
{
  int ilev, id, i, npe, ios, count, len, ncell[MAXLEV], npart;
  bool fail;
  char path[SZNAM], tail[SZNAM], ptail[SZNAM], output[MAXPE][SZNAM], partfile[MAXPE][SZNAM];
  float ocpu, ncpu, *inbufr, *outbufr, *pinbufr, *poutbufr;
  const float clocks = CLOCKS_PER_SEC;
  struct cell *cptr;
  struct part *pptr;
  struct part_list *start, *end, *plist, *tplist;
  ocpu = clock()/clocks;
  if (nom){
    /* Create ins and outs arrays */
    ins = new int[limpe];
    outs = new int[limpe];
  }
  if (mype == 0){
    /* Output file is path/name.nnn.dsb */
    usrpath(path);
    if (getenv(path) == NULL){
      cout << path << "not defined" << endl;
      fail = true;
    }
    else{
      /* Get output file name */
      strcpy(head, getenv(path));
      strcat(head, "/");
      strcat(head, name);
      strcpy(output[0], head);
      strcat(output[0], ".h.dsb");
      /*  Open output file on unit (write only, create, user +rwx)
	  do not overwrite existing files
      */
      unit = open(output[0], O_WRONLY|O_CREAT|O_EXCL, S_IRWXU);
      if (unit < 0){
	cout << "Error : unable to open file (exists?)" << endl;
	fail = true;
      }
      else
	fail = false;
    }
  }
  /* All processors return on failure */
  broadcast(&fail);
  if (fail)
    return(false);
  /* Calculate number of cells on each level (assume defragged) */
  for (ilev = 0; ilev < nlevs; ilev++){
    if (lev[ilev]->fcell != NULL)
      ncell[ilev] = lev[ilev]->lcell->icell
	- lev[ilev]->fcell->icell + 1;
    else
      ncell[ilev] = 0;
  }
  /* count number of particles here BM */
  npart = 0;
  if (particle) {
    /* Particle information on each proc */
    for (npe = 0; npe < limpe; npe++){
      if (mype == npe){
	// make part list on this processor 
	mk_partlist(start, end);
        /* Count particles */
        npart = 0;
        plist = start;
        while (plist != NULL){
	  npart ++;
	  plist = plist->nxt;
        }
      }
    }
  }
  /* Calculate maxcom */
  maxcom = 0;
  maxpcom = 0;
  for (ilev = 0; ilev < nlevs; ilev++)
    maxcom = max(maxcom, ncell[ilev]);
  maxpcom = max(maxpcom, npart);
  /* Global maxcom */
  maxcom = (int)global_max((float)maxcom);
  maxpcom = (int)global_max((float)maxpcom);
  if (mype == 0){
    /* Master writes globals to unit */
    write_global();
    /* Master closes unit */
    if (close(unit) != 0){
      cout << "failed to close file " << output[0] << endl;
      fail = true;
    }
  }
  /* All processors return on failure */
  broadcast(&fail);
  if (fail)
    return(false);
  /* Set up communications arrays */
  inbufr = new float[(nq + nd)*maxcom+3*nlevs+SAFETY];
  outbufr = new float[(nq + nd)*maxcom+3*nlevs+SAFETY];
  /* Zero ins/outs for safety */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Set file names for each processor */
  for (npe = 0; npe < limpe; npe++){
    strcpy(output[npe], head);
    if (limpe < 10)
      sprintf(tail,".00%d.dsb", npe);
    else if (limpe < 100)
      sprintf(tail,".0%d.dsb", npe);
    else
      sprintf(tail,".%d.dsb", npe);
    strcat(output[npe], tail);
    /* particle filename for each processor BM */
    if (particle && writepart) {
      strcpy(partfile[npe], head);
      if (limpe < 10)
        sprintf(ptail,".00%d.p.dsb", npe);
      else if (limpe < 100)
        sprintf(ptail,".0%d.p.dsb", npe);
      else
        sprintf(ptail,".%d.p.dsb", npe);
      strcat(partfile[npe], ptail);
    }
  }
  /* Write model information for each processor */
  for (npe = 0; npe < limpe; npe++){
    /* Each processor in turn sends information to master */
    outs[0] = 0;
    if (mype == npe){
      count = -1;
      /* No of cells and processor boundaries on each level */
      for (ilev = 0; ilev < nlevs; ilev++){
	count++;
	outbufr[count] = lev[ilev]->xl[nd-1];
	count++;
	outbufr[count] = lev[ilev]->xr[nd-1];
	count++;
	outbufr[count] = (float)ncell[ilev];
      }
      outs[0] = count + 1;
    }
    ios = 0;
    parc_xdatar(&limpe, outs, outbufr, ins, inbufr, &ios);
    if (ios != 0)
      cout << "Communications failure in model write!" << endl;
    else{
      /* Master writes cell information from processor npe */
      if (mype == 0){
	/* Open file */
	unit = open(output[npe], O_WRONLY|O_CREAT|O_EXCL, S_IRWXU);
	if (unit < 0){
	  cout << "unable to open file " << output[npe] << endl;
	  fail = true;
	}
      }
      /* Broadcast failure */
      broadcast(&fail);
      if (fail)
	break;
      /* Calculate total length */
      len = 3*nlevs;
      for (ilev= 0; ilev < nlevs; ilev++)
	len += ncell[ilev]*(nd + nq);
      len *= sizeof(float);
      if (mype == 0){
	safe_write(unit, &len, sizeof(int));
	if (npe == 0)
	  /* Master writes its own data from outbufr */
	  safe_write(unit, outbufr, 3*nlevs*sizeof(float));
	else
	  /* Master writes sybling data from inbufr */
	  safe_write(unit, inbufr, 3*nlevs*sizeof(float));
	if (close(unit) != 0){
	  cout << "failed to close file " << output[npe] << endl;
	  fail = true;
	}
      }
      /* Broadcast failure */
      broadcast(&fail);
      if (fail)
	break;
    }
  }
  /* Cell information on each level */
  for (ilev = 0; ilev < nlevs; ilev++){
    for (npe = 0; npe < limpe; npe++){
      /* Each processor in turn sends information to master */
      outs[0] = 0;
      if (mype == npe){
	count = -1;
	cptr = lev[ilev]->fcell;
	while (cptr != NULL){
	  /* Cell coordinates */
	  for (id = 1; id <= nd; id++){
	    count++;
	    outbufr[count] = cptr->xc[id-1];
	  }
	  /* Cell solution */
	  for (i = 0; i < nq; i++){
	    count++;
	    outbufr[count] = cptr->qa[i];
	  }
	  cptr = cptr->nxt;
	}
	outs[0] = count + 1;
      }
      ios = 0;
      parc_xdatar(&limpe, outs, outbufr, ins, inbufr, &ios);
      if (ios != 0)
	cout << "Communications failure in model write!" << endl;
      else{
	/* Reset parc buffer size */
	parc_reset();
	/* Master writes cell information from processor npe */
	if (mype == 0){
	  /* Open file */
	  unit = open(output[npe], O_WRONLY|O_APPEND, S_IRWXU);
	  if (unit < 0){
	    cout << "unable to open file " << output[npe] << endl;
	    fail = true;
	  }
	}
	/* Broadcast failure */
	broadcast(&fail);
	if (fail)
	  break;
	if (mype == 0){
	  if (npe == 0){
	    /* Master writes its own data from outbufr */
	    len = outs[npe]*sizeof(float);
	    safe_write(unit, outbufr, len);
	  }
	  else{
	    /* Master writes sybling data from inbufr */
	    len = ins[npe]*sizeof(float);
	    safe_write(unit, inbufr, len);
	  }
	  if (close(unit) != 0){
	    cout << "failed to close file " << output[npe] << endl;
	    fail = true;
	  }
	}
	/* Broadcast failure */
	broadcast(&fail);
	if (fail)
	  break;
      }
    }
  }
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  /* Return on failure */
  if (fail)
    return(false);
  if (particle  && writepart){
    /* Set up communications arrays */
    pinbufr = new float[(2*nd+4+MAXPROP)*maxpcom+2+SAFETY];	// NOTE:  CHANGE NUMBER OF PRPERTIES LIKE EN, MASS, etc.. HERE.
    poutbufr = new float[(2*nd+4+MAXPROP)*maxpcom+2+SAFETY];
    /* Zero ins/outs for safety */
    for (npe = 0; npe < limpe; npe++){
      ins[npe] = 0;
      outs[npe] = 0;
    }
    /* Particle information on each proc */
    for (npe = 0; npe < limpe; npe++){
      /* Each processor in turn sends information to master */
      outs[0] = 0;
      if (mype == npe){
	// make part list on this processor 
	mk_partlist(start, end);
	count = -1;
        /* Count particles */
        npart = 0;
        plist = start;
        while (plist != NULL){
	  npart ++;
	  plist = plist->nxt;
        }
	if (npart > maxpcom)
	  cout << "npart on cpu " << mype << " greater than maxpcom " << maxpcom << endl;
        float fnpart, fmaxprop;
        fnpart = (float)npart;
        fmaxprop = (float)MAXPROP;
	count++;
	poutbufr[count] = fnpart;
	count++;
	poutbufr[count] = fmaxprop;

        plist = start;
        while (plist != NULL){
	  pptr = plist->pptr;
	  /* Particle coordinates and velocities */
	  for (id = 1; id <= nd; id++){
	    count++;
	    poutbufr[count] = pptr->xp[id-1];
	    count++;
	    poutbufr[count] = pptr->vp[id-1];
	  }
	  /* Particle level, time, mass and thermal energy */
	  count++;
	  poutbufr[count] = (float)pptr->ilev;
	  count++;
	  poutbufr[count] = pptr->t;
	  count++;
	  poutbufr[count] = pptr->mass;
	  count++;
	  poutbufr[count] = pptr->en;
	  /* Other particle properties */
	  for (i = 0; i < MAXPROP; i++){
	    count++;
	    poutbufr[count] = pptr->props[i];
	  }
	  plist = plist->nxt;
        }
        outs[0] = count + 1;
      }
      ios = 0;
      parc_xdatar(&limpe, outs, poutbufr, ins, pinbufr, &ios);
      if (ios != 0)
	cout << "Communications failure in model write!" << endl;
      else{
	/* Reset parc buffer size */
	parc_reset();
	/* Master writes cell information from processor npe */
	if (mype == 0){
	  /* Open file */
          unit = open(partfile[npe], O_WRONLY|O_CREAT|O_EXCL, S_IRWXU);
	  if (unit < 0){
	    cout << "unable to open file " << partfile[npe] << endl;
	    fail = true;
	  }
	  else
	    fail = false;
	}
	/* Broadcast failure */
	broadcast(&fail);
	if (fail)
	  break;
	if (mype == 0){
	  if (npe == 0){
	    /* Master writes its own data from outbufr */
	    len = outs[npe]*sizeof(float);
	    safe_write(unit, poutbufr, len);
	  }
	  else{
	    /* Master writes sybling data from inbufr */
	    len = ins[npe]*sizeof(float);
	    safe_write(unit, pinbufr, len);
	  }
	  if (close(unit) != 0){
	    cout << "failed to close file " << partfile[npe] << endl;
	    fail = true;
	  }
          else
            fail = false;
	}
	/* Broadcast failure */
	broadcast(&fail);
	if (fail)
	  break;
      }
    }
    /* Delete communications arrays */
    delete [] pinbufr;
    delete [] poutbufr;
    if (fail)
      return(false);
    /* Delete other processor particles and list */
    plist = start;
    while (plist != NULL){
      /* Delete particle */
      //delete plist->pptr;
      /* Delete list element */
      tplist = plist;
      plist = plist->nxt;
      delete tplist;
    }
  }
  /* Broadcast failure */
  broadcast(&fail);
  if (fail)
    return(fail);
  if (mype == 0){
    ncpu = clock()/clocks;
    cout << "wrote " << head << endl;
    cout.setf(ios::scientific);
    cout << "that took " << setprecision(4)
	   << ncpu - ocpu << " secs" << endl;
  }
  return(true);
}
//end writem.c

//file write_global.c
/*
+ write_global - write globals to file.
     Description :
        This writes the global variables to unit.
*/
void write_global()
{
  int len, id;
  char bcond_c[MAXD][SZITEM], vname[SZITEM];
  /* Header is "start" */
  strcpy(vname, "start");
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  write_var("limpe", limpe);
  write_var("title", title);
  write_var("geom", geom);
  write_var("cas", cas);
  write_var("nd", nd);
  write_var("nq", nq);
  write_var("nal", nal);
  write_var("scad", scad);
  write_var("nlevs", nlevs);
  write_var("nrcell", nrcell);
  write_var("order", order);
  write_var("particle", particle);  
  write_var("levc", levc);
  write_var("levb", levb);
  write_var("nc", nd, nc);
  write_var("bcl", nd, bcl);
  write_var("bcr", nd, bcr);
  strcpy(vname, "coord");
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  memset(safe_vname, 0, sizeof(safe_vname));	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = MAXD*SZITEM;
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, coord, len);
  strcpy(vname, "bcond");
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  memset(safe_vname, 0, sizeof(safe_vname));	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  for (id = 1; id <= nd; id++)
    strcpy(bcond_c[id-1], bcond[id-1].c_str());
  len = MAXD*SZITEM;
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, bcond_c, len);
  strcpy(vname, "tnupdat");  
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  memset(safe_vname, 0, sizeof(safe_vname));	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = sizeof(int);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, &tnupdat, len);
  write_var("admin", admin);
  write_var("debug", debug);
  write_var("dump", dump);
  write_var("prof", prof);
  write_var("diffbad", diffbad);
  write_var("remhole", remhole);
  write_var("correct", correct);
  write_var("defrag", defrag);
  write_var("regrid", regrid);
  write_var("repart", repart);
  write_var("timpart", timpart);
  write_var("ndiffstep", ndiffstep);
  write_var("rtol", rtol);
  write_var("derefcrit", derefcrit);
  write_var("refcrit", refcrit);
  write_var("cn", cn);
  write_var("cnf", cnf);
  write_var("cnv", cnv);
  write_var("cng", cng);
  write_var("grav", grav);
  write_var("boundpot", boundpot);
  write_var("nrel", nrel);
  write_var("normrat", normrat);
  write_var("maxmug", maxmug);
  write_var("mugi", mugi);  
  write_var("mcnv", mcnv);
  write_var("fdq", fdq);  
  write_var("t", t);
  write_var("tcpu", tcpu);
  write_var("target", target);
  write_var("xl", nd, xl);
  write_var("xr", nd, xr);
  write_var("maxcom", maxcom);
  write_var("superbee", superbee);
  write_var("mc", mc);
  write_var("vanleer", vanleer);
  write_var("minmod", minmod);
  /* User parameters */
  usrwr();
  /* Write end of globals */
  strcpy(vname, "end");
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  memset(safe_vname, 0, sizeof(safe_vname));	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  return;
}
//end write_global.c

//file write_varr0.c
/*
+ write_varr0 - write a float to file.
     Description :
        This writes the real variable vname to the model file.
*/
void write_var(const char vname[], float value)
{
  int len;
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = sizeof(float);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, &value, len);
  return;
}
//end write_varr0.c

//file write_varr1.c
/*
+ write_varr1 - write a float array to file.
     Description :
        This writes the real variable array vname to the model file.
*/
void write_var(const char vname[], int size, float value[])
{
  int len;
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = size*sizeof(float);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, value, len);
  return;
}
//end write_varr1.c

//file write_vari0.c
/*
+ write_vari0 - write an integer to file.
     Description :
        This writes the integer variable vname to the model file.
*/
void write_var(const char vname[], int value)
{
  int len;
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = sizeof(int);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, &value, len);
  return;
}
//end write_vari0.c

//file write_vari1.c
/*
+ write_vari1 - write an integer array to file.
     Description :
        This writes the integer array vname to the model file.
*/
void write_var(const char vname[], int size, int value[])
{
  int len;
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = size*sizeof(int);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, value, len);
  return;
}
//end write_vari1.c

//file write_varb.c
/*
+ write_varb - write a bool to file.
     Description :
        This writes the boolean variable vname to the model file.
*/
void write_var(const char vname[], bool value)
{
  int len;
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = sizeof(bool);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, &value, len);
  return;
}
//end write_varb.c

//file write_vars.c
/*
+ write_vars - write a string to file.
     Description :
        This writes the string vname to the model file.
*/
void write_var(const char vname[], char value[])
{
  int len;
  //write(unit, vname, SZITEM);
  // Guarantee string will be SZITEM in length  BM
  char safe_vname[SZITEM] = {0};	// pad char array with zeros
  strncpy(safe_vname, vname, SZITEM - 1);
  safe_write(unit, safe_vname, SZITEM);
  len = (strlen(value) + 1)*sizeof(char);
  safe_write(unit, &len, sizeof(int));
  safe_write(unit, value, len);
  return;
}
//end write_vars.c

//file readm.c
/*
+ readm - read model from file.
     Description :
        This reads  the information that is necessary  to recreate the
        model from a file.  The  global data is in filename.h.dsb, and
        the data  for each  processor in files  filename.npe.dsb where
        npe is  the processor number.  It returns  true if successful,
        false if not.
*/
bool readm(char name[])
{
  int ilev, id, i, npe, npe2, count, len, lenf, nprop, npart, ipart;
  bool fail;
  char path[SZNAM], ptail[SZNAM], partfile[MAXPE][SZNAM];
  float *inbufr, temp;
  float ocpu, ncpu;
  const float clocks = CLOCKS_PER_SEC;
  struct part *pptr;
  struct part_list *start, *pmov, *plist, *nplist, *tplist, *prev;
  ocpu = clock()/clocks;
  nprop = 0;
  npart = 0;
  start = NULL;
  /* Set workwt */
  workwt = 1.0;
  if (mype == 0){
    /* Input file is path/name.nnn.dsb */
    usrpath(path);
    if (getenv(path) == NULL){
      cout << path << " not defined" << endl;
      fail = true;
    }
    else{
      fail = false;
      strcpy(head, getenv(path));
      strcat(head, "/");
      strcat(head, name);
      strcpy(input, head);
      strcat(input, ".h.dsb");
      /* Open input file on unit (read only) */
      unit = open(input, O_RDONLY);
      if (unit < 0){
	cout << "unable to open file (header) " << input << endl;
	fail = true;
      }
    }
  }
  /* All processors return on failure */
  broadcast(&fail);
  if (fail)
    return(false);
  /* Clear previous model */
  if (!nom)
    clear();
  /* Read and broadcast globals */
  fail = !read_global();
  /* All processors return on failure */
  broadcast(&fail);
  if (fail)
    return(false);
  if (mype == 0){
    /* Master closes file name.h.dsb */
    if (close(unit) != 0){
      cout << "failed to close file " << input << "!" << endl;
      fail = true;
    }
    else
      fail = false;
  }
  /* All processors return on failure */
  broadcast(&fail);
  if (fail)
    return(false);
  /* Make levels and plevels*/
  mk_level();
  mk_plevel();
  /* Set initial memory size */
  carrtop = CINCR;
  jarrtop = JINCR;
  /* Make memory */
  mk_memory();
  /* Create work arrays */
  wka = new float[nc[nd-1]*ipow(nlevs-1)+SAFETY];
  wkda = new float[nc[nd-1]*ipow(nlevs-1)+SAFETY];
  /* Distribute cell data to all processors */
  fail = !distrib();
  /* All processors return on failure */
  broadcast(&fail);
  if (fail)
    return(false);
  if (particle && writepart){
    /* Read particle data */
    /* Call de_frag with particle= false to connect cells on master */
    particle = false;
    de_frag();
    particle = true;
    for (npe = 0; npe < olimpe; npe++){
      if (mype == 0){
      /* Set file names for each processor */
        strcpy(partfile[npe], head);
        if (olimpe < 10)
          sprintf(ptail,".00%d.p.dsb", npe);
        else if (olimpe < 100)
          sprintf(ptail,".0%d.p.dsb", npe);
        else
          sprintf(ptail,".%d.p.dsb", npe);
        strcat(partfile[npe], ptail);

        unit = open(partfile[npe], O_RDONLY);
        if (unit < 0){
	  cout << "unable to open file (part data) " << partfile[npe] << " npe = " << npe << endl;
	  fail = true;
        }
        else
	  fail = false;
	if (fail)
	  break;

        /* Get number of particles and nprop */
        lenf = sizeof(float);
        safe_read(unit, &temp, lenf);
        npart = (int)temp;
        safe_read(unit, &temp, lenf);
        nprop = (int)temp;
        if (nprop >  MAXPROP){
	  cout << "Error : nprop = " << nprop << " > MAXPROP = " 
            << MAXPROP << " Change MAXPROP in struct.h" << endl;
	  fail = true;
	  if (close(unit) != 0)
	    cout << "failed to close file " << partfile << endl;
        }
        else
	  fail = false;
      }
      broadcast(&fail);
      if (fail)
        return(false);
      /* Initialise particle move list */
      pmov = NULL;
      if (mype == 0){
        /* Assemble particle list */
        count = -1;
        start = NULL;
        prev = NULL;
        for (ipart = 1; ipart <= npart; ipart++){
	  /* Fix (BM) - Load each particle to inbufr individually */
	  len = (2*nd+4+nprop);		// BM	// NOTE:  CHANGE NUMBER OF PRPERTIES LIKE EN, MASS, etc.. HERE.
          inbufr = new float[len+SAFETY];		// BM
	  len = lenf*len;
          safe_read(unit, inbufr, len);		//BM
          count = -1;				// BM

	  plist = new part_list;
	  pptr = new part;
	  /* Particle pointers to NULL */
	  pptr->prev = NULL;
	  pptr->nxt = NULL;
	  pptr->cptr = NULL;
	  plist->pptr = pptr;
          /* original */
	  //plist->nxt = start;
	  //start = plist;
	  /* fix for particle list order on "read" command BM */
	  if (prev != NULL)
	  prev->nxt = plist;
          if (ipart == 1)
 	    start = plist;	// set the first one */

	  /* Particle coordinates and velocities */
	  for (id = 1; id <= nd; id++){
	    count++;
	    pptr->xp[id-1] = inbufr[count]; 
	    count++;
	    pptr->vp[id-1] = inbufr[count];
	  }
	  /* Particle level, time, mass and thermal energy */
	  count++;
	  pptr->ilev = (int)inbufr[count]; 
	  count++;
	  pptr->t = inbufr[count]; 
	  count++;
	  pptr->mass = inbufr[count]; 
	  count++;
	  pptr->en = inbufr[count]; 
	  /* Other particle properties */
	  //for (i = 0; i < MAXPROP; i++){	// note nprop instead of maxprop
	  for (i = 0; i < nprop; i++){
	    count++;
	    pptr->props[i] = inbufr[count]; 
	  }
          prev = plist;	// BM
          delete [] inbufr;	// BM
        }
        if (npart > 0)
          plist->nxt = NULL;	// BM bug fix
	/* Close unit - BM */
        if (close(unit) != 0) {
	   cout << "failed to close file " << partfile[npe] << endl;
	   fail = true;
        }
      }
      /* Call de_frag with particle= false to connect cells on master */
      //particle = false;
      //de_frag();
      //particle = true;
      if (limpe > 1){
        pmov = NULL;
        if (mype == 0){
	  /* Master prunes particle list and creates move list */
	  plist = start;
	  prev = NULL;
	  while (plist != NULL){
	    nplist = plist->nxt;
	    pptr = plist->pptr;
	    /* Particle level */
	    ilev = pptr->ilev;
	    /* Find destination processor */
	    for (npe2 = 0; npe2 < limpe; npe2++){
	      if ((pptr->xp[nd-1] >= levpbd[ilev][npe2])
	  	&& (pptr->xp[nd-1] < levpbd[ilev][npe2+1]))
	        break;
	    }
	    if (npe2 != mype){
	      /* Add to move list */
	      tplist = new part_list;
	      tplist->pproc = npe2;
	      tplist->pptr = pptr;
	      tplist->nxt = pmov;
	      pmov = tplist;
	      /* Delete original list element*/
	      if (prev == NULL)
	        start = nplist;
	      else
	        prev->nxt = nplist;
	      delete plist;
	    }
	    else
	      prev = plist;
	    plist = nplist;
	  }
        }
        /* Move particles */
        pmov = mov_part(false, pmov);
      }
      if (mype == 0) {
       /* Assign particles on master */
        assign_part(start);
      }
      else {
        /* Assign particles on other processors */
        assign_part(pmov);
      }
    }
  }
  /* Clear model and return on failure */
  broadcast(&fail);
  if (fail){
    clear();
    return(false);
  }
  /* Plotting levels */
  plmin = 0;
  plmax = nlevs - 1;
  if (mype == 0){
    ncpu = clock()/clocks;
    cout << "read " << name << endl;
    cout.setf(ios::scientific);
    cout << "that took " << setprecision(4)
	 << ncpu - ocpu << " secs" << endl;
  }
  return(true);
}
//end readm.c

//file distrib.c
/*
+ distrib - distribute cell data.
     Description :
        The  master reads  the data  from the  files  head.npd.dsb and
        distributes  it  among  the  processors. It  returns  true  on
        success and false on failure. re_distrib is used, which starts
        with  a uniform  partition,  distributes each  level and  then
        calls re_part.
*/
bool distrib()
{
  int ilev, tnlevs;
  bool fail;
  float *inbufr, *outbufr;
  /* Create ins and outs arrays */
  ins = new int[limpe];
  outs = new int[limpe];
  /* Create initial communications arrays */
  bufsiz = SAFETY;
  inbufr = new float[bufsiz];
  outbufr =new float[bufsiz];
  /* Start with uniform partition */
  if (!init_part())
    return(false);
  /* Master reads and distributes cell information */
  tnlevs = nlevs;
  for (ilev = 0; ilev < tnlevs; ilev++){
    /* Set temporary level */
    nlevs = min(ilev+1, tnlevs);
    fail = !re_distrib(tnlevs, ilev, inbufr, outbufr);
    /* All processors break on failure */
    broadcast(&fail);
    if (fail)
      break;
    /* Repartition */
    re_part();
  }
  /* Final repartition */
  re_part();
  /* Reset parc buffer size */
  parc_reset();
  /* Delete communications arrays */
  delete [] inbufr;
  delete [] outbufr;
  return(!fail);
}
//end distrib.c

//file re_distrib.c
/*
+ re_distrib - re_distribute cell data.
     Description :
        This  re_distributes the  cell data  for level  ilev  when the
        number  of processors  is different  from the  number  used to
        create the  model. The  master reads the  data from  the files
        head.npd.dsb  and  distributes it  among  the processors.   It
        returns true on success and false on failure.
*/
bool re_distrib(int tnlevs, int ilev, float *&inbufr, float *&outbufr)
{
  int onpe, lnpe, hnpe, npe, tilev, lncell[MAXLEV],
    ndata, cdloc, ncell, icell, i, lenf, count, ios;
  float temp, dx, levpbd[MAXPE];
  bool fail = false;
  /* Zero ins/outs (safety) */
  for (npe = 0; npe < limpe; npe++){
    ins[npe] = 0;
    outs[npe] = 0;
  }
  /* Data sizes */
  lenf = sizeof(float);
    //float len = 3*tnlevs*lenf;
  ndata = nd + nq;
  /* Mesh spacing */
  dx = lev[ilev]->dxbar[nd-1];
  /* Calculate processor boundaries on this level */
  levpbd[0] = xl[nd-1];
  for (npe = 1; npe <  limpe; npe++)
    levpbd[npe] = xl[nd-1] + 
      dx*(int)((pbd[npe] - xl[nd-1] + 0.5*dx)/dx);
  levpbd[limpe] = xr[nd-1];
  /* Run through processors */
  for (npe = limpe - 1; npe >= 0; npe--){
    /* Master finds file numbers relevant to this processor */
    if (mype == 0)
      fail = !find_data_files(npe, ilev, lnpe, hnpe);
    /* All processors break on failure */
    broadcast(&fail);
    if (fail)
      break;
    if (mype == 0){
      /* Initialize (cell number, count, processor boundaries) */
      ncell = 0;
      count = 2;
      /* Master reads data from relevant data files */
      for (onpe = lnpe; onpe <= hnpe; onpe++){
	/* Open data file */
	open_input(onpe);
	lseek(unit, (off_t)sizeof(int), SEEK_SET);
	/* Old number of cells on each level */
	for (tilev = 0; tilev < tnlevs; tilev++){
	  lseek(unit, (off_t)(sizeof(int)+(3*tilev+2)*lenf), SEEK_SET);
	  safe_read(unit, &temp, lenf);
	  lncell[tilev] = (int)temp;
	}
	/* Check buffer size */
	if (count + ndata*lncell[ilev]+1 > bufsiz){
	  bufsiz = count + ndata*lncell[ilev] + 1;
	  outbufr = (float *)realloc((void *)outbufr,
	    (size_t)((bufsiz+SAFETY)*sizeof(float)));
	  inbufr = (float *)realloc((void *)inbufr,
	    (size_t)((bufsiz+SAFETY)*sizeof(float)));
	}
	/* Location of cell data for this level */
	cdloc = 3*tnlevs;
	for (tilev = 0; tilev < ilev; tilev++)
	  cdloc += ndata*lncell[tilev];
	/* Read cell data */
	lseek(unit, (off_t)(sizeof(int)+cdloc*lenf), SEEK_SET);
	safe_read(unit, inbufr, lncell[ilev]*ndata*lenf);
	/* Close unit */
	if (close(unit) != 0){
	  cout << "failed to close file " << input << endl;
	  fail = true;
	  break;
	}
	/* Master puts relevant cell data in outbufr */
	for (icell = 0; icell < lncell[ilev]; icell++){
	  temp = inbufr[icell*ndata+nd-1];
	  if ((temp > levpbd[npe]) && (temp < levpbd[npe+1])){
	    /* Cell belongs on processor npe */
	    ncell++;
	    for (i = 0; i < ndata; i++){
	      count++;
	      outbufr[count] = inbufr[icell*ndata+i];
	    }
	  }
	}
      }
      /* Processor boundaries */
      outbufr[0] = levpbd[npe];
      outbufr[1] = levpbd[npe+1];
      /* Number of cells on processor npe */
      outbufr[2] = ncell;
      /* Set number of outs */
      outs[npe] = count + 1;
    }
    /* All processors break on failure */
    broadcast(&fail);
    if (fail)
      break;
    /* Broadcast buffer size */
    broadcast(&bufsiz);
    /* Fix buffer sizes for siblings */
    if (mype != 0){
      outbufr = (float *)realloc((void *)outbufr,
        (size_t)((bufsiz+SAFETY)*sizeof(float)));
      inbufr = (float *)realloc((void *)inbufr,
       (size_t)((bufsiz+SAFETY)*sizeof(float)));
    }
    if (npe == 0){
      /* Master just uses outbufr */
      if (mype == npe)
	rec_cell_data(ilev, outbufr);
    }
    else{
      /* Send cell information to processor npe */
      ios = 0;
      parc_xdatar(&limpe, outs, outbufr, ins, inbufr, &ios);
      /* Reset parc buffer size */
      parc_reset();
      outs[npe] = 0;
      if (mype == npe)
	rec_cell_data(ilev, inbufr);
    }
  }
  return(!fail);
}
//end re_distrib.c

//file open_input.c
/*
+ open_input  - open input file from processor npe.
     Description :
        The  master opens  the  file head.npe.dsb,  where  npe is  the
        number of the processor to  which the data applies. It returns
        true if successful, false if not.
*/
bool open_input(int npe)
{
  char tail[SZNAM];
  /* Open file head.npe.dsb */
  strcpy(input, head);
  if (olimpe < 10)
    sprintf(tail,".00%d.dsb", npe);
  else if (olimpe < 100)
    sprintf(tail,".0%d.dsb", npe);
  else
    sprintf(tail,".%d.dsb", npe);
  strcat(input, tail);
  /* Open input file on unit (read only) */
  unit = open(input, O_RDONLY);
  if (unit < 0){
    cout << "unable to open file (cell data) " << input << endl;
    return(false);
  }
  else
    return(true);
}
//end file open_input.c

//file find_data_files.c
/*
+ find_data_files - find data files for a processor.
     Description :
        The master finds the data files that contain the cell data for
        processor    npe.    The    files    are   head.lnpe.dsb    to
        head.hnpe.dsb. This returns true on success, false on failure.
*/
bool find_data_files(int npe, int ilev, int &lnpe, int &hnpe)
{
  int onpe;
  float txl[MAXPE], txr[MAXPE];
  bool fail;
  /* Read left/right boundaries */
  fail = false;
  for (onpe = 0; onpe < olimpe; onpe++){
    /* Open data file */
    if (!open_input(onpe)){
      fail = true;
      break;
    }
    /* Read left/right boundaries for this level */
    lseek(unit, (off_t)(sizeof(int)+3*ilev*sizeof(float)),
	  SEEK_SET);
    safe_read(unit, &txl[onpe], sizeof(float));
    safe_read(unit, &txr[onpe], sizeof(float));
    /* Close unit */
    if (close(unit) != 0){
      cout << "failed to close file " << input << endl;
      fail = true;
      break;
    }
  }
  /* Return on failure */
  if (fail)
    return(false);
  /* Find lowest file number */
  for (onpe = 0; onpe < olimpe; onpe++){
    if (txr[onpe] > pbd[npe]){
      lnpe = onpe;
      break;
    }
  } 
  /* Find highest file number */
  for (onpe = olimpe-1; onpe >= 0; onpe--){
    if (txl[onpe] < pbd[npe+1]){
      hnpe = onpe;
      break;
    }
  }
  return(true); 
}
//end find_data_files.c

//file rec_cell_data.c
/*
+ rec_cell_data - receive model cell data.
     Description :
        The cell data  in inbufr is used to create  and fill the cells
        with the solution for level ilev.
*/
void rec_cell_data(int ilev, float *&inbufr)
{
  int count, icell, fcell, lcell, id, i, ncell;
  struct cell *cptr;
  /* Receive processor boundaries and no of cells on each level */
  count = -1;
  count++;
  lev[ilev]->xl[nd-1] = inbufr[count];
  count++;
  lev[ilev]->xr[nd-1] = inbufr[count];
  count++;
  ncell = (int)inbufr[count];
  /* Receive cell information */
  if (ilev == 0)
    celltop = 0;
  if (ncell == 0){
    lev[ilev]->fcell = NULL;
    lev[ilev]->lcell = NULL;
  }
  else{
    fcell = celltop + 1;
    for (icell = 1; icell <= ncell; icell++){
      /* Make cell */
      cptr = new_cell();
      /* Cell coordinates */
      for (id = 1; id <= nd; id++){
	count++;
	cptr->xc[id-1] = inbufr[count];
      }
      /* Cell volume */
      cptr->vol = cell_vol(ilev, cptr->xc);
      /* Solution */
      for (i = 0; i < nq; i++){
	count++;
	cptr->qa[i] = inbufr[count];
      }
    }
    lcell = celltop;
    /* Make cell list */
    for (icell = fcell; icell < lcell; icell++)
      carr[icell]->nxt = carr[icell+1];
    carr[lcell]->nxt = NULL;
    for (icell = fcell + 1; icell <= lcell; icell++)
      carr[icell]->prev = carr[icell-1];
    carr[fcell]->prev = NULL;
    lev[ilev]->fcell = carr[fcell];
    lev[ilev]->lcell = carr[lcell];
  }
  /* Recalculate max no of cells in nd direction on each level */
  lev[ilev]->nc[nd-1] = (int)((lev[ilev]->xr[nd-1] 
    - lev[ilev]->xl[nd-1]+ 0.1*lev[ilev]->dxbar[nd-1])/
    lev[ilev]->dxbar[nd-1]);
  return;
}
//end rec_cell_data.c
 
//file read_global.c
/*
+ read_global - read globals from file.
     Description :
        This  reads  the global  variables  from  unit and  broadcasts
        them. Returns false if there is no start header.
 */
bool read_global()
{
  int id, len;
  bool fail;
  char bcond_c[MAXD][SZITEM], vname[SZITEM];
  /* Check for start header */
  if (mype == 0){
    safe_read(unit, vname, SZITEM);
    if (strcmp(vname, "start") != 0){
      cout << "No start header " << endl;
      fail = true;
    }
    else
      fail = false;
  }
  /* Return on failure */
  broadcast(&fail);
  if (fail)
    return(false);

  /* Master reads global data from unit and broadcasts */
  while (strcmp(vname, "end") != 0){
    if (mype == 0)
      safe_read(unit, vname, SZITEM);
    broadcast(vname);
    if (strcmp(vname, "limpe") == 0)
      read_var(&olimpe);
    else if (strcmp(vname, "title") == 0)
      read_var(title);
    else if (strcmp(vname, "geom") == 0)
      read_var(geom);
    else if (strcmp(vname, "cas") == 0)
      read_var(cas);
    else if (strcmp(vname, "nd") == 0)
      read_var(&nd);
    else if (strcmp(vname, "nq") == 0)
      read_var(&nq);
    else if (strcmp(vname, "nal") == 0)
      read_var(&nal);
    else if (strcmp(vname, "scad") == 0)
      read_var(&scad);
    else if (strcmp(vname, "nlevs") == 0)
      read_var(&nlevs);
    else if (strcmp(vname, "nrcell") == 0)
      read_var(&nrcell);
    else if (strcmp(vname, "order") == 0)
      read_var(&order);
    else if (strcmp(vname, "particle") == 0)
      read_var(&particle);    
    else if (strcmp(vname, "levc") == 0)
      read_var(&levc);
    else if (strcmp(vname, "levb") == 0)
      read_var(&levb);
    else if (strcmp(vname, "nc") == 0)
      read_var(nc);
    else if (strcmp(vname, "bcl") == 0)
      read_var(bcl);
    else if (strcmp(vname, "bcr") == 0)
      read_var(bcr);
    else if (strcmp(vname, "coord") == 0){
      if (mype == 0){
	safe_read(unit, &len, sizeof(int));
	safe_read(unit, coord, len);
      }
      for (id = 1; id <= nd; id++){
	broadcast(coord[id-1]);
      } 
    }
    else if (strncmp(vname, "bcond", SZITEM) == 0){
      if (mype == 0){
	safe_read(unit, &len, sizeof(int));
	safe_read(unit, bcond_c, len);
      }
      for (id = 1; id <= nd; id++){
	broadcast(bcond_c[id-1]);
	bcond[id-1] = bcond_c[id-1];
      }
    }
    else if (strcmp(vname, "tnupdat") == 0){
      if (mype == 0){
	safe_read(unit, &len, sizeof(int));
	safe_read(unit, &tnupdat, len);
      }
    }
    else if (strcmp(vname, "admin") == 0)
      read_var(&admin);
    else if (strcmp(vname, "debug") == 0)
      read_var(&debug);
    else if (strcmp(vname, "dump") == 0)
      read_var(&dump);
    else if (strcmp(vname, "prof") == 0)
      read_var(&prof);
    else if (strcmp(vname, "diffbad") == 0)
      read_var(&diffbad);
    else if (strcmp(vname, "remhole") == 0)
      read_var(&remhole);
    else if (strcmp(vname, "correct") == 0)
      read_var(&correct);
    else if (strcmp(vname, "defrag") == 0)
      read_var(&defrag);
    else if (strcmp(vname, "regrid") == 0)
      read_var(&regrid);
    else if (strcmp(vname, "repart") == 0)
      read_var(&repart);
    else if (strcmp(vname, "timpart") == 0)
      read_var(&timpart);
    else if (strcmp(vname, "ndiffstep") == 0)
      read_var(&ndiffstep);
    else if (strcmp(vname, "rtol") == 0)
      read_var(&rtol);
    else if (strcmp(vname, "derefcrit") == 0)
      read_var(&derefcrit);
    else if (strcmp(vname, "refcrit") == 0)
      read_var(&refcrit);
    else if (strcmp(vname, "cn") == 0)
      read_var(&cn);
    else if (strcmp(vname, "cnf") == 0)
      read_var(&cnf);
    else if (strcmp(vname, "cnv") == 0)
      read_var(&cnv);
    else if (strcmp(vname, "cng") == 0)
      read_var(&cng);
    else if (strcmp(vname, "grav") == 0)
      read_var(&grav);
    else if (strcmp(vname, "boundpot") == 0)
      read_var(&boundpot);
    else if (strcmp(vname, "nrel") == 0)
      read_var(&nrel);
    else if (strcmp(vname, "normrat") == 0)
      read_var(&normrat);
    else if (strcmp(vname, "maxmug") == 0)
      read_var(&maxmug);
    else if (strcmp(vname, "mugi") == 0)
      read_var(&mugi);
    else if (strcmp(vname, "mcnv") == 0)
      read_var(&mcnv);
    else if (strcmp(vname, "fdq") == 0)
      read_var(&fdq);    
    else if (strcmp(vname, "t") == 0)
      read_var(&t);
    else if (strcmp(vname, "tcpu") == 0)
      read_var(&tcpu);
    else if (strcmp(vname, "target") == 0)
      read_var(&target);
    else if (strcmp(vname, "xl") == 0)
      read_var(xl);
    else if (strcmp(vname, "xr") == 0)
      read_var(xr);
    else if (strcmp(vname, "maxcom") == 0)
      read_var(&maxcom);
    else if (strcmp(vname, "superbee") == 0)
      read_var(&superbee);
    else if (strcmp(vname, "mc") == 0)
      read_var(&mc);
    else if (strcmp(vname, "vanleer") == 0)
      read_var(&vanleer);
    else if (strcmp(vname, "minmod") == 0)
      read_var(&minmod);
    /* User parameters */
    else if (usrrd(vname));
    else if (strcmp(vname, "end") != 0){
      if (mype == 0){
	cout << "Unknown parameter " << vname << '\n';
	safe_read(unit, &len, sizeof(int));
	lseek(unit, (off_t)len, SEEK_CUR);
      }
    }
  }
  /* Checks */
  if (nq > MAXQ){
    cout << "Error: nq = " << nq << " > MAXQ = " << MAXQ
	 << " Change MAXQ in struct.h" << endl;
    return(false);
  }
  else if (nd > MAXD){
    cout << "Error: nde = " << nd << " > MAXD = " << MAXD
	 << " Change MAXD in struct.h" << endl;
    return(false);
  }
  else
    return(true);
}
//end read_global.c

//file read_varf.c
/*
+ read_varf - reads a float or float array from file.
     Description :
        This reads a  real or real array from  file and broadcasts it.
 */
void read_var(float *value)
{
  int len, nvalues;
  if (mype == 0){
    safe_read(unit, &len, sizeof(int));
    safe_read(unit, value, len);
  }
  /* Broadcast */
  nvalues = (int)(len/sizeof(float));
  broadcast(&nvalues);
  if (nvalues < 2)
    broadcast(value);
  else
    broadcast(value, nvalues);
  return;
}
//end read_varf.c

//file read_vari.c
/*
+ read_vari - reads an integer or integer array from file.
     Description :
        This  reads  an  integer   or  integer  array  from  file  and
        broadcasts it.
*/
void read_var(int *value)
{
  int len, nvalues;
  if (mype == 0){
    safe_read(unit, &len, sizeof(int));
    safe_read(unit, value, len);
  }
  /* Broadcast */
  nvalues = (int)(len/sizeof(int));
  broadcast(&nvalues);
  if (nvalues < 2)
    broadcast(value);
  else
    broadcast(value, nvalues);
  return;
}
//end read_vari.c

//file read_varb.c
/*
+ read_varb - reads a boolean from file.
     Description :
        This reads a boolean from file and  broadcasts it.

*/
void read_var(bool *value)
{
  int len;
  if (mype == 0){
    safe_read(unit, &len, sizeof(int));
    safe_read(unit, value, len);
  }
  /* Broadcast */
  broadcast(value);
  return;
}
//end read_varb.c

//file read_vars.c
/*
+ read_vars - reads a string from file.
     Description :
        This reads a string from file and broadcasts it.
*/
void read_var(char *value)
{
  int len;
  if (mype == 0){
    safe_read(unit, &len, sizeof(int));
    safe_read(unit, value, len);
  }
  broadcast(value);
  return;
}
//end read_vars.c


// A safe read function that verifies the result automatically BM
void safe_read(int fd, void *buf, int count) {
    if (read(fd, buf, count) != (int)count) {
        fprintf(stderr, "CRITICAL ERROR: read failed!\n");
        exit(1);
    }
    return;
}

// A safe write function that verifies the result automatically BM
void safe_write(int fd, const void *buf, int count) {
    if (write(fd, buf, count) != (int)count) {
        fprintf(stderr, "CRITICAL ERROR: write failed! Disk full or file error.\n");
        exit(1);
    }
    return;
}




