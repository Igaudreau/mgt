/*
          MG  Copyright (C) 2023  Samuel A.E.G. Falle
     
    This program comes with ABSOLUTELY NO WARRANTY; for details type
    `$w'.
    
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `$c' for details.
*/
#define MAIN
#include <stdio.h>
#include <iostream>
#include <iostream>
#include<iomanip>
#include <cstring>
#include <stdio.h>
#include <new>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"

using namespace std;
extern void activ(char *prompt);
extern void deactiv();
extern void plact();
extern void pldac();
extern void deactiv();
extern bool writem(char filename[]);
extern bool readm(char filename[]);
extern void clear();
extern void help(char topic[]);
extern void peek_cell(int ilev, float xc[]);
extern void de_frag();
extern void mem_calc(bool all);
extern bool get_cmd(char key[]);
extern bool get_i(int &value);
extern bool get_r(float &value);
extern bool get_n(char key[]);
extern bool gen_cmd(char key[]);
extern float global_max(float value);
extern bool mvcmd(char key[]);
extern bool usrcmd(char key[]);
extern void usrquit();
extern void usrdef();
extern bool usroff();
extern bool usrsum();
extern void usrprompt(char prompt[]);
extern void report();
extern bool start();
extern void usrprt();
extern bool plcmd(char key[]);
extern void inc_lev(bool doref);
extern void dec_lev();
extern void rem_g0();
extern bool newdom(float nxl[], float nxr[]);
extern bool march(int nsteps);
extern void merge_array(float buf[], int nvalues[]);
extern void clear_mug();
extern void mk_partlist(struct part_list *&start, 
  struct part_list *&end);
extern struct part_list *mov_part(bool all, struct part_list *start);
extern void re_grid(int ilev);
extern void set_alev();
void control();
bool command(char key[]);

//file main.c
/*
+ main - main programme
   Description :
     This is the main program for  the MG program.  It simply sets the
     parc   arguments,  registers  and   calls  the   control  program
     control().
*/
int main(int argc, char *argv[])
{
  int tmype, tlimpe, ios;
  /* Set parc arguments */
  parc_args(argc, argv);
  /* Register */
  ios = 0;
  parc_register(&tmype, &tlimpe, &ios);
  /* Copy of mype and limpe */
  mype = tmype;
  limpe = tlimpe;
  /* Check no of processors */
  if (limpe > MAXPE){
    if (mype == 0){
    cout << "Error:  no processors = " << limpe << " > MAXPE = " << 
      MAXPE << " edit struct.h" <<'\n';
    }
  }
  else
  /* Call control routine */
      control();
  parc_exit();
  return(0);
}
//end main.c

//file control.c
/*
+ control - control routine
   Description :
     This is  the top  level control routine  for the MG  program.  It
     makes  calls  to  principal  application  routines  to  activate,
     deactivate  and parse  commands.  Each  invocation is  assumed to
     start afresh, without any remnant context.
*/
void control()
{
  char prompt[SZNAM], key[SZNAM];
  bool more;
  if (mype == 0) {
    /* Print banner */
    cout <<"      MG  Copyright (C) 2023  Samuel A.E.G. Falle \n";
    cout << '\n' ;    
    cout << "This program comes with ABSOLUTELY NO WARRANTY; ";
    cout << "for details type `$w'. \n";
    cout << "This is free software, and you are welcome to" <<
      " redistribute it \n";
    cout << "under certain conditions; for details type `$c'. \n";    
    cout << '\n' ;
    cout << '\n' ;
    cout << " This version has been modified and maintained by B. Maxwell (bmaxwell@uottawa.ca)" << '\n';
    cout << " Department of Mechanical Engineering, University of Ottawa" << '\n';
    cout << '\n' ;
  }
  /* Activate command interpreter, set command prompt */
  usrprompt(prompt);
  activ(prompt);
  /* No model */
  nom = true;
  /* Initialise plotting */
  plact();
  /* Defaults */
  rtol = 0.05;
  levc = 1;
  levb = 1;
  debug = false;
  admin = false;
  dump = false;
  prof = false;
  remhole = true;
  diffbad = true;
  ndiffstep = 5;
  refcrit = 0.3;
  derefcrit = 0.2;
  correct = true;
  defrag = true;
  regrid = true;
  repart = true;
  timpart = true;
  workwt = 1.0;
  fast = true;
  superbee = false;
  mc = false;
  minmod = false;
  vanleer = false;
  cputime = big;
  normrat = 0.25;
  nrel = 5;
  nomug = true;
  mugprt = false;
  mugi = false;
  mcnv = 0.9;
  fdq = 0.01;
  maxmug = 5;
  grav = false;
  boundpot = false;
  /* Default scalars */
  nal = 1;
  scad = 1;
  /* Default no turbulence */
  keps = false;
  /* User defaults */
  usrdef();
  /* Command loop */
  more = true;
  while (more){
    /* Get next command */
    if (get_cmd(key)){
      /* Quit and return to caller */
      if (strcmp(key, "Q") == 0){
	/* Clear model memory */
	if (!nom)
	  clear();
	more = false;
      }
      /* Other commands */
      else{
	/*Generic-specific commands (system, $LOG etc) */
	if (gen_cmd(key));
	/* Plotting commands */
	else if (plcmd(key));
	/* Muvi commands */
	else if (mvcmd(key));
	/* PMG_M commands proper */
      	else if (command(key));
	/* User commands */
	else if (usrcmd(key));
	else{
	  /* Command not found */
	  if (mype == 0)
	    cout << "unknown command:  " << key << '\n';
	  report();
	}
      }
    }
  }
  /* Deactivate plotting */
  pldac();
  /* Deactivate command interpreter */
  deactiv();
  return;
}
//end control.c

//file command.c
/*
+ command - model command interface
    Description :
     This checks to see if it  knows the command string in key, and if
     necessary performs the associated actions, returning true to stop
     other command menus from seeing it.
*/
bool command(char key[])
{
  int nsteps, i, ilev, id, npe, nvals[MAXPE], count, nchar;
  float xc[MAXD], nxl[MAXD], nxr[MAXD], avspeed, *temp;
  bool cmd, found, fail;
  found = true;
  cmd = true;
  /* [no]admin - plot every step */
  if (strcmp(key, "ADMIN") == 0)
    admin = true;
  else if (strcmp(key, "NOADMIN") == 0) 
    admin = false;
  /* [no]debug - switch debug mode (print every step) */
  else if (strcmp(key, "DEBUG") == 0) 
    debug = true;
  else if (strcmp(key, "NODEBUG") == 0) 
    debug = false;
  /* [no]mugprt - switch mugprt mode (print every mug cycle) */
  else if (strcmp(key, "MUGPRT") == 0) 
    mugprt = true;
  else if (strcmp(key, "NOMUGPRT") == 0) 
    mugprt = false;
  /* [no]dump - switch dump mode (print every step) */
  else if (strcmp(key, "DUMP") == 0) 
    dump = true;
  else if (strcmp(key, "NODUMP") == 0) 
    dump = false;
  /* [no]profile - switch on verbose timing information */
  else if (strcmp(key, "PROF") == 0) 
    prof = true;
  else if (strcmp(key, "NOPROF") == 0) 
    prof = false;
  else if (strcmp(key, "FAST") == 0) 
    fast = true;
  else if (strcmp(key, "SLOW") == 0) 
    fast = false;
  else if (strcmp(key, "DIFFBAD") == 0){
    if (!get_i(ndiffstep)){
      cout << "Error : reading ndiffstep" << '\n';
      report();
    }
    else
     diffbad = true;
  }
  else if (strcmp(key, "NODIFFBAD") == 0)
    diffbad = false;
  else if (strcmp(key, "REMHOLE") == 0)
    remhole = true;
  else if (strcmp(key, "NOREMHOLE") == 0)
    remhole = false;
  else if (strcmp(key, "CORRECT") == 0)
    correct = true;
  else if (strcmp(key, "NOCORRECT") == 0)
    correct = false;
  else if (strcmp(key, "DEFRAG") == 0)
    defrag = true;
  else if (strcmp(key, "NODEFRAG") == 0)
    defrag = false;
  else if (strcmp(key, "REGRID") == 0)
    regrid = true;
  else if (strcmp(key, "NOREGRID") == 0)
    regrid = false;
  else if (strcmp(key, "REPART") == 0)
    repart = true;
  else if (strcmp(key, "NOREPART") == 0)
    repart = false;
  else if (strcmp(key, "TIMPART") == 0)
    timpart = true;
  else if (strcmp(key, "NOTIMPART") == 0)
    timpart = false;
  else if (strcmp(key, "SUPERBEE") == 0) {
    superbee = true;
    mc = false;
    vanleer = false;
    minmod = false;
  }
  else if (strcmp(key, "MINMOD") == 0) {
    superbee = false;
    mc = false;
    vanleer = false;
    minmod = true;
  }
  else if (strcmp(key, "MC") == 0) {
    superbee = false;
    mc = true;
    vanleer = false;
    minmod = false;
  }
  else if (strcmp(key, "VANLEER") == 0) {
    superbee = false;
    mc = false;
    vanleer = true;
    minmod = false;
  }
  else if (strcmp(key, "VANALBADA") == 0) {
    superbee = false;
    mc = false;
    vanleer = false;
    minmod = false;
  }
  else if (strcmp(key, "INCLEV") == 0){
    if (nlevs == MAXLEV){
      cout << "Error : nlevs > MAXLEV" << '\n';
      report();
    }
    else
      inc_lev(true);
  }
  else if (strcmp(key, "MGIMP") == 0){
    if (nlevs < 3){
      mugi = false;
      if (mype == 0)
	cout << "Error: need at least 3 levels for mutligrid implicit \n";
      report();
    }
    else if (order == 1){
      mugi = false;
      if (mype == 0)
	cout << "Error: need at least 3 levels for mutligrid implicit \n";
      report();
    }
    else {
      mugi = true;
      nomug = true;
    }
  }
  else if (strcmp(key, "NOMGIMP") == 0){
    mugi = false;
    clear_mug();
  }  
  else if (strcmp(key, "GRAV") == 0){
    grav = true;
    if (usroff()){
      clear_mug();
      if (!nom){
	/* Zero potential */
	for (i = 1; i <= celltop; i++)
	  carr[i]->qa[nq-1] = 0.0;
      }
    }
    else{
      grav = false;
      cout << "Error : cannot set gravity" << '\n';
    }
  }
  else if (strcmp(key, "NOGRAV") == 0){
    grav = false;
    usroff();
    clear_mug();
  }
  else if (strcmp(key, "BOUNDPOT") == 0)
    boundpot = true;
  else if (strcmp(key, "NOBOUNDPOT") == 0)
    boundpot = false;
  else if (strcmp(key, "DECLEV") == 0){
    if (nlevs == 1){
      cout << "Error : only one level" << '\n';
      report();
    }
    else
      dec_lev();
  }
  else if (strcmp(key, "REMG0") == 0){
    if (nlevs == 1){
      cout << "Error : only one level" << '\n';
      report();
    }
    else
      rem_g0();
  }
  else if (strcmp(key, "NEWDOM") == 0){
// TODO something not working right when trying to extend in y in parallel
    found = true;
    for (id = 1; id <= nd; id++){
      if (!get_r(nxl[id-1])){
	if (mype == 0)
	  cout << "Error : reading left coord " << id  << '\n';
	found = false;
	break;
      }
    }
    for (id = 1; id <= nd; id++){
      if (!get_r(nxr[id-1])){
	if (mype == 0)
	  cout << "Error : reading right coord " << id  << '\n';
	found = false;
	break;
      }
    }
    if (found){
      if (nom){
	cout << "no model" << '\n';
	report();
      }
      else{
	/* Checks */
	fail = false;
        int extid = -1;
	for (id = 1; id <= nd; id++){
	  if (nxl[id-1] > nxr[id-1]){
	    fail = true;
	    if (mype == 0)
	      cout << "Error : left boundary > right boundary for dir "
		   << id << '\n';
	    fail = (bool)global_max((float)fail);
	    break;
	  }
          /* first check if more than one extension in left BM*/
          if (nxl[id - 1] != xl[id - 1] && extid > 0)	{	// this left coord has already been changed
            fail = true;
	    if (mype == 0)
              cout << "Error : Can only extend in one direction at a time!" << '\n';
	    fail = (bool)global_max((float)fail);
            break;
          }
          else if (nxl[id - 1] != xl[id - 1] && extid < 0) {	// this left coord has just been changed
            extid = id;		// save this id for comparing right
          }
          /* now check if more than one extension in right BM*/
          if (nxr[id - 1] != xr[id - 1] && extid < 0)	{	// no change in left or right yet
            extid = id;
          }
          else if (nxr[id - 1] != xr[id - 1] && id != extid) {	// wrong direction has changed
            fail = true;
	    if (mype == 0)
              cout << "Error : Can only extend in one direction at a time!" << '\n';
	    fail = (bool)global_max((float)fail);
            break;
          } 
	}
	if (fail)
	  report();
	else{
	  /* Do each boundary individually */
	  for (id = 1; id <= nd; id++){
	    /* Left boundaries*/
	    if (!(bool)global_max((float)newdom(nxl, xr))){
	      if (mype == 0)
		cout << "Error : failure in newdom " << '\n';
	      break;
	      report();
	    }
	    /* Right boundaries*/
	    if (!(bool)global_max((float)newdom(nxl, nxr))){
	      if (mype == 0)
		cout << "Error : failure in newdom " << '\n';
	      break;
	      report();
	    }       
	  }
	}
      }
    }
  }
  else if (strcmp(key, "BOUNDL") == 0){
    if (!get_i(id)){
      if (mype == 0)
	cout << "Error : reading DIRECTION" << '\n';
      report();
    }
    else if (!get_cmd(key)){
      if (mype == 0)
	cout << "Error : reading boundary condition" << '\n';
      report();
    }
    else if ((id < 1) || (id > nd)){
	if (mype == 0)
	  cout << "Error : invalid DIRECTION" << '\n';
	report();
    }
    else{
      if (strcmp(key, "FREE") == 0){
	bcl[id-1] = bcfree;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "FREE")
	    bcond[id-1] = "FREE";
	  else
	    bcond[id-1] = "FREE-" + bcond[id-1];
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(nchar+1);
	  if (bcond[id-1] == "FREE")
	    bcond[id-1] = "FREE";
	  else
	    bcond[id-1] = "FREE-" + bcond[id-1];
	}
      }
      else if (strcmp(key, "SYMMETRIC") == 0){
	bcl[id-1] = bcref;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "SYMMETRIC")
	    bcond[id-1] = "SYMMETRIC";
	  else
	    bcond[id-1] = "SYMMETRIC-" + bcond[id-1];
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(nchar+1);
	  if (bcond[id-1] == "SYMMETRIC")
	    bcond[id-1] = "SYMMETRIC";
	  else
	    bcond[id-1] = "SYMMETRIC-" + bcond[id-1];
	}
      }
      else if (strcmp(key, "FIX") == 0){
	bcl[id-1] = bcfix;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "FIX")
	    bcond[id-1] = "FIX";
	  else
	    bcond[id-1] = "FIX-" + bcond[id-1];
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(nchar+1);
	  if (bcond[id-1] == "FIX")
	    bcond[id-1] = "FIX";
	  else
	    bcond[id-1] = "FIX-" + bcond[id-1];
	}
      }
      else if (strcmp(key, "WALL") == 0){
	bcl[id-1] = bcwall;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "WALL")
	    bcond[id-1] = "WALL";
	  else
	    bcond[id-1] = "WALL-" + bcond[id-1];
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(nchar+1);
	  if (bcond[id-1] == "WALL")
	    bcond[id-1] = "WALL";
	  else
	    bcond[id-1] = "WALL-" + bcond[id-1];
	}
      }
      else if (strcmp(key, "ISOTHERMAL") == 0){
	bcl[id-1] = bcisotherm;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "ISOTHERMAL")
	    bcond[id-1] = "ISOTHERMAL";
	  else
	    bcond[id-1] = "ISOTHERMAL-" + bcond[id-1];
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(nchar+1);
	  if (bcond[id-1] == "ISOTHERMAL")
	    bcond[id-1] = "ISOTHERMAL";
	  else
	    bcond[id-1] = "ISOTHERMAL-" + bcond[id-1];
	}
      }
      else if (strcmp(key, "PRESSURE") == 0){
	bcl[id-1] = bcpressure;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "PRESSURE")
	    bcond[id-1] = "PRESSURE";
	  else
	    bcond[id-1] = "PRESSURE-" + bcond[id-1];
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(nchar+1);
	  if (bcond[id-1] == "PRESSURE")
	    bcond[id-1] = "PRESSURE";
	  else
	    bcond[id-1] = "PRESSURE-" + bcond[id-1];
	}
      }
      else if (strcmp(key, "PERIODIC") == 0){
	bcl[id-1] = bcper;
	bcr[id-1] = bcper;
	bcond[id-1] = "PERIODIC";
	de_frag();
      }
    }
  }
  else if (strcmp(key, "BOUNDR") == 0){
    if (!get_i(id)){
      if (mype == 0)
	cout << "Error : reading DIRECTION" << '\n';
      report();
    }
    else if (!get_cmd(key)){
      if (mype == 0)
	cout << "Error : reading boundary condition" << '\n';
      report();
    }
    else if ((id < 1) || (id > nd)){
	if (mype == 0)
	  cout << "Error : invalid DIRECTION" << '\n';
	report();
    }
    else{
      if (strcmp(key, "FREE") == 0){
	bcr[id-1] = bcfree;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "FREE")
	    bcond[id-1] = "FREE";
	  else
	    bcond[id-1] = bcond[id-1] + "-FREE";
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(0, nchar);
	  if (bcond[id-1] == "FREE")
	    bcond[id-1] = "FREE";
	  else
	    bcond[id-1] = bcond[id-1] + "-FREE";
	}
      }
      else if (strcmp(key, "SYMMETRIC") == 0){
	bcr[id-1] = bcref;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "SYMMETRIC")
	    bcond[id-1] = "SYMMETRIC";
	  else
	    bcond[id-1] = bcond[id-1] + "-SYMMETRIC";
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(0, nchar);
	  if (bcond[id-1] == "SYMMETRIC")
	    bcond[id-1] = "SYMMETRIC";
	  else
	    bcond[id-1] = bcond[id-1] + "-SYMMETRIC";
	}
      }
      else if (strcmp(key, "FIX") == 0){
	bcr[id-1] = bcfix;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "FIX")
	    bcond[id-1] = "FIX";
	  else
	    bcond[id-1] = bcond[id-1] + "-FIX";
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(0, nchar);
	  if (bcond[id-1] == "FIX")
	    bcond[id-1] = "FIX";
	  else
	    bcond[id-1] = bcond[id-1] + "-FIX"; 
	}
      }
      else if (strcmp(key, "WALL") == 0){
	bcr[id-1] = bcwall;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "WALL")
	    bcond[id-1] = "WALL";
	  else
	    bcond[id-1] = bcond[id-1] + "-WALL";
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(0, nchar);
	  if (bcond[id-1] == "WALL")
	    bcond[id-1] = "WALL";
	  else
	    bcond[id-1] = bcond[id-1] + "-WALL";
	}
      }
      else if (strcmp(key, "ISOTHERMAL") == 0){
	bcr[id-1] = bcisotherm;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "ISOTHERMAL")
	    bcond[id-1] = "ISOTHERMAL";
	  else
	    bcond[id-1] = bcond[id-1] + "-ISOTHERMAL";
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(0, nchar);
	  if (bcond[id-1] == "ISOTHERMAL")
	    bcond[id-1] = "ISOTHERMAL";
	  else
	    bcond[id-1] = bcond[id-1] + "-ISOTHERMAL";
	}
      }
      else if (strcmp(key, "PRESSURE") == 0){
	bcr[id-1] = bcpressure;
	nchar = bcond[id-1].find("-");
	if (nchar <= 0){
	  if (bcond[id-1] == "PRESSURE")
	    bcond[id-1] = "PRESSURE";
	  else
	    bcond[id-1] = bcond[id-1] + "-PRESSURE";
	}
	else{
	  bcond[id-1] = bcond[id-1].substr(0, nchar);
	  if (bcond[id-1] == "PRESSURE")
	    bcond[id-1] = "PRESSURE";
	  else
	    bcond[id-1] = bcond[id-1] + "-PRESSURE";
	}
      }
      else if (strcmp(key, "PERIODIC") == 0){
	bcl[id-1] = bcper;
	bcr[id-1] = bcper;
	bcond[id-1] = "PERIODIC";
	de_frag();
      }
    }
  }
  else if (strcmp(key, "ORDER") == 0){
    if (!get_i(order)){
      if (mype == 0)
	cout << "Error : reading ORDER" << '\n';
      report();
    }
  }
  else if (strcmp(key, "LEVC") == 0){
    if (!get_i(levc)){
      if (mype == 0)
	cout << "Error : reading LEVC" << '\n';
      report();
    }
    levc = max(1, levc);
    if (!nom){
      levc = max(1, min(nlevs-1,levc));      
      /* Regrid */
      for (ilev = 0; ilev <= levc; ilev++)
	re_grid(ilev);
    }
  }
  else if (strcmp(key, "LEVB") == 0){
    if (!get_i(levb)){
      if (mype == 0)
	cout << "Error : reading LEVB" << '\n';
      report();
    }
  }
  else if (strcmp(key, "START") == 0){
    /* Clear previous model */
    if (!nom)
      clear();
    nom = !start();
    if (nom){
      if (mype == 0)
	cout << "Unable to create model" << '\n';
      report();
    }
    else
      set_alev();
  }
  else if (strcmp(key, "WRITE") == 0){
    /* Write model to file */
    if (nom){
      if (mype == 0)
	cout << "no model" << '\n';
      report();
    }
    else{
      if (!get_n(key)){
	if (mype == 0)
	  cout << "Error : reading filename" << '\n';
	report();
      }
      else{
	if (writem(key));
	else{
	  if (mype == 0)
	    cout << "Error : failed to write model" << '\n';
	  report();
	}
      }
    }
  }
  else if (strcmp(key, "READ") == 0){
    /* Read model from file */
    if (!get_n(key)){
      if (mype == 0)
	cout << "Error : reading filename" << '\n';
      report();
    }
    else{
      usrquit();
      if (readm(key)){
	nom = false;
	set_alev();
      }
      else{
	if (mype == 0)
	  cout << "Error : failed to read model" << '\n';
	report();
      }
    }
  }
  else if (strcmp(key, "HELP") == 0){
    /* Help */
    if (!get_cmd(key)){
      if (mype == 0)
	cout << "Error : reading topic" << '\n';
      report();
    }
    else
      help(key);
  }
  else if (strcmp(key, "MARCH") == 0){
    if (!get_i(nsteps)){
      if (mype == 0)
	cout << "Error : reading no of steps" << '\n';
      report();
    }
    else{
      if (nom){
	if (mype == 0)
	  cout << "no model" << '\n';
	report();
      }
      else
	march(nsteps);
    }
  }
  else if (strcmp(key, "SHOW") == 0){
    if (nom){
      if (mype == 0)
	cout << "no model" << '\n';
      report();
    }
    else if (mype == 0){
      cout << title << '\n';
      cout.setf(ios::scientific);
      cout << setprecision(4);
      for (id = 1; id <= nd; id++){
	cout << setw(1) << right << coord[id-1] << " " << setw(5)
	     << nc[id-1] << " " <<  setw(12)
	     << xl[id-1] << " " << xr[id-1]
	     << " " << setw(20) << left << bcond[id-1] << '\n';
      }
      cout << "case " << cas << '\n';
      cout << "geom " << geom << '\n';
      cout << "finest level " << nlevs - 1 << " finest active level " <<
	alev <<  '\n';
      cout << "order " << order << '\n';
      cout << "cn " << cn << " cnv " << cnv << " cng " << cng << '\n';
      cout << "cnf " << cnf << '\n';
      cout << "target " << target << '\n';
      cout << setprecision(8);
      cout  << "time " << t << '\n';
      cout << setprecision(4);
      cout  << "cputime " << cputime << '\n';
    }
  }
  else if (strcmp(key, "PARTITION") == 0){
    if (nom){
      if (mype == 0)
	cout << "no model" << '\n';
      report();
    }
    else{
      /* Get partition info from other processors */
      for (npe = 0; npe < limpe; npe++)
	nvals[npe] = 2*nlevs;
      count = -1;
      /* Create temp array */
      temp = new float[2*nlevs*limpe];
      for (ilev = 0; ilev < nlevs; ilev++){
	count++;
	temp[count] = lev[ilev]->xl[nd-1];
	count++;
	temp[count] = lev[ilev]->xr[nd-1];
      }
      merge_array(temp, nvals);
      if (mype == 0){
	cout.setf(ios::scientific);
	cout << setprecision(4);
	for (npe = 0; npe < limpe; npe++){
	  for (ilev = 0; ilev < nlevs; ilev++){
	    cout << "Proc " << npe << " level " << ilev << "  left = "
	      << setw(11) << temp[2*npe*nlevs+2*ilev] << 
              "  right = " << temp[2*npe*nlevs+2*ilev+1]
	    << '\n';
	  }
	}
      }
      /* Delete temp array */
      delete [] temp;
    }
  }
  else if (strcmp(key, "WORK") == 0){
    if (nom){
      if (mype == 0)
	cout << "no model" << '\n';
      report();
    }
    else{
      if (mype == 0){
	for (npe = 0; npe < limpe; npe++)
	  cout << "Proc " << npe << " work " << procwork[npe] << '\n';
      }
    }
  }
  else if (strcmp(key, "USERS") == 0){
    if (mype == 0){
      if (nom){
	cout << "no model" << '\n';
	report();
      }
      else
	usrprt();
    }
  }
  else if (strcmp(key, "CONTROLS") == 0){
    if (nom){
      cout << "no model" << '\n';
      report();
    }
    else{
      if (mype == 0){
	cout << "rtol " << rtol << " refcrit " << refcrit <<
	  " derefcrit " << derefcrit << " levc " << levc << '\n';
	if (diffbad)
	  cout << "diffbad ndiffstep " << ndiffstep << '\n';
	else
	  cout << "no diffbad " << '\n';
	if (remhole)
	  cout << "remhole" << '\n';
	else
	  cout << "no remhole" << '\n';
	if (correct)
	  cout << "correct" << '\n';
	else
	  cout << "no correct" << '\n';
	if (regrid)
	  cout << "regrid" << '\n';
	else
	  cout << "no regrid" << '\n';
	if (repart)
	  cout << "repart" << '\n';
	else
	  cout << "no repart" << '\n';
	if (timpart)
	  cout << "timpart" << '\n';
	else
	  cout << "no timpart" << '\n';
	if (debug)
	  cout << "debug" << '\n';
	else
	  cout << "no debug" << '\n';
	if (dump)
	  cout << "dump" << '\n';
	else
	  cout << "no dump" << '\n';
      }
    }
  }
  else if (strcmp(key, "SUM") == 0)
    usrsum();
  else if (strcmp(key, "PSIZE") == 0)
    parc_size();
  else if (strcmp(key, "PRESET") == 0)
    parc_reset();
  else if (strcmp(key, "NUPRESET") == 0)
    tnupdat = 0;
  else if (strcmp(key, "STATS") == 0){
    if (nom){
      if (mype == 0)
	cout << "no model " << '\n';
    }
    else{
      mem_calc(true);
      if (mype == 0){
	if (tcpu < tiny)
	  avspeed = 0.0;
	else
	  avspeed = (float)(tnupdat/tcpu);
	cout << "cpu time " << tcpu << " cell updates "
	     << (double)tnupdat << " average speed " << avspeed << '\n';
      }
    }
  }
  else if (strcmp(key, "GSTATS") == 0){
    if (nom){
      if (mype == 0)
	cout << "no model " << '\n';
    }
    else{
      mem_calc(false);
      if (mype == 0){
	if (tcpu < tiny)
	  avspeed = 0.0;
	else
	  avspeed = (float)(tnupdat/tcpu);
	cout << "cpu time " << tcpu << " cell updates "
	     << (double)tnupdat << " average speed " << avspeed << '\n';
      }
    }
  }
  else if (strcmp(key, "CN") == 0){
    if (!get_r(cn)){
      if (mype == 0)
	cout << "Error : reading CN" << '\n';
      report();
    }
  }
  else if (strcmp(key, "MCNV") == 0){
    if (!get_r(mcnv)){
      if (mype == 0)
	cout << "Error : reading MCNV" << '\n';
      report();
    }
  }
  else if (strcmp(key, "FDQ") == 0){
    if (!get_r(fdq)){
      if (mype == 0)
	cout << "Error : reading FDQ" << '\n';
      report();
    }
  }
  else if (strcmp(key, "MAXMUG") == 0){
    if (!get_i(maxmug)){
      if (mype == 0)
	cout << "Error : reading MAXMUG" << '\n';
      report();
    }
  }    
  else if (strcmp(key, "CNV") == 0){
    if (!get_r(cnv)){
      if (mype == 0)
	cout << "Error : reading CNV" << '\n';
      report();
    }
  }
  else if (strcmp(key, "CNG") == 0){
    if (!get_r(cng)){
      if (mype == 0)
	cout << "Error : reading CNG" << '\n';
      report();
    }
  }
  else if (strcmp(key, "CNF") == 0){
    if (!get_r(cnf)){
      if (mype == 0)      
	cout << "Error : reading CNF" << '\n';
      report();
    }
  }
  else if (strcmp(key, "NORMRAT") == 0){
    if (!get_r(normrat)){
      if (mype == 0)      
	cout << "Error : reading NORMRAT" << '\n';
      report();
    }
  }
  else if (strcmp(key, "NDIFFSTEP") == 0){
    if (!get_i(ndiffstep)){
      if (mype == 0)      
	cout << "Error : reading no of diff steps" << '\n';
      report();
    }
  }
  else if (strcmp(key, "REFCRIT") == 0){
    if (!get_r(refcrit)){
      if (mype == 0)            
	cout << "Error : reading refinement criterion" << '\n';
      report();
    }
  }
  else if (strcmp(key, "DEREFCRIT") == 0){
    if (!get_r(derefcrit)){
      if (mype == 0)            
	cout << "Error : reading derefinement criterion" << '\n';
      report();
    }
  }
  else if (strcmp(key, "NREL") == 0){
    if (!get_i(nrel)){
      if (mype == 0)            
	cout << "Error : reading NREL" << '\n';
      report();
    }
  }
  else if (strcmp(key, "TARGET") == 0){
    if (!get_r(target)){
      if (mype == 0)
	cout << "Error : reading TARGET" << '\n';
      report();
    }
  }
  else if (strcmp(key, "CPUTIME") == 0){
    if (!get_r(cputime)){
      if (mype == 0)
	cout << "Error : reading CPUTIME" << '\n';
      report();
    }
  }
  else if (strcmp(key, "RTOL") == 0){
    if (!get_r(rtol)){
      if (mype == 0)            
	cout << "Error : reading RTOL" << '\n';
      report();
    }
  }
  else if (strcmp(key, "LEV") == 0){
    if (!get_i(plmin)){
      if (mype == 0)            
	cout << "Error : reading plot level" << '\n';
      report();
    }
    else{
      if (plmin < 0){
	if (mype == 0)	
	  cout << "Error : plot level < 0" << '\n';
	report();
	plmin = 0;
      }
      else if (plmin > nlevs -1){
	if (mype == 0)      	
	  cout << "Error : plot level  > max level" << '\n';
	report();
	plmin = nlevs - 1;
      }
      else
	plmax = plmin;
    }
  }
  else if (strcmp(key, "LEVS") == 0){
    if (!get_i(plmin)){
      if (mype == 0)            
	cout << "Error : reading plot level" << '\n';
      report();
    }
    else if (!get_i(plmax)){
      if (mype == 0)            
	cout << "Error : reading plot level" << '\n';
      report();
    }
    else{
      if (plmin < 0){
	cout << "Error : plot level < 0" << '\n';
	report();
	plmin = 0;
      }
      else if (plmax > nlevs -1){
	if (mype == 0)      	
	  cout << "Error : plot level  > max level" << '\n';
	report();
	plmax = nlevs - 1;
      }
    }
  }
  else if (strcmp(key, "PEEKC") == 0){
    /* Show cell info */
    if (!get_i(ilev)){
      if (mype == 0)            
	cout << "Error : reading LEVEL" << '\n';
      report();
    }
    else{
      found = true;
      for (id = 1; id <= nd; id++){
	if (!get_r(xc[id-1])){
	  if (mype == 0)
	    cout << "Error : reading coord " << id  << '\n';
	  found = false;
	  report();
	}
      }
    }
    if (found){
      if (nom){
	if (mype == 0)	
	  cout << "no model" << '\n';
	report();
      }
      else{
	if (ilev <= nlevs -1)
	  peek_cell(ilev, xc);
	else{
	  if (mype == 0)      	  
	    cout << "Finest level is " << nlevs -1 << '\n';
	  report();
	}
      }
    }
  }
  else
    /* Command not known here, pass it back to caller */
    cmd = false;
  return(cmd);
}
//end command.c
