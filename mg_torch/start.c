#undef MAIN
#include <iostream>
#include <cstring>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "para.h"
#include "funcs.h"
#include "user.h"			// specific to this MG version.
extern bool get_i(int &value);
extern bool get_r(float &value);
extern bool get_n(char key[]);
extern bool get_cmd(char key[]);
extern bool usrqis(float xc[], float qa[]);
extern bool usroff();
extern void set_proc_bound();
extern void re_part();
extern void refine(int ilev);
extern void proj_init(int ilev, int flev, struct cell *cptr);
extern void set_bound_flag(int ilev);
extern void mk_level();
extern void mk_plevel();
extern void mk_locarr();
extern void mk_memory();
extern void i_to_it(int ndim, int nt[], int i, int it[]);
extern float cell_vol(int ilev, float xc[]);
extern void report();
extern void clear();
extern void create_halo(int ilev);
extern void connect_halo(int ilev, struct cell_list *shalo);
extern void connect_cell(int ilev, int id, struct join *&start, 
  struct join *&end, struct cell *cptr);
extern void connect(int ilev);
extern void broadcast(bool *value);
extern bool check_grid();
extern void merge(int ilev);
extern float global_max(float value);
bool start();
bool create();
bool init_part();
bool create_g0();

//file start.c
/*
+ start - read start command structure (simple)
     Description :
        This reads a command structure of the form:
       start title
          levels nlevs
          geom nd
             coord nx xl xr boundaries
             ..
          case demo
       end
*/
bool start()
{
  char tstring[SZNAM], key[SZNAM];
  int id, nchar; 
  bool found;
  string token;
  if (!get_n(key)){
    if (mype == 0)
      cout << "Error:  reading title." << '\n';
    report();
    return(false);
  }
  else
    strcpy(title, key);
  if (!get_cmd(key)){
    if (mype == 0)
      cout << "Error:  levels missing." << '\n';
    report();
    return(false);
  }
  strcpy(tstring, key);
  if (strcmp(tstring, "LEVELS") != 0){
    if (mype == 0)
      cout << "Error:  wanted level command" << '\n';
    report();
    return(false);
  }
  else{
    if (!get_i(nlevs)){
      if (mype == 0)
	cout << "Error:  reading no of levels." << '\n';
      report();
      return(false);
    }
  }
  /* Check levc and MAXLEV */
  levc = min(max(nlevs-1, 1),levc);      
  if (nlevs > MAXLEV){
    if (mype == 0)
      cout << "Error:  nlevs = " << nlevs << " > MAXLEV = " << MAXLEV <<
        " edit struct.h" << '\n';
    report();
    return(false);
  }
  if (!get_cmd(key)){
    if (mype == 0)
      cout << "Error:  reading geom." << '\n';
    report();
    return(false);
  }
  else{
    strcpy(geom,key);
    if ((strcmp(geom, "XYZ") != 0)
	&& (strcmp(geom, "RZP") != 0) 
	&& (strcmp(geom, "RTP") != 0)){
      if (mype == 0) {
        cout << "Error:  unknown geometry " << geom << '\n';
      }
      report();
      return(false);
    }
    else{
      if (!get_i(nd)){
	if (mype == 0)
	  cout << "Error:  reading nd." << '\n';
	report();
	return(false);
      }
      else{
	if (nd < 1){
	  if (mype == 0)
	    cout << "Error:  nd must be > 0 " << nd << '\n';
	  report();
	  return(false);
	}
	else if (nd > MAXD){
	  if (mype == 0)
	    cout << "Error:  nd= " << nd << " > MAXD= " << MAXD << 
              " Change MAXD in struct.h " << '\n';
	  report();
	  return(false);
	}
	else{
	  for (id = 1; id <= nd; id++){
	    if (!get_cmd(key)){
	      if (mype == 0)
		cout << "Error:  reading coord " << id << '\n';
	      report();
	      return(false);
	    }
	    else{
	      strcpy(coord[id-1], key);
	      if (!get_i(nc[id-1])){
		if (mype == 0)
		  cout << "Error:  reading nc for direction " <<
		    id << '\n';
		report();
		return(false);
	      }
	      else if (!get_r(xl[id-1])){
		if (mype == 0)
		  cout << "Error:  reading xl for direction " <<
		    id << '\n';
		report();
		return(false);
	      }
	      else if (!get_r(xr[id-1])){
		if (mype == 0)
		  cout << "Error:  reading xr for direction " <<
		    id << '\n';
		report();
		return(false);
	      }
	      else if (xl[id-1] >= xr[id-1]){
		if (mype == 0)
		  cout << "Error:  illegal coordinate range" <<
		    " for id " << id << '\n';
		report();
		return(false);
	      }
	      else{
		if (!get_cmd(key)){
		  if (mype == 0)
		    cout << "Error:  reading boundary " << id << '\n';
		  report();
		  return(false);
		}
		else{
		  found = true;
		  strcpy(tstring, key);
		  token = tstring;
		  bcond[id-1] = token;
		  if (strcmp(tstring, "FREE") == 0){
		    bcl[id-1] = bcfree;
		    bcr[id-1] = bcfree;
		  }
		  else if (strcmp(tstring,  "SYMMETRIC") == 0){
		    bcl[id-1] = bcref;
		    bcr[id-1] = bcref;
		  }
		  else if (strcmp(tstring,  "WALL") == 0){
		    bcl[id-1] = bcwall;
		    bcr[id-1] = bcwall;
		  }
		  else if (strcmp(tstring,  "FIX") == 0){
		    bcl[id-1] = bcfix;
		    bcr[id-1] = bcfix;
		  }
		  else if (strcmp(tstring,  "PERIODIC") == 0){
		    bcl[id-1] = bcper;
		    bcr[id-1] = bcper;
		  }
		  else{
		    nchar = token.find("-");
		    if (token.substr(0,nchar) == "FREE")
		      bcl[id-1] = bcfree;
		    else if (token.substr(0,nchar) == "SYMMETRIC")
		      bcl[id-1] = bcref;
		    else if (token.substr(0,nchar) == "WALL")
		      bcl[id-1] = bcwall;
		    else if (token.substr(0,nchar) == "FIX")
		      bcl[id-1] = bcfix;
		    else
		      found = false;
		    if (found){
		      if (token.substr(nchar+1) == "FREE")
			bcr[id-1] = bcfree;
		      else if (token.substr(nchar+1) == "SYMMETRIC")
			bcr[id-1] = bcref;
		      else if (token.substr(nchar+1) == "WALL")
			bcr[id-1] = bcwall;
		      else if (token.substr(nchar+1) == "FIX")
			bcr[id-1] = bcfix;
		      else
			found = false;
		    }
		    if(!found){
		      if (mype == 0)
			cout << "Error:  unknown boundary type." 
			     << '\n';
		      report();
		      return(false);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  strcpy(cas, "NONE");
  /* Other commands */
  while(get_cmd(key)){
    /* Check for scalars */
    if (strcmp(key, "SCALARS") == 0){
      if (!get_i(nal)){
	if (mype == 0)
	  cout << "Error:  reading no of scalars." << '\n';
	report();
	return(false);
      }
      scad = nal;
      // TODO copy what's done in mg_cti to automatically get the scalars
      if (nal < niso){
	if (mype == 0)
	  cout << "Error:  nal=" << nal << " niso=" << niso <<
	    '\n';
	report();
	return(false);
      }
    }
    else if (strcmp(key, "ADVECT") == 0){
      if (!get_i(scad)){
	if (mype == 0)
	  cout << "Error:  reading no of advecting scalars." << '\n';
	report();
	return(false);
      }
      if (scad > nal){
	if (mype == 0)
	  cout << "Error: no of advecting scalars > no of scalars."
	       << '\n';
	report();
	return(false);
      }
    }
    /* Case */
    else if (strcmp(key, "CASE") == 0){
      if (!get_cmd(key)){
	if (mype == 0)
	  cout << "Error: no case name" << '\n';
      }
      else
	strcpy(cas, key);
    }
    /* End command */
    else if (strcmp(key, "END") == 0)
      break;
    else{
      /* Unknown command */
      if (mype == 0)
	cout << "Error: unknown command or missing END" << '\n';
      return(false);
    }
  }
  /* Check for case */
  if (strcmp(cas, "NONE") == 0){
    if (mype == 0)
      cout << "Error: no case command" << '\n';
    report();
    return(false);
  }
  /* Offsets */
  if (!usroff()){
    report();
    return(false);
  }
  /* Set defaults */
  t = 0.0;
  nupdat = 0;
  tnupdat = 0;
  tcpu = 0.0;
  target = big;
  cn = 0.4;
  cnv = cn;
  cng = 1.0;
  if (nlevs == 1)
    cnf = 1.0;
  else
    cnf = 0.9;
  order = 2;
  plmin = 0;
  plmax = nlevs - 1;
  nrel = 5;
  /* Create model */
  return(create());
}
//end start.c

//file create.c
/*
+ create - create start model
     Description :
        This creates the model using the information from start. 
*/   
bool create()
{
  int ilev, tnlevs;
  float temp;
  (void)temp;
  struct cell *cptr;
  bool fail;
  /* Make level structures */
  mk_level();
  mk_plevel();
  /* Create ins and outs arrays */
  ins = new int[limpe];
  outs = new int[limpe];
  /* Create work arrays */
  wka = new float[nc[nd-1]*ipow(nlevs-1)+SAFETY];
  wkda = new float[nc[nd-1]*ipow(nlevs-1)+SAFETY];
  /* carr, jarr, locarr to NULL (safety) */
  carr = NULL;
  jarr = NULL;
  locarr = NULL;
  /* Initial partitioning */
  if (!init_part())
    return(false);
  /* Create coarsest level  */
  if (!create_g0()){
    clear();
    return(false);
  }
  tnlevs = nlevs;
  fail = false;
  if (nlevs != 1){
    /* Multiple levels --refine */
    for (ilev = 0; ilev < tnlevs - 1; ilev++){
      if (mype == 0)
	cout << "refining level " << ilev << '\n';
      /* Set temporary level */
      nlevs = min(ilev+2, tnlevs);
      refine(ilev);
      fail = !check_grid();
      fail = (bool)global_max((float)fail);
      if (fail)
	break;
      /* Fill ilev+1 cells with initial solution */
      cptr = lev[ilev+1]->fcell;
      while (cptr != NULL){
	usrqis(cptr->xc, cptr->qa);
	cptr = cptr->nxt;
      }
      if (fast)
	/* Fast initialisation - just use merge */
	merge(ilev);
      else{
	/* Slow initialisation - project initial solution onto ilev */
	cptr = lev[ilev]->fcell;
	while (cptr != NULL){
	  proj_init(ilev, tnlevs-1, cptr);
	  cptr = cptr->nxt;
	}
      }
      if (limpe > 1)
        re_part();
    }
  }
  /* Check for failure */
  if (fail){
    /* Clear model and return */
    clear();
    if (mype==0)
      cout << "Error: double refinements - suggest increase diffbad" 
	   << " temporarily." << '\n';
    return(false);
  }
  if (fast){
    for (ilev = nlevs - 2; ilev >= 0; ilev--)
      merge(ilev);
  }
  return(true);
}
//end create.c

//file init_part.c
/*
+ init_part - calculate initial partition.
     Description :
        This calculates  processor boundaries for  a uniform partition
        and  then  hands  over  to  set_proc_bound to  get  the  final
        partition. It also initialises workwt.
*/   
bool init_part()
{
  int npe;
  float dxnd; 
  /* Initialise workwt */
  workwt = 1.0;
  /* Initial guess for processor boundaries */
  dxnd = (xr[nd-1] - xl[nd-1])/(float)(limpe);
  for (npe = 0; npe < limpe; npe++)
    pbd[npe] = xl[nd-1] + npe*dxnd;
  pbd[limpe] = xr[nd-1];
  /* Set processor boundaries */
  set_proc_bound();
  return(true);
}
//end init_part.c

//file create_g0.c
/*
+ create_g0 - create uniform grid on level 0.
     Description :
        Creates a uniform  grid on level 0 using  the information from
        start. This may be empty.
*/   
bool create_g0()
{
  int id, icell, ncell, it[MAXD];
  float dx[MAXD], xc[MAXD], qa[MAXQ];
  /* No of cells and mesh spacing*/
  ncell = 1;
  for (id = 1; id <= nd; id++){
    ncell = ncell*lev[0]->nc[id-1];
    dx[id-1] = lev[0]->dxbar[id-1];
  }
  /* Set initial memory size */
  carrtop = ncell;
  if (nlevs > 1)
    carrtop += ipow(nd)*ncell;
  jarrtop = 2*carrtop*nd;
  /* Make memory */
  mk_memory();
  /* Start of free space in carr  and free cell list */
  celltop = ncell;
  freecell = NULL;
  if (ncell > 0){
    /* Fill array of cell pointers, set cell geometry */
    for (icell = 1; icell <= ncell; icell++){
      i_to_it(nd, nc, icell, it);
      for (id = 1; id <=nd; id++)
	carr[icell]->xc[id-1] = lev[0]->xl[id-1] 
	  + (it[id-1] -0.5)*dx[id-1];
      carr[icell]->vol = cell_vol(0, carr[icell]->xc);
      carr[icell]->empty = false;
    }
    lev[0]->fcell = carr[1];
    lev[0]->lcell = carr[ncell];
    /* Create linked list of cells */
    for (icell = 1; icell < ncell; icell++)
      carr[icell]->nxt = carr[icell+1];
    carr[ncell]->nxt = NULL;
    for (icell = ncell; icell > 1; icell--)
      carr[icell]->prev = carr[icell-1];
    carr[1]->prev = NULL;
  }
  /* Check usrqis to see if case exists */
  for (id = 1; id <= nd; id++)
    xc[id-1] = 0.0;
  if (!usrqis(xc, qa)){
    if (mype == 0)
      cout << "Error : CASE " << cas << " not found \n"; 
    return(false);
  }
  /* Fill cells with initial solution */
  for (icell = 1; icell <= ncell; icell++)
    usrqis(carr[icell]->xc, carr[icell]->qa);
  /* Create level 0 halo */
  create_halo(0);
  /* Make location array */
  mk_locarr();
  /* Connect level 0 */
  connect(0);
  /* Connect level 0 left halo cells */
  connect_halo(0, plev[0]->lhalo);
  /* Connect level 0 right halo cells */
  connect_halo(0, plev[0]->rhalo);
  /* Set level 0 boundary flags */
  set_bound_flag(0);
  return(true);
}
//end create_g0.c

