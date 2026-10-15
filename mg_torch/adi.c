#undef MAIN

void implicit_vflux(float dt, float tfrac, int ilev, int norder, float rate);


/*
////////////////////////////////////////////////////////////////

Implicit diffusion routines
The plan:  Do explicit on ilev=0
Then use course grid values as boundary conditions for implicit calculations on different processors
*/
void implicit_vflux(float dt, float tfrac, int ilev, int norder, float rate){


  /* Empty for now */

  return;
}




