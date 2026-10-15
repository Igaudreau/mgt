extern bool vanleer, mc, superbee, minmod;

inline float min(float a, float b)
{
 return (a < b ? a : b);
}

inline float max(float a, float b)
{
 return (a > b ? a : b);
}

inline int min(int a, int b)
{
 return (a < b ? a : b);
}

inline int max(int a, int b)
{
 return (a > b ? a : b);
}

inline float sgn(float a)
{
 return (a < 0 ? -1 : 1);
}

inline int ipow (int j)
{
  return(1 << j);
}

inline float mabs(float a)
{
  return(a > 0.0 ? a : -a);
}

inline float av(float a, float b)
{
  if (vanleer)
    return(b*a < tiny ? 0 : (a*b + mabs(a*b))/(a + b));	// Van Leer
  else if (minmod)
    return(b*a < tiny ? 0 : sgn(a)*max(0.0,min(mabs(a),b*sgn(a))));	//Minmod
  else if (superbee)
    return(b*a < tiny ? 0 : max(0.0, max(min(2.0*a,b),min(a,2.0*b))));	// Superbee
  else if (mc)
    return(b*a < tiny ? 0 : max((float)0.0,min((float)2.0*a,min((float)0.5*(a+b),(float)2.0*b))));	// MC
  else
    return(b*a < tiny ? 0 : a*b*(a + b)/ (a*a + b*b));		// Original (same as Van Albada 1)
}


