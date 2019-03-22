#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif


void bolometric_luminosity (amrex::Real time, amrex::Real& lum_bol);


#ifdef DO_PROBLEM_POST_TIMESTEP
void problem_post_timestep();
#endif

static amrex::Real lum_bol;
