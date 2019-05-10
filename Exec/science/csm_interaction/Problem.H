#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif


void bolometric_luminosity (amrex::Real time, amrex::Real& lum_bol);


#ifdef DO_PROBLEM_POST_TIMESTEP
void problem_post_timestep();
void cd_shock_radius(amrex::Real time, amrex::Real cd_radius);

#endif

static amrex::Real lum_bol;
static amrex::Real radius_cd;
