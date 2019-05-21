#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif


void bolometric_luminosity (amrex::Real time, amrex::Real& lum_bol);
void fs_lum(amrex::Real time, amrex::Real& lum_fs);
void rs_lum(amrex::Real time, amrex::Real& lum_rs);


#ifdef DO_PROBLEM_POST_TIMESTEP
void problem_post_timestep();
void cd_shock_radius(amrex::Real time, amrex::Real& radius_cd);
void csm_edge_radius(amrex::Real time, amrex::Real& radius_csm);

#endif

static amrex::Real lum_bol;
static amrex::Real radius_cd;
static amrex::Real radius_csm;
static amrex::Real radius_rs;
static amrex::Real lum_fs;
static amrex::Real lum_rs;
