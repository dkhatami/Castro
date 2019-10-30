#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif


void bolometric_luminosity (amrex::Real time, amrex::Real& lum_bol);
void fs_lum(amrex::Real time, amrex::Real& lum_fs);
void rs_lum(amrex::Real time, amrex::Real& lum_rs);


#ifdef DO_PROBLEM_POST_TIMESTEP
void problem_post_timestep();
void cd_shock_radius(amrex::Real time, amrex::Real& radius_cd);
void reverse_shock_radius(amrex::Real time, amrex::Real& radius_rs);
void csm_edge_radius(amrex::Real time, amrex::Real& radius_csm);

void integ_optical_depth(amrex::Real time);


void photospheric_radius(amrex::Real time, amrex::Real& radius_phot);

void cd_optical_depth(amrex::Real time, amrex::Real& tau_cd);
void fs_optical_depth(amrex::Real time, amrex::Real& tau_fs);
void rs_optical_depth(amrex::Real time, amrex::Real& tau_rs);



#endif

static amrex::Real lum_bol;
static amrex::Real radius_cd;
static amrex::Real radius_csm;
static amrex::Real radius_rs;
static amrex::Real lum_fs;
static amrex::Real lum_rs;
static amrex::Real radius_phot;
static amrex::Real tau_cd;
static amrex::Real tau_fs;
static amrex::Real tau_rs;
