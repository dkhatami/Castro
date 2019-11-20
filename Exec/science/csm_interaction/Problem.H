#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif


void bolometric_luminosity (amrex::Real time, amrex::Real& lum_bol);
void fs_lum(amrex::Real time, amrex::Real& lum_fs,int& ind_fs);
void rs_lum(amrex::Real time, amrex::Real& lum_rs);



#ifdef DO_PROBLEM_POST_TIMESTEP
void problem_post_timestep();
void get_shock_lums(amrex::Real time,int& ind_rs,int& ind_cd,int& ind_fs,
  amrex::Real& lum_rs_lab,amrex::Real& lum_rs_com,amrex::Real& lum_cd_lab,
  amrex::Real& lum_cd_com, amrex::Real& lum_fs_lab, amrex::Real& lum_fs_com,amrex::Real& lum_bol);
void cd_shock_radius(amrex::Real time, amrex::Real& radius_cd,int& ind_cd);
void reverse_shock_radius(amrex::Real time, amrex::Real& radius_rs, int& ind_rs);
void csm_edge_radius(amrex::Real time, amrex::Real& radius_csm,int& ind_csm);
void forward_shock_radius(amrex::Real time, amrex::Real& radius_fs,int& ind_fs);

void integ_optical_depth(amrex::Real time,int& ind_tauphot,int& ind_taucol,
                         amrex::Real& tau_fs,amrex::Real& tau_rs,amrex::Real& tau_cd,
                          int& ind_fs, int& ind_rs, int& ind_cd);

void rad_temperatures(amrex::Real time, int& ind_rs, int& ind_cd, int& ind_fs, int& ind_phot,
                    amrex::Real& temp_rs, amrex::Real& temp_cd, amrex::Real& temp_fs, amrex::Real& temp_phot);

void gas_temperatures(amrex::Real time, int& ind_rs, int& ind_cd, int& ind_fs, int& ind_phot,
                    amrex::Real& temp_rs, amrex::Real& temp_cd, amrex::Real& temp_fs, amrex::Real& temp_phot);

void shock_velocities(amrex::Real time, int& ind_rs, int& ind_cd, int& ind_fs,
                    amrex::Real& vel_rs, amrex::Real& vel_cd, amrex::Real& vel_fs);


void comoving_shock_luminosities(amrex::Real time, int& ind_rs, int& ind_cd, int& ind_fs,
                        amrex::Real& lum_rs, amrex::Real& lum_cd, amrex::Real& lum_fs);

void get_photospheric_velocity(amrex::Real time, int& ind_tauphot,
                          amrex::Real& vel_phot);


#endif

static amrex::Real lum_bol;
static amrex::Real radius_cd;
static amrex::Real radius_csm;
static amrex::Real radius_rs;
static amrex::Real radius_fs;
static amrex::Real velocity_rs;
static amrex::Real velocity_cd;
static amrex::Real velocity_fs;
static amrex::Real velocity_phot;
static amrex::Real radius_phot;
static amrex::Real tau_cd;
static amrex::Real tau_fs;
static amrex::Real tau_rs;
static amrex::Real trad_rs;
static amrex::Real trad_cd;
static amrex::Real trad_fs;
static amrex::Real tgas_rs;
static amrex::Real tgas_cd;
static amrex::Real tgas_fs;
static amrex::Real trad_phot;
static amrex::Real tgas_phot;
static amrex::Real temp_eff;
static amrex::Real temp_col;
static amrex::Real radius_col;
static amrex::Real lum_rs;
static amrex::Real lum_cd;
static amrex::Real lum_fs_lab;
static amrex::Real lum_fs_com;
static amrex::Real lum_rs_lab;
static amrex::Real lum_rs_com;
static amrex::Real lum_cd_lab;
static amrex::Real lum_cd_com;
static amrex::Real lum_fs_est;

static int ind_tauphot;
static int ind_taucol;
static int ind_fs;
static int ind_rs;
static int ind_cd;
static int ind_csm;
