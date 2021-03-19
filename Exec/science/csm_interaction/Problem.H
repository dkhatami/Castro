

#ifndef DO_PROBLEM_POST_INIT
#define DO_PROBLEM_POST_INIT
#endif

#ifndef DO_PROBLEM_POST_RESTART
#define DO_PROBLEM_POST_RESTART
#endif


#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif

void bolometric_lum(amrex::Real time, amrex::Real dt,int ind_phot,int ind_surf,int ind_col, amrex::Real& lum_phot,amrex::Real& lum_edge,
					amrex::Real& lum_phot_lab, amrex::Real& lum_phot_com,
					amrex::Real& lum_edge_lab, amrex::Real& lum_edge_com,
					amrex::Real& lum_surf_lab, amrex::Real& lum_surf_com,
					amrex::Real& T_col,amrex::Real& v_phot, amrex::Real& rho_phot,
					amrex::Real& v_col, amrex::Real& rho_col);


void shock_properties(amrex::Real time, amrex::Real dt, int& ind_rs, int& ind_sh, int& ind_fs,
					  amrex::Real& radius_rs, amrex::Real& radius_sh, amrex::Real& radius_fs,
					  amrex::Real& vel_rs, amrex::Real& vel_sh, amrex::Real& vel_fs);

void integ_optical_depth(amrex::Real time, int &ind_tau_phot, int &ind_tau_surf,int &ind_tau_col, int ind_rs, int ind_sh, int ind_fs, amrex::Real& tau_rs, amrex::Real& tau_sh, amrex::Real& tau_fs);


void sum_energy(amrex::Real time, amrex::Real dt, int ind_rs, int ind_fs,
				amrex::Real& Ek_csm, amrex::Real& Eg_csm, amrex::Real& Er_csm,
				amrex::Real& Ek_ej,  amrex::Real& Eg_ej,  amrex::Real& Er_ej);

void problem_post_init();



void problem_post_restart();

void problem_post_timestep();


static amrex::Real lum_phot;
static amrex::Real lum_edge;
static amrex::Real T_col;
static amrex::Real rho_chrom;
static amrex::Real v_chrom;
static amrex::Real lum_phot_com;
static amrex::Real lum_phot_lab;
static amrex::Real lum_edge_com;
static amrex::Real lum_edge_lab;
static amrex::Real lum_surf_com;
static amrex::Real lum_surf_lab;

static amrex::Real r_phot;
static amrex::Real r_col;
static amrex::Real T_eff;
static amrex::Real rho_phot;
static amrex::Real v_phot;
static int ind_tau_phot;
static  int ind_tau_surf;
static int ind_tau_col;

static amrex::Real radius_shock;
static amrex::Real vel_shock;
static amrex::Real tau_shock;
static int ind_shock;
static int ind_rs;
static int ind_fs;
static int ind_sh;
static amrex::Real radius_rs;
static amrex::Real vel_rs;
static amrex::Real radius_fs;
static amrex::Real vel_fs;
static amrex::Real radius_sh;
static amrex::Real vel_sh;
static amrex::Real tau_fs;
static amrex::Real tau_rs;
static amrex::Real tau_sh;

static amrex::Real Ek_csm;
static amrex::Real Eg_csm;
static amrex::Real Er_csm;
static amrex::Real Ek_ej;
static amrex::Real Eg_ej;
static amrex::Real Er_ej;
