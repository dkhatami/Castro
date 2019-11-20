#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
/* problem-specific stuff goes here */


  void problem_checkpoint(int * int_dir_name, int * len);

  void problem_restart(int * int_dir_name, int * len);

  void bolometric_lum(BL_FORT_FAB_ARG_3D(rad_state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      amrex::Real* lum_bol);

  void rad_temp_rs(BL_FORT_FAB_ARG_3D(rad_state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_rs, amrex::Real* temp );

  void rad_temp_cd(BL_FORT_FAB_ARG_3D(rad_state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_cd, amrex::Real* temp );

  void rad_temp_fs(BL_FORT_FAB_ARG_3D(rad_state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_fs, amrex::Real* temp );

  void rad_temp_phot(BL_FORT_FAB_ARG_3D(rad_state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_phot, amrex::Real* temp );

  void gas_temp_rs(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_rs, amrex::Real* temp );

  void gas_temp_cd(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_cd, amrex::Real* temp );

  void gas_temp_fs(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_fs, amrex::Real* temp );

  void gas_temp_phot(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_phot, amrex::Real* temp );

  void vel_rs(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_rs, amrex::Real* vel );

  void vel_cd(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_cd, amrex::Real* vel );

  void vel_fs(BL_FORT_FAB_ARG_3D(state),
      const int* lo, const int* hi,
      const amrex::Real* dx, const amrex::Real* time,
      int* ind_fs, amrex::Real* vel );

  void cdshock(BL_FORT_FAB_ARG_3D(cd_mask),
                const int* lo, const int* hi,
                const amrex::Real* dx, const amrex::Real* time,
                amrex::Real* r_cd,int* i_cd);

  void lum_fs_shock(BL_FORT_FAB_ARG_3D(state),
              const int* lo, const int* hi,
              const amrex::Real* dx, const amrex::Real* time,
              amrex::Real* lum_fs,int* ind_fs);

  void lum_rs_shock(BL_FORT_FAB_ARG_3D(state),
              const int* lo, const int* hi,
              const amrex::Real* dx, const amrex::Real* time,
              amrex::Real* lum_rs);

  void csm_edge(BL_FORT_FAB_ARG_3D(csm_mask),
                const int* lo, const int* hi,
                const amrex::Real* dx, const amrex::Real* time,
                amrex::Real* r_csm,int* i_csm);

  void rs_radius(BL_FORT_FAB_ARG_3D(state),
                 const int* lo, const int* hi,
                 const amrex::Real* dx, const amrex::Real* time,
                 amrex::Real* r_rs,int* i_rs);

  void fs_radius(BL_FORT_FAB_ARG_3D(state),
                const int* lo, const int* hi,
                const amrex::Real* dx, const amrex::Real* time,
                amrex::Real* r_fs,int* i_fs);

  void phot_velocity(BL_FORT_FAB_ARG_3D(state),
                const int* lo, const int* hi,
                const amrex::Real* dx, const amrex::Real* time,
                int* i_phot,amrex::Real* v_phot);
}
