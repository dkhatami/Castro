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

  void cdshock(BL_FORT_FAB_ARG_3D(cd_mask),
                const int* lo, const int* hi,
                const amrex::Real* dx, const amrex::Real* time,
                amrex::Real* r_cd);

  void lum_fs_shock(BL_FORT_FAB_ARG_3D(cd_mask),
              BL_FORT_FAB_ARG_3D(state),
              const int* lo, const int* hi,
              const amrex::Real* dx, const amrex::Real* time,
              amrex::Real* lum_fs);

  void lum_rs_shock(BL_FORT_FAB_ARG_3D(state),
              const int* lo, const int* hi,
              const amrex::Real* dx, const amrex::Real* time,
              amrex::Real* lum_rs);

  void csm_edge(BL_FORT_FAB_ARG_3D(csm_mask),
                const int* lo, const int* hi,
                const amrex::Real* dx, const amrex::Real* time,
                amrex::Real* r_csm);

  void rs_radius(BL_FORT_FAB_ARG_3D(state),
                 const int* lo, const int* hi,
                 const amrex::Real* dx, const amrex::Real* time,
                 amrex::Real* r_rs);
}
