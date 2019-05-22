/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

#include <iomanip>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>

using namespace amrex;

Real Castro::lum_bol = 0.0;
Real Castro::radius_cd  = 0.0;
Real Castro::radius_csm = 0.0;
Real Castro::radius_rs = 0.0;
Real Castro::lum_fs = 0.0;
Real Castro::lum_rs = 0.0;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    bolometric_luminosity(time,lum_bol);
    fs_lum(time,lum_fs);
    rs_lum(time,lum_rs);
    cd_shock_radius(time,radius_cd);
    csm_edge_radius(time,radius_csm);
    reverse_shock_radius(time,radius_rs);


    ParallelDescriptor::ReduceRealMax(lum_bol);
    ParallelDescriptor::ReduceRealMax(radius_cd);
    ParallelDescriptor::ReduceRealMax(radius_csm);
    ParallelDescriptor::ReduceRealMax(lum_fs);
    ParallelDescriptor::ReduceRealMax(lum_rs);
    ParallelDescriptor::ReduceRealMax(radius_rs);


    if(amrex::ParallelDescriptor::IOProcessor()){

      int dataprecision = 16; // Number of digits after the decimal point, for float data

      int datwidth      = 25; // Floating point data in scientific notation
      int fixwidth      = 25; // Floating point data not in scientific notation
      int intwidth      = 12; // Integer data

      std::ostream& log = parent->DataLog(1);

      log << std::fixed;

      log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;
      log << std::scientific;
      log << std::setw(datwidth) << std::setprecision(dataprecision) << lum_bol;
      log << std::setw(datwidth) << std::setprecision(dataprecision) << lum_fs;
      log << std::setw(datwidth) << std::setprecision(dataprecision) << lum_rs;
      log << std::setw(datwidth) << std::setprecision(dataprecision) << radius_cd;
      log << std::setw(datwidth) << std::setprecision(dataprecision) << radius_csm;
      log << std::setw(datwidth) << std::setprecision(dataprecision) << radius_rs;

      log << std::endl;

      std::cout << "Bolometric Luminosity    =   " << lum_bol << "      erg/s    " << "\n";
      std::cout << "Forward Shock Luminosity    =   " << lum_fs << "      erg/s    " << "\n";
      std::cout << "Reverse Shock Luminosity    =   " << lum_rs << "      erg/s    " << "\n";
      std::cout << "Forward Shock Radius    =   " << radius_cd << "      cm    " << "\n";
      std::cout << "Reverse Shock Radius    =   " << radius_rs << "      cm    " << "\n";
      std::cout << "Outer CSM Radius  =  " << radius_csm << "    cm     " << "\n";

    }

    lum_bol = 0.;
    lum_fs = 0.;
    radius_cd = 0.;
    radius_csm = 0.;
    lum_rs = 0.;
    radius_rs = 0.;

}
#endif

void
Castro::bolometric_luminosity(Real time, Real& lum)
{
  const Real* dx = geom.CellSize();


  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    //ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    //Castro& c_lev = getLevel(lev);

    MultiFab& Er = getLevel(lev).get_new_data(Rad_Type);

    for(MFIter mfi(Er,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        bolometric_lum(BL_TO_FORTRAN_3D(Er[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &lum);

            }
  }
}

void
Castro::fs_lum(Real time, Real& lum)
{
  const Real* dx = geom.CellSize();

  for(int lev = 0; lev <= parent -> finestLevel(); lev++)
  {
    Castro& c_lev = getLevel(lev);
    MultiFab& Es = c_lev.get_new_data(State_Type);
    auto mfcdmask = c_lev.derive("cd_mask",time,0);

    for(MFIter mfi(*mfcdmask,true);mfi.isValid();++mfi)
    {
      const Box& box = mfi.tilebox();
      const int* lo = box.loVect();
      const int * hi = box.hiVect();
      FArrayBox& fabcdmask = (*mfcdmask)[mfi];
      lum_fs_shock(BL_TO_FORTRAN_3D(fabcdmask),
      BL_TO_FORTRAN_3D(Es[mfi]),
      ARLIM_3D(lo),ARLIM_3D(hi),
      ZFILL(dx),&time,&lum);
    }

  }
}

void
Castro::rs_lum(Real time, Real& lum)
{
  const Real* dx = geom.CellSize();

  for(int lev = 0; lev <= parent -> finestLevel(); lev++)
  {
    Castro& c_lev = getLevel(lev);
    MultiFab& Es = c_lev.get_new_data(State_Type);
	auto mfcdmask = c_lev.derive("cd_mask",time,0);
    for(MFIter mfi(*mfcdmask,true);mfi.isValid();++mfi)
    {
      const Box& box = mfi.tilebox();
      const int* lo = box.loVect();
      const int * hi = box.hiVect();
      FArrayBox& fabcdmask = (*mfcdmask)[mfi];
      lum_rs_shock(BL_TO_FORTRAN_3D(Es[mfi]),
      ARLIM_3D(lo),ARLIM_3D(hi),
      ZFILL(dx),&time,&lum);
    }

  }
}


void
Castro::cd_shock_radius(Real time, Real& radius_cd)
{
  const Real* dx = geom.CellSize();


  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    Castro& c_lev = getLevel(lev);
	auto mfcdmask = c_lev.derive("cd_mask",time,0);


    for(MFIter mfi(*mfcdmask,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
		FArrayBox& fabcdmask = (*mfcdmask)[mfi];
        cdshock(BL_TO_FORTRAN_3D(fabcdmask),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &radius_cd);

            }
  }
}


void
Castro::csm_edge_radius(Real time, Real& radius_csm)
{
  const Real* dx = geom.CellSize();


  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    Castro& c_lev = getLevel(lev);
	auto mfcsmmask = c_lev.derive("csm_mask",time,0);


    for(MFIter mfi(*mfcsmmask,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
		FArrayBox& fabcsmmask = (*mfcsmmask)[mfi];
        csm_edge(BL_TO_FORTRAN_3D(fabcsmmask),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &radius_csm);

            }
  }
}


void
Castro::reverse_shock_radius(Real time, Real& radius_rs)
{
  const Real* dx = geom.CellSize();


  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    Castro& c_lev = getLevel(lev);
	  auto mfcsmmask = c_lev.derive("csm_mask",time,0);
    MultiFab& Es = c_lev.get_new_data(State_Type);


    for(MFIter mfi(*mfcsmmask,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
        csm_edge(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &radius_rs);

            }
  }
}
