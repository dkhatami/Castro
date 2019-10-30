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
    ParallelDescriptor::ReduceRealMin(radius_cd);
    ParallelDescriptor::ReduceRealMax(radius_csm);
    ParallelDescriptor::ReduceRealMax(lum_fs);
    ParallelDescriptor::ReduceRealMax(lum_rs);
    ParallelDescriptor::ReduceRealMin(radius_rs);

    integ_optical_depth(time);


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
    radius_cd = 1e99;
    radius_csm = 0.;
    lum_rs = 0.;
    radius_rs = 1e99;

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
        rs_radius(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &radius_rs);

            }
  }
}

void
Castro::integ_optical_depth(Real time){


  const Real* dx = geom.CellSize();

  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);


    Castro& c_lev = getLevel(lev);

    BoxArray ba = c_lev.BoxArray();

    long num_grids = ba.size();

    struct GidLo{
      long gid;
      int lo;
      Real exclusive_sum;


    GidLo(long id, int smallend) : gid(id), lo(smallend), exclusive_sum(0.0) {};

    bool operator<(const GidLo& other) const
    {
      return lo > other.lo;
    }


  };


  Vector<GidLo> grids;

  for(long i=0; i < num_grids; ++i){
    grids.emplace_back(i,ba[i].smallEnd(0),1);
  }

  std::sort(grids.begin(),grids.end());

  Vector<long>grid_inv(num_grids);
  for(long i = 0; i < num_grids; ++i) grid_inv[grids[i].gid] = i;

  Vector<Real> prefix_sums(num_grids,0.0);



	  auto mfdtau = c_lev.derive("dtau",time,0);
    auto mftaur = c_lev.derive("taur",time,0);

    for(MFIter mfi(*mfdtau); mfi.isValid(); ++mfi)
    {
        long gidi = mfi.index();
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
        auto dtau_fab = mfdtau.array(mfi);
        auto dtau_psum_fab = mftaur.array(mfi);


        Real fab_prefix_sum = 0.0;
        for (int k = hi.z; k >= lo.z; --k) {
        for (int j = hi.y; j >= lo.y; --j) {
        for (int i = hi.x; i >= lo.x; --i) {

          fab_prefix_sum += dtau_fab(i,j,k);
          dtau_psum_fab(i,j,k) = fab_prefix_sum;

        }
      }
    }

    prefix_sums[gidi] = fab_prefix_sum;

            }


            ParallelDescriptor::ReduceRealSum(pefix_sums.dataPtr(),num_grids);


        for(long i = 1; i < num_grids; ++i){
          grids[i].exclusive_sum = grids[i-1].exclusive_sum + prefix_sums[grid_ind[i-1]];
        }


        for(MFIter mfi(*mftaur);mfi.isValid();++mfi)
        {
          long gidi = mfi.index();
          auto dtau_psum_fab = mftaur[mfi];
          dtau_psum_fab.plus(grids[grid_inv[gidi]].exclusive_sum,0,1);
        }
  }

}
