/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

#include <iomanip>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>

using namespace amrex;

Real Castro::lum_bol = 0.0;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    bolometric_luminosity(time,lum_bol);

    ParallelDescriptor::ReduceRealMax(lum_bol);


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

      log << std::endl;

      std::cout << "Luminosity    =   " << lum_bol << "      erg/s    " << "\n";

    }
    lum_bol = 0.;
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
