/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Radiation.H"
#include "Castro_F.H"
#include "Problem_F.H"

#include <iomanip>
#include <cmath>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>

using namespace amrex;

Real Castro::lum_bol = 0.0;
Real Castro::radius_cd  = 0.0;
Real Castro::radius_csm = 0.0;
Real Castro::radius_rs = 0.0;
Real Castro::radius_fs = 0.0;
Real Castro::velocity_phot = 0.0;
Real Castro::velocity_fs = 0.0;
Real Castro::velocity_cd = 0.0;
Real Castro::velocity_rs = 0.0;
Real Castro::temp_eff = 0.0;
Real Castro::temp_col = 0.0;
Real Castro::trad_rs = 0.0;
Real Castro::trad_cd = 0.0;
Real Castro::trad_fs = 0.0;
Real Castro::trad_phot = 0.0;
Real Castro::tgas_rs = 0.0;
Real Castro::tgas_cd = 0.0;
Real Castro::tgas_fs = 0.0;
Real Castro::tgas_phot = 0.0;
Real Castro::lum_fs_lab = 0.0;
Real Castro::lum_fs_com = 0.0;
Real Castro::lum_rs_lab = 0.0;
Real Castro::lum_rs_com = 0.0;
Real Castro::lum_cd_lab = 0.0;
Real Castro::lum_cd_com = 0.0;
Real Castro::lum_rs = 0.0;
Real Castro::tau_fs = 0.0;
Real Castro::tau_rs = 0.0;
Real Castro::tau_cd = 0.0;
int Castro::ind_tauphot = 0;
int Castro::ind_taucol = 0;
Real Castro::radius_phot = 0.0;
Real Castro::radius_col = 0.0;
int Castro::ind_fs = 0;
int Castro::ind_rs = 0;
int Castro::ind_cd = 0;
int Castro::ind_csm = 0;
Real Castro::lum_fs_est = 0.0;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep()
{

    if (level != 0) return;



    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

  //  bolometric_luminosity(time,lum_bol);
    // fs_lum(time,lum_fs);
    // rs_lum(time,lum_rs);
    cd_shock_radius(time,radius_cd,ind_cd);
    csm_edge_radius(time,radius_csm,ind_csm);
    reverse_shock_radius(time,radius_rs,ind_rs);
    forward_shock_radius(time,radius_fs,ind_fs);


//    ParallelDescriptor::ReduceRealMax(lum_bol);
    ParallelDescriptor::ReduceRealMin(radius_cd);
    ParallelDescriptor::ReduceRealMax(radius_csm);
  //  ParallelDescriptor::ReduceRealMax(lum_fs);
  //  ParallelDescriptor::ReduceRealMax(lum_rs);
    ParallelDescriptor::ReduceRealMin(radius_rs);
    ParallelDescriptor::ReduceRealMin(radius_fs);

     ParallelDescriptor::ReduceIntMin(ind_cd);
     ParallelDescriptor::ReduceIntMax(ind_csm);
     ParallelDescriptor::ReduceIntMin(ind_rs);
     ParallelDescriptor::ReduceIntMin(ind_fs);

    // ParallelDescriptor::ReduceIntMax(ind_cd);
    // ParallelDescriptor::ReduceIntMax(ind_cd);


    shock_velocities(time,ind_rs,ind_cd,ind_fs,velocity_rs,velocity_cd,velocity_fs);
    ParallelDescriptor::ReduceRealMax(velocity_rs);
    ParallelDescriptor::ReduceRealMax(velocity_cd);
    ParallelDescriptor::ReduceRealMax(velocity_fs);


    integ_optical_depth(time,ind_tauphot,ind_taucol,tau_fs,tau_rs,tau_cd,ind_fs,ind_rs,ind_cd);
    ParallelDescriptor::ReduceIntMin(ind_tauphot);
    ParallelDescriptor::ReduceIntMin(ind_taucol);
    ParallelDescriptor::ReduceRealMax(tau_fs);
    ParallelDescriptor::ReduceRealMax(tau_rs);
    ParallelDescriptor::ReduceRealMax(tau_cd);

    get_shock_lums(time,ind_rs,ind_cd,ind_fs,lum_rs_lab,lum_rs_com,lum_cd_lab,lum_cd_com,lum_fs_lab,lum_fs_com,lum_bol);
    ParallelDescriptor::ReduceRealSum(lum_rs_lab);
    ParallelDescriptor::ReduceRealSum(lum_rs_com);
    ParallelDescriptor::ReduceRealSum(lum_cd_lab);
    ParallelDescriptor::ReduceRealSum(lum_cd_com);
    ParallelDescriptor::ReduceRealSum(lum_fs_lab);
    ParallelDescriptor::ReduceRealSum(lum_fs_com);
    ParallelDescriptor::ReduceRealSum(lum_bol);

    get_photospheric_velocity(time,ind_tauphot,velocity_phot);
    ParallelDescriptor::ReduceRealMax(velocity_phot);

    rad_temperatures(time,ind_rs,ind_cd,ind_fs,ind_tauphot,trad_rs,trad_cd,trad_fs,trad_phot);
    ParallelDescriptor::ReduceRealMax(trad_rs);
    ParallelDescriptor::ReduceRealMax(trad_cd);
    ParallelDescriptor::ReduceRealMax(trad_fs);
    ParallelDescriptor::ReduceRealMax(trad_phot);

    gas_temperatures(time,ind_rs,ind_cd,ind_fs,ind_tauphot,tgas_rs,tgas_cd,tgas_fs,tgas_phot);
    ParallelDescriptor::ReduceRealMax(tgas_rs);
    ParallelDescriptor::ReduceRealMax(tgas_cd);
    ParallelDescriptor::ReduceRealMax(tgas_fs);
    ParallelDescriptor::ReduceRealMax(tgas_phot);

    fs_lum(time,lum_fs_est,ind_fs);
    ParallelDescriptor::ReduceRealMax(lum_fs_est);



    radius_phot = ind_tauphot*geom.CellSize(0);
    radius_col = ind_taucol*geom.CellSize(0);
    temp_eff = pow(lum_bol/(4.*3.14159265359*radius_phot*radius_phot*5.67e-5),0.25);
    temp_col = pow(lum_bol/(4.*3.14159265359*radius_col*radius_col*5.67e-5),0.25);







    if(amrex::ParallelDescriptor::IOProcessor()){

      int dataprecision = 16; // Number of digits after the decimal point, for float data

      int datwidth      = 25; // Floating point data in scientific notation
      int fixwidth      = 25; // Floating point data not in scientific notation
      int intwidth      = 12; // Integer data

      // logs = {luminosities,radii,velocities,taus}
      std::ostream& log_lc = parent->DataLog(1);
      std::ostream& log_lums = parent->DataLog(2);
      std::ostream& log_radii = parent->DataLog(3);
      std::ostream& log_vels = parent->DataLog(4);
      std::ostream& log_temps = parent->DataLog(5);
      std::ostream& log_taus = parent-> DataLog(6);


      if(time == dt)
      {
        log_lc << "Time (s)" << std::setw(fixwidth) << "L_bol" << std::setw(fixwidth) << "R_phot" << std::setw(fixwidth) << "T_eff" << std::setw(fixwidth) << "T_col";
        log_lums << "Time (s)"  << std::setw(fixwidth) << "L_bol" << std::setw(fixwidth) <<  "L_com_fs" << std::setw(fixwidth) << "L_com_cd" << std::setw(fixwidth) << "L_com_fs";
        log_lums << std::setw(fixwidth) << "L_lab_rs" << std::setw(fixwidth) << "L_lab_cd" << std::setw(fixwidth) << "L_lab_fs" << std::setw(fixwidth) << "L_est_fs";
        log_radii << "Time (s)" << std::setw(fixwidth)<< "R_phot" << std::setw(fixwidth) << "R_col"<< std::setw(fixwidth)  << "R_rs" << std::setw(fixwidth) << "R_cd" << std::setw(fixwidth) << "R_fs"<< std::setw(fixwidth)  << "R_csm";
        log_vels << "Time (s)"<< std::setw(fixwidth)  << "v_phot" << std::setw(fixwidth) << "v_rs" << std::setw(fixwidth) << "v_cd"<< std::setw(fixwidth)  << "v_fs";
        log_temps << "Time (s)" << std::setw(fixwidth) << "T_eff" << std::setw(fixwidth) << "T_col" << std::setw(fixwidth) << "Tg_phot" << std::setw(fixwidth) << "Tr_phot" << std::setw(fixwidth) << "Tg_fs" << std::setw(fixwidth) << "Tg_cd" << std::setw(fixwidth) << "Tg_rs";
        log_temps << std::setw(fixwidth) << "Tr_fs" << std::setw(fixwidth) << "Tr_cd" << std::setw(fixwidth) << "Tr_rs";
        log_taus << "Time (s)" << std::setw(fixwidth) << "tau_rs"<< std::setw(fixwidth)  << "tau_cd" << std::setw(fixwidth) << "tau_fs";
        log_lc << std::endl;
        log_lums << std::endl;
        log_radii << std::endl;
        log_vels << std::endl;
        log_temps << std::endl;
        log_taus << std::endl;

      }

      log_lc << std::fixed;
      log_lums << std::fixed;
      log_radii << std::fixed;
      log_vels << std::fixed;
      log_temps << std::fixed;
      log_taus << std::fixed;


      log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << time;
      log_lc << std::scientific;
      log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << lum_bol;
      log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << radius_phot;
      log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << temp_eff;
      log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << temp_col;

      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << time;
      log_lums << std::scientific;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_bol;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_rs_com;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_cd_com;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_fs_com;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_rs_lab;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_cd_lab;
      log_lums << std::setw(datwidth) << std::setprecision(dataprecision) << lum_fs_lab;
      log_lums  << std::setw(datwidth) << std::setprecision(dataprecision) << lum_fs_est;

      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << time;
      log_radii << std::scientific;
      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << radius_phot;
      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << radius_col;
      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << radius_rs;
      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << radius_cd;
      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << radius_fs;
      log_radii << std::setw(datwidth) << std::setprecision(dataprecision) << radius_csm;


      log_vels << std::setw(datwidth) << std::setprecision(dataprecision) << time;
      log_vels << std::scientific;
      log_vels << std::setw(datwidth) << std::setprecision(dataprecision) << velocity_phot;
      log_vels << std::setw(datwidth) << std::setprecision(dataprecision) << velocity_rs;
      log_vels << std::setw(datwidth) << std::setprecision(dataprecision) << velocity_cd;
      log_vels << std::setw(datwidth) << std::setprecision(dataprecision) << velocity_fs;


      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << time;
      log_temps << std::scientific;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << temp_eff;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << temp_col;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << tgas_phot;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << trad_phot;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << tgas_fs;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << tgas_cd;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << tgas_rs;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << trad_fs;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << trad_cd;
      log_temps << std::setw(datwidth) << std::setprecision(dataprecision) << trad_rs;


      log_taus << std::setw(datwidth) << std::setprecision(dataprecision) << time;
      log_taus << std::scientific;
      log_taus << std::setw(datwidth) << std::setprecision(dataprecision) << tau_rs;
      log_taus << std::setw(datwidth) << std::setprecision(dataprecision) << tau_cd;
      log_taus << std::setw(datwidth) << std::setprecision(dataprecision) << tau_fs;

    //  log << std::fixed;

    //   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;
    //   log << std::scientific;
    //   log << std::setw(datwidth) << std::setprecision(dataprecision) << lum_bol;
    // //  log << std::setw(datwidth) << std::setprecision(dataprecision) << lum_fs;
    //   log << std::setw(datwidth) << std::setprecision(dataprecision) << lum_rs;
    //   log << std::setw(datwidth) << std::setprecision(dataprecision) << radius_cd;
    //   log << std::setw(datwidth) << std::setprecision(dataprecision) << radius_csm;
    //   log << std::setw(datwidth) << std::setprecision(dataprecision) << radius_rs;

      log_lc << std::endl;
      log_lums << std::endl;
      log_radii << std::endl;
      log_vels << std::endl;
      log_temps << std::endl;
      log_taus << std::endl;


      std::cout << "Bolometric Luminosity    =   " << lum_bol << "      erg/s    " << "\n";
    //  std::cout << "Forward Shock Luminosity    =   " << lum_fs << "      erg/s    " << "\n";
      std::cout << "Reverse Shock Luminosity    =   " << lum_rs << "      erg/s    " << "\n";
      std::cout << "Forward Shock Radius    =   " << radius_fs << "      cm    " << "\n";
      std::cout << "Forward Shock Radius IND   =   " << ind_fs*geom.CellSize(0) << "      cm    " << "\n";
      std::cout << "Reverse Shock Radius    =   " << radius_rs << "      cm    " << "\n";
      std::cout << "Reverse Shock IND    =   " << ind_rs*geom.CellSize(0) << "      cm    " << "\n";
      std::cout << "CD Shock Radius    =   " << radius_cd << "      cm    " << "\n";
      std::cout << "CD Shock IND    =   " << ind_cd*geom.CellSize(0) << "      cm    " << "\n";
      std::cout << "Outer CSM Radius  =  " << radius_csm << "    cm     " << "\n";
      std::cout <<"Photospheric Radius = " << ind_tauphot*geom.CellSize(0) << "cm \n";
      std::cout <<"Photospheric Velocity = " << velocity_phot << "cm/s \n";
      std::cout <<"Chromospheric Radius = " << ind_taucol*geom.CellSize(0) << "cm \n";
      std::cout <<"Forward Shock Tau = " << tau_fs<< "\n";
      std::cout <<"Reverse Shock Tau = " << tau_rs<< "\n";
      std::cout <<"CD Shock Tau = " << tau_cd<< "\n";
      std::cout <<"Reverse Shock Temperature = " << trad_rs << "    K    " << "\n";
      std::cout <<"CD Shock Temperature = " << trad_cd << "    K    " << "\n";
      std::cout <<"Forward Shock Temperature = " << trad_fs << "    K    " << "\n";
      std::cout <<"Reverse Shock Velocity = " << velocity_rs << "    cm/s    " << "\n";
      std::cout <<"CD Shock Velocity = " << velocity_cd << "    cm/s    " << "\n";
      std::cout <<"Forward Shock Velocity = " << velocity_fs << "    cm/s    " << "\n";
      std::cout << "Reverse Shock Lab Luminosity = " << lum_rs_lab << "\n";
      std::cout << "Reverse Shock Comoving Luminosity = " << lum_rs_com << "\n";
      std::cout << "CD Shock Lab Luminosity = " << lum_cd_lab << "\n";
      std::cout << "CD Shock Comoving Luminosity = " << lum_cd_com << "\n";
      std::cout << "Forward Shock Lab Luminosity = " << lum_fs_lab << "\n";
      std::cout << "Forward Shock Comoving Luminosity = " << lum_fs_com << "\n";
      std::cout << "Forward Shock EST Luminosity = " << lum_fs_est << "\n";



    }

    lum_bol = 0.;
    lum_fs_lab = 0.;
    lum_fs_com = 0.;
    lum_rs_lab = 0.;
    lum_rs_com = 0.;
    lum_cd_lab = 0.;
    lum_cd_com = 0.;
    radius_cd = 1e99;
    radius_csm = 0.;
    lum_rs = 0.;
    ind_cd = geom.Domain().bigEnd(0);
    ind_rs = geom.Domain().bigEnd(0);
    ind_fs = geom.Domain().bigEnd(0);
    ind_csm = 0;
    radius_rs = 1e99;
    radius_fs = 1e99;
    ind_tauphot = 0;
    ind_taucol = 0;
    tau_fs = 0.0;
    tau_rs = 0.0;
    tau_cd = 0.0;
    trad_rs = 0.0;
    trad_fs = 0.0;
    trad_cd = 0.0;
    tgas_rs = 0.0;
    tgas_cd = 0.0;
    tgas_fs = 0.0;
    trad_phot = 0.0;
    tgas_phot = 0.0;
    velocity_rs = 0.0;
    velocity_cd = 0.0;
    velocity_fs = 0.0;
    velocity_phot = 0.0;
    lum_fs_est = 0;

}
#endif


void
Castro::get_shock_lums(Real time,int& ind_rs, int& ind_cd, int& ind_fs,
   Real& lum_lab_rs,Real& lum_com_rs,Real& lum_lab_cd, Real& lum_com_cd,
   Real& lum_lab_fs, Real& lum_com_fs, Real& lum_bol)
{

  for (int lev = 0; lev <= parent->finestLevel(); lev++) {
    ca_set_amr_info(lev,-1,-1,-1.0,-1.0);
    Castro& c_lev = getLevel(lev);
  // std::vector<std::string>::iterator it_lab = std::find(radiation->plotvar_names.begin(),radiation->plotvar_names.end(),"Frlabx");
  // std::vector<std::string>::iterator it_com = std::find(radiation->plotvar_names.begin(),radiation->plotvar_names.end(),"Frcomx");
  //
  // int idx_lab = std::distance(radiation->plotvar_names.begin(),it_lab);
  // int idx_com = std::distance(radiation->plotvar_names.begin(),it_com);
  //
  // // std::cout << idx_lab << "\n";
  // // std::cout << idx_com << "\n";
  int idx_lab = radiation->icomp_lab_Fr;
  int idx_com = radiation->icomp_com_Fr;

  auto& flx_lab_mfab = radiation->plotvar[lev];

//  auto& flx_com_mfab = radiation->plotvar[lev];

  for(MFIter mfi(*flx_lab_mfab,true);mfi.isValid();++mfi)
  {
    const Box& box  = mfi.validbox();
    const auto lo   = amrex::lbound(box);
    const auto hi   = amrex::ubound(box);
    auto flx_fab = (*flx_lab_mfab).array(mfi);
    if((ind_rs >= lo.x) && (ind_rs <= hi.x )) {
      lum_lab_rs = 4.*3.14159265359*pow(ind_rs*geom.CellSize(0),2)*(flx_fab(ind_rs,0,0,idx_lab));
      lum_com_rs = 4.*3.14159265359*pow(ind_rs*geom.CellSize(0),2)*(flx_fab(ind_rs,0,0,idx_com));}
    if((ind_cd >= lo.x) && (ind_cd <= hi.x )) {
      lum_lab_cd = 4.*3.14159265359*pow(ind_cd*geom.CellSize(0),2)*(flx_fab(ind_cd,0,0,idx_lab));
      lum_com_cd = 4.*3.14159265359*pow(ind_cd*geom.CellSize(0),2)*(flx_fab(ind_cd,0,0,idx_com));}
    if((ind_fs >= lo.x) && (ind_fs <= hi.x )) {
      lum_lab_fs = 4.*3.14159265359*pow(ind_fs*geom.CellSize(0),2)*(flx_fab(ind_fs,0,0,idx_lab));
      lum_com_fs = 4.*3.14159265359*pow(ind_fs*geom.CellSize(0),2)*(flx_fab(ind_fs,0,0,idx_com));}
    if(hi.x == geom.Domain().bigEnd(0)) lum_bol = 4.*3.14159265359*pow(geom.ProbHi(0),2)*(flx_fab(geom.Domain().bigEnd(0),0,0,idx_lab));
  }

  // for(MFIter mfi(*flx_com_mfab,true);mfi.isValid();++mfi)
  // {
  //   const Box& box  = mfi.validbox();
  //   const auto lo   = amrex::lbound(box);
  //   const auto hi   = amrex::ubound(box);
  //   auto flx_fab = (*flx_com_mfab).array(mfi);
  //   if((ind_fs >= lo.x) && (ind_fs <= hi.x )) lum_com_fs = 4.*3.14159265359*pow(ind_fs*geom.CellSize(0),2)*abs(flx_fab(ind_fs,0,0));
  // }
  //std::cout << "plotvar ind = " << idx << "\n";
}
}

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
Castro::fs_lum(Real time, Real& lum,int& i_fs)
{
  const Real* dx = geom.CellSize();

  for(int lev = 0; lev <= parent -> finestLevel(); lev++)
  {
    Castro& c_lev = getLevel(lev);
    MultiFab& Es = c_lev.get_new_data(State_Type);

    for(MFIter mfi(Es);mfi.isValid();++mfi)
    {
      const Box& box = mfi.tilebox();
      const int* lo = box.loVect();
      const int * hi = box.hiVect();
      lum_fs_shock(BL_TO_FORTRAN_3D(Es[mfi]),
      ARLIM_3D(lo),ARLIM_3D(hi),
      ZFILL(dx),&time,&lum,&i_fs);
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
Castro::cd_shock_radius(Real time, Real& r_cd,int& i_cd)
{
  const Real* dx = geom.CellSize();


  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    Castro& c_lev = getLevel(lev);
	auto mfcdmask = c_lev.derive("cd_mask",time,0);
  //MultiFab& Es = c_lev.get_new_data(State_Type);
  // r_cd = 0.0;
  // i_cd = 0;


    for(MFIter mfi(*mfcdmask); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
        FArrayBox& fabcdmask = (*mfcdmask)[mfi];
        cdshock(BL_TO_FORTRAN_3D(fabcdmask),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &r_cd,&i_cd);

            }
  }
}


void
Castro::csm_edge_radius(Real time, Real& r_csm,int& i_csm)
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
                  &r_csm,&i_csm);

            }
  }
}


void
Castro::reverse_shock_radius(Real time, Real& r_rs,int& i_rs)
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
                  &r_rs,&i_rs);

            }
  }
}


void
Castro::forward_shock_radius(Real time, Real& r_fs,int& i_fs)
{
  const Real* dx = geom.CellSize();


  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    Castro& c_lev = getLevel(lev);
//	  auto mfcsmmask = c_lev.derive("csm_mask",time,0);
    MultiFab& Es = c_lev.get_new_data(State_Type);


    for(MFIter mfi(Es); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
        fs_radius(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &r_fs,&i_fs);

            }
  }
}

void
Castro::rad_temperatures(Real time, int& ind_rs, int& ind_cd, int& ind_fs, int& ind_phot, Real& temp_rs, Real& temp_cd, Real& temp_fs, Real& temp_phot)
{

  const Real* dx = geom.CellSize();

  for(int lev = 0; lev <= parent->finestLevel(); lev++){

    MultiFab& Er = getLevel(lev).get_new_data(Rad_Type);

    for(MFIter mfi(Er,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        rad_temp_rs(BL_TO_FORTRAN_3D(Er[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_rs,&temp_rs);
        rad_temp_cd(BL_TO_FORTRAN_3D(Er[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_cd,&temp_cd);
        rad_temp_fs(BL_TO_FORTRAN_3D(Er[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_fs,&temp_fs);
        rad_temp_phot(BL_TO_FORTRAN_3D(Er[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_phot,&temp_phot);

            }

  }
}

void
Castro::gas_temperatures(Real time, int& ind_rs, int& ind_cd, int& ind_fs, int& ind_phot, Real& temp_rs, Real& temp_cd, Real& temp_fs, Real& temp_phot)
{

  const Real* dx = geom.CellSize();

  for(int lev = 0; lev <= parent->finestLevel(); lev++){

    MultiFab& Es = getLevel(lev).get_new_data(State_Type);

    for(MFIter mfi(Es,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        gas_temp_rs(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_rs,&temp_rs);
        gas_temp_cd(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_cd,&temp_cd);
        gas_temp_fs(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_fs,&temp_fs);
        gas_temp_phot(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_phot,&temp_phot);
            }

  }
}



void
Castro::shock_velocities(Real time, int& ind_rs, int& ind_cd, int& ind_fs, Real& veloc_rs, Real& veloc_cd, Real& veloc_fs)
{

  const Real* dx = geom.CellSize();

  for(int lev = 0; lev <= parent->finestLevel(); lev++){

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

    Castro& c_lev = getLevel(lev);
//	  auto mfcsmmask = c_lev.derive("csm_mask",time,0);
    MultiFab& Es = c_lev.get_new_data(State_Type);

    for(MFIter mfi(Es,true); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        vel_rs(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_rs,&veloc_rs);
        vel_cd(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_cd,&veloc_cd);
        vel_fs(BL_TO_FORTRAN_3D(Es[mfi]),
                  ARLIM_3D(lo),ARLIM_3D(hi),
                  ZFILL(dx),&time,
                  &ind_fs,&veloc_fs);

            }

  }
}


void
Castro::get_photospheric_velocity(Real time,int& ind_tauphot,Real& vel_phot){


  vel_phot = 0.0;
  const Real* dx = geom.CellSize();
  for(int lev = 0; lev <= parent->finestLevel();lev++){
    ca_set_amr_info(lev,-1,-1,-1.0,-1.0);

    Castro& c_lev = getLevel(lev);
    MultiFab& Es = c_lev.get_new_data(State_Type);

    for(MFIter mfi(Es,true);mfi.isValid();++mfi){
      const Box& box  = mfi.tilebox();
      const int* lo   = box.loVect();
      const int* hi   = box.hiVect();

      phot_velocity(BL_TO_FORTRAN_3D(Es[mfi]),
                ARLIM_3D(lo),ARLIM_3D(hi),
                ZFILL(dx),&time,
                &ind_tauphot,&vel_phot);
    }
  }
}

void
Castro::integ_optical_depth(Real time,int& ind_tauphot,int& ind_taucol,Real& tau_fs, Real& tau_rs, Real& tau_cd,
                            int& i_fs,int& i_rs,int& i_cd){




  const Real* dx = geom.CellSize();

  // tau_fs = 0.0;
  // tau_rs = 0.0;
  // tau_cd = 0.0;


  // this loops through levels, but this so far is only implemented for single-level (ref=0)
  for (int lev = 0; lev <= parent->finestLevel(); lev++) {

    ca_set_amr_info(lev, -1, -1, -1.0, -1.0);


    Castro& c_lev = getLevel(lev);

    // need to work off of the boxarray
    BoxArray ba = c_lev.boxArray();

    DistributionMapping dm(ba);

    MultiFab tau_r(ba,dm,1,0);

    tau_r = 0.0;

    // how many grids has the domain been divided into?
    long num_grids = ba.size();


    // struct that holds the global grid ID (across all ranks), its lo index, and the sum across the grid
    // reverse_order is not implemented here, default assumes integrating from hi to lo
    struct GidLo{
      long gid;
      int lo;
      int reverse;
      Real exclusive_sum;


    GidLo(long id, int smallend, int reverse_order) : gid(id), lo(smallend), reverse(reverse_order), exclusive_sum(0.0) {};

    // this comparison is in reverse order (sorts from hi to lo)
    bool operator<(const GidLo& other) const
    {
      return lo > other.lo;
    }


  };


  Vector<GidLo> grids;

  // put the grids into the vec
  for(long i=0; i < num_grids; ++i){
    grids.emplace_back(i,ba[i].smallEnd(0),1);
  //  std::cout << "presort: i = " << i << "  j = " << grids[i].gid << "  smallend = " << ba[grids[i].gid].smallEnd(0) << "\n";
  }

  // sort grid id's based on hi to lo op
  std::sort(grids.begin(),grids.end());


  // need to make an inverse mapping of the ith position of the jth grid, where j is sorted by hi-lo
  Vector<long>grid_inv(num_grids);
  for(long i = 0; i < num_grids; ++i) {
  //  grid_inv[i] = grids[i].gid;
    grid_inv[grids[i].gid] = i;
  //  std::cout << "i = " << i << "  j = " << grids[i].gid << "  smallend = " << ba[grids[i].gid].smallEnd(0) << "\n";
  }

  Vector<Real> prefix_sums(num_grids,0.0);



  // get multifab for tau in each zone (a derived quantity)
	  auto mfdtau = c_lev.derive("dtau",time,0);
  //  auto mftaur = c_lev.derive("tau_r",time,0);

    for(MFIter mfi(*mfdtau,true); mfi.isValid(); ++mfi)
    {
        long gidi = mfi.index();
        const Box& box  = mfi.validbox();
        const auto lo   = amrex::lbound(box);
        const auto hi   = amrex::ubound(box);
        auto dtau_fab = (*mfdtau).array(mfi);
        auto dtau_psum_fab = tau_r.array(mfi);

        // calculate the prefix sums for each grid
        Real fab_prefix_sum = 0.0;
        for (int kk = hi.z; kk >= lo.z; --kk) {
        for (int jj = hi.y; jj >= lo.y; --jj) {
        for (int ii = hi.x; ii >= lo.x; --ii) {
          fab_prefix_sum += dtau_fab(ii,jj,kk);
          dtau_psum_fab(ii,jj,kk) = fab_prefix_sum;

        }
      }
    }

    //  std::cout << "grid " << gidi << "  prefix sum = " << fab_prefix_sum << "\n";


    // need to store prefix sums for the AlltoAll operation
    prefix_sums[gidi] = fab_prefix_sum;

            }

            // communicate partial sums across ranks
            ParallelDescriptor::ReduceRealSum(prefix_sums.dataPtr(),num_grids);


      //      std::cout << "grid " << 0 << "  sum = " << grids[0].exclusive_sum << "\n";
        for(long i = 1; i < num_grids; ++i){
          grids[i].exclusive_sum = grids[i-1].exclusive_sum + prefix_sums[grids[i-1].gid];
        //  std::cout << "grid i = " << i << "  j =  " << grids[i-1].gid << "  psum = " << prefix_sums[grids[i-1].gid] << "\n";
        }


        for(MFIter mfi(tau_r);mfi.isValid();++mfi)
        {
          long gidi = mfi.index();
          auto& dtau_psum_fab = tau_r[mfi];
          dtau_psum_fab.plus(grids[grid_inv[gidi]].exclusive_sum,0,1);
        }

        Real radius = -1.0;
        Real tempval = 0.0;

      //  auto mftemp = c_lev.derive("Temp",time,0);



      ind_tauphot = 1e40;
      ind_taucol = 1e40;
      // tau_fs = 0;
      // tau_rs = 0;
      // tau_cd = 0;

      for(MFIter mfi(tau_r); mfi.isValid(); ++mfi)
      {
          long gidi = mfi.index();
          const Box& box  = mfi.validbox();
          const auto lo   = amrex::lbound(box);
          const auto hi   = amrex::ubound(box);
          auto tau_arr = tau_r.array(mfi);

          Real tmp_ind = -1;
          Real tmp_tau = -1;

          for (int kk = hi.z; kk >= lo.z; --kk) {
          for (int jj = hi.y; jj >= lo.y; --jj) {
          for (int ii = hi.x; ii >= lo.x; --ii) {

            if(ii==i_fs) tau_fs = tau_arr(ii,jj,kk);
            if(ii==i_rs) tau_rs = tau_arr(ii,jj,kk);
            if(ii==i_cd) tau_cd = tau_arr(ii,jj,kk);
        //   std::cout << tau_arr(ii,jj,kk) << "\n";
              if((tau_arr(ii,jj,kk)<31.623) && ii<ind_taucol) ind_taucol = ii;
              if((tau_arr(ii,jj,kk)<2.0/3.0) && ii<ind_tauphot) ind_tauphot = ii;
                //tmp_tau = tau_arr(ii,jj,kk);}


             }
           }
        }

      //  std::cout << "tau= " << tmp_tau << "ind = " << tmp_ind << "\n";

      }








  }

}
