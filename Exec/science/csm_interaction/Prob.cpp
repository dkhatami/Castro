
#include <Castro.H>
#include <Castro_F.H>
#include <Radiation.H>

#include <fstream>

using namespace amrex;


Real Castro::lum_edge = 0.0;
Real Castro::lum_phot = 0.0;
Real Castro::r_phot = 0.0;
Real Castro::r_col = 0.0;
Real Castro::T_eff = 0.0;
Real Castro::T_col = 0.0;
Real Castro::v_chrom = 0.0;
Real Castro::rho_chrom = 0.0;
Real Castro::rho_phot = 0.0;
Real Castro::v_phot = 0.0;
Real Castro::lum_phot_com = 0.0;
Real Castro::lum_phot_lab = 0.0;
Real Castro::lum_edge_com = 0.0;
Real Castro::lum_edge_lab = 0.0;
Real Castro::lum_surf_com = 0.0;
Real Castro::lum_surf_lab = 0.0;

int Castro::ind_tau_phot = 0;
int Castro::ind_tau_surf = 0;
int Castro::ind_tau_col = 0;
int Castro::ind_rs = 0;
int Castro::ind_sh = 0;
int Castro::ind_fs = 0;
Real Castro::radius_rs = 0.0;
Real Castro::radius_sh = 0.0;
Real Castro::radius_fs = 0.0;
Real Castro::vel_rs = 0.0;
Real Castro::vel_sh = 0.0;
Real Castro::vel_fs = 0.0;
Real Castro::tau_fs = 0.0;
Real Castro::tau_sh = 0.0;
Real Castro::tau_rs = 0.0;
Real Castro::Ek_csm = 0.0;
Real Castro::Eg_csm = 0.0;
Real Castro::Er_csm = 0.0;
Real Castro::Ek_ej = 0.0;
Real Castro::Eg_ej = 0.0;
Real Castro::Er_ej = 0.0;
void
Castro::problem_post_timestep()
{
	BL_PROFILE("Castro::problem_post_timestep()");

	using namespace problem;

	if(level != 0) return;

	int finest_level = parent->finestLevel();
	Real time = state[State_Type].curTime();
	Real dt = parent->dtLevel(0);


	shock_properties(time,dt,ind_rs,ind_sh,ind_fs,radius_rs,radius_sh,radius_fs,vel_rs,vel_sh,vel_fs);

	sum_energy(time,dt,ind_rs,ind_fs,Ek_csm,Eg_csm,Er_csm,Ek_ej,Eg_ej,Er_ej);


	integ_optical_depth(time,ind_tau_phot,ind_tau_surf,ind_tau_col,ind_rs,ind_sh,ind_fs,tau_rs,tau_sh,tau_fs);



	ParallelDescriptor::ReduceIntMin(ind_tau_phot);
	ParallelDescriptor::ReduceIntMin(ind_tau_surf);
	ParallelDescriptor::ReduceIntMin(ind_tau_col);
	ParallelDescriptor::ReduceRealMax(tau_rs);
	ParallelDescriptor::ReduceRealMax(tau_sh);
	ParallelDescriptor::ReduceRealMax(tau_fs);

	bolometric_lum(time,dt,ind_tau_phot,ind_tau_surf,ind_tau_col,lum_phot,lum_edge,
					lum_phot_lab, lum_phot_com, lum_edge_lab, lum_edge_com,
					lum_surf_lab, lum_surf_com,T_col,v_phot,rho_phot,v_chrom,rho_chrom);

	r_phot = ind_tau_phot*geom.CellSize(0);
	r_col = ind_tau_col*geom.CellSize(0);
	T_eff = pow(lum_phot/(4.*3.14159265359*r_phot*r_phot*5.67e-5),0.25);

	if(amrex::ParallelDescriptor::IOProcessor()){
		int dataprecision = 8;
		int datwidth = 25;
		int fixwidth = 25;
		int intwidth = 12;

		std::ostream& log_lc = parent->DataLog(1);
		std::ostream& log_egy = parent->DataLog(2);
		std::ostream& log_shock = parent->DataLog(3);
		std::ostream& log_lum = parent->DataLog(4);
		std::ostream& log_eps = parent->DataLog(5);

		if(time == dt)
		{
			log_lc << "Time (s)" << std::setw(fixwidth) << "L_bol" << std::setw(fixwidth) << "R_phot" << std::setw(fixwidth) << "T_eff" << std::setw(fixwidth) << "T_col" << std::setw(fixwidth);
			log_lc << std::endl;

			log_egy << "Time (s)" << std::setw(fixwidth) << "Ej Kinetic Energy" << std::setw(fixwidth) << "Csm Kinetic Energy" << std::setw(fixwidth)
					<< "Ej Internal Energy" << std::setw(fixwidth) << "Csm Internal Energy" << std::setw(fixwidth)
				    << "Ej Radiation Energy" << std::setw(fixwidth) << "Csm Radiation Energy" << std::setw(fixwidth);

			log_egy << std::endl;


			log_shock << "Time (s)" << std::setw(fixwidth)
					<< "RS  Radius (cm) " <<  std::setw(fixwidth) << "RS Velocity (cm/s) " << std::setw(fixwidth) << "RS Tau " << std::setw(fixwidth) 
					<< "FS  Radius (cm) " <<  std::setw(fixwidth) << "FS Velocity (cm/s) " << std::setw(fixwidth) << "FS Tau " << std::setw(fixwidth)
					<< "Shock Radius (cm)" << std::setw(fixwidth) << "Shock Velocity (cm/s) " << std::setw(fixwidth) << "Shock Tau " << std::setw(fixwidth);

			log_shock << std::endl;

			log_lum << "Time (s)" << std::setw(fixwidth) << "L_phot_lab" << std::setw(fixwidth) << "L_phot_com" << std::setw(fixwidth) << "L_edge_lab" << std::setw(fixwidth)
			        << "L_edge_com" << std::setw(fixwidth) << "L_surf_lab" << std::setw(fixwidth) << "L_surf_com" << std::setw(fixwidth);
			log_lum << std::endl;

			log_eps << "Time (s)" << std::setw(fixwidth) << "T_eff" << std::setw(fixwidth) << "T_col" << std::setw(fixwidth) << "R_phot" << 
			   std::setw(fixwidth) << "R_col" << std::setw(fixwidth) <<	"v_phot" << std::setw(fixwidth) << "rho_phot" <<
				std::setw(fixwidth) << "v_col" << std::setw(fixwidth) << "rho_col" << std::setw(fixwidth) << std::endl;
		}

		log_lc << std::fixed;

		log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << time;
		log_lc << std::scientific;
		log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << lum_phot;
		log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << r_phot;
		log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << T_eff;
		log_lc << std::setw(datwidth) << std::setprecision(dataprecision) << T_col;

		log_lc << std::endl;

		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << time;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << radius_rs;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << vel_rs;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << tau_rs;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << radius_fs;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << vel_fs;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << tau_fs;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << radius_sh;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << vel_sh;
		log_shock << std::setw(datwidth) << std::setprecision(dataprecision) << tau_sh;

		log_shock << std::endl;

		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << time;
		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << Ek_ej;
		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << Ek_csm;
		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << Eg_ej;
		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << Eg_csm;
		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << Er_ej;
		log_egy << std::setw(datwidth) << std::setprecision(dataprecision) << Er_csm;

		log_egy << std::endl;

		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << time;
		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << lum_phot_lab;
		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << lum_phot_com;
		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << lum_edge_lab;
		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << lum_edge_com;
		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << lum_surf_lab;
		log_lum << std::setw(datwidth) << std::setprecision(dataprecision) << lum_surf_com;

		log_lum << std::endl;

		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << time;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << T_eff;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << T_col;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << r_phot;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << r_col;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << v_phot;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << rho_phot;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << v_chrom;
		log_eps << std::setw(datwidth) << std::setprecision(dataprecision) << rho_chrom;

		log_eps << std::endl;


		std::cout << "Bolometric Luminosity    =   " << lum_phot << "     erg/s     " << "\n";
		std::cout << "Edge Luminosity          =   " << lum_edge << "     erg/s     " << "\n";
		std::cout << "Photospheric Radius      =   " << r_phot << "     cm       " << "\n";
		std::cout << "Effective Temperature    =   " << T_eff << "        K       " << "\n";
		std::cout << "Shock Radius             =   " << radius_fs << "        cm      " << "\n";
		std::cout << "Shock Velocity           =   " << vel_fs/1.0e8_rt    << "      (1e3 km/s)    " << "\n";
		std::cout << "Shock Tau                =   " << tau_fs << "\n";

	}

	lum_edge = 0.0_rt;
	lum_phot = 0.0_rt;
	lum_phot_lab = 0.0_rt;
	lum_phot_com = 0.0_rt;
	lum_edge_lab = 0.0_rt;
	lum_edge_com = 0.0_rt;
	lum_surf_lab = 0.0_rt;
	lum_surf_com = 0.0_rt;
	r_phot = 0.0_rt;
	T_eff = 0.0_rt;
	T_col = 0.0_rt;
	rho_chrom = 0.0_rt;
	v_chrom = 0.0_rt;
	rho_phot = 0.0_rt;
	v_phot = 0.0_rt;
	ind_tau_phot = 0;
	ind_tau_surf = 0;
	ind_tau_col = 0;
	ind_fs = -1;
	ind_sh = -1;
	ind_rs = -1;
	radius_rs = 0.0;
	radius_sh = 0.0;
	radius_fs = 0.0;
	vel_rs = 0.0;
	vel_sh = 0.0;
	vel_fs = 0.0;
	tau_rs = 0.0;
	tau_sh = 0.0;
	tau_fs = 0.0;
	Ek_csm = 0.0;
	Eg_csm = 0.0;
	Er_csm = 0.0;
	Ek_ej = 0.0;
	Eg_ej = 0.0;
	Er_ej = 0.0;


}

void
Castro::sum_energy(Real time, Real dt, int ind_rs, int ind_fs,
				   Real& Ek_csm, Real& Eg_csm, Real& Er_csm,
				   Real& Ek_ej,  Real& Eg_ej,  Real& Er_ej)
{
	for(int lev = 0; lev <= parent->finestLevel();lev++){

		Castro& c_lev = getLevel(lev);

		GeometryData geomdata = c_lev.geom.data();

		const auto dx = c_lev.geom.CellSizeArray();
		const auto problo = c_lev.geom.ProbLoArray();
		const auto probhi = c_lev.geom.ProbHiArray();


		auto mf_rhoE = c_lev.derive("eint_E",time,0);
		auto mf_rho  = c_lev.derive("density",time,0);
		auto mf_rhoK = c_lev.derive("kineng",time,0);
		auto mf_Er = c_lev.derive("rad",time,0);


		auto mf_csm = c_lev.derive("csm_mask",time,0);
		auto mf_ej = c_lev.derive("ej_mask",time,0);

		Ek_csm = 0.0;
		Eg_csm = 0.0;
		Er_csm = 0.0;
		Ek_ej = 0.0;
		Eg_ej = 0.0;
		Er_ej = 0.0;



		for(MFIter mfi(*mf_csm,true);mfi.isValid();++mfi){

			const Box& box = mfi.validbox();
			const auto lo = amrex::lbound(box);
			const auto hi = amrex::ubound(box);
			auto vol = c_lev.volume[mfi].array();


			auto csm_mask = (*mf_csm)[mfi].array();
			auto ej_mask = (*mf_ej)[mfi].array();


			for(int kk = hi.z; kk >= lo.z; --kk){
				for(int jj = hi.y; jj >= lo.y; --jj){
					for(int ii = hi.x; ii >= lo.x; --ii){

						Real ek_tmp = (*mf_rhoK)[mfi].array()(ii,jj,kk);
						Real eg_tmp = (*mf_rhoE)[mfi].array()(ii,jj,kk)*(*mf_rho)[mfi].array()(ii,jj,kk);
						Real er_tmp = (*mf_Er)[mfi].array()(ii,jj,kk);





						Ek_csm += ek_tmp*csm_mask(ii,jj,kk)*vol(ii,jj,kk);
						Eg_csm += eg_tmp*csm_mask(ii,jj,kk)*vol(ii,jj,kk);
						Er_csm += er_tmp*csm_mask(ii,jj,kk)*vol(ii,jj,kk);

						Ek_ej += ek_tmp*ej_mask(ii,jj,kk)*vol(ii,jj,kk);
						Eg_ej += eg_tmp*ej_mask(ii,jj,kk)*vol(ii,jj,kk);
						Er_ej += er_tmp*ej_mask(ii,jj,kk)*vol(ii,jj,kk);


					}
				}
			}



		}

		ParallelDescriptor::ReduceRealSum(Ek_csm);
		ParallelDescriptor::ReduceRealSum(Eg_csm);
		ParallelDescriptor::ReduceRealSum(Er_csm);
		ParallelDescriptor::ReduceRealSum(Ek_ej);
		ParallelDescriptor::ReduceRealSum(Eg_ej);
		ParallelDescriptor::ReduceRealSum(Er_ej);










	}



}


void
Castro::bolometric_lum(Real time, Real dt,int ind_phot,int ind_surf,int ind_col, Real& lum_phot, Real& lum_edge,
						Real& lum_phot_lab, Real& lum_phot_com,
						Real& lum_edge_lab, Real& lum_edge_com,
						Real& lum_surf_lab, Real& lum_surf_com,
						Real& T_col,Real& v_phot,Real& rho_phot,
						Real& v_col, Real& rho_col)
{
	BL_PROFILE("Castro::bolometric_lum()");

	using namespace problem;

	BL_ASSERT(level == 0 || (!parent->subCycle() && level == parent->finestLevel()));

	ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,
			  ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,ReduceOpSum,
			  ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
	ReduceData<Real,Real,Real,Real,Real,Real,Real,Real,Real,Real,Real,Real,Real> reduce_data(reduce_op);

	using ReduceTuple = typename decltype(reduce_data)::Type;

	for (int lev = 0; lev <= parent->finestLevel(); lev++){
		
		Castro& c_lev = getLevel(lev);

		GeometryData geomdata = c_lev.geom.data();

		const auto dx = c_lev.geom.CellSizeArray();
		const auto problo = c_lev.geom.ProbLoArray();
		const auto probhi = c_lev.geom.ProbHiArray();

		int idx_lab = radiation->icomp_lab_Fr;
		int idx_com = radiation->icomp_com_Fr;

		auto& flx_lab_mfab = radiation->plotvar[lev];
		//auto rad_mfab = c_lev.derive("rad",time,0);
		auto temp_mfab = c_lev.derive("Temp",time,0);
		auto rho_mfab = c_lev.derive("density",time,0);
		auto vel_mfab = c_lev.derive("magvel",time,0);


#ifdef _OPENMP
#pragma omp parallel
#endif
		for (MFIter mfi(*flx_lab_mfab,TilingIfNotGPU()); mfi.isValid(); ++mfi){

			const Box& box = mfi.tilebox();
			const auto hi = amrex::ubound(box);
			auto flx_fab = (*flx_lab_mfab)[mfi].array();
			auto vel_fab = (*vel_mfab)[mfi].array();
			auto rho_fab = (*rho_mfab)[mfi].array();
			//auto rad_fab = (*rad_mfab)[mfi].array();
			auto temp_fab = (*temp_mfab)[mfi].array();

			reduce_op.eval(box, reduce_data,
			[=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
			{

				Real lbol = 0.0_rt;
				Real lphot = 0.0_rt;
				Real Tc = 0.0_rt;
				Real vp = 0.0_rt;
				Real rhop = 0.0_rt;
				Real vc = 0.0_rt;
				Real rhoc = 0.0_rt;
				Real lphot_c = 0.0_rt;
				Real lphot_l = 0.0_rt;
				Real ledge_c = 0.0_rt;
				Real ledge_l = 0.0_rt;
				Real lsurf_c = 0.0_rt;
				Real lsurf_l = 0.0_rt;
				Real lc_tmp = 0.0_rt;
				if(i == geom.Domain().bigEnd(0))
				{
				ledge_l = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_lab));
				ledge_c = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_com));

				}
				if(i==ind_phot) lphot = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_com));
				if(i==ind_phot){
				lphot_l = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_lab));
				lphot_c = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_com));
				vp = vel_fab(i,0,0);
				rhop = rho_fab(i,0,0);
				}
				if(i==ind_surf){
				lsurf_l = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_lab));
				lsurf_c = 4.0_rt*3.1415926359*pow(i*geom.CellSize(0),2)*std::abs(flx_fab(i,0,0,idx_com));
				}
				if(i==ind_col){
					//Real er0 = rad_fab(i,0,0);
					//Tc = pow(er0/7.5646e-15_rt,0.25);
					Tc = temp_fab(i,0,0);
					vc = vel_fab(i,0,0);
					rhoc = rho_fab(i,0,0);
				}



					






				return {lbol,lphot,lphot_l,lphot_c,ledge_l,ledge_c,lsurf_l,lsurf_c,Tc,vp,rhop,vc,rhoc};
			});

		}

	}

	ReduceTuple hv = reduce_data.value();

	lum_edge = amrex::get<0>(hv);
	lum_phot = amrex::get<1>(hv);
	lum_phot_lab = amrex::get<2>(hv);
	lum_phot_com = amrex::get<3>(hv);
	lum_edge_lab = amrex::get<4>(hv);
	lum_edge_com = amrex::get<5>(hv);
	lum_surf_lab = amrex::get<6>(hv);
	lum_surf_com = amrex::get<7>(hv);
	T_col = amrex::get<8>(hv);
	v_phot = amrex::get<9>(hv);
	rho_phot = amrex::get<10>(hv);
	v_col = amrex::get<11>(hv);
	rho_col = amrex::get<12>(hv);

	const int nfoo_sum = 13;
	Real foo_sum[nfoo_sum] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	foo_sum[0] = lum_edge;
	foo_sum[1] = lum_phot;
	foo_sum[2] = lum_phot_lab;
	foo_sum[3] = lum_phot_com;
	foo_sum[4] = lum_edge_lab;
	foo_sum[5] = lum_edge_com;
	foo_sum[6] = lum_surf_lab;
	foo_sum[7] = lum_surf_com;
	foo_sum[8] = T_col;
	foo_sum[9] = v_phot;
	foo_sum[10] = rho_phot;
	foo_sum[11] = v_col;
	foo_sum[12] = rho_col;

	amrex::ParallelDescriptor::ReduceRealSum(foo_sum,nfoo_sum);

	lum_edge = foo_sum[0];
	lum_phot = foo_sum[1];
	lum_phot_lab = foo_sum[2];
	lum_phot_com = foo_sum[3];
	lum_edge_lab = foo_sum[4];
	lum_edge_com = foo_sum[5];
	lum_surf_lab = foo_sum[6];
	lum_surf_com = foo_sum[7];
	T_col = foo_sum[8];
	v_phot = foo_sum[9];
	rho_phot = foo_sum[10];
	v_col = foo_sum[11];
	rho_col = foo_sum[12];

}

void
Castro::shock_properties(Real time, Real dt, int& ind_rs, int& ind_sh, int& ind_fs,
						Real& radius_rs, Real& radius_sh, Real& radius_fs,
						Real& vel_rs, Real& vel_sh, Real& vel_fs)
{
	for(int lev = 0; lev <= parent->finestLevel(); lev++)
	{
		ca_set_amr_info(lev,-1,-1,-1.0,-1.0);
		Castro& c_lev = getLevel(lev);

		GeometryData geomdata = c_lev.geom.data();
		const auto dx = c_lev.geom.CellSizeArray();
		const auto problo = c_lev.geom.ProbLoArray();
		const auto probhi = c_lev.geom.ProbHiArray();

		auto mfej_mask = c_lev.derive("ej_mask",time,0);
		auto mffs_mask = c_lev.derive("fs_mask",time,0);
		auto mfcsm_mask = c_lev.derive("csm_mask",time,0);

		ind_rs = -1;
		ind_sh = 1e40;
		ind_fs = 1e40;

		for(MFIter mfi(*mfcsm_mask,true);mfi.isValid();++mfi){

			const Box& box = mfi.validbox();
			const auto lo = amrex::lbound(box);
			const auto hi = amrex::ubound(box);
			auto csm_mask_fab = (*mfcsm_mask)[mfi].array();
			auto ej_mask_fab = (*mfej_mask)[mfi].array();
			auto fs_mask_fab = (*mffs_mask)[mfi].array();


			for(int kk = hi.z; kk >= lo.z; --kk){
				for(int jj = hi.y; jj >= lo.y; --jj){
					for(int ii = hi.x; ii >= lo.x; --ii){


						if(csm_mask_fab(ii,jj,kk) > 0 && ii < ind_fs) ind_fs = ii;
						if(ej_mask_fab(ii,jj,kk) > 0 && ii > ind_rs) ind_rs = ii;
						if(fs_mask_fab(ii,jj,kk) > 0 && ii < ind_sh) ind_sh = ii;


					}
				}
			}
		}

		ParallelDescriptor::ReduceIntMax(ind_rs);
		ParallelDescriptor::ReduceIntMin(ind_sh);
		ParallelDescriptor::ReduceIntMin(ind_fs);


		auto mfvel = c_lev.derive("magvel",time,0);


		for(MFIter mfi(*mfvel,true);mfi.isValid();++mfi){

			radius_fs = 0.0_rt;
			vel_fs = 0.0_rt;
			radius_rs = 0.0_rt;
			vel_rs = 0.0_rt;
			radius_sh = 0.0_rt;
			vel_sh = 0.0_rt;


			const Box& box = mfi.validbox();
			const auto lo = amrex::lbound(box);
			const auto hi = amrex::ubound(box);
			auto  vel_fab = (*mfvel)[mfi].array();


			for(int kk = hi.z; kk >= lo.z; --kk){
				for(int jj = hi.y; jj >= lo.y; --jj){
					for(int ii = hi.x; ii >= lo.x; --ii){

						if(ii==ind_rs){
							vel_rs = vel_fab(ii,jj,kk);
							radius_rs = ind_rs*dx[0];

						}
						else if(ii==ind_sh){
							vel_sh = vel_fab(ii,jj,kk);
							radius_sh = ind_sh*dx[0];
						}
						else if(ii==ind_fs){
							vel_fs = vel_fab(ii,jj,kk);
							radius_fs = ind_fs*dx[0];

						}


					}
				}
			}


		}

		ParallelDescriptor::ReduceRealMax(vel_rs);
		ParallelDescriptor::ReduceRealMax(radius_rs);

		ParallelDescriptor::ReduceRealMax(vel_sh);
		ParallelDescriptor::ReduceRealMax(radius_sh);
		ParallelDescriptor::ReduceRealMax(vel_fs);
		ParallelDescriptor::ReduceRealMax(radius_fs);





	}






}

void
Castro::integ_optical_depth(Real time, int& ind_tau_phot,int& ind_tau_surf,int& ind_tau_col, int ind_rs, int ind_sh, int ind_fs, Real& tau_rs, Real& tau_sh, Real& tau_fs)
{
	for(int lev = 0; lev <= parent->finestLevel(); lev++)
	{
		ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

		Castro& c_lev = getLevel(lev);

		BoxArray ba = c_lev.boxArray();

		DistributionMapping dm(ba);

		MultiFab tau_r(ba,dm,2,0);
		
		tau_r = 0.0;

		long num_grids = ba.size();

		struct GidLo{
			long gid;
			int lo;
			int reverse;
			Real exclusive_sum;
			Real exclusive_sum_a;

		GidLo(long id, int smallend, int reverse_order) : gid(id), lo(smallend), reverse(reverse_order), exclusive_sum(0.0), exclusive_sum_a(0.0) {};

		bool operator<(const GidLo& other) const
		{
			return lo > other.lo;
		}
	};

	Vector<GidLo> grids;

	for(long i=0; i < num_grids; ++i)
	{
		grids.emplace_back(i,ba[i].smallEnd(0),1);
	}

		std::sort(grids.begin(),grids.end());

		Vector<long> grid_inv(num_grids);
		for(long i = 0; i < num_grids; ++i)
		{
			grid_inv[grids[i].gid] = i;
		}

		Vector<Real> prefix_sums(num_grids,0.0);
		Vector<Real> prefix_sums_a(num_grids,0.0);

		auto mfdtau = c_lev.derive("dtau",time,0);
		auto mfdtau_a = c_lev.derive("dtau_abs",time,0);

		for(MFIter mfi(*mfdtau,true); mfi.isValid(); ++mfi)
		{
			long gidi = mfi.index();
			const Box& box = mfi.validbox();
			const auto lo = amrex::lbound(box);
			const auto hi = amrex::ubound(box);
			auto dtau_fab = (*mfdtau)[mfi].array();
			auto dtau_psum_fab = tau_r[mfi].array();

			auto dtau_a_fab = (*mfdtau_a)[mfi].array();

			Real fab_prefix_sum = 0.0;
			Real fab_a_prefix_sum = 0.0;
			for(int kk = hi.z; kk >= lo.z; --kk){
				for(int jj = hi.y; jj >= lo.y; --jj){
					for(int ii = hi.x; ii >= lo.x; --ii){
						fab_prefix_sum += dtau_fab(ii,jj,kk);
						fab_a_prefix_sum += dtau_a_fab(ii,jj,kk);
						dtau_psum_fab(ii,jj,kk,0) = fab_prefix_sum;
						dtau_psum_fab(ii,jj,kk,1) = fab_a_prefix_sum;
					}
				}
			}

			prefix_sums[gidi] = fab_prefix_sum;
			prefix_sums_a[gidi] = fab_a_prefix_sum;

		}

		ParallelDescriptor::ReduceRealSum(prefix_sums.dataPtr(),num_grids);
		ParallelDescriptor::ReduceRealSum(prefix_sums_a.dataPtr(),num_grids);

		for(long i = 1; i < num_grids; ++i)
		{
			grids[i].exclusive_sum = grids[i-1].exclusive_sum + prefix_sums[grids[i-1].gid];
			grids[i].exclusive_sum_a = grids[i-1].exclusive_sum_a + prefix_sums_a[grids[i-1].gid];
		}

		for(MFIter mfi(tau_r); mfi.isValid();++mfi)
		{
			long gidi = mfi.index();
			auto& dtau_psum_fab = tau_r[mfi];
			dtau_psum_fab.plus(grids[grid_inv[gidi]].exclusive_sum,0,1);
			dtau_psum_fab.plus(grids[grid_inv[gidi]].exclusive_sum_a,1,1);
		}

		Real radius = -1.0;

		ind_tau_phot = 1e40;
		ind_tau_surf = 1e40;
		ind_tau_col = 1e40;
		tau_fs = 0;
		tau_sh = 0;
		tau_rs = 0;






		for(MFIter mfi(tau_r); mfi.isValid();++mfi)
		{
			long gidi = mfi.index();
			const Box& box = mfi.validbox();
			const auto lo = amrex::lbound(box);
			const auto hi = amrex::ubound(box);
			auto tau_arr = tau_r[mfi].array();

			Real tmp_ind = -1;
			Real tmp_tau = -1;

			for (int kk = hi.z; kk >= lo.z; --kk) {
			for (int jj = hi.y; jj >= lo.y; --jj) {
			for (int ii = hi.x; ii >= lo.x; --ii) {

				if((tau_arr(ii,jj,kk,0)<2.0/3.0) && ii<ind_tau_phot) ind_tau_phot = ii;
				if((tau_arr(ii,jj,kk,0)<0.01_rt) && ii<ind_tau_surf) ind_tau_surf = ii;
				if(sqrt(3.0*tau_arr(ii,jj,kk,0)*tau_arr(ii,jj,kk,1))<2.0/3.0  && ii<ind_tau_col) ind_tau_col = ii;
				if(ii==ind_rs) tau_rs = tau_arr(ii,jj,kk);
				if(ii==ind_sh) tau_sh = tau_arr(ii,jj,kk);
				if(ii==ind_fs) tau_fs = tau_arr(ii,jj,kk);

			}
			}
			}


		}

	}






}	


void
Castro::problem_post_restart()
{

}

void
Castro::problem_post_init()
{

}
