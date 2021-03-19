#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <ambient.H>
#include <math.h>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
									Array4<Real> const& state,
									const GeometryData& geomdata)

{
	constexpr int nsub = 1;
	constexpr Real Pi = 3.1415926535897932_rt;
	constexpr Real Msun = 1.989e33_rt;
	constexpr Real clight = 2.998e10_rt;
	constexpr Real yr_to_sec = 3.154e7_rt;
	constexpr Real km_to_cm = 1.e5_rt;
	const Real* dx = geomdata.CellSize();
	const Real* problo = geomdata.ProbLo();
	const Real* probhi = geomdata.ProbHi();
	Real n = problem::n;
	Real d = problem::d;
	Real s = problem::s;
	Real p = problem::p;


	Real xcl = problo[0]+dx[0] * static_cast<Real>(i);
	Real dx_sub = dx[0] / static_cast<Real>(nsub);

	Real R_bo = (1.+problem::delt+problem::xi_bo*problem::delt)*problem::r_0;
	Real R_t = problem::f_r*problem::r_0;
	Real M_ej = problem::m_0*Msun;
	Real M_csm = M_ej*problem::eta;
	Real R_csm = problem::r_0;
	Real R_out = (1.+problem::delt)*problem::r_0;
	Real R_w = problem::R_w;
	Real rho_0 = M_csm/std::pow(R_csm,3.0_rt);
	Real v_t = problem::beta*clight;
	Real kap = 0.4_rt;
	Real rho_a = min(problem::tau_a/(kap*(probhi[0]-R_w)),problem::f_a*rho_0);

	Real T_ej = 0.5e0_rt*M_ej*v_t*v_t*pow(problem::r_0,-3)/7.5646e-15_rt;
	T_ej = pow(T_ej,0.25);
	T_ej *= problem::r_star/problem::r_0;

	Real D_wind = problem::Mdot*Msun/yr_to_sec/(4.0_rt*Pi*problem::v_w*km_to_cm);

	problem::filter_rhomax = 30.0_rt*D_wind*pow(R_bo,-2.0_rt);

	Real vtot = 0.0_rt;
	Real rho = 0.0_rt;
	Real T = 0.0_rt;
	Real vel = 0.0_rt;

	for(int ii = 0; ii < nsub; ii++){

		Real xx = xcl + (static_cast<Real>(ii)+0.5_rt)*dx_sub;

		Real vsub = xx * xx;
		vtot += vsub;
		Real rho0_ej, rho0_csm, rho_csm, rho_ww, rho_bo, rho_cw, rho_wa, rho_ca, rho_cc;
		Real rho_sub, vel_sub, T_sub;
		if (problem::pl_ej > 0)
			rho0_ej = (n-3.0_rt)*(3.0_rt-d)/(4.0_rt*Pi*(n-d))*M_ej/std::pow(R_t,3);
		else rho0_ej = 3.0_rt/(4.0_rt*Pi)*M_ej/std::pow(R_csm,3);

		rho0_csm = (3.0_rt-s)/(4.0_rt*Pi)/(std::pow(1.0_rt+problem::delt,3-s)-1.0_rt)*M_csm/std::pow(R_csm,3);

		rho_csm = rho0_csm*std::pow(xx/R_csm,-1.0_rt*s);
		rho_ww = D_wind*std::pow(xx,-2.0_rt);
		rho_bo = 1.0_rt*rho0_csm*std::pow(R_out/R_csm,-1.0_rt*s)*std::pow(xx/R_out,-1.0_rt*p);
		rho_cc = rho_csm;// + 0.5_rt*(rho_bo-rho_csm)*(1.0_rt+tanh((xx-R_out)/(problem::h_csm*problem::delt*R_out)));
		rho_cw = rho_csm;

		if (R_w >= 0.0_rt){
			Real h_c = problem::h_csm;
			Real dR_t = problem::h_csm*R_csm;
			rho_cw = rho_bo + 0.5_rt*(rho_ww-rho_bo)*(1.0_rt+tanh((xx-R_bo)/dR_t));

			if(R_w > 0.0_rt)
				rho_wa = rho_ww + 0.5_rt*(rho_a-rho_ww)*(1.0_rt+tanh((xx-(R_out+R_w))/(h_c*(R_out+R_w))));
		}
		else
		{
			Real dR_t = problem::h_csm*R_csm;
			rho_ca = rho_bo + 0.5_rt*(rho_a-rho_bo)*(1.0_rt+tanh((xx-R_bo)/dR_t));
		}

		if (xx < R_csm)
		{
			T_sub = T_ej;
			vel_sub = v_t*(xx/R_t);

			rho_sub = rho0_ej;

			if (problem::pl_ej > 0)
			{
				if (xx < R_t) rho_sub *= std::pow(xx/R_t,-1.0_rt*d);
				else rho_sub *= std::pow(xx/R_t,-1.0_rt*n);
			}

		}
		else if (xx < R_bo)
		{
			vel_sub = problem::v_csm*km_to_cm;
			T_sub = problem::T_csm;


			if (xx < R_out) rho_sub = rho_cc;
			else rho_sub = rho_cw;
		}
		else
		{
			T_sub = problem::T_0;

			if ((xx < R_w) or (R_w==0.0_rt))
			{
				rho_sub = rho_cw;
				vel_sub = problem::v_w*km_to_cm;
			}
			else
			{
				vel_sub = 0.0_rt;

				if (R_w < 0.0_rt) rho_sub = rho_ca;
				else rho_sub = rho_wa;
			}
		}

		rho += rho_sub*vsub;
		T += T_sub*vsub;
		vel += vel_sub*vsub;




	}

	rho /= vtot;
	T /= vtot;
	vel /= vtot;


	state(i,j,k,URHO) = rho;
	state(i,j,k,UTEMP) = T;
	state(i,j,k,UMX) = state(i,j,k,URHO)*vel;
	state(i,j,k,UMY) = 0.0_rt;
	state(i,j,k,UMZ) = 0.0_rt;
	
	for (int m = 0; m < NumSpec; m++){
		state(i,j,k,UFS+m) =0.0_rt;
	}
	for(int m = 0; m < NumAux; m++){
		state(i,j,k,UFX+m) = 0.0_rt;
	}
	//state(i,j,k,UFS) = state(i,j,k,URHO);


	if (xcl < R_csm)
	{
		state(i,j,k,UFX+2) = 1.0_rt*rho;
		state(i,j,k,UFX+3) = 0.0_rt;
		state(i,j,k,UFX+4) = 0.0_rt;
	}
	else if (xcl < R_bo)
	{
		state(i,j,k,UFX+2) = 0.0_rt;
		state(i,j,k,UFX+3) = 1.0_rt*rho;
		state(i,j,k,UFX+4) = 0.0_rt;
	}
	else{
		state(i,j,k,UFX+2) = 0.0_rt;
		state(i,j,k,UFX+3) = 0.0_rt;
		state(i,j,k,UFX+4) = 1.0_rt*rho;
	}

	state(i,j,k,UFS) = problem::x_H*rho;
	state(i,j,k,UFS+1) = problem::x_He*rho;
	state(i,j,k,UFS+2) = problem::x_C*rho;
	state(i,j,k,UFS+3) = problem::x_O*rho;



	eos_t eos_state;

	eos_state.rho = rho;
	eos_state.T = T;

	for(int m = 0; m < NumSpec; m++) eos_state.xn[m] = state(i,j,k,UFS+m)/rho;
	for(int m = 2; m < NumAux; m++) eos_state.aux[m] = state(i,j,k,UFX+m)/rho;
	eos(eos_input_rt,eos_state);

	state(i,j,k,UFX) = eos_state.aux[0]*rho;
	state(i,j,k,UFX+1) = eos_state.aux[1]*rho;
	

	state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
	state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * state(i,j,k,UMX) * state(i,j,k,UMX) / state(i,j,k,URHO);



	




}
#endif
