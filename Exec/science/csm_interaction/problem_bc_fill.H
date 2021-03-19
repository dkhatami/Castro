#ifndef problem_bc_fill_H
#define problem_bc_fill_H

///
/// Problem-specific boundary condition fill
///
/// @param i         x-index
/// @param j         y-index
/// @param k         z-index
/// @param state     simulation state (Fab)
/// @param time      simulation time
/// @param bcs       boundary conditions
/// @param geomdata  geometry data
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     Array1D<BCRec, 0, NUM_STATE-1> bcs,
                     GeometryData geomdata) {

	const int* domlo = geomdata.Domain().loVect();
	const int* domhi = geomdata.Domain().hiVect();
	const Real* problo = geomdata.ProbLo();
	const Real* dx = geomdata.CellSize();

	if((bcs(URHO).lo(0) == EXT_DIR) && i < domlo[0]){

		Real t_0 = problem::r_0/problem::beta/2.998e10_rt;
		Real n = problem::n;
		Real d = problem::d;
		Real rho0_ej;
		Real R_t = problem::f_r*problem::r_0;
		Real T_ej = 0.5e0_rt*problem::m_0*1.989e33_rt*pow(problem::beta*2.998e10,2)*pow(problem::r_0,-3)/7.5646e-15_rt;
		T_ej = pow(T_ej,0.25);
		T_ej *= problem::r_star/problem::r_0;
		if(problem::pl_ej > 0) rho0_ej = (n-3.0_rt)*(3.0_rt-d)/(4.0_rt*3.14159_rt*(n-d))*problem::m_0*1.989e33_rt/pow(R_t,3);
		else rho0_ej = 3.0_rt/(4.0_rt*3.14159_rt)*problem::m_0*1.989e33_rt/pow(problem::r_0,3);

		state(i,j,k,URHO) = 5.e-1_rt*rho0_ej*pow(problo[0]/R_t,-d)*pow(1+time/t_0,-3);
		state(i,j,k,UMX) =  state(i,j,k,URHO)*problo[0]/(t_0+time);
		state(i,j,k,UMY) = 0.0_rt;
		state(i,j,k,UMZ) = 0.0_rt;
		state(i,j,k,UTEMP) = T_ej/(1.0_rt+time/t_0);
		state(i,j,k,UEINT) = state(i,j,k,URHO)*1.38e-16_rt*state(i,j,k,UTEMP)/1.67e-24_rt;
		state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5e0_rt*state(i,j,k,UMX)*state(i,j,k,UMX)/state(i,j,k,URHO);
		
		state(i,j,k,UFS) = problem::x_H*state(i,j,k,URHO);
		state(i,j,k,UFS+1) = problem::x_He*state(i,j,k,URHO);
		state(i,j,k,UFX) = 0.5_rt*state(i,j,k,URHO);
		state(i,j,k,UFX+1) = 2.0_rt*state(i,j,k,URHO);
		state(i,j,k,UFX+2) = 1.0_rt*state(i,j,k,URHO);
		state(i,j,k,UFX+3) = 0.0_rt;
		state(i,j,k,UFX+4) = 0.0_rt;

	}



}

#endif
