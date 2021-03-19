#ifndef problem_initialize_rad_data_H
#define problem_initialize_rad_data_H

#include <fundamental_constants.H>
#include <blackbody.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_rad_data (int i, int j, int k,
								  Array4<Real> const& rad,
								  GpuArray<Real, NGROUPS+1> const& xnu,
								  GpuArray<Real, NGROUPS> const& nugroup,
								  GpuArray<Real, NGROUPS> const& dnugroup,
								  const GeometryData& geomdata)
{
	constexpr int nsub = 16;

	const Real* dx = geomdata.CellSize();
	const Real* problo = geomdata.ProbLo();

	Real R_ej = problem::r_0;
	Real R_csm = (1.0_rt+problem::delt+problem::xi_bo)*problem::r_0;
	Real T_ej = 0.5_rt*problem::m_0*1.989e33_rt*pow(problem::beta*2.998e10_rt,2)*pow(problem::r_0,-3)/C::a_rad;
	T_ej *= problem::r_star/problem::r_0;
	T_ej = pow(T_ej,0.25);

	Real dx_sub = dx[0] / static_cast<Real>(nsub);

	Real xcl = problo[0] + dx[0]*static_cast<Real>(i);

	Real vtot = 0.0_rt;
	Real T = 0.0_rt;

	for(int ii = 0; ii < nsub; ii++){
		Real xx = xcl + (static_cast<Real>(ii)+0.5_rt)*dx_sub;
		Real vsub = xx * xx;
		vtot += vsub;

		if (xx < R_ej) T += vsub*T_ej;
		else if (xx < R_csm) T += vsub*problem::T_csm;
		else T += vsub*problem::T_0;

	}

	T /= vtot;

	if (NGROUPS == 1) {
		rad(i,j,k,0) = C::a_rad*std::pow(T,4);
	}
	else{
		for(int g = 0; g < NGROUPS; g++){
			rad(i,j,k,g) = BGroup(T,xnu[g],xnu[g+1]);
		}
	}



}
#endif
