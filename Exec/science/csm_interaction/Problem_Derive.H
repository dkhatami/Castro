
void ca_derfs_mask
	(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
	 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
	 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);


void ca_dercsm_mask
	(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
	 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
	 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derej_mask
	(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
	 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
	 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derdtau
	(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
	 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
	 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derdtauabs
	(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
	 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
	 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);



