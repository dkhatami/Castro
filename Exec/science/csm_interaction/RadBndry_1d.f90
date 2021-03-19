! THIS IS ONLY A TEMPLATE!!!

! If radiation.[lo|hi]_bcflag is set to 1, this subroutine will be used to set
! boundary conditions for radiation diffusion.
!
! For LO_DIRICHLET, this should set Dirichlet value of radiation energy density.
! For LO_NEUMANN, this should set inward radiation flux.
! For LO_MARSHAK & LO_SANCHEZ_POMRANING, this should set incident radiation flux.

subroutine rbndry(  &
     bf, b_l1, b_h1, &
     &   d_l1, d_h1, &
     dx, xlo, t, dir, face) 

  use rad_params_module, only : ngroups
  use probdata_module, only: M_ni, r_0, beta

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: b_l1, b_h1
  integer, intent(in) :: d_l1, d_h1  ! computational domain index
  integer, intent(in) :: dir         ! always 0 (i.e., x-direction) in 1D
  integer, intent(in) :: face        ! 0: low     1: high
  real(rt)        , intent(out) :: bf(b_l1:b_h1,0:ngroups-1)
  real(rt)        , intent(in) :: dx(1), xlo(1), t

  real(rt)  ::  L_ni, eps_ni, eps_co, t_ni, t_co,t_0,tt,day

  day = 24.0*3600.0

  eps_ni = 3.9e10_rt !nickel decay specific energy
  eps_co = 6.8e9_rt ! cobalt decay specific energy
  t_ni = 8.8*day 
  t_co = 111.3*day
  t_0 = r_0/(beta*2.998e10_rt) ! time at which simulation starts
  ! we assume anything before was lost to adiabatic expansion

  tt = t_0+t

  L_ni = M_ni*1.989e33_rt*((eps_ni-eps_co)*exp(-tt/t_ni)+eps_co*exp(-tt/t_co))

  if (face .eq. 0) then

      bf = L_ni/(4.0_rt*3.14159_rt*xlo(1)*xlo(1)) !convert to inner boundary flux

    

  endif

end subroutine rbndry
