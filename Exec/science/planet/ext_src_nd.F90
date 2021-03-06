subroutine ca_ext_src(lo, hi, &
                      old_state, os_lo, os_hi, &
                      new_state, ns_lo, ns_hi, &
                      src, src_lo, src_hi, &
                      problo, dx, time, dt) bind(C, name='ca_ext_src')
  ! Compute the external sources for all the conservative equations.
  !
  ! This is called twice in the evolution:
  !
  ! first, for the predictor, it is called with (old, old) states
  !
  ! This is also used in the first pass of the conservative update
  ! (adding dt * S there)
  !
  ! Next we correct the source terms in the conservative update to
  ! time-center them.  Here we call ext_src(old, new), and then
  ! in time_center_source_terms we subtract off 1/2 of the first S
  ! and add 1/2 of the new S.
  !
  ! Therefore, to get a properly time-centered source, generally
  ! speaking, you always want to use the "new" state here.  That
  ! will be the time n state in the first call and the n+1 in the
  ! second call.


  use amrex_constants_module, only : ZERO, TWO, HALF
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, NVAR, UTEMP, UFS, UFX, UEINT
  use prob_params_module, only: center, probhi
  use actual_network, only: nspec, naux
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_re
  use eos_module, only: eos

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: os_lo(3), os_hi(3)
  integer,  intent(in   ) :: ns_lo(3), ns_hi(3)
  integer,  intent(in   ) :: src_lo(3), src_hi(3)
  real(rt), intent(in   ) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
  real(rt), intent(in   ) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
  real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
  real(rt), intent(in   ) :: problo(3), dx(3)
  real(rt), intent(in   ), value :: time, dt

  real(rt) :: source(3)
  integer  :: i, j, k
  real(rt) :: ang_vel, y, beta, p_orb_vel, p_radius,reset_center

  p_orb_vel = 1.d-3
  p_radius = 1.d10
  ang_vel = ZERO
  beta = p_orb_vel * TWO / p_radius
  reset_center = (problo(2) + probhi(2)) * HALF

  source(:) = ZERO

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + HALF) - reset_center
        do i = lo(1), hi(1)
           ang_vel = 0.5 * beta * y

#if AMREX_SPACEDIM == 3
           source(1) = + TWO * new_state(i,j,k,UMY) * ang_vel
           source(2) = - TWO * new_state(i,j,k,UMX) * ang_vel
           source(3) = ZERO
#endif

           src(i,j,k,UMX) = src(i,j,k,UMX) + source(1)
           src(i,j,k,UMY) = src(i,j,k,UMY) + source(2)
           src(i,j,k,UMZ) = src(i,j,k,UMZ) + source(3)

        end do
     end do
  end do

end subroutine ca_ext_src
