module rad_util_module

  use castro_error_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  use rad_params_module, only : ngroups

  implicit none

contains

  subroutine compute_ptot_ctot(lam, q, cg, ptot, ctot, gamc_tot)

    use meth_params_module, only : QPRES, QRHO, comoving, QRAD, QPTOT, NQ
    use fluxlimiter_module, only : Edd_factor

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , intent(in) :: lam(0:ngroups-1)
    real(rt)        , intent(in) :: q(NQ)
    real(rt)        , intent(in) :: cg
    real(rt)        , intent(out) :: ptot
    real(rt)        , intent(out) :: ctot
    real(rt)        , intent(out) :: gamc_tot

    integer :: g

    real(rt)         :: csrad2, Eddf, gamr, prad

    csrad2 = ZERO
    prad = ZERO

    do g = 0, ngroups-1
       if (comoving) then
          Eddf = Edd_factor(lam(g))
          gamr = (THREE - Eddf)/TWO
       else
          gamr = lam(g) + ONE
       end if

       prad = prad + lam(g)*q(QRAD+g)
       csrad2 = csrad2 + gamr * (lam(g)*q(QRAD+g)) / q(QRHO)
    end do

    ptot = q(QPRES) + prad

    ctot = cg**2 + csrad2
    gamc_tot = ctot * q(QRHO) / ptot

    ctot = sqrt(ctot)

  end subroutine compute_ptot_ctot

  function FLDlambda(r, limiter) result (lambda)

    use amrex_fort_module, only : rt => amrex_real
    real(rt)         :: r
    integer :: limiter

    real(rt)         :: lambda

    if (limiter .eq. 0) then
       ! no limiter
       lambda = 1.e0_rt/3.e0_rt

    else if (limiter < 10) then
       ! approximate LP
       lambda = (2.e0_rt + r) / (6.e0_rt + r * (3.e0_rt + r))

    else if (limiter < 20) then
       ! Bruenn
       lambda = 1.e0_rt / (3.e0_rt + r)

    else if (limiter < 30) then
       ! Larsen's square root
       lambda = 1.e0_rt / sqrt(9.e0_rt + r**2)

    else if (limiter < 40) then 
       ! Minerbo
       if (r .lt. 1.5e0_rt) then
          lambda = 2.e0_rt/(3.e0_rt + sqrt(9.e0_rt+12.e0_rt*r**2))
       else 
          lambda = 1.e0_rt/(1.e0_rt+r+sqrt(1.e0_rt+2.e0_rt*r))
       end if

    else
       print *, "limiter = ", limiter
       call castro_error("Unknown limiter type")
    endif
  end function FLDlambda

end module rad_util_module
