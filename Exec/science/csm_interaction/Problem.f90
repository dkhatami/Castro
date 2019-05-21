! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len) bind(C, name="problem_checkpoint")

  ! called by the IO processor during checkpoint

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo



end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len) bind(C, name="problem_restart")

  ! called by ALL processors during restart

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_restart



subroutine bolometric_lum(rad_state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          lum) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use fundamental_constants_module, only: c_light
    use rad_params_module, only : ngroups
    use meth_params_module, only: NVAR, QRAD
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: rad_state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:ngroups-1)
    real(rt)        :: lum

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             r = position(i,j,k)
             xx = abs(probhi(1)-r(1))

             if (xx .lt. 1.1e0_rt*dx) then

               lum = four*M_PI*r(1)*r(1)*c_light*rad_state(i,j,k,0)

             else
               lum = -1.e0_rt

             endif

          enddo
       enddo
    enddo




end subroutine bolometric_lum


subroutine lum_fs(cd_mask,cd_lo,cd_hi,state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          lum) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module
    use fundamental_constants_module, only: c_light
    use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFS
    use prob_params_module, only : probhi
    use probdata_module, only : kap
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision,  intent(in) :: cd_mask(cd_lo(1):cd_hi(1),cd_lo(2):cd_hi(2),cd_lo(3):cd_hi(3))
    real(rt)        :: lum

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx
    real(rt)         :: lfac


    lfac = 9.e0_rt*M_PI/8.e0_rt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             r = position(i,j,k)


             if (cd_mask(i,j,k) > ZERO) then
              lum = lfac*r(1)*r(1)*state(i+1,j,k,UMX)**3/state(i+1,j,k,URHO)**2
              lum = abs(lum)
            endif


          enddo
       enddo
    enddo




end subroutine lum_fs


subroutine lum_rs(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          lum) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use fundamental_constants_module, only: c_light
    use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFS
    use prob_params_module, only : probhi
    use probdata_module, only : kap
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: lum

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx
    real(rt)         :: lfac


    lfac = 9.e0_rt*M_PI/8.e0_rt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             r = position(i,j,k)

             if (state(i,j,k,USHK) > ZERO .and. state(i,j,k,UFS)/state(i,j,k,URHO) > HALF) then
              lum = lfac*r(1)*r(1)*(state(i-1,j,k,UMX)/state(i-1,j,k,URHO)-state(i+1,j,k,UMX)/state(i+1,j,k,URHO))**3*state(i,j,k,URHO)
              lum = abs(lum)
            endif


          enddo
       enddo
    enddo




end subroutine lum_rs


subroutine cdshock(cd_mask, r_lo, r_hi, lo, hi, dx, time, r_cd) bind(C,name='cdshock')

  use amrex_constants_module, only: ZERO
  use castro_util_module, only : position

  implicit none


  integer         ,  intent(in) :: r_lo(3), r_hi(3)

  double precision,  intent(in) :: cd_mask(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))

  integer         ,  intent(in) :: lo(3), hi(3)

  double precision,  intent(in) :: dx, time

  double precision,  intent(inout) :: r_cd

  integer :: i,j,k

  double precision :: r(3)

  do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        r = position(i,j,k)

        if (cd_mask(i,j,k) > ZERO) then
          r_cd = r(1)
        endif

      enddo
    enddo
  enddo



end subroutine cdshock


subroutine csm_edge(csm_mask, r_lo, r_hi, lo, hi, dx, time, r_csm) bind(C,name='csm_edge')

  use amrex_constants_module, only: HALF
  use castro_util_module, only : position

  implicit none


  integer         ,  intent(in) :: r_lo(3), r_hi(3)

  double precision,  intent(in) :: csm_mask(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))

  integer         ,  intent(in) :: lo(3), hi(3)

  double precision,  intent(in) :: dx, time

  double precision,  intent(inout) :: r_csm

  integer :: i,j,k

  double precision :: r(3)

  do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        r = position(i,j,k)

        if (csm_mask(i,j,k) > HALF) then
          r_csm = r(1)
        endif

      enddo
    enddo
  enddo



end subroutine csm_edge
