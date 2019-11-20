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



subroutine rad_temp_rs(rad_state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use fundamental_constants_module, only: c_light, a_rad
    use rad_params_module, only : ngroups
    use meth_params_module, only: NVAR, QRAD
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: rad_state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:ngroups-1)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = (rad_state(i,j,k,0)/a_rad)**0.25
              endif



          enddo
       enddo
    enddo




end subroutine rad_temp_rs



subroutine rad_temp_cd(rad_state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use fundamental_constants_module, only: c_light, a_rad
    use rad_params_module, only : ngroups
    use meth_params_module, only: NVAR, QRAD
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: rad_state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:ngroups-1)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = (rad_state(i,j,k,0)/a_rad)**0.25
              endif



          enddo
       enddo
    enddo




end subroutine rad_temp_cd



subroutine rad_temp_fs(rad_state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use fundamental_constants_module, only: c_light, a_rad
    use rad_params_module, only : ngroups
    use meth_params_module, only: NVAR, QRAD
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: rad_state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:ngroups-1)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = (rad_state(i,j,k,0)/a_rad)**0.25
              endif



          enddo
       enddo
    enddo




end subroutine rad_temp_fs


subroutine rad_temp_phot(rad_state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use fundamental_constants_module, only: c_light, a_rad
    use rad_params_module, only : ngroups
    use meth_params_module, only: NVAR, QRAD
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: rad_state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:ngroups-1)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = (rad_state(i,j,k,0)/a_rad)**0.25
              endif



          enddo
       enddo
    enddo




end subroutine rad_temp_phot



subroutine gas_temp_rs(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UTEMP
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = state(i,j,k,UTEMP)
              endif



          enddo
       enddo
    enddo




end subroutine gas_temp_rs

subroutine gas_temp_cd(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UTEMP
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = state(i,j,k,UTEMP)
              endif



          enddo
       enddo
    enddo




end subroutine gas_temp_cd

subroutine gas_temp_fs(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UTEMP
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = state(i,j,k,UTEMP)
              endif



          enddo
       enddo
    enddo




end subroutine gas_temp_fs

subroutine gas_temp_phot(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, temp) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UTEMP
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: temp

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                temp = state(i,j,k,UTEMP)
              endif



          enddo
       enddo
    enddo




end subroutine gas_temp_phot


subroutine vel_rs(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, vel) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UMX
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: vel

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                vel = state(i,j,k,UMX)/state(i,j,k,URHO)
              endif



          enddo
       enddo
    enddo




end subroutine vel_rs



subroutine vel_cd(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, vel) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UMX
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: vel

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                vel = state(i,j,k,UMX)/state(i,j,k,URHO)
              endif



          enddo
       enddo
    enddo




end subroutine vel_cd


subroutine vel_fs(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, vel) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UMX
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: vel

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                vel = state(i,j,k,UMX)/state(i,j,k,URHO)
              endif



          enddo
       enddo
    enddo




end subroutine vel_fs


subroutine phot_velocity(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          ind, vel) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module, only: four, M_PI
    use meth_params_module, only: NVAR, QRAD, URHO, UMX
    use prob_params_module, only : probhi
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer, intent(in) :: ind
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: vel

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if(i .eq. ind) then
                vel = state(i,j,k,UMX)/state(i,j,k,URHO)
              endif



          enddo
       enddo
    enddo




end subroutine phot_velocity


subroutine lum_fs_shock(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          lum,i_fs) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module
    use fundamental_constants_module, only: c_light
    use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFX
    use prob_params_module, only : probhi
    use probdata_module, only : kap
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3), cd_lo(3), cd_hi(3)
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        :: lum

    integer :: i_fs

    ! Local variables

    integer :: i, j, k
    real(rt)         :: r(3), xx, dx
    real(rt)         :: lfac


    lfac = 9.e0_rt*M_PI/8.e0_rt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1)+2, hi(1)-2

             r = position(i,j,k)


             if (i .eq. i_fs) then
              lum = lfac*r(1)*r(1)*state(i,j,k,UMX)**3/state(i,j,k,URHO)**2
            endif


          enddo
       enddo
    enddo




end subroutine lum_fs_shock


subroutine lum_rs_shock(state,s_lo,s_hi, &
                          lo, hi, dx, time, &
                          lum) bind(C)
    use castro_util_module, only : position
    use amrex_constants_module
    use fundamental_constants_module, only: c_light
    use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFX
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
          do i = lo(1)+2, hi(1)-2

             r = position(i,j,k)

             if (state(i,j,k,USHK) > ZERO .and. state(i,j,k,UFX)/state(i,j,k,URHO) > HALF) then
              lum = lfac*r(1)*r(1)*(state(i-2,j,k,UMX)/state(i-2,j,k,URHO)-state(i+2,j,k,UMX)/state(i+2,j,k,URHO))**3*state(i,j,k,URHO)
              lum = abs(lum)
            endif


          enddo
       enddo
    enddo




end subroutine lum_rs_shock


subroutine cdshock(cd_mask, s_lo, s_hi, lo, hi, dx, time, r_cd,i_cd) bind(C,name='cdshock')

  use amrex_constants_module, only: ZERO,HALF
  use castro_util_module, only : position
  !use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFX
  use amrex_fort_module, only : rt => amrex_real
  implicit none


  integer         ,  intent(in) :: s_lo(3), s_hi(3)

  double precision,  intent(in) :: cd_mask(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

  integer         ,  intent(in) :: lo(3), hi(3)

  double precision,  intent(in) :: dx, time

  double precision,  intent(inout) :: r_cd

  integer         ,  intent(inout) :: i_cd

  integer :: i,j,k

  double precision :: r(3)

  !i_cd = 0

  do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        r = position(i,j,k)

        if (cd_mask(i,j,k) > ZERO) then
          if(r(1) < r_cd)then
          r_cd = r(1)
          i_cd = i
        endif
        endif

      enddo
    enddo
  enddo



end subroutine cdshock


subroutine csm_edge(csm_mask, r_lo, r_hi, lo, hi, dx, time, r_csm,i_csm) bind(C,name='csm_edge')

  use amrex_constants_module, only: HALF
  use castro_util_module, only : position

  implicit none


  integer         ,  intent(in) :: r_lo(3), r_hi(3)

  double precision,  intent(in) :: csm_mask(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))

  integer         ,  intent(in) :: lo(3), hi(3)

  double precision,  intent(in) :: dx, time

  double precision,  intent(inout) :: r_csm

  integer         ,  intent(inout) :: i_csm

  integer :: i,j,k

  double precision :: r(3)


  do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        r = position(i,j,k)

        if (csm_mask(i,j,k) > HALF) then
          if(r(1) .gt. r_csm) then
          r_csm = r(1)
          i_csm = i
        endif
        endif

      enddo
    enddo
  enddo



end subroutine csm_edge



subroutine rs_radius(state, s_lo, s_hi, lo, hi, dx, time, r_rs,i_rs) bind(C,name='rs_radius')

  use amrex_constants_module, only: ZERO
  use castro_util_module, only : position
  use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFX
  use amrex_fort_module, only : rt => amrex_real
  implicit none


  integer         ,  intent(in) :: s_lo(3), s_hi(3)

  double precision,  intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

  integer         ,  intent(in) :: lo(3), hi(3)

  double precision,  intent(in) :: dx, time

  double precision,  intent(inout) :: r_rs

  integer         ,  intent(inout) :: i_rs

  integer :: i,j,k

  double precision :: r(3)



  do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        r = position(i,j,k)

        if (state(i,j,k,USHK) > ZERO .and. state(i,j,k,UFX)/state(i,j,k,URHO) > 0.9e0_rt) then
          if(r(1) < r_rs) then
          r_rs = r(1)
          i_rs = i
        endif
        endif

      enddo
    enddo
  enddo



end subroutine rs_radius


subroutine fs_radius(state, s_lo, s_hi, lo, hi, dx, time, r_fs,i_fs) bind(C,name='fs_radius')

  use amrex_constants_module, only: ZERO
  use castro_util_module, only : position
  use meth_params_module, only: NVAR, QRAD, URHO, UMX, USHK, UFX, QTEMP
  use probdata_module, only: rho_a
  use amrex_fort_module, only : rt => amrex_real
  implicit none


  integer         ,  intent(in) :: s_lo(3), s_hi(3)

  double precision,  intent(in) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

  integer         ,  intent(in) :: lo(3), hi(3)

  double precision,  intent(in) :: dx, time

  double precision,  intent(inout) :: r_fs

  integer         ,  intent(inout) :: i_fs

  integer :: i,j,k

  double precision :: r(3)


  do k = lo(3),hi(3)
    do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        r = position(i,j,k)

        if (state(i,j,k,USHK) > ZERO .and. state(i,j,k,UFX+1)/state(i,j,k,URHO) > 0.9e0_rt) then
        !  if(state(i,j,k,UMX)/state(i,j,k,URHO) .gt. 1.0e5_rt) then
          if(r(1) < r_fs) then
          r_fs = r(1)
          i_fs = i
      !  endif
        endif
        endif

      enddo
    enddo
  enddo



end subroutine fs_radius
