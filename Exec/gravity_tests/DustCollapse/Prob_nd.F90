subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, ONE
  use castro_error_module, only: castro_error
  use probdata_module, only: rho_0, T_0, X_0, p_0, rho_ambient, T_ambient, &
                             r_old, r_old_s, r_0, smooth_delta, r_offset, offset_smooth_delta, &
                             center_x, center_y, center_z, nsub
  use eos_type_module, only: eos_t, eos_input_rp
  use eos_module, only: eos
  use network, only: nspec
  use meth_params_module, only: small_temp
  use prob_params_module, only: center

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i
  type (eos_t) :: eos_state

  namelist /fortin/ &
       rho_0, r_0, r_old, p_0, rho_ambient, smooth_delta, &
       center_x, center_y, center_z, r_offset, offset_smooth_delta, nsub

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 127
  character :: probin*(maxlen)

  if (namlen > maxlen) then
     call castro_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  rho_0 = 1.e9_rt
  r_0 = 6.5e8_rt
  r_old = r_0
  p_0 = 1.e10_rt
  rho_ambient = 1.e0_rt
  smooth_delta = 1.e-5_rt
  r_offset = ZERO
  offset_smooth_delta = 1.e0_rt
  nsub = 5

  ! Read namelists in probin file
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  r_old_s = r_old

  ! in 3-d we center the sphere at (center_x, center_y, center_z)

  ! in 2-d we are going to enforce that the lower left corner of the
  ! domain is 0.0 (i.e., we only model a quadrant)

  ! in 1-d, we enforce that the center is the origin (since we are
  ! spherical)

  center(:) = ZERO

#if AMREX_SPACEDIM == 1
  center(1) = ZERO
#else
  center(1) = center_x
  center(2) = center_y
#if AMREX_SPACEDIM == 3
  center(3) = center_z
#endif
#endif

  if (problo(1) /= ZERO) call castro_error("ERROR: xmin should be 0!")
  if (problo(2) /= ZERO) call castro_error("ERROR: ymin should be 0!")
  if (problo(3) /= ZERO) call castro_error("ERROR: zmin should be 0!")

  ! set the composition to be uniform
  allocate(X_0(nspec))

  X_0(:) = ZERO
  X_0(1) = ONE

  ! get the ambient temperature and sphere temperature, T_0

  eos_state % rho = rho_0
  eos_state % p   = p_0
  eos_state % xn  = x_0
  eos_state % T   = small_temp ! Initial guess for the EOS

  call eos(eos_input_rp, eos_state)

  T_0 = eos_state % T

  eos_state % rho = rho_ambient

  call eos(eos_input_rp, eos_state)

  T_ambient = eos_state % T

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, state_lo, state_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module, only: ZERO, HALF, ONE
  use probdata_module, only: rho_0, X_0, p_0, rho_ambient, r_0, smooth_delta, r_offset, offset_smooth_delta, nsub
  use eos_type_module, only : eos_t, eos_input_rp
  use eos_module, only: eos
  use network, only: nspec
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, small_temp
  use prob_params_module, only: problo, center, dg
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: state_lo(3), state_hi(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
  real(rt), intent(in   ) :: time, delta(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3)

  real(rt) :: xl, yl, zl, xx, yy, zz
  real(rt) :: dist
  real(rt) :: pres, eint, temp, avg_rho, rho_n
  real(rt) :: volinv
  real(rt) :: dx_sub, dy_sub, dz_sub
  integer  :: i, j, k, ii, jj, kk, n

  type (eos_t) :: eos_state

#if AMREX_SPACEDIM == 1
  volinv = ONE/dble(nsub)
#elif AMREX_SPACEDIM == 2
  volinv = ONE/dble(nsub*nsub)
#else
  volinv = ONE/dble(nsub*nsub*nsub)
#endif

  dx_sub = delta(1)/dble(nsub)
  dy_sub = delta(2)/dble(nsub)
  dz_sub = delta(3)/dble(nsub)

  do k = lo(3), hi(3)
     zl = problo(1) + dble(k) * delta(3)

     do j = lo(2), hi(2)
        yl = problo(2) + dble(j) * delta(2)

        do i = lo(1), hi(1)
           xl = problo(3) + dble(i) * delta(1)

           avg_rho = ZERO

           do kk = 0, dg(3)*(nsub-1)
              zz = zl + (dble(kk) + HALF) * dz_sub

              do jj = 0, dg(2)*(nsub-1)
                 yy = yl + (dble(jj) + HALF) * dy_sub

                 do ii = 0, nsub-1
                    xx = xl + (dble(ii) + HALF) * dx_sub

                    dist = sqrt((xx-center(1))**2 + (yy-center(2))**2 + (zz-center(3))**2)

                    ! use a tanh profile to smooth the transition between rho_0
                    ! and rho_ambient
                    rho_n = rho_0 - HALF * (rho_0 - rho_ambient) * (ONE + tanh((dist - r_0)/smooth_delta))

                    ! allow for the center to be empty
                    if (r_offset > ZERO) then
                       rho_n = rho_n - HALF * (rho_n - rho_ambient) * (ONE + tanh((r_offset - dist) / offset_smooth_delta))
                    end if

                    avg_rho = avg_rho + rho_n

                 enddo
              enddo
           enddo

           state(i,j,k,URHO) = avg_rho * volinv

           eos_state % rho = state(i,j,k,URHO)
           eos_state % p   = p_0
           eos_state % T   = small_temp ! Initial guess for the EOS
           eos_state % xn  = X_0

           call eos(eos_input_rp, eos_state)

           temp = eos_state % T
           eint = eos_state % e

           state(i,j,k,UTEMP) = temp
           state(i,j,k,UMX) = ZERO
           state(i,j,k,UMY) = ZERO
           state(i,j,k,UMZ) = ZERO
           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eint
           state(i,j,k,UEINT) = state(i,j,k,URHO) * eint
           state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO) * X_0(1:nspec)

        enddo
     enddo
  enddo

end subroutine ca_initdata
