subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use probdata_module
  use actual_eos_module, only : gamma_const
  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer untin, i

  namelist /fortin/ const_c_v, c_v_exp_n

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen=127
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call castro_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  const_c_v = 8.0744926545150113e-24_rt
  c_v_exp_n = -3.0_rt

  xmin = problo(1)
  xmax = probhi(1)

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

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

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only : nspec

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), nscal)

  integer :: i, j, k
  real(rt) :: c_v, eint

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           state(i,j,k,URHO) = 1.0_rt
           state(i,j,k,UTEMP) = 0.0_rt
           state(i,j,k,UMX:UMZ) = 0.0_rt

           ! set the composition to be all in the first species
           state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
           state(i,j,k,UFS) = state(i,j,k,URHO)

           state(i,j,k,UEINT) = 0.0_rt
           state(i,j,k,UEDEN) = 0.0_rt

        end do
     end do
  end do

end subroutine ca_initdata


! :::
! ::: -----------------------------------------------------------
! :::

subroutine ca_initrad(level, time, lo, hi, nrad, &
                      rad_state, rad_state_lo, rad_state_hi, &
                      delta, xlo, xhi)

  use probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nrad
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: rad_state_lo(3), rad_state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: rad_state(rad_state_lo(1):rad_state_hi(1), &
                                       rad_state_lo(2):rad_state_hi(2), &
                                       rad_state_lo(3):rad_state_hi(3), 0:nrad-1)

  ! local variables
  integer :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rad_state(i,j,k,:) = 0.0
        end do
     end do
  end do

end subroutine ca_initrad
