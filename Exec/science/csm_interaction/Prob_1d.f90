subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)
    
      use probdata_module
      use prob_params_module, only : center
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer :: init, namlen
      integer :: name(namlen)
      real(rt)      :: problo(1), probhi(1)

      integer :: untin,i

      namelist /fortin/ M_ej, E_ej, d, n, v_ej, T_ej, &
                        R_c, s, T_w, Mdot_w, v_w, &
                        rho_a, T_a, &
                        R_sh, M_sh, N_sh, f_sh, v_sh,T_sh, &
                        R_x0, R_x1, M_x, &
                        filter_rhomax, filter_timemax

      integer, parameter :: maxlen = 256
      character :: probin*(maxlen)


      do i = 1, namlen
        probin(i:i) = char(name(i))
      end do

      !namelist defaults
       
      M_ej = 5.e0_rt ! ejecta mass in solar masses
      E_ej = 1.e0_rt ! ejecta energy in 1e51 ergs
      d = 1.e0_rt  ! inner ejecta profile
      n = 10.e0_rt  ! outer ejecta profile
      T_ej = 1.e3_rt ! ejecta temperature
      v_ej = 1.e9_rt ! ejecta velocity in cm/s
      Mdot_w = 1.e-6_rt ! wind mass-loss in Msun/yr
      v_w = 1.e7_rt ! wind velocity in cm/s
      T_w = 2.e3_rt ! wind temperature
      rho_a = 1.e-18_rt ! ambient density
      T_a = 1.e3_rt ! ambient temperature




      R_c = 2.e14_rt ! ejecta->CSM radius
      s = 2.e0_rt ! CSM density profile

      R_sh = 1.e15_rt !shell radius
      M_sh = 1.e0_rt !shell mass in Msun
      N_sh = 3.e0_rt !num sigma (just needs to be >1)
      f_sh = 1.e-2_rt !shell thickness
      v_sh = 1.e7_rt !shell velocity in cm/s
      T_sh = 2.e3_rt ! shell temperature

      R_x0 = 1.e15_rt
      R_x1 = 2.e15_rt
      M_x = 0.e0_rt


      filter_rhomax = 1.e-99_rt
      filter_timemax = 0.e0_rt


      open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      center(1) = 0.e0_rt
      

end subroutine amrex_probinit


subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_h1, &
                       delta,xlo,xhi)

    use probdata_module
    use eos_module
    use network, only : nspec,naux
    use eos_type_module
    use amrex_constants_module, only: M_PI, FOUR3RD
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMZ, UEDEN, UEINT, UFS, UFX
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: level, nscal
    integer, parameter :: nsub = 1
    integer :: lo(1), hi(1)
    integer :: state_l1, state_h1
    real(rt)        :: xlo(1), xhi(1), time, delta(1)
    real(rt)        :: state(state_l1:state_h1,NVAR)

    real(rt)        :: xx, xxs, vol, dx_sub, vol_sub
    real(rt)        :: R_s, eta_r, eta_v, v_t, M_ch, rho_w
    real(rt)        :: rho0, rho1, rhow, rhoej, rho_sh, rho_x
    real(rt)        :: rho, vel, T, rho_sub, vel_sub, T_sub
    real(rt)		:: h_c, h_ej
    integer :: i, ii

    type(eos_t)  :: eos_state
    

    dx_sub = delta(1)/dble(nsub)

    eta_r = (n-3.e0_rt)*(3.e0_rt-d)/(n-d)/(4.e0_rt*M_PI)
    eta_v = sqrt(2.e0_rt*(5.e0_rt-d)*(n-5.e0_rt) &
                /((3.e0_rt-d)*(n-3.e0_rt)))

    M_ch = M_ej / 1.4e0_rt
    v_t = 6.e8_rt*eta_v*sqrt(E_ej/M_ch)
    R_s = v_t*R_c/v_ej


    h_ej = 0.1e0_rt*R_s
    h_c = 0.02e0_rt*R_x1
    rho_w = Mdot_w*2.e33_rt/(3.154e7_rt)/(4.*M_PI*R_c**2*v_w)


    rho_sh = M_sh*2.e33_rt/(4.e0_rt*M_PI*(f_sh*R_sh)**3)
    rho_sh = rho_sh/(sqrt(M_PI)/2.e0_rt*(1.e0_rt+2.e0_rt/f_sh**2)-N_sh*exp(-N_sh*N_sh))

    rho_x = M_x*2.e33_rt/(4.e0_rt*M_PI*R_x0*R_x0*(R_x1-R_x0))

    do i = lo(1), hi(1)

     xx = xlo(1)+delta(1)*(dble(i-lo(1)+0.5e0_rt))

     state(i,UMX:UMZ) = 0.0e0_rt

     vol  = 0.e0_rt
     rho     = 0.e0_rt
     T       = 0.e0_rt
     vel       = 0.e0_rt

     do ii=0,nsub-1

      xxs = xx + (dble(ii)+0.5e0_rt)*dx_sub

      vol_sub = xxs*xxs
      vol = vol + vol_sub


      rho0 = eta_r*M_ej*2.e33_rt/(R_s)**3*(R_s/xxs)**d
      rho1 = eta_r*M_ej*2.e33_rt/(R_s)**3*(R_s/xxs)**n
      rhow = rho_w*(R_c/xxs)**s

      rhoej = rho0 + 0.5e0_rt*(rho1-rho0)*(1.e0_rt+tanh((xxs-R_s)/h_ej))

      if (xxs .lt. R_c) then
          rho_sub = rhoej
          vel_sub = v_ej*(xxs/R_c)
          T_sub = T_ej
      else if (xxs .gt. R_x0) then
          rho_sub = rho_x + 0.5e0_rt*(rho_a-rho_x)*(1.e0_rt+tanh((xxs-R_x1)/h_c))
          T_sub = T_sh
          vel_sub = v_sh
      else
          rho_sub = rho_a
          T_sub = T_a
          vel_sub = 1.e5_rt
      endif
     ! rho_sub = rhoej + 0.5e0_rt*(rhow-rhoej)*(1.e0_rt+tanh((xxs-R_c)/h_c))

     ! rho_sub = rho_sub !+ rho_sh*exp(-(xxs/R_sh-1)**2/(f_sh*f_sh))




      rho = rho + rho_sub*vol_sub
      T = T + T_sub*vol_sub
      vel = vel + vel_sub*vol_sub

     enddo

     rho = rho / vol
     T = T / vol
     vel = vel / vol
     rho = MAX(rho,rho_a)
     rho = MIN(rho,1.e-7_rt)
     state(i,URHO) = rho
     state(i,UTEMP) = T
     state(i,UMX) = state(i,URHO)*vel

     eos_state % rho = state(i,URHO)
     eos_state % T   = state(i,UTEMP)

     call eos(eos_input_rt, eos_state)

     state(i,UEINT) = state(i,URHO) * eos_state % e 

     state(i,UEDEN) = state(i,UEINT) + 0.5e0_rt*state(i,UMX)**2/state(i,URHO)
     !state(i,UFS:UFS-1+nspec) = 0.e0_rt
     state(i,UFS) = state(i,URHO)
     !state(i,UFX) = state(i,URHO)
     !state(i,UFX+1) = state(i,URHO)



    enddo


end subroutine ca_initdata


subroutine ca_initrad(level,time,lo,hi,nrad, &
                      rad_state,rad_state_l1, &
                      rad_state_h1, &
                      delta,xlo,xhi)

    use probdata_module
    use meth_params_module, ONLY : UTEMP
    use fundamental_constants_module, only: a_rad
    use amrex_fort_module, only : rt => amrex_real
    use rad_params_module, only : xnu
    use blackbody_module, only : BGroup
    implicit none
        
    integer :: level, nrad
    integer :: lo(1), hi(1)
    integer :: rad_state_l1, rad_state_h1
    real(rt)    :: xlo(1), xhi(1), time, delta(1)
    real(rt)    :: rad_state(rad_state_l1:rad_state_h1, 0:nrad-1)
    integer :: i, igroup
    real(rt) :: xx, T

    do i = lo(1),hi(1)

    	xx = xlo(1)+delta(1)*(dble(i-lo(1)+0.5e0_rt))
      T = T_ej

      if (nrad .eq. 1) then
        rad_state(i,0) = a_rad*T**4
      else
        do igroup=0,nrad-1
          rad_state(i,igroup) = BGroup(T,xnu(igroup),xnu(igroup+1))
        enddo
      endif

    enddo






end subroutine ca_initrad




