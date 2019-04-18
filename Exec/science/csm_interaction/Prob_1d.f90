subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)

      use probdata_module
      use prob_params_module, only : center
      use amrex_constants_module, only: M_PI, FOUR3RD
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer :: init, namlen
      integer :: name(namlen)
      real(rt)      :: problo(1), probhi(1)

      integer :: untin,i

      namelist /fortin/ m_0, r_0, v_0, &
                        eta, beta, delt, tau, &
                        f_a, TT_0, &
                        M_csm, M_ej, dR_csm, kap, &
                        t_0, rho_0, E_0, h_csm, &
                        filter_rhomax, filter_timemax


      integer, parameter :: maxlen = 256
      character :: probin*(maxlen)


      do i = 1, namlen
        probin(i:i) = char(name(i))
      end do

      !namelist defaults

      m_0 = 1.98e33_rt ! characteristic mass in g (M_csm)
      r_0 = 1.e14_rt   ! characteristic radius in cm (R_csm)
      !v_0 = 1.e9_rt   ! characteristic velocity in cm/s (v_ej)

      eta = 1.e-2_rt ! ratio M_csm/M_ej
      beta = 1.e-1_rt ! = v_ej/c
      delt = 1.e0_rt ! ratio dR/R_csm
      tau = 1.e2_rt ! csm optical depth

      f_a = 1.e-6_rt ! ratio of ambient to csm density
      TT_0 = 1.e2_rt ! ambient temperature

      h_csm = 1.e-1_rt ! smoothing length ratio h/R_csm

      filter_rhomax = 1.e-99_rt
      filter_timemax = 0.e0_rt


      open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      v_0 = beta*2.998e10_rt
      M_csm = m_0
      M_ej = m_0/eta

      t_0 = r_0/v_0 !characteristic time
      rho_0 = M_ej/(FOUR3RD*M_PI*r_0*r_0*r_0) !characteristic density
      E_0 = m_0*v_0*v_0 ! characteristic energy

      dR_csm = delt*r_0
      !kap = tau*r_0*r_0/m_0
      kap = FOUR3RD*M_PI*tau*r_0*r_0/m_0/delt*((1.e0_rt+delt)**3.e0_rt-1.e0_rt)

      t_d = sqrt(kap*m_0/(v_0*2.998e10_rt))



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
!    real(rt)        :: R_s, eta_r, eta_v, v_t, M_ch, rho_w
    real(rt)        :: rho_ej, rho_csm, rho_o, rho_a
    real(rt)        :: rho, vel, T, rho_sub, vel_sub, T_sub
!    real(rt)		:: h_c, h_ej
    integer :: i, ii

    type(eos_t)  :: eos_state


    dx_sub = delta(1)/dble(nsub)

    ! eta_r = (n-3.e0_rt)*(3.e0_rt-d)/(n-d)/(4.e0_rt*M_PI)
    ! eta_v = sqrt(2.e0_rt*(5.e0_rt-d)*(n-5.e0_rt) &
    !             /((3.e0_rt-d)*(n-3.e0_rt)))
    !
    ! M_ch = M_ej / 1.4e0_rt
    ! v_t = 6.e8_rt*eta_v*sqrt(E_ej/M_ch)
    ! R_s = v_t*R_c/v_ej
    !
    !
    ! h_ej = 0.1e0_rt*R_s
    ! h_c = 0.02e0_rt*R_x1
    ! rho_w = Mdot_w*2.e33_rt/(3.154e7_rt)/(4.*M_PI*R_c**2*v_w)
    !
    !
    ! rho_sh = M_sh*2.e33_rt/(4.e0_rt*M_PI*(f_sh*R_sh)**3)
    ! rho_sh = rho_sh/(sqrt(M_PI)/2.e0_rt*(1.e0_rt+2.e0_rt/f_sh**2)-N_sh*exp(-N_sh*N_sh))
    !
    ! rho_x = M_x*2.e33_rt/(4.e0_rt*M_PI*R_x0*R_x0*(R_x1-R_x0))

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


      ! rho0 = eta_r*M_ej*2.e33_rt/(R_s)**3*(R_s/xxs)**d
      ! rho1 = eta_r*M_ej*2.e33_rt/(R_s)**3*(R_s/xxs)**n
      ! rhow = rho_w*(R_c/xxs)**s

      ! rhoej = rho0 + 0.5e0_rt*(rho1-rho0)*(1.e0_rt+tanh((xxs-R_s)/h_ej))

      rho_csm = M_csm/(FOUR3RD*M_PI)/((r_0+dR_csm)**3-r_0**3)
      rho_a = f_a*rho_csm
      rho_o = rho_csm + 0.5e0_rt*(rho_a-rho_csm)*(1.e0_rt+tanh((xxs-(r_0+dR_csm))/(h_csm*r_0)))

      if (xxs .lt. r_0) then
          rho_sub = rho_0
          vel_sub = v_0*(xxs/r_0)
          T_sub = TT_0
      else
          rho_sub = rho_o
          T_sub = TT_0
          vel_sub = 0.e0_rt
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
      T = TT_0

      if (nrad .eq. 1) then
        rad_state(i,0) = a_rad*T**4
      else
        do igroup=0,nrad-1
          rad_state(i,igroup) = BGroup(T,xnu(igroup),xnu(igroup+1))
        enddo
      endif

    enddo






end subroutine ca_initrad
