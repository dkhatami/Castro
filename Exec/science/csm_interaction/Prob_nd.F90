subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(C,name="amrex_probinit")

      use probdata_module
      !use prob_params_module, only : center
      use amrex_constants_module, only: M_PI, FOUR3RD
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer :: init, namlen
      integer :: name(namlen)
      real(rt)      :: problo(3), probhi(3)

      integer :: untin,i

      namelist /fortin/ m_0, r_0, v_0, &
                        eta, beta, delt, tau, &
                        f_a, T_init,TT_0, &
                        M_csm, M_ej, dR_csm, kap, &
                        t_0, rho_0, E_0, h_csm, &
                        filter_rhomax, filter_timemax, v_max, &
                        n, d, s, f_r, tau_a, pl_ej, use_Trec, Trec, dT_rec, &
                        dR_w, rho_w, eps, p, fm_bo, ft_bo, xi_bo


      integer, parameter :: maxlen = 256
      character :: probin*(maxlen)




      do i = 1, namlen
        probin(i:i) = char(name(i))
      end do

      !namelist defaults

      m_0 = 1.98e33_rt ! characteristic mass in g (M_csm)
      r_0 = 1.e14_rt   ! characteristic radius in cm (R_csm)
      !v_0 = 1.e9_rt   ! characteristic velocity in cm/s (v_ej)

      kap = 0.2e0_rt ! characteristic opacity (cm^2/g)
      eps = 1.e-3_rt ! planck mean (thermalization frac)
      M_ej = 3.e0_rt*1.989e33_rt ! ejecta mass (g)

      use_Trec = -1.e0_rt ! whether or not to use a recombination opacity
      Trec = 5.5e3_rt ! recombination temperature
      dT_rec = 0.1e0_rt*Trec ! smoothing length

      eta = 1.e-2_rt ! ratio M_csm/M_ej
      beta = 3.333333e-2_rt ! = v_ej/c
      delt = 1.e0_rt ! ratio dR/R_csm
      tau = 1.e1_rt ! csm optical depth

      f_a = 1.e-4_rt ! ratio of maximum possible ambient to csm density
                     ! this is only set if the optical depth param is too large
                     ! to make sure the ambient density doesn't contribute to
                     ! late-time LC in cases of small kappa
      T_init = 1.e1_rt ! initial csm/ejecta temperature
      TT_0 = 1.e4_rt ! ambient temperature

      rho_w = 1.e-15_rt !interface wind density
      dR_w = 2e15_rt !wind width

      h_csm = 1.e-1_rt ! smoothing length ratio h/R_csm

      d = 0.e0_rt ! inner ejecta density profile power law
      n = 10.e0_rt ! outer ejecta density profile power law
      pl_ej = -1.e0_rt ! whether or not to use power law ejecta prescription
      s = 0.e0_rt ! CSM density profile power law

      f_r = 0.1e0_rt ! ejecta transition velocity coordinate

      tau_a = 1.e-3_rt ! ambient density optical depth

      p = 10.e0_rt ! breakout layer power law index
      fm_bo = -1.e0_rt ! fractional mass of breakout layer
      ft_bo = -1.e0_rt ! fractional optical depth of breakout layer

      filter_rhomax = 1.e-99_rt
      filter_timemax = 0.e0_rt

      xi_bo = 0.e0_rt

      v_max = 3.e9_rt


      open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      v_0 = beta*2.998e10_rt
      M_csm = eta*M_ej
      m_0 = M_csm

      !r_0 = sqrt(kap*eta*M_ej/tau)*sqrt((3.e0_rt-s)/(4.e0_rt*M_PI*abs(1.e0_rt-s)))
      !r_0 = r_0*sqrt(abs((1.e0_rt+delt)**(1.e0_rt-s)-1.e0_rt)/abs((1.e0_rt+delt)**(3.e0_rt-s)-1.e0_rt))
      t_0 = r_0/v_0 !characteristic time
      rho_0 = M_ej/(FOUR3RD*M_PI*r_0*r_0*r_0) !characteristic density
      E_0 = m_0*v_0*v_0 ! characteristic energy

      dR_csm = delt*r_0
      !kap = tau*r_0*r_0/m_0
      !kap = FOUR3RD*M_PI*tau*r_0*r_0/m_0/delt*((1.e0_rt+delt)**3.e0_rt-1.e0_rt)
      !kap = tau*r_0*r_0/m_0*4.e0_rt*M_PI*(1.e0_rt-s)/(3.e0_rt-s)
      !kap = kap*((1.e0_rt+delt)**(3.e0_rt-s)-1.e0_rt)/((1.e0_rt+delt)**(1.e0_rt-s)-1.e0_rt)

      t_d = sqrt(kap*m_0/(v_0*2.998e10_rt))

      rho_a = min(tau_a/(kap*probhi(1)),f_a*rho_0)



    !  center(1) = 0.e0_rt


end subroutine amrex_probinit


subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)

    use probdata_module
    use eos_module, only : eos
    use network, only : nspec,naux
    use eos_type_module, only: eos_t, eos_input_rt
    use amrex_constants_module, only: M_PI, FOUR3RD
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMZ, UEDEN, UEINT, UFS, UFX
    use amrex_fort_module, only : rt => amrex_real
    ! use prob_params_module, only: problo
    implicit none

    integer :: level, nscal
    integer, parameter :: nsub = 1
    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    real(rt)        :: xlo(3), xhi(3), time, delta(3)
    real(rt)        :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),nscal)

    real(rt)        :: xx, xxs, vol, dx_sub, vol_sub
!    real(rt)        :: R_s, eta_r, eta_v, v_t, M_ch, rho_w
    real(rt)        :: rho_ej, rho_csm, rho_o, R_t, rho_wa, rho_cw,rho_ca,rho_ww
    real(rt)        :: rho, vel, T, rho_sub, vel_sub, T_sub
    real(rt)        :: rho0_ej, rho0_csm, xi, R_bo, rho_bo
!    real(rt)		:: h_c, h_ej
    integer :: i,j,k, ii

    type(eos_t)  :: eos_state


    dx_sub = delta(1)/dble(nsub)

    R_t = f_r*r_0




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

    xi = xi_bo

!    if (fm_bo .gt. 0) then
!        xi = 1.e0_rt+(1.e0_rt-p/3.e0_rt)*(1.e0_rt-(1.e0_rt+delt)**(-3))*fm_bo
!        xi = xi**(1.e0_rt/(3.e0_rt-p))
!    else if (ft_bo .gt. 0) then
!        xi = 1.e0_rt - (p-1.e0_rt)/(1.e0_rt+1.e0_rt/delt)*ft_bo
!        xi = xi**(1.e0_rt/(1.e0_rt-p))
!    endif

    R_bo = r_0+dR_csm+xi*dR_csm



    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           xx = xlo(1)+delta(1)*(dble(i-lo(1)+0.5e0_rt))

           state(i,j,k,UMX:UMZ) = 0.0e0_rt


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

            if (pl_ej .gt. 0 ) then
              rho0_ej = (n-3.e0_rt)*(3.e0_rt-d)/(4.e0_rt*M_PI*(n-d))*M_ej/(f_r*r_0)**3
            else
              rho0_ej = 3.e0_rt/(4.e0_rt*M_PI)*M_ej/r_0**3
            endif

            rho0_csm = (3.e0_rt-s)/(4.e0_rt*M_PI)/((1.e0_rt+delt)**(3.e0_rt-s)-1.e0_rt)*M_csm/r_0**3



            rho_csm = rho0_csm*(xxs/r_0)**(-s)
            rho_ww = rho_w*(xxs/(r_0+dR_csm))**(-2.e0_rt)
            rho_bo = rho0_csm*((r_0+dR_csm)/r_0)**(-s)*(xxs/(r_0+dR_csm))**(-p)
            if (dR_w .gt. 0.e0_rt) then
                rho_cw = rho_csm + 0.5e0_rt*(rho_ww-rho_csm)*(1.e0_rt+tanh((xxs-(r_0+dR_csm))/(h_csm*r_0)))
                rho_wa = rho_ww + 0.5e0_rt*(rho_a-rho_ww)*(1.e0_rt+tanh((xxs-(r_0+dR_csm+dR_w))/(h_csm*(r_0+dR_csm+dR_w))))
                else
                    rho_ca = rho_csm + 0.5e0_rt*(rho_a-rho_csm)*(1.e0_rt+tanh((xxs-(r_0+dR_csm))/(h_csm*r_0)))
                endif
            if (xxs .lt. R_t) then
                rho_sub = rho0_ej!*(xxs/R_t)**(-d)
                if (pl_ej .gt. 0) then
                  rho_sub = rho_sub*(xxs/R_t)**(-d)
                endif
                vel_sub = v_0*(xxs/r_0)
                T_sub = T_init
            else if (xxs .lt. r_0) then
                rho_sub = rho0_ej!*(xxs/R_t)**(-n)
                if (pl_ej .gt. 0) then
                  rho_sub = rho_sub*(xxs/R_t)**(-n)
                endif
                vel_sub = v_0*(xxs/r_0)
                T_sub = T_init
            else if (xxs .lt. r_0+dR_csm) then
            !    if (dR_w .gt. 0.e0_rt) then
            !        rho_sub = rho_cw
            !    else
            !        rho_sub = rho_ca
            !    endif
                rho_sub = rho_csm
                vel_sub = 0.e0_rt
                T_sub = T_init
            else if (xxs .lt. R_bo) then
                rho_sub =  rho_bo
                vel_sub = 0.e0_rt
                T_sub = T_init
            else
              !if (dR_w .gt. 0.e0_rt) then
              !    rho_sub = rho_wa
              !else
              !    rho_sub = rho_ca
              !endif
                rho_sub = rho_a
                vel_sub = 0.e0_rt
                T_sub =TT_0
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
           state(i,j,k,URHO) = rho
           state(i,j,k,UTEMP) = T
           state(i,j,k,UMX) = state(i,j,k,URHO)*vel
           state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO)


          if (xx .lt. r_0) then
             state(i,j,k,UFX) = 1.e0_rt*rho
             state(i,j,k,UFX+1)= 0.e0_rt
             state(i,j,k,UFX+2) = 0.e0_rt
         else if (xx .lt. R_bo) then
           state(i,j,k,UFX) = 0.e0_rt
           state(i,j,k,UFX+1) = 1.e0_rt*rho
           state(i,j,k,UFX+2) = 0.e0_rt
         else
           state(i,j,k,UFX) = 0.e0_rt
           state(i,j,k,UFX+1) = 0.e0_rt
           state(i,j,k,UFX+2) = 1.e0_rt*rho
         endif

           eos_state % rho = state(i,j,k,URHO)
           eos_state % T   = state(i,j,k,UTEMP)
           eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
           eos_state % aux = state(i,j,k,UFX:UFX+naux-1) / state(i,j,k,URHO)

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e

           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5e0_rt*state(i,j,k,UMX)**2/state(i,j,k,URHO)
           !state(i,UFS:UFS-1+nspec) = 0.e0_rt

           !state(i,UFS) = state(i,URHO)
           !state(i,UFX) = state(i,URHO)
           !state(i,UFX+1) = state(i,URHO)


         enddo
       enddo
    enddo


end subroutine ca_initdata


subroutine ca_initrad(level,time,lo,hi,nrad, &
                      rad_state,rad_state_lo, &
                      rad_state_hi, &
                      delta,xlo,xhi)

    use probdata_module
    use meth_params_module, ONLY : UTEMP
    use fundamental_constants_module, only: a_rad
    use amrex_fort_module, only : rt => amrex_real
    use rad_params_module, only : xnu
    use blackbody_module, only : BGroup
    implicit none

    integer :: level, nrad
    integer :: lo(3), hi(3)
    integer :: rad_state_lo(3), rad_state_hi(3)
    real(rt)    :: xlo(3), xhi(3), time, delta(3)
    real(rt)    :: rad_state(rad_state_lo(1):rad_state_hi(1),rad_state_lo(2):rad_state_hi(2),rad_state_lo(3):rad_state_hi(3), 0:nrad-1)
    integer :: i,j,k, igroup
    real(rt) :: xx, T

    do k = lo(3),hi(3)
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)

        	xx = xlo(1)+delta(1)*(dble(i-lo(1)+0.5e0_rt))

          ! if (xx .lt. r_0+dR_csm+xi_bo*dR_csm) then
          !    T = T_init
          !  else
          !    T = TT_0
          !  endif
          T = T_init

          if (nrad .eq. 1) then
            rad_state(i,j,k,0) = a_rad*T**4
          else
            do igroup=0,nrad-1
              rad_state(i,j,k,igroup) = BGroup(T,xnu(igroup),xnu(igroup+1))
            enddo
          endif

        enddo
      enddo
    enddo






end subroutine ca_initrad
