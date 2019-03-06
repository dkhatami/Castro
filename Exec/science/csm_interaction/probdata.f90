module probdata_module

      use amrex_fort_module, only : rt=>amrex_real
      real(rt)          , save :: M_ej, E_ej, d, n, v_ej, T_ej, R_c
      real(rt)			, save :: Mdot_w, v_w, T_w, rho_a, T_a, s
      real(rt)			, save :: R_sh, M_sh, N_sh, f_sh, v_sh, T_sh
      real(rt)          , save :: M_x, R_x0, R_x1
      real(rt)          , save :: filter_rhomax, filter_timemax
end module probdata_module

