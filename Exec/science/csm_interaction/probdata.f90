module probdata_module

      use amrex_fort_module, only : rt=>amrex_real

      real(rt)          , save :: filter_rhomax, filter_timemax

      real(rt)        , save :: m_0, r_0, v_0
      real(rt)        , save :: n, d, s, f_r
      real(rt)        , save :: pl_ej
      real(rt)        , save :: eta, beta, delt, tau
      real(rt)        , save :: tau_a
      real(rt)        , save :: M_csm, M_ej, dR_csm, kap
      real(rt)        , save :: f_a, TT_0, h_csm
      real(rt)        , save :: t_0, rho_0, E_0, t_d,v_max
      real(rt)        , save :: use_Trec, Trec, dT_rec
      real(rt)        , save :: rho_w, dR_w, eps
      real(rt)        , save :: p, fm_bo, ft_bo, xi_bo
end module probdata_module
