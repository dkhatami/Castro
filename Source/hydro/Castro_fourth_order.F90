! advection routines in support of method of lines integration
!
module fourth_order_hydro

  implicit none

contains

  subroutine ca_fourth_single_stage(lo, hi, time, domlo, domhi, &
                                    uin, uin_lo, uin_hi, &
                                    uout, uout_lo, uout_hi, &
                                    q, q_lo, q_hi, &
                                    q_bar, q_bar_lo, q_bar_hi, &
                                    qaux, qa_lo, qa_hi, &
                                    qaux_bar, qa_bar_lo, qa_bar_hi, &
#ifdef DIFFUSION
                                    T_cc, Tcc_lo, Tcc_hi, &
#endif
                                    srcU, srU_lo, srU_hi, &
                                    update, updt_lo, updt_hi, &
                                    dx, dt, &
                                    flx, flx_lo, flx_hi, &
#if AMREX_SPACEDIM >= 2
                                    fly, fly_lo, fly_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                    flz, flz_lo, flz_hi, &
#endif
                                    area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                                    area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                    area3, area3_lo, area3_hi, &
#endif
#if AMREX_SPACEDIM < 3
                                    pradial, p_lo, p_hi, &
                                    dloga, dloga_lo, dloga_hi, &
#endif
                                    vol, vol_lo, vol_hi, &
                                    verbose) bind(C, name="ca_fourth_single_stage")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : NQ, NVAR, NGDNV, NQAUX, GDPRES, &
                                   UTEMP, UEINT, USHK, GDU, GDV, GDW, UMX, &
                                   use_flattening, QPRES, NQAUX, &
                                   QTEMP, QFS, QFX, QREINT, QRHO, QGAME, QGC, &
                                   first_order_hydro, difmag, hybrid_riemann, &
                                   limit_fluxes_on_small_dens, ppm_temp_fix, &
#ifdef DIFFUSION
                                   diffuse_temp, &
#endif
                                   do_hydro
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens, ca_shock, &
                                      normalize_species_fluxes, avisc

    use castro_error_module
    use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
    use flatten_module, only: ca_uflatten
    use riemann_module, only: riemann_state
    use riemann_util_module, only: compute_flux_q
    use fourth_order
    use amrex_fort_module, only : rt => amrex_real
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
    use riemann_util_module, only : ca_store_godunov_state
#endif
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux
    use prob_params_module, only : dg, coord_type

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  domlo(3), domhi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: q_bar_lo(3), q_bar_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: qa_bar_lo(3), qa_bar_hi(3)
#ifdef DIFFUSION
    integer, intent(in) :: Tcc_lo(3), Tcc_hi(3)
#endif
    integer, intent(in) :: srU_lo(3), srU_hi(3)
    integer, intent(in) :: updt_lo(3), updt_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: fly_lo(3), fly_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flz_lo(3), flz_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: p_lo(3), p_hi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    integer, intent(in) :: vol_lo(3), vol_hi(3)

    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: q_bar(q_bar_lo(1):q_bar_hi(1), q_bar_lo(2):q_bar_hi(2), q_bar_lo(3):q_bar_hi(3), NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: qaux_bar(qa_bar_lo(1):qa_bar_hi(1), qa_bar_lo(2):qa_bar_hi(2), qa_bar_lo(3):qa_bar_hi(3), NQAUX)
#ifdef DIFFUSION
    real(rt), intent(inout) :: T_cc(Tcc_lo(1):Tcc_hi(1), Tcc_lo(2):Tcc_hi(2), Tcc_lo(3):Tcc_hi(3), 1)
#endif
    real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1), flx_lo(2):flx_hi(2), flx_lo(3):flx_hi(3), NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: fly(fly_lo(1):fly_hi(1), fly_lo(2):fly_hi(2), fly_lo(3):fly_hi(3), NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: flz(flz_lo(1):flz_hi(1), flz_lo(2):flz_hi(2), flz_lo(3):flz_hi(3), NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
#if AMREX_SPACEDIM <= 2
    real(rt), intent(inout) :: pradial(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1), dloga_lo(2):dloga_hi(2), dloga_lo(3):dloga_hi(3))
#endif
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(in) :: dx(3), dt, time

#ifndef RADIATION
    ! Automatic arrays for workspace
    real(rt), pointer :: flatn(:,:,:)
    real(rt), pointer :: avisx(:,:,:), avisy(:,:,:), avisz(:,:,:)

    ! Edge-centered primitive variables (Riemann state)
    real(rt), pointer :: qx_avg(:,:,:,:)
    real(rt), pointer :: qy_avg(:,:,:,:)
    real(rt), pointer :: qz_avg(:,:,:,:)

    real(rt), pointer :: qx(:,:,:,:)
    real(rt), pointer :: qy(:,:,:,:)
    real(rt), pointer :: qz(:,:,:,:)

#ifdef HYBRID_MOMENTUM
    real(rt), pointer :: qgdnvx(:,:,:,:)
    real(rt), pointer :: qgdnvy(:,:,:,:)
    real(rt), pointer :: qgdnvz(:,:,:,:)
#endif

    ! Temporaries (for now)
    real(rt), pointer :: flx_avg(:,:,:,:)
    real(rt), pointer :: fly_avg(:,:,:,:)
    real(rt), pointer :: flz_avg(:,:,:,:)

    real(rt), pointer :: shk(:,:,:)

    real(rt), pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
    real(rt), pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

    integer :: It_lo(3), It_hi(3)
    integer :: st_lo(3), st_hi(3)
    integer :: shk_lo(3), shk_hi(3)

    real(rt) :: lap
    integer :: i, j, k, n, m
    integer :: is_avg

    type (eos_t) :: eos_state

    ! artifical viscosity strength
    real(rt), parameter :: alpha = 0.3_rt
    real(rt) :: avisc_coeff

    ! to do 4th order for axisymmetry, we need to derive the transformations between
    ! averages and cell-centers with the correct volume terms in the integral.
#ifndef AMREX_USE_CUDA
    if (coord_type > 0) then
       call castro_error("Error: fourth order not implemented for axisymmetric")
    endif
#endif

    It_lo = lo(:) - dg(:)
    It_hi = hi(:) + dg(:)

    shk_lo(:) = lo(:) - dg(:)
    shk_hi(:) = hi(:) + dg(:)

    call bl_allocate(avisx, lo, hi+dg)
#if BL_SPACEDIM >= 2
    call bl_allocate(avisy, lo, hi+dg)
#endif
#if BL_SPACEDIM == 3
    call bl_allocate(avisz, lo, hi+dg)
#endif

    call bl_allocate(qx_avg, q_lo, q_hi, NQ)
    call bl_allocate(flx_avg, q_lo, q_hi, NVAR)
#if AMREX_SPACEDIM >= 2
    call bl_allocate(qx, q_lo, q_hi, NQ)

    call bl_allocate(qy_avg, q_lo, q_hi, NQ)
    call bl_allocate(qy, q_lo, q_hi, NQ)
    call bl_allocate(fly_avg, q_lo, q_hi, NVAR)
#endif
#if AMREX_SPACEDIM == 3
    call bl_allocate(qz_avg, q_lo, q_hi, NQ)
    call bl_allocate(qz, q_lo, q_hi, NQ)
    call bl_allocate(flz_avg, q_lo, q_hi, NVAR)
#endif

    call bl_allocate(qxm, q_lo, q_hi, NQ)
    call bl_allocate(qxp, q_lo, q_hi, NQ)

#if AMREX_SPACEDIM >= 2
    call bl_allocate(qym, q_lo, q_hi, NQ)
    call bl_allocate(qyp, q_lo, q_hi, NQ)
#endif
#if AMREX_SPACEDIM == 3
    call bl_allocate(qzm, q_lo, q_hi, NQ)
    call bl_allocate(qzp, q_lo, q_hi, NQ)
#endif

    call bl_allocate(shk, shk_lo, shk_hi)

#ifdef SHOCK_VAR
     call ca_shock(lo-dg, hi+dg, &
                   q_bar, q_bar_lo, q_bar_hi, &
                   shk, shk_lo, shk_hi, &
                   dx)

     ! We'll update the shock data for future use in the burning step.
     ! For the update, we are starting from USHK == 0 (set at the
     ! beginning of the timestep) and we need to divide by dt since
     ! we'll be multiplying that for the update calculation.

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              update(i,j,k,USHK) = shk(i,j,k) / dt
           enddo
        enddo
     enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:,:) = ZERO
    end if
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call ca_shock(lo-dg, hi+dg, &
                     q_bar, q_bar_lo, q_bar_hi, &
                     shk, shk_lo, shk_hi, &
                     dx)
    else
       shk(:,:,:) = ZERO
    end if
#endif

    ! Compute flattening coefficient for slope calculations -- we do
    ! this with q_bar, since we need all of the ghost cells
    call bl_allocate(flatn, q_bar_lo, q_bar_hi)

    if (use_flattening == 1) then
       call ca_uflatten(lo-dg, hi+dg, &
                        q_bar, q_bar_lo, q_bar_hi, &
                        flatn, q_bar_lo, q_bar_hi, QPRES)
    else
       flatn = ONE
    end if

    ! do the reconstruction here -- get the interface states

    ! x-interfaces
    call states(1, &
                q, q_lo, q_hi, &
                flatn, q_bar_lo, q_bar_hi, &
                qxm, qxp, q_lo, q_hi, &
                lo, hi, &
                domlo, domhi)

#if AMREX_SPACEDIM >= 2
    ! y-interfaces
    call states(2, &
                q, q_lo, q_hi, &
                flatn, q_bar_lo, q_bar_hi, &
                qym, qyp, q_lo, q_hi, &
                lo, hi, &
                domlo, domhi)
#endif

#if AMREX_SPACEDIM == 3
    ! z-interfaces
    call states(3, &
                q, q_lo, q_hi, &
                flatn, q_bar_lo, q_bar_hi, &
                qzm, qzp, q_lo, q_hi, &
                lo, hi, &
                domlo, domhi)
#endif

    ! this is where we would implement ppm_temp_fix


    ! solve the Riemann problems -- we just require the interface state
    ! at this point

    ! get <q> and F(<q>) on the x interfaces
    call riemann_state(qxm, q_lo, q_hi, &
                       qxp, q_lo, q_hi, 1, 1, &
                       qx_avg, q_lo, q_hi, &
                       qaux, qa_lo, qa_hi, &
                       1, &
                       [lo(1), lo(2)-dg(2), lo(3)-dg(3)], &
                       [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                       domlo, domhi)

    call compute_flux_q([lo(1), lo(2)-dg(2), lo(3)-dg(3)], &
                        [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                        qx_avg, q_lo, q_hi, &
                        flx_avg, q_lo, q_hi, &
                        1)

    if (do_hydro == 0) then
       flx_avg(:,:,:,:) = ZERO
    end if

#ifdef DIFFUSION
    if (diffuse_temp == 1) then
       is_avg = 1
       call add_diffusive_flux([lo(1), lo(2)-dg(2), lo(3)-dg(3)], &
                               [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                               q, q_lo, q_hi, NQ, QTEMP, &
                               qx_avg, q_lo, q_hi, &
                               flx_avg, q_lo, q_hi, &
                               dx, 1, is_avg)
    end if
#endif

#if AMREX_SPACEDIM >= 2
    ! get <q> and F(<q>) on the y interfaces
    call riemann_state(qym, q_lo, q_hi, &
                       qyp, q_lo, q_hi, 1, 1, &
                       qy_avg, q_lo, q_hi, &
                       qaux, qa_lo, qa_hi, &
                       2, &
                       [lo(1)-1, lo(2), lo(3)-dg(3)], &
                       [hi(1)+1, hi(2)+1, hi(3)+dg(3)], &
                       domlo, domhi)

    call compute_flux_q([lo(1)-1, lo(2), lo(3)-dg(3)], &
                        [hi(1)+1, hi(2)+1, hi(3)+dg(3)], &
                        qy_avg, q_lo, q_hi, &
                        fly_avg, q_lo, q_hi, &
                        2)

    if (do_hydro == 0) then
       fly_avg(:,:,:,:) = ZERO
    end if

#ifdef DIFFUSION
    if (diffuse_temp == 1) then
       is_avg = 1
       call add_diffusive_flux([lo(1)-1, lo(2), lo(3)-dg(3)], &
                               [hi(1)+1, hi(2)+1, hi(3)+dg(3)], &
                               q, q_lo, q_hi, NQ, QTEMP, &
                               qy_avg, q_lo, q_hi, &
                               fly_avg, q_lo, q_hi, &
                               dx, 2, is_avg)
    end if
#endif
#endif

#if AMREX_SPACEDIM == 3
    ! get <q> and F(<q>) on the z interfaces
    call riemann_state(qzm, q_lo, q_hi, &
                       qzp, q_lo, q_hi, 1, 1, &
                       qz_avg, q_lo, q_hi, &
                       qaux, qa_lo, qa_hi, &
                       3, &
                       [lo(1)-1, lo(2)-1, lo(3)], &
                       [hi(1)+1, hi(2)+1, hi(3)+1], &
                       domlo, domhi)

    call compute_flux_q([lo(1)-1, lo(2)-1, lo(3)], &
                        [hi(1)+1, hi(2)+1, hi(3)+1], &
                        qz_avg, q_lo, q_hi, &
                        flz_avg, q_lo, q_hi, &
                        3)

    if (do_hydro == 0) then
       flz_avg(:,:,:,:) = ZERO
    end if

#ifdef DIFFUSION
    if (diffuse_temp == 1) then
       is_avg = 1
       call add_diffusive_flux([lo(1)-1, lo(2)-1, lo(3)], &
                               [hi(1)+1, hi(2)+1, hi(3)+1], &
                               q, q_lo, q_hi, NQ, QTEMP, &
                               qz_avg, q_lo, q_hi, &
                               flz_avg, q_lo, q_hi, &
                               dx, 3, is_avg)
    end if
#endif
#endif


    call bl_deallocate(flatn)

    call bl_deallocate(qxm)
    call bl_deallocate(qxp)

#if AMREX_SPACEDIM >= 2
    call bl_deallocate(qym)
    call bl_deallocate(qyp)
#endif

#if AMREX_SPACEDIM == 3
    call bl_deallocate(qzm)
    call bl_deallocate(qzp)
#endif

    call bl_deallocate(shk)


    ! we now have the face-average interface states and fluxes evaluated with these

    ! Note: for 1-d, we are done


    ! construct the face-center interface states

#if AMREX_SPACEDIM >= 2
    ! x-interfaces
    do n = 1, NQ
       if (n == QGAME .or. n == QGC .or. n == QTEMP) cycle

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                lap = transx_laplacian(i, j, k, n, &
                                       qx_avg, q_lo, q_hi, NQ, &
                                       domlo, domhi)

                qx(i,j,k,n) = qx_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap

             end do
          end do
       end do
    end do

    ! y-interfaces
    do n = 1, NQ
       if (n == QGAME .or. n == QGC .or. n == QTEMP) cycle

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                lap = transy_laplacian(i, j, k, n, &
                                       qy_avg, q_lo, q_hi, NQ, &
                                       domlo, domhi)

                qy(i,j,k,n) = qy_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap

             end do
          end do
       end do
    end do

#if AMREX_SPACEDIM == 3
    ! z-interfaces
    do n = 1, NQ
       if (n == QGAME .or. n == QGC .or. n == QTEMP) cycle

       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lap = transz_laplacian(i, j, k, n, &
                                       qz_avg, q_lo, q_hi, NQ, &
                                       domlo, domhi)

                qz(i,j,k,n) = qz_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap

             end do
          end do
       end do
    end do

#endif

    ! compute face-centered fluxes
    ! these will be stored in flx, fly, flz
    call compute_flux_q([lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], &
                         qx, q_lo, q_hi, &
                         flx, flx_lo, flx_hi, &
                         1)

    if (do_hydro == 0) then
       flx(flx_lo(1):flx_hi(1), flx_lo(2):flx_hi(2), flx_lo(3):flx_hi(3), :) = ZERO
    end if

#ifdef DIFFUSION
    if (diffuse_temp == 1) then
       is_avg = 0
       call add_diffusive_flux([lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], &
                                T_cc, Tcc_lo, Tcc_hi, 1, 1, &
                                qx, q_lo, q_hi, &
                                flx, flx_lo, flx_hi, &
                                dx, 1, is_avg)
    end if
#endif

#if AMREX_SPACEDIM >= 2
    call compute_flux_q([lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], &
                        qy, q_lo, q_hi, &
                        fly, fly_lo, fly_hi, &
                        2)

    if (do_hydro == 0) then
       fly(fly_lo(1):fly_hi(1), fly_lo(2):fly_hi(2), fly_lo(3):fly_hi(3), :) = ZERO
    end if

#ifdef DIFFUSION
    if (diffuse_temp == 1) then
       is_avg = 0
       call add_diffusive_flux([lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], &
                               T_cc, Tcc_lo, Tcc_hi, 1, 1, &
                               qy, q_lo, q_hi, &
                               fly, fly_lo, fly_hi, &
                               dx, 2, is_avg)
    end if
#endif

#endif

#if AMREX_SPACEDIM == 3
    call compute_flux_q([lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1], &
                        qz, q_lo, q_hi, &
                        flz, flz_lo, flz_hi, &
                        3)

    if (do_hydro == 0) then
       flz(flz_lo(1):flz_hi(1), flz_lo(2):flz_hi(2), flz_lo(3):flz_hi(3), :) = ZERO
    end if

#ifdef DIFFUSION
    if (diffuse_temp == 1) then
       is_avg = 0
       call add_diffusive_flux([lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1], &
                               T_cc, Tcc_lo, Tcc_hi, 1, 1, &
                               qz, q_lo, q_hi, &
                               flz, flz_lo, flz_hi, &
                               dx, 3, is_avg)
    end if
#endif

#endif

    call bl_deallocate(qx)
#if AMREX_SPACEDIM >= 2
    call bl_deallocate(qy)
#endif
#if AMREX_SPACEDIM == 3
    call bl_deallocate(qz)
#endif

    ! compute the final fluxes (as an average over the interface), this
    ! requires a transverse correction.  Note, we don't need to do anything
    ! to get the average of the Godunov states over the interface--this is
    ! essentially what qx_avg already is

    ! x-interfaces
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1

                lap = transx_laplacian(i, j, k, n, &
                                       flx_avg, q_lo, q_hi, NVAR, &
                                       domlo, domhi)

                flx(i,j,k,n) = flx(i,j,k,n) + 1.0_rt/24.0_rt * lap

             end do
          end do
       end do
    end do

    ! y-interfaces
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)

                lap = transy_laplacian(i, j, k, n, &
                                       fly_avg, q_lo, q_hi, NVAR, &
                                       domlo, domhi)

                fly(i,j,k,n) = fly(i,j,k,n) + 1.0_rt/24.0_rt * lap

             end do
          end do
       end do
    end do


#if AMREX_SPACEDIM == 3
    ! z-interfaces
    do n = 1, NVAR
       do k = lo(3), hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                lap = transz_laplacian(i, j, k, n, &
                                       flz_avg, q_lo, q_hi, NVAR, &
                                       domlo, domhi)

                flz(i,j,k,n) = flz(i,j,k,n) + 1.0_rt/24.0_rt * lap

             end do
          end do
       end do
    end do
#endif

#else
    ! for 1-d, we just copy flx_avg -> flx, since there is no face averaging
    flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = flx_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
#endif

    if (do_hydro == 1) then

       ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
       call avisc(lo, hi, &
                  q_bar, q_bar_lo, q_bar_hi, &
                  qaux_bar, qa_bar_lo, qa_bar_hi, &
                  dx, avisx, lo, hi+dg, 1)

#if BL_SPACEDIM >= 2
       call avisc(lo, hi, &
                  q_bar, q_bar_lo, q_bar_hi, &
                  qaux_bar, qa_bar_lo, qa_bar_hi, &
                  dx, avisy, lo, hi+dg, 2)
#endif

#if BL_SPACEDIM == 3
       call avisc(lo, hi, &
                  q_bar, q_bar_lo, q_bar_hi, &
                  qaux_bar, qa_bar_lo, qa_bar_hi, &
                  dx, avisz, lo, hi+dg, 3)
#endif

       ! avisc_coefficient is the coefficent we use.  The McCorquodale &
       ! Colella paper suggest alpha = 0.3, but our other hydro solvers use
       ! a coefficient on the divergence that defaults to 0.1, so we
       ! normalize to that value, to allow for adjustments
       avisc_coeff = alpha * (difmag / 0.1_rt)

       do n = 1, NVAR

          if ( n == UTEMP ) then
             flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
             fly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
             flz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

#ifdef SHOCK_VAR
          else if ( n == USHK ) then
             flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
             fly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
             flz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif
#endif

          else

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)+1

                      flx(i,j,k,n) = flx(i,j,k,n) + &
                           avisc_coeff * avisx(i,j,k) * (uin(i,j,k,n) - uin(i-1,j,k,n))
                   end do
                end do
             end do
#if AMREX_SPACEDIM >= 2
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)+1
                   do i = lo(1), hi(1)

                      fly(i,j,k,n) = fly(i,j,k,n) + &
                           avisc_coeff * avisy(i,j,k) * (uin(i,j,k,n) - uin(i,j-1,k,n))
                   end do
                end do
             end do
#endif
#if AMREX_SPACEDIM == 3
             do k = lo(3), hi(3)+1
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

                      flz(i,j,k,n) = flz(i,j,k,n) + &
                           avisc_coeff * avisz(i,j,k) * (uin(i,j,k,n) - uin(i,j,k-1,n))
                   end do
                end do
             end do
#endif
          end if

       end do

       call normalize_species_fluxes(flx_lo, flx_hi, flx, flx_lo, flx_hi)
#if AMREX_SPACEDIM >= 2
       call normalize_species_fluxes(fly_lo, fly_hi, fly, fly_lo, fly_hi)
#endif
#if AMREX_SPACEDIM == 3
       call normalize_species_fluxes(flz_lo, flz_hi, flz, flz_lo, flz_hi)
#endif

    endif

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 1
                update(i,j,k,n) = update(i,j,k,n) + &
                     (flx(i,j,k,n) * area1(i,j,k) - flx(i+1,j,k,n) * area1(i+1,j,k) ) / vol(i,j,k)

#elif AMREX_SPACEDIM == 2
                update(i,j,k,n) = update(i,j,k,n) + &
                     (flx(i,j,k,n) * area1(i,j,k) - flx(i+1,j,k,n) * area1(i+1,j,k) + &
                      fly(i,j,k,n) * area2(i,j,k) - fly(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)

#else
                update(i,j,k,n) = update(i,j,k,n) + &
                     (flx(i,j,k,n) * area1(i,j,k) - flx(i+1,j,k,n) * area1(i+1,j,k) + &
                      fly(i,j,k,n) * area2(i,j,k) - fly(i,j+1,k,n) * area2(i,j+1,k) + &
                      flz(i,j,k,n) * area3(i,j,k) - flz(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)
#endif

#if AMREX_SPACEDIM == 1
                if (do_hydro == 1) then
                   if (n == UMX) then
                      update(i,j,k,UMX) = update(i,j,k,UMX) - &
                           ( qx_avg(i+1,j,k,QPRES) - qx_avg(i,j,k,QPRES) ) / dx(1)
                   end if
                endif
#endif

#if AMREX_SPACEDIM == 2
                if (do_hydro == 1) then
                   if (n == UMX) then
                      ! add the pressure source term for axisymmetry
                      if (coord_type > 0) then
                         update(i,j,k,n) = update(i,j,k,n) - (qx_avg(i+1,j,k,QPRES) - qx_avg(i,j,k,QPRES))/ dx(1)
                      end if
                   end if
                endif
#endif

                update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

             end do
          end do
       end do
    end do

#if AMREX_SPACEDIM == 3
#ifdef HYBRID_MOMENTUM
    call bl_allocate(qgdnvx, q_lo, q_hi, NGDNV)
    call bl_allocate(qgdnvy, q_lo, q_hi, NGDNV)
    call bl_allocate(qgdnvz, q_lo, q_hi, NGDNV)

    call ca_store_godunov_state(lo, hi+dg, &
                                qx_avg, q_lo, q_hi, &
                                qgdnvx, q_lo, q_hi)

    call ca_store_godunov_state(lo, hi+dg, &
                                qy_avg, q_lo, q_hi, &
                                qgdnvy, q_lo, q_hi)

    call ca_store_godunov_state(lo, hi+dg, &
                                qz_avg, q_lo, q_hi, &
                                qgdnvz, q_lo, q_hi)

    call add_hybrid_advection_source(lo, hi, dt, &
                                     update, uout_lo, uout_hi, &
                                     qgdnvx, flx_lo, flx_hi, &
                                     qgdnvy, fly_lo, fly_hi, &
                                     qgdnvz, flz_lo, flz_hi)
    call bl_deallocate(qgdnvx)
    call bl_deallocate(qgdnvy)
    call bl_deallocate(qgdnvz)
#endif
#endif



    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                flx(i,j,k,n) = dt * flx(i,j,k,n) * area1(i,j,k)

#if AMREX_SPACEDIM == 1
                if (coord_type .eq. 0 .and. n == UMX) then
                   flx(i,j,k,n) = flx(i,j,k,n) + &
                        dt * area1(i,j,k) * qx_avg(i,j,k,QPRES)
                endif
#endif

             end do
          end do
       end do
    end do

#if AMREX_SPACEDIM >= 2
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                fly(i,j,k,n) = dt * fly(i,j,k,n) * area2(i,j,k)
             end do
          end do
       end do
    end do
#endif

#if AMREX_SPACEDIM == 3
    do n = 1, NVAR
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flz(i,j,k,n) = dt * flz(i,j,k,n) * area3(i,j,k)
             end do
          end do
       end do
    end do
#endif

#if AMREX_SPACEDIM < 3
    if (coord_type > 0) then
       pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = &
            qx_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),QPRES) * dt
    end if
#endif

    call bl_deallocate(avisx)
#if BL_SPACEDIM >= 2
    call bl_deallocate(avisy)
#endif
#if BL_SPACEDIM == 3
    call bl_deallocate(avisz)
#endif

    call bl_deallocate(qx_avg)
    call bl_deallocate(flx_avg)

#if AMREX_SPACEDIM >= 2
    call bl_deallocate(qy_avg)
    call bl_deallocate(fly_avg)
#endif
#if AMREX_SPACEDIM == 3
    call bl_deallocate(qz_avg)
    call bl_deallocate(flz_avg)
#endif
#else
#ifndef AMREX_USE_CUDA
   ! RADIATION check
    call castro_error("ERROR: ca_fourth_single_stage does not support radiation")
#endif
#endif
  end subroutine ca_fourth_single_stage

#ifdef DIFFUSION
  subroutine add_diffusive_flux(lo, hi, &
                                q, q_lo, q_hi, ncomp, temp_comp, &
                                qint, qi_lo, qi_hi, &
                                F, F_lo, F_hi, &
                                dx, idir, is_avg)
    ! add the diffusive flux to the energy fluxes
    !

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NQ, NVAR, &
                                   URHO, &
                                   UEDEN, UEINT, UTEMP, &
                                   QRHO, QREINT, QFS
    use eos_type_module, only : eos_t, eos_input_re
    use conductivity_module, only : conducteos
    use network, only : nspec
    use castro_error_module, only : castro_error

    integer, intent(in) :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qi_lo(3), qi_hi(3)
    integer, intent(in) :: F_lo(3), F_hi(3)
    integer, intent(in) :: ncomp, temp_comp

    real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), ncomp)
    real(rt), intent(in) :: qint(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), NQ)
    real(rt), intent(out) :: F(F_lo(1):F_hi(1), F_lo(2):F_hi(2), F_lo(3):F_hi(3), NVAR)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    integer, intent(in) :: is_avg

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: dTdx

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state % rho = qint(i,j,k,QRHO)
             eos_state % T = q(i,j,k,temp_comp)   ! initial guess
             eos_state % e = qint(i,j,k,QREINT) / qint(i,j,k,QRHO)
             eos_state % xn(:) = qint(i,j,k,QFS:QFS-1+nspec)

             call conducteos(eos_input_re, eos_state)

             if (idir == 1) then

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i+1,j,k,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i-1,j,k,temp_comp) + q(i-2,j,k,temp_comp))/(24.0_rt * dx(1))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i+1,j,k,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i-1,j,k,temp_comp) + q(i-2,j,k,temp_comp))/(12.0_rt * dx(1))
                end if

             else if (idir == 2) then

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i,j+1,k,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i,j-1,k,temp_comp) + q(i,j-2,k,temp_comp))/(24.0_rt * dx(2))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i,j+1,k,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i,j-1,k,temp_comp) + q(i,j-2,k,temp_comp))/(12.0_rt * dx(2))
                end if

             else

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i,j,k+1,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i,j,k-1,temp_comp) + q(i,j,k-2,temp_comp))/(24.0_rt * dx(3))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i,j,k+1,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i,j,k-1,temp_comp) + q(i,j,k-2,temp_comp))/(12.0_rt * dx(3))
                end if

             endif

             F(i,j,k,UEINT) = F(i,j,k,UEINT) - eos_state % conductivity * dTdx
             F(i,j,k,UEDEN) = F(i,j,k,UEDEN) - eos_state % conductivity * dTdx

          end do
       end do
    end do

  end subroutine add_diffusive_flux
#endif

end module fourth_order_hydro
