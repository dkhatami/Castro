subroutine ca_dercd_mask(mask,mask_lo,mask_hi,ncomp_mask, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo,hi,domlo,domhi, &
                         dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_dercd_mask')


    use amrex_constants_module
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: mask_lo(3), mask_hi(3), ncomp_mask ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 7 (Density, 3xMomentum, Eint, Shk,Adv)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    double precision :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3),ncomp_mask)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no


   integer :: i,j,k

   mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO


   do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

    !       if(u(i,j,k,6) > ZERO) then

                 if(u(i,j,k,7)/u(i,j,k,1) >= 0.999e0_rt) then
                   mask(i,j,k,1) = ONE
          !   endif
           endif
         enddo
      enddo
    enddo



end subroutine ca_dercd_mask


subroutine ca_dercsm_mask(mask,mask_lo,mask_hi,ncomp_mask, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo,hi,domlo,domhi, &
                         dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_dercsm_mask')


    use amrex_constants_module

    implicit none

    integer          :: mask_lo(3), mask_hi(3), ncomp_mask ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 3 (Density, Aux2, Aux3)
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    double precision :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3),ncomp_mask)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no


   integer :: i,j,k

   mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO


   do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

           if(u(i,j,k,2) > ZERO) then
                 mask(i,j,k,1) = ONE
           endif
         enddo
      enddo
    enddo



end subroutine ca_dercsm_mask


subroutine ca_derdtau(dtau,dtau_lo,dtau_hi,ncomp_dtau, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo,hi,domlo,domhi, &
                         dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derdtau')


    use amrex_constants_module
    use opacity_table_module
    use probdata_module, only : kap

    implicit none

    integer          :: dtau_lo(3), dtau_hi(3), ncomp_dtau ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 2
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    double precision :: dtau(dtau_lo(1):dtau_hi(1),dtau_lo(2):dtau_hi(2),dtau_lo(3):dtau_hi(3),ncomp_dtau)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no


   integer :: i,j,k

   double precision :: kp, kr, temp, rho, nu

   dtau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO
   nu = ZERO


   do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          rho  = u(i,j,k,1)
          temp = u(i,j,k,2)

           call get_opacities(kp,kr,rho,temp,rho,nu,.false.,.true.)

           dtau(i,j,k,1) = kr*dx(1)

         enddo
      enddo
    enddo



end subroutine ca_derdtau


subroutine ca_dertau_r(dtau,dtau_lo,dtau_hi,ncomp_dtau, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo,hi,domlo,domhi, &
                         dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_dertau_r')


    use amrex_constants_module
    use probdata_module, only : kap

    implicit none

    integer          :: dtau_lo(3), dtau_hi(3), ncomp_dtau ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 1
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    double precision :: dtau(dtau_lo(1):dtau_hi(1),dtau_lo(2):dtau_hi(2),dtau_lo(3):dtau_hi(3),ncomp_dtau)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    double precision :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no


   integer :: i,j,k

   dtau(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO


   do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

           dtau(i,j,k,1) = 0.0

         enddo
      enddo
    enddo



end subroutine ca_dertau_r
