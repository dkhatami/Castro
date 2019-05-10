subroutine ca_dercd_mask(mask,mask_lo,mask_hi,ncomp_mask, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo,hi,domlo,domhi, &
                         dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_dercd_mask')


    use amrex_constants_module

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

           if(u(i,j,k,6) > ZERO) then

                 if(u(i,j,k,7)/u(i,j,k,1) > HALF) then
                   mask(i,j,k,1) = ONE
             endif
           endif
         enddo
      enddo
    enddo



end subroutine ca_dercd_mask
