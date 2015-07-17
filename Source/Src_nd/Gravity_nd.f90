
! ::
! :: ----------------------------------------------------------
! ::

      ! Returns the gravitational constant, G

      subroutine get_grav_const(Gconst_out)

         use fundamental_constants_module, only: Gconst

         double precision :: Gconst_out

         Gconst_out = Gconst

      end subroutine get_grav_const

! ::
! :: ----------------------------------------------------------
! ::

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Given a radial mass distribution, this computes the gravitational
      ! acceleration as a function of radius by computing the mass enclosed
      ! in successive spherical shells.
      ! Inputs: mass(r), dr, numpts_1d
      ! Outputs: grav(r)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_integrate_grav (mass,den,grav,max_radius,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst
      use bl_constants_module

      implicit none
      integer          :: numpts_1d
      double precision :: mass(0:numpts_1d-1)
      double precision ::  den(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: max_radius,dr

      integer          :: i
      double precision :: rc,rlo,rhi,halfdr
      double precision :: mass_encl
      double precision :: vol_inner_shell, vol_outer_shell
      double precision :: vol_lower_shell, vol_upper_shell
      double precision :: vol_total_im1, vol_total_i

      double precision, parameter ::  fourthirdspi = 4.d0 * M_PI / 3.d0

      halfdr = 0.5d0 * dr

      mass_encl = 0.d0
      do i = 0,numpts_1d-1
         rlo = (dble(i)      ) * dr
         rc  = (dble(i)+0.5d0) * dr
         rhi = (dble(i)+1.0d0) * dr

         if (i.eq.0) then

            ! The mass at (i) is distributed into these two regions
            vol_outer_shell = fourthirdspi * rc**3 
            vol_upper_shell = fourthirdspi * (rhi**3 - rc**3)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = vol_outer_shell * mass(i)  / vol_total_i

         else if (rc .lt. max_radius) then

            ! The mass at (i-1) is distributed into these two shells
            vol_lower_shell = vol_outer_shell   ! This copies from the previous i
            vol_inner_shell = vol_upper_shell   ! This copies from the previous i
            vol_total_im1   = vol_total_i       ! This copies from the previous i

            ! The mass at (i)   is distributed into these two shells
            vol_outer_shell = fourthirdspi * halfdr * ( rc**2 + rlo*rc + rlo**2)
            vol_upper_shell = fourthirdspi * halfdr * ( rc**2 + rhi*rc + rhi**2)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = mass_encl + (vol_inner_shell/vol_total_im1) * mass(i-1) + & 
                                    (vol_outer_shell/vol_total_i  ) * mass(i  ) 

         else 

            ! The mass at (i-1) is distributed into these two shells
            vol_lower_shell = vol_outer_shell   ! This copies from the previous i
            vol_inner_shell = vol_upper_shell   ! This copies from the previous i
            vol_total_im1   = vol_total_i       ! This copies from the previous i

            ! The mass at (i)   is distributed into these two shells
            vol_outer_shell = fourthirdspi * halfdr * ( rc**2 + rlo*rc + rlo**2)
            vol_upper_shell = fourthirdspi * halfdr * ( rc**2 + rhi*rc + rhi**2)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = mass_encl + vol_inner_shell*den(i-1) + vol_outer_shell*den(i  )
         end if

         grav(i) = -Gconst * mass_encl / rc**2
!        print *,'GRAV MASS ',rc, mass_encl

      enddo

      end subroutine ca_integrate_grav

! ::
! :: ----------------------------------------------------------
! ::

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Integrates radial mass elements of a spherically symmetric          
      ! mass distribution to calculate both the gravitational acceleration,  
      ! grav, and the gravitational potential, phi. Here the mass variable   
      ! gives the mass contained in each radial shell.                      
      !                                                                     
      ! The convention in Castro for Poisson's equation is                  
      !                                                                     
      !     laplacian(phi) = -4*pi*G*rho                                    
      !
      ! The gravitational acceleration is then
      !
      !     g(r) = -G*M(r) / r**2
      !
      ! where M(r) is the mass interior to radius r.
      !
      ! The strategy for calculating the potential is to calculate the potential
      ! at the boundary assuming all the mass is enclosed:
      !
      !     phi(R) = G * M / R 
      !
      ! Then, the potential in all other zones can be found using
      !
      !     d(phi)/dr = g    ==>    phi(r < R) = phi(R) - int(g * dr)
      !
      ! Inputs: mass, grav, dr, numpts_1d
      ! Outputs: phi
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_integrate_phi (mass,grav,phi,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst

      implicit none
      integer          :: numpts_1d
      double precision :: mass(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: phi(0:numpts_1d-1)
      double precision :: dr
      double precision :: gravBC, phiBC

      integer          :: i
      double precision :: mass_encl,rhi

      mass_encl = 0.d0
      grav(0)   = 0.d0
      do i = 1,numpts_1d-1
         rhi = dble(i) * dr
         mass_encl = mass_encl + mass(i-1)
         grav(i) = -Gconst * mass_encl / rhi**2
      enddo

      mass_encl = mass_encl + mass(numpts_1d-1)
      phiBC = Gconst * mass_encl / (numpts_1d*dr)
      gravBC = -Gconst * mass_encl / (numpts_1d*dr)**2
      phi(numpts_1d-1) = phiBC - gravBC * dr
       
      do i = numpts_1d-2,0,-1
        phi(i) = phi(i+1) - grav(i+1) * dr
      enddo
      
      end subroutine ca_integrate_phi

! ::
! :: ----------------------------------------------------------
! ::

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Same as ca_integrate_grav above, but includes general relativistic effects.
      ! Inputs: rho, mass, pressure, dr, numpts_1d
      ! Outputs: grav
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ca_integrate_gr_grav (rho,mass,pres,grav,dr,numpts_1d)

      use fundamental_constants_module, only : Gconst
      use bl_constants_module

      implicit none
      integer          :: numpts_1d
      double precision ::  rho(0:numpts_1d-1)
      double precision :: mass(0:numpts_1d-1)
      double precision :: pres(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: dr

      integer          :: i
      double precision :: mass_encl,rc,rlo,rhi,halfdr
      double precision :: ga, gb, gc, P,R
      double precision :: vol_inner_shell, vol_outer_shell
      double precision :: vol_lower_shell, vol_upper_shell
      double precision :: vol_total_im1, vol_total_i

      double precision, parameter ::  fourpi       = 4.d0 * M_PI
      double precision, parameter ::  fourthirdspi = 4.d0 * M_PI / 3.d0
      double precision, parameter ::  sqvc         = 29979245800.d0**2

      halfdr = 0.5d0 * dr

      mass_encl = 0.d0
      do i = 0,numpts_1d-1
         rlo = (dble(i)      ) * dr
         rc  = (dble(i)+0.5d0) * dr
         rhi = (dble(i)+1.0d0) * dr

         if (i.eq.0) then

            ! The mass at (i) is distributed into these two regions
            vol_outer_shell = fourthirdspi * rc**3 
            vol_upper_shell = fourthirdspi * (rhi**3 - rc**3)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = vol_outer_shell * mass(i)  / vol_total_i

         else

            ! The mass at (i-1) is distributed into these two shells
            vol_lower_shell = vol_outer_shell   ! This copies from the previous i
            vol_inner_shell = vol_upper_shell   ! This copies from the previous i
            vol_total_im1   = vol_total_i       ! This copies from the previous i

            ! The mass at (i)   is distributed into these two shells
            vol_outer_shell = fourthirdspi * halfdr * ( rc**2 + rlo*rc + rlo**2)
            vol_upper_shell = fourthirdspi * halfdr * ( rc**2 + rhi*rc + rhi**2)
            vol_total_i     = vol_outer_shell + vol_upper_shell

            mass_encl = mass_encl + vol_inner_shell / vol_total_im1 * mass(i-1) + & 
                                    vol_outer_shell / vol_total_i   * mass(i  ) 
         end if
         grav(i) = -Gconst * mass_encl / rc**2

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This adds the post-Newtonian correction
!!       Added by Ken Chen, 6/9 2010
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!       Tolman-Oppenheimer-Volkoff(TOV) case

         if (rho(i) .gt. 0.d0) then
            P =  pres(i)
            R =  rho(i)
            ga = (1.d0 + P/(R*sqvc))
            gb = (1.d0 + fourpi * rc**3 * P / (mass_encl*sqvc))
            gc = 1.d0 / (1.d0 - 2.d0 * Gconst * mass_encl / (rc*sqvc))

            grav(i) = grav(i)*ga*gb*gc
         end if

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This ends the post-Newtonian correction
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

      end subroutine ca_integrate_gr_grav


! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_const_grav_phi (lo,hi,phi,p_lo,p_hi,dx,const_grav)

      use bl_constants_module
      use prob_params_module, only: problo, k3d, j2d
        
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: p_lo(3), p_hi(3)
      double precision :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
      double precision :: dx(3)
      double precision :: const_grav(3)
      
      integer          :: i, j, k
      double precision :: loc(3)

      integer :: ng = 1 ! We only need one ghost cell for phi

      ! The zero-point for the potential doesn't matter for
      ! constant gravity, so we'll just arbitrarily set it
      ! to zero out at the low end of the domain on the relevant axis.

      do k = lo(3)-ng*k3d, hi(3)+ng*k3d
         loc(3) = problo(3) + (k + HALF) * dx(3)
         do j = lo(2)-ng*j2d, hi(2)+ng*j2d
            loc(2) = problo(2) + (j + HALF) * dx(2)
            do i = lo(1)-ng, hi(1)+ng
               loc(1) = problo(1) + (i + HALF) * dx(1)

               phi(i,j,k) = - sum(loc * const_grav)
            enddo
         enddo
      enddo
      
      end subroutine ca_const_grav_phi


      
! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_ec_grad_phi (lo,hi,phi,p_lo,p_hi,&
                                 gphi,g_lo,g_hi,dx,idir)

      use bl_constants_module
      use prob_params_module, only: problo, k3d, j2d, dim
        
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: p_lo(3), p_hi(3)
      integer          :: g_lo(3), g_hi(3)
      double precision ::  phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
      double precision :: gphi(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3))
      double precision :: dx(3)
      integer          :: idir
      
      integer          :: i, j, k

      do k = lo(3), hi(3)+1*k3d
         do j = lo(2), hi(2)+1*j2d
            do i = lo(1), hi(1)+1

               if (idir .eq. 0) then
                  gphi(i,j,k) = -(phi(i,j,k) - phi(i-1,j,k)) / dx(1)
               else if (idir .eq. 1 .and. dim .ge. 2) then
                  gphi(i,j,k) = -(phi(i,j,k) - phi(i,j-1,k)) / dx(2)
               else if (idir .eq. 2 .and. dim .eq. 3) then
                  gphi(i,j,k) = -(phi(i,j,k) - phi(i,j,k-1)) / dx(3)
               endif
               
            enddo
         enddo
      enddo
      
      end subroutine ca_ec_grad_phi
      
