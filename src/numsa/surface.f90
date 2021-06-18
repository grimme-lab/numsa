! This file is part of numsa.
! SPDX-Identifier: LGPL-3.0-or-later
!
! numsa is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! numsa is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with numsa.  If not, see <https://www.gnu.org/licenses/>.

!> Integration of surface area
module numsa_surface
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use numsa_search, only : list_bisection
   use numsa_lebedev, only : get_angular_grid, grid_size
   implicit none
   private

   public :: surface_integrator, new_surface_integrator


   type :: surface_integrator
      !> number of atoms
      integer :: nat
      !> atom types
      integer, allocatable :: at(:)
      !> number of pairs
      integer :: ntpair
      !> number of angular grid points
      integer :: nang 
      !> angular grid
      real(wp), allocatable :: ang_grid(:, :)
      real(wp), allocatable :: ang_weight(:)
      !> cut-off radius for the nearest neighbour list
      real(wp) :: srcut
      !> number of neighbors for surface computation
      integer, allocatable :: nnsas(:)
      !> neighbors of an atom in surface computation
      integer, allocatable :: nnlists(:, :)
      !> all pairs indeces array
      integer, allocatable :: ppind(:, :)
      !> Atom specific surface data
      real(wp), allocatable :: vdwsa(:)
      real(wp), allocatable :: wrp(:)
      real(wp), allocatable :: trj2(:, :)
      real(wp) :: ah0, ah1, ah3
   contains
      procedure :: get_surface
   end type surface_integrator

   !> Smoothing dielectric function parameters
   real(wp), parameter :: w = 0.3_wp*aatoau

   !> real space cut-offs
   real(wp), parameter :: tolsesp=1.e-6_wp


   real(wp), parameter :: default_offset = 2.0_wp * aatoau
   real(wp), parameter :: default_smoothing = 0.3 * aatoau

contains


!> Initialize data straucture
subroutine new_surface_integrator(self, num, rad, probe, nang, offset, smoothing)
   !> Instance of the surface integrator
   type(surface_integrator), intent(out) :: self
   !> Atomic numbers
   integer, intent(in) :: num(:)
   !> Van-der-Waals Radii
   real(wp), intent(in) :: rad(:)
   !> Probe radius of the solvent
   real(wp), intent(in) :: probe
   !> Number of angular grid points for integration
   integer, intent(in) :: nang
   !> Offset for surface integration cutoff
   real(wp), intent(in), optional :: offset
   !> Smooting function parameter
   real(wp), intent(in), optional :: smoothing

   integer :: iat, izp, jat, ij, iang, ierr
   real(wp) :: r, w, w3

   self%nat = size(num)
   self%at = num
   self%ntpair = self%nat*(self%nat-1)/2
   allocate(self%ppind(2, self%ntpair))
   allocate(self%nnsas(self%nat))
   allocate(self%nnlists(self%nat, self%nat))
   ij = 0
   do iat = 1, self%nat
      do jat = 1, iat-1
         ij = ij+1
         self%ppind(1, ij) = iat
         self%ppind(2, ij) = jat
      enddo
   enddo

   allocate(self%vdwsa(self%nat))
   allocate(self%trj2(2, self%nat))
   allocate(self%wrp(self%nat))
   if (present(smoothing)) then
      w = smoothing
   else
      w = default_smoothing
   end if
   w3 = w*(w*w)
   self%ah0 = 0.5_wp
   self%ah1 = 3._wp/(4.0_wp*w)
   self%ah3 = -1._wp/(4.0_wp*w3)
   do iat = 1, self%nat
      izp = num(iat)
      self%vdwsa(iat) = rad(izp) + probe
      self%trj2(1, iat) = (self%vdwsa(iat)-w)**2
      self%trj2(2, iat) = (self%vdwsa(iat)+w)**2
      r=self%vdwsa(iat)+w
      self%wrp(iat)=(0.25_wp/w+ &
         &            3.0_wp*self%ah3*(0.2_wp*r*r-0.5_wp*r*self%vdwsa(iat)+ &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
      r=self%vdwsa(iat)-w
      self%wrp(iat)=self%wrp(iat)-(0.25/w+ &
         &    3.0_wp*self%ah3*(0.2_wp*r*r-0.5_wp*r*self%vdwsa(iat)+ &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
   end do
   self%srcut = 2 * (w + maxval(self%vdwsa))
   if (present(offset)) then
      self%srcut = self%srcut + offset
   else
      self%srcut = self%srcut + default_offset
   end if

   iang = list_bisection(grid_size, nang)
   allocate(self%ang_grid(3, grid_size(iang)))
   allocate(self%ang_weight(grid_size(iang)))
   call get_angular_grid(iang, self%ang_grid, self%ang_weight, ierr)
   self%ang_weight(:) = self%ang_weight * 4*pi

end subroutine new_surface_integrator


subroutine get_surface(self, num, xyz, surface, dsdr)
   !> Instance of the surface integrator
   class(surface_integrator), intent(inout) :: self
   !> Atomic numbers
   integer, intent(in) :: num(:)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Surface area for each sphere
   real(wp), intent(out) :: surface(:)
   !> Derivative of surface area w.r.t. to coordinate displacements
   real(wp), intent(out) :: dsdr(:, :, :)

   ! initialize the neighbor list
   call update_nnlist(self%nat, self%ntpair, self%ppind, xyz, &
      & self%srcut, self%nnsas, self%nnlists)

   ! compute solvent accessible surface and its derivatives
   call compute_numsa(self%nat, self%nnsas, self%nnlists, xyz, self%vdwsa, &
      & self%wrp, self%trj2, self%ah0, self%ah1, self%ah3, self%ang_weight, self%ang_grid, &
      & surface, dsdr)

end subroutine get_surface


subroutine update_nnlist(nat, ntpair, ppind, xyz, srcut, nnsas, nnlists)
!$ use omp_lib

   integer, intent(in) :: nat
   integer, intent(in) :: ntpair
   integer, intent(in) :: ppind(:, :)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: srcut
   integer, intent(out) :: nnsas(:)
   integer, intent(out) :: nnlists(:, :)

   integer kk, i1, i2
   real(wp) rcutn2, srcut2
   real(wp) x, y, z, dr2
   integer ip, ip2, thrid, nproc
   integer, allocatable :: npid(:)
   integer, allocatable :: nntmp(:)
   integer, allocatable :: nnls(:, :)

   nproc=1
!$ nproc=omp_get_max_threads()

   allocate(nnls(nat, nat))
   allocate(nntmp(nat), npid(nproc))
   npid = 0

   srcut2 = srcut*srcut

   nnsas=0
   nnlists=0
!$omp parallel default(none) &
!$omp&         shared ( xyz,srcut2,ntpair,ppind,nat,nnlists,nnsas) &
!$omp&         private( i1,i2,x,y,z,dr2,ip,ip2,thrid,nntmp,nnls ) &
!$omp&         shared ( npid )
   ip=0
   ip2=0
   nntmp=0
   nnls=0
   thrid=1
!$ thrid=omp_get_thread_num() + 1
!$omp do
   do kk=1,ntpair
      i1=ppind(1,kk)
      i2=ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      if(dr2.lt.srcut2) then
         nntmp(i1) = nntmp(i1) + 1
         nntmp(i2) = nntmp(i2) + 1
         nnls(nntmp(i1),i1)=i2
         nnls(nntmp(i2),i2)=i1
      endif
   enddo
!$omp end do
   npid(thrid)=ip
!$omp critical
   do i1=1,nat
      do i2=1,nntmp(i1)
         nnlists(nnsas(i1)+i2,i1)=nnls(i2,i1)
      enddo
      nnsas(i1)=nnsas(i1)+nntmp(i1)
   enddo
!$omp end critical
!$omp end parallel

end subroutine update_nnlist


subroutine compute_numsa(nat, nnsas, nnlists, xyz, vdwsa, &
      & wrp, trj2, ah0, ah1, ah3, ang_weight, ang_grid, surface, dsdrt)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Number of neighbours to consider
   integer, intent(in) :: nnsas(:)
   !> Neighbourlist
   integer, intent(in) :: nnlists(:, :)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Van-der-Waals radii including probe radius of solvent
   real(wp), intent(in) :: vdwsa(:)
   !> Radial weights for the numerical integration
   real(wp), intent(in) :: wrp(:)
   !> Radial smoothing function
   real(wp), intent(in) :: trj2(:, :)
   real(wp),intent(in)  :: ah0,ah1,ah3
   !> Angular weights for the numerical integration
   real(wp), intent(in) :: ang_weight(:)
   !> Angular grid for each atom
   real(wp), intent(in) :: ang_grid(:, :)
   !> Surface area for each atom, including surface tension
   real(wp), intent(out) :: surface(:)
   !> Derivative of surface area w.r.t. cartesian coordinates
   real(wp), intent(out) :: dsdrt(:, :, :)

   integer :: iat, jat, ip, jj, nnj, nnk, nni, nno
   real(wp) :: rsas, sasai, xyza(3), xyzp(3), sasap, wr, wsa, drjj(3)
   real(wp), allocatable :: grds(:, :), grads(:, :)
   integer, allocatable :: grdi(:)

   surface(:) = 0.0_wp
   dsdrt(:, :, :) = 0.0_wp

   ! allocate space for the gradient storage
   allocate(grads(3,nat), source = 0.0_wp)
   allocate(grds(3,maxval(nnsas)))
   allocate(grdi(maxval(nnsas)))

   !$omp parallel do default(none) shared(surface, dsdrt, ah0, ah1, ah3) &
   !$omp shared(nat, vdwsa, nnsas, xyz, wrp, ang_grid, ang_weight, nnlists, trj2) &
   !$omp private(iat, rsas, nno, grads, sasai, xyza, wr, ip, xyzp, wsa, &
   !$omp& sasap, jj, nni, nnj, grdi, grds, drjj)
   do iat = 1, nat

      rsas = vdwsa(iat)
      nno = nnsas(iat)

      ! initialize storage
      grads = 0.0_wp
      sasai = 0.0_wp

      ! atomic position
      xyza(:) = xyz(:,iat)
      ! radial atomic weight
      wr = wrp(iat)

      ! loop over grid points
      do ip = 1, size(ang_grid, 2)
         ! grid point position
         xyzp(:) = xyza(:) + rsas*ang_grid(:,ip)
         ! atomic surface function at the grid point
         call compute_w_sp(nat, nnlists(:nno, iat), trj2, vdwsa, xyz, nno, xyzp, &
            & ah0, ah1, ah3, sasap, grds, nni, grdi)

         if(sasap.gt.tolsesp) then
            ! numerical quadrature weight
            wsa = ang_weight(ip)*wr*sasap
            ! accumulate the surface area
            sasai = sasai + wsa
            ! accumulate the surface gradient
            do jj = 1, nni
               nnj = grdi(jj)
               drjj(:) = wsa*grds(:,jj)
               grads(:,iat) = grads(:,iat)+drjj(:)
               grads(:,nnj) = grads(:,nnj)-drjj(:)
            end do
         end if
      end do

      surface(iat) = sasai
      dsdrt(:,:,iat) = grads

   end do

end subroutine compute_numsa


pure subroutine compute_w_sp(nat,nnlists,trj2,vdwsa,xyza,nno,xyzp,ah0,ah1,ah3, &
      & sasap,grds,nni,grdi)
   integer, intent(in)  :: nat
   integer, intent(in)  :: nnlists(nno)
   integer, intent(in)  :: nno
   integer, intent(out) :: nni
   real(wp),intent(in)  :: xyza(3,nat)
   real(wp),intent(in)  :: xyzp(3)
   real(wp),intent(in)  :: ah0,ah1,ah3
   real(wp),intent(out) :: sasap
   real(wp),intent(out) :: grds(3,nno)
   integer, intent(out) :: grdi(nno)
   real(wp),intent(in)  :: trj2(2,nat)
   real(wp),intent(in)  :: vdwsa(nat)

   integer  :: i,ia
   real(wp) :: tj(3),tj2,sqtj
   real(wp) :: uj,uj3,ah3uj2
   real(wp) :: sasaij,dsasaij

   ! initialize storage
   nni=0
   sasap=1.0_wp
   do i = 1, nno
      ia = nnlists(i)
      ! compute the distance to the atom
      tj(:)=xyzp(:)-xyza(:,ia)
      tj2=dot_product(tj, tj)
      ! if within the outer cut-off compute
      if(tj2.lt.trj2(2,ia)) then
         if(tj2.le.trj2(1,ia)) then
            sasap=0.0_wp
            return
         else
            sqtj = sqrt(tj2)
            uj = sqtj - vdwsa(ia)
            ah3uj2 = ah3*uj*uj
            dsasaij = ah1+3.0_wp*ah3uj2
            sasaij =  ah0+(ah1+ah3uj2)*uj

            ! accumulate the molecular surface
            sasap = sasap*sasaij
            ! compute the gradient wrt the neighbor
            dsasaij = dsasaij/(sasaij*sqtj)
            nni=nni+1
            grdi(nni)=ia
            grds(:,nni) = dsasaij*tj(:)
         endif
         ! check if the point is already buried
!        if(sasap.lt.tolsesp) return
      endif
   enddo

end subroutine compute_w_sp

end module numsa_surface
