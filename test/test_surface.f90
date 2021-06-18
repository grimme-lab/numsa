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

module test_surface
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use mstore, only : get_structure
   use numsa
   implicit none
   private

   public :: collect_surface

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_surface(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("test-1", test_mb01), &
      new_unittest("test-2", test_mb02), &
      new_unittest("test-3", test_mb03) &
      ]

end subroutine collect_surface


subroutine test_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110
   real(wp), parameter :: ref(16) = 4*pi * [&
      & 1.57762257884252E+1_wp, 7.44023682724103E+0_wp, 5.78326352983607E+0_wp, &
      & 2.96273887964889E+0_wp, 7.96228449837119E+0_wp, 6.94475560013532E+0_wp, &
      & 1.39709090557345E+0_wp, 4.61011360744476E+0_wp, 7.81217049844597E-4_wp, &
      & 8.37602325360240E+0_wp, 5.27091803232173E+0_wp, 1.15343158333901E+1_wp, &
      & 2.65268867545965E+0_wp, 4.61347534035499E+0_wp, 5.32574658431686E-1_wp, &
      & 3.87132091420996E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface
   end if

end subroutine test_mb01


subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 1.2_wp * aatoau
   integer, parameter :: nang = 230
   real(wp), parameter :: ref(16) = 4*pi * [&
      & 2.27659153452265E+0_wp, 5.97577009756815E+0_wp, 6.41298531564733E+0_wp, &
      & 6.55734382952262E+0_wp, 5.15770714137722E+0_wp, 1.57234609721010E+0_wp, &
      & 3.90432810917537E+0_wp, 4.21140518928471E+0_wp, 7.27814416367210E+0_wp, &
      & 1.10051618677536E+0_wp, 7.17814405439378E+0_wp, 9.04904257798528E+0_wp, &
      & 7.82899127218094E+0_wp, 4.74222917138838E+0_wp, 2.36038948903930E-1_wp, &
      & 1.15287630553696E+1_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface
   end if

end subroutine test_mb02


subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp), parameter :: probe = 0.2_wp * aatoau
   integer, parameter :: nang = 111
   real(wp), parameter :: ref(16) = 4*pi * [&
      & 3.92672838745645E+0_wp, 4.31618518252200E+0_wp, 2.05344924554658E+0_wp, &
      & 2.60133065338031E+0_wp, 1.01849685835522E+0_wp, 7.52651997890806E+0_wp, &
      & 2.73374455863003E+0_wp, 2.19905522537305E+0_wp, 2.18761474083904E+0_wp, &
      & 2.27442777032111E+0_wp, 6.36072952076570E+0_wp, 1.00473034014832E+1_wp, &
      & 4.28139813443965E+0_wp, 3.31270794097678E+0_wp, 7.92536272758048E+0_wp, &
      & 1.87822532886867E+0_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_cosmo(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface
   end if

end subroutine test_mb03

end module test_surface
