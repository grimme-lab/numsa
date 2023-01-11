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
      & 1.57762308784207E+1_wp, 7.44023828473916E+0_wp, 5.78326430494107E+0_wp, &
      & 2.96273853143212E+0_wp, 7.96228620441620E+0_wp, 6.94475617383051E+0_wp, &
      & 1.39709036162651E+0_wp, 4.61011426631971E+0_wp, 7.81213684924535E-4_wp, &
      & 8.37602534092320E+0_wp, 5.27091918116814E+0_wp, 1.15343190424705E+1_wp, &
      & 2.65268998376027E+0_wp, 4.61347672833828E+0_wp, 5.32574605544628E-1_wp, &
      & 3.87132219837450E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface/(4*pi)
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
      & 2.27659103959728E+0_wp, 5.97577118950122E+0_wp, 6.41298824748133E+0_wp, &
      & 6.55734631230962E+0_wp, 5.15770822160020E+0_wp, 1.57234572258234E+0_wp, &
      & 3.90432768274610E+0_wp, 4.21140479966022E+0_wp, 7.27814781080846E+0_wp, &
      & 1.10051545734758E+0_wp, 7.17814855288648E+0_wp, 9.04904553399203E+0_wp, &
      & 7.82899298504117E+0_wp, 4.74222890087015E+0_wp, 2.36038439880367E-1_wp, &
      & 1.15287664294478E+1_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface/(4*pi)
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
      & 3.92672956420918E+0_wp, 4.31618535417720E+0_wp, 2.05344888225025E+0_wp, &
      & 2.60133027199439E+0_wp, 1.01849623754110E+0_wp, 7.52652188411287E+0_wp, &
      & 2.73374452715462E+0_wp, 2.19905511071582E+0_wp, 2.18761463996078E+0_wp, &
      & 2.27442772580975E+0_wp, 6.36073079991843E+0_wp, 1.00473063578636E+1_wp, &
      & 4.28139986348257E+0_wp, 3.31270864378470E+0_wp, 7.92536426567952E+0_wp, &
      & 1.87822502816429E+0_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_cosmo(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   if (any(abs(surface - ref) > thr2)) then
      call test_failed(error, "Surface area values do not match")
      print '(es20.14e1)', surface/(4*pi)
   end if

end subroutine test_mb03

end module test_surface
