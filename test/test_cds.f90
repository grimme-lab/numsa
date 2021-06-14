
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

module test_cds
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   use mstore, only : get_structure
   use numsa
   use smd
   implicit none
   private

   public :: collect_cds

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_cds(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("test-1", test_mb01), &
      new_unittest("test-2", test_mb02), &
      new_unittest("test-3", test_mb03) &
      ]

end subroutine collect_cds


subroutine test_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(smd_param) :: param
   type(smd_surft) :: surft
   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), cds(:)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110
   real(wp), parameter :: ref_h2o(16) = [&
      & 0.00000000000000E+0_wp, 1.01444613143073E+2_wp, 0.00000000000000E+0_wp, &
      & 4.03957436394365E+1_wp, 8.51287137325784E+1_wp, 9.46889274612189E+1_wp, &
      & 1.90487681686013E+1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 1.14203681733362E+2_wp, 7.18668307358171E+1_wp, 3.17179991434577E+1_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp] 
   real(wp), parameter :: ref_methanol(16) = [&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.32489573687157E+1_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 5.04187432438299E+1_wp, 9.50094225504006E-3_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.04541278420455E+2_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.47702335383448E+0_wp, &
      & 0.00000000000000E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   Call t_cds(mol,surface,"h2o",cds)

   if (any(abs(cds - ref_h2o) > thr2)) then
      call test_failed(error, "CDS Part of the Energy for H2O does not match")
      print '(es20.14e1)', cds
   end if
   deallocate(cds)

   Call t_cds(mol,surface,"methanol",cds)

   if (any(abs(cds - ref_methanol) > thr2)) then
      call test_failed(error, "CDS Part of the Energy for DMSO does not match")
      print '(es20.14e1)', cds
   end if
   deallocate(cds)

end subroutine test_mb01

subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(smd_param) :: param
   type(smd_surft) :: surft

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), cds(:)
   real(wp), parameter :: probe = 1.2_wp * aatoau
   integer, parameter :: nang = 230
   real(wp), parameter :: ref_h2o(16) = [& 
      & 3.10404027273529E+1_wp,-1.52278198113227E+1_wp, 0.00000000000000E+0_wp, & 
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.14383016646053E+1_wp, & 
      & 5.32339310986410E+1_wp, 5.74208026082603E+1_wp, 0.00000000000000E+0_wp, & 
      & 1.50050920981926E+1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 8.37036100582293E+1_wp, 6.46583724278679E+1_wp, 3.21829538686014E+0_wp, &
      &-2.93782932723676E+1_wp] 
   real(wp), parameter :: ref_dmso(16) = [&
      & 0.00000000000000E+0_wp,-8.19938175886519E+1_wp, 0.00000000000000E+0_wp, &
      &-1.18391529570206E+2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-5.43123934832201E+1_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      &-1.58186690510307E+2_wp] 
   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   Call t_cds(mol,surface,"h2o",cds)

   if (any(abs(cds - ref_h2o) > thr2)) then
      call test_failed(error, "CDS Part of the Energy for H2O does not match")
      print '(es20.14e1)', cds
   end if
   deallocate(cds)

   Call t_cds(mol,surface,"dmso",cds)

   if (any(abs(cds - ref_dmso) > thr2)) then
      call test_failed(error, "CDS Part of the Energy for DMSO does not match")
      print '(es20.14e1)', cds
   end if
   deallocate(cds)
end subroutine test_mb02


subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(smd_param) :: param
   type(smd_surft) :: surft
   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), cds(:)
   real(wp), parameter :: probe = 0.2_wp * aatoau
   integer, parameter :: nang = 111
   real(wp), parameter :: ref_h2o(16) = [& 
      & 1.77491382598271E+1_wp,-1.08611521278540E+1_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.77106099318681E+1_wp, &
      & 2.29077445894260E+1_wp, 2.27586923792508E+1_wp, 0.00000000000000E+0_wp, &
      & 1.71835788013756E+1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 4.64829064939494E+1_wp, 2.64570923166574E+1_wp, 4.60570209042278E+0_wp, &
      &-1.81487909333211E+1_wp]
   real(wp), parameter :: ref_dmso(16) = [&
      & 0.00000000000000E+0_wp,-5.31330629474944E+1_wp, 0.00000000000000E+0_wp, &
      &-2.50614157045337E+1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.98591793242355E+1_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      &-8.87843977995727E+1_wp]
      
   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   Call t_cds(mol,surface,"h2o",cds)

   if (any(abs(cds - ref_h2o) > thr2)) then
      call test_failed(error, "CDS Part of the Energy for H2O does not match")
      print '(es20.14e1)', cds
   end if
   deallocate(cds)


   Call t_cds(mol,surface,"acetonitrile",cds)

   if (any(abs(cds - ref_dmso) > thr2)) then
      call test_failed(error, "CDS Part of the Energy for Acetonitrile does not match")
      print '(es20.14e1)', cds
   end if

   deallocate(cds)

end subroutine test_mb03

subroutine t_cds(mol,surface,solvent,cds)

   type(smd_param) :: param
   type(structure_type) :: mol
   type(smd_surft) :: surft
   character(len=*) :: solvent
   real(wp), allocatable :: cds(:), surface(:)

   Call init_smd(param,solvent,"test/")
   Call calc_surft(mol%xyz,mol%id,mol%sym,param,surft)
   Call calc_cds(surft,surface,mol%sym,mol%id,cds)

end subroutine t_cds

end module test_cds