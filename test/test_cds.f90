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
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp) :: cds
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110
   real(wp), parameter :: ref_h2o = 7.01825666242813E+3_wp
   real(wp), parameter :: ref_methanol = 3.16981152608071E+3_wp
   call get_structure(mol, "MB16-43", "01")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_d3(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   Call t_cds(mol,surface,"h2o",cds)

   if (abs(cds - ref_h2o) > thr2) then
      call test_failed(error, "CDS Part of the Energy for H2O does not match")
      print '(es20.14e1)', cds
   end if

   Call t_cds(mol,surface,"methanol",cds)

   if (abs(cds - ref_methanol) > thr2) then
      call test_failed(error, "CDS Part of the Energy for methanol does not match")
      print '(es20.14e1)', cds
   end if

end subroutine test_mb01

subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(smd_param) :: param
   type(smd_surft) :: surft

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp) :: cds
   real(wp), parameter :: probe = 1.2_wp * aatoau
   integer, parameter :: nang = 230
   real(wp), parameter :: ref_h2o = 3.58282976254987E+3_wp
   real(wp), parameter :: ref_dmso = 1.27231102159210E+3_wp
   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   Call t_cds(mol,surface,"h2o",cds)

   if (abs(cds - ref_h2o) > thr2) then
      call test_failed(error, "CDS Part of the Energy for H2O does not match")
      print '(es20.14e1)', cds
   end if


   Call t_cds(mol,surface,"dmso",cds)

   if (abs(cds - ref_dmso) > thr2) then
      call test_failed(error, "CDS Part of the Energy for DMSO does not match")
      print '(es20.14e1)', cds
   end if
end subroutine test_mb02


subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(smd_param) :: param
   type(smd_surft) :: surft
   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :)
   real(wp) :: cds
   real(wp), parameter :: probe = 0.2_wp * aatoau
   integer, parameter :: nang = 111
      
   real(wp), parameter :: ref_h2o = 1.84531408914288E+3_wp
   real(wp), parameter :: ref_acetonitrile = 2.30829715019427E+1_wp
   call get_structure(mol, "MB16-43", "02")

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))
   rad = get_vdw_rad_bondi(mol%num)

   call new_surface_integrator(sasa, mol%id, rad, probe, nang)
   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   Call t_cds(mol,surface,"h2o",cds)

   if (abs(cds - ref_h2o) > thr2) then
      call test_failed(error, "CDS Part of the Energy for H2O does not match")
      print '(es20.14e1)', cds
   end if


   Call t_cds(mol,surface,"acetonitrile",cds)

   if (abs(cds - ref_acetonitrile) > thr2) then
      call test_failed(error, "CDS Part of the Energy for Acetonitrile does not match")
      print '(es20.14e1)', cds
   end if


end subroutine test_mb03

subroutine t_cds(mol,surface,solvent,cds_total)

   type(smd_param) :: param
   type(structure_type) :: mol
   type(smd_surft) :: surft
   character(len=*) :: solvent
   real(wp), allocatable :: cds(:), surface(:)
   real(wp) :: cds_sm, cds_total

   Call init_smd(param,solvent,"test/")
   Call calc_surft(mol%xyz,mol%id,mol%sym,param,surft)
   Call calc_cds(surft,surface,cds,cds_sm)

   cds_total=sum(cds)+cds_sm

end subroutine t_cds

end module test_cds
