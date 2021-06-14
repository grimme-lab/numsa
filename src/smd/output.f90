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

!> Routines for displaying results
module smd_output
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   implicit none
   private

   public :: ascii_cds

contains

subroutine ascii_cds(unit, mol, cds,cds_sm)
   !> Unit for output
   integer, intent(in) :: unit
   !> Molecular structure data
   class(structure_type), intent(in) :: mol
   !> Surface area for each atom
   real(wp), intent(in) :: cds(:), cds_sm

   integer :: iat, isp

   write(unit, '(a)')
   write(unit, '(a,":")') "SMD CDS Part in kcal/mol"
   write(unit, '(28("-"))')
   write(unit, '(a6,1x,a4,5x,*(1x,a10))') "#", "Z", "cds"
   write(unit, '(28("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i6,1x,i4,1x,a4,*(1x,f10.4))') &
         & iat, mol%num(isp), mol%sym(isp), cds(iat)/1000.0_wp
   end do
   write(unit, '(28("-"))') 
   write(unit, '(1x, a, t15, f13.4)') &
      "Solvent Contribution:",  (cds_sm)/1000.0_wp
   write(unit, '(1x, a, t15, f13.4)') &
      "total",  (sum(cds)+cds_sm)/1000.0_wp
   write(unit, '(28("-"))')

end subroutine ascii_cds

end module smd_output
