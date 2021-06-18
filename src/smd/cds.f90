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

!> Routines for calculating the CDS Part of the Solvation Free Energy
module smd_cds
   use smd_sigma, only: smd_surft
   use mctc_env, only: wp
   use mctc_io_symbols, only: to_number
   use mctc_io_convert, only: autoaa
   implicit none
   private
   
   public :: calc_cds

   ! interface calc_cds
   !    module procedure :: calc_cds_symbol
   !    module procedure :: calc_cds_number
   ! end interface

contains

   subroutine calc_cds(surft,surface,cds,cds_sm)
      !>smd Surface Tensions per Atom and for the Solvent
      type(smd_surft), intent(in) :: surft
      !>SASA Surface per Atom
      real(wp), intent(in) :: surface(:)
      !>CDS Part of the Energy per Atom
      real(wp),allocatable, intent(out) :: cds(:)
      !> Solvent Contribution to CDS
      real(wp), intent(out) :: cds_sm
      
      !> Laufvariable
      integer :: i
      !> Number of Atoms
      integer :: nat

      nat=size(surface)
      allocate(cds(nat))

      cds=0.0_wp
      cds_sm=surft%sm*sum(surface)
      do i=1,nat
         cds(i)=surft%sk(i)*surface(i)
      end do

   end subroutine calc_cds
         
   ! subroutine calc_cds_number(surft,surface,numbers,cds)
   !    !>smd Surface Tensions per Atom and for the Solvent
   !    type(smd_surft) :: surft
   !    !>SASA Surface per Atom
   !    real(wp) :: surface(:)
   !    !>Identifiers per Atom as Numbers
   !    integer :: numbers(:)
   !    !>CDS Part of the Energy
   !    real(wp) :: cds

   !    !>CDS Part of the Energy per Atom
   !    real(wp) :: atom_cds(size(surface))
   !    !> Laufvariable
   !    integer :: i

   !    cds=0.0_wp
   !    atom_cds=0.0_wp

   !    do concurrent (i=1:size(surface))
   !       atom_cds(i)=surft%sk(numbers(i))*surface(i)*autoaa**2
   !    end do

   !    cds=sum(atom_cds)+surft%sm*sum(surface)
   ! end subroutine calc_cds_number


end module smd_cds
