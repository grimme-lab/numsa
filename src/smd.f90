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

!> Public API of the numsa library
module smd
   use smd_init, only: init_smd, smd_param
   use smd_sigma, only: calc_surft, smd_surft
   use smd_cds, only: calc_cds
   implicit none
   public

end module smd
