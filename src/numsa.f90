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
module numsa
   use numsa_data, only : get_vdw_rad_d3, get_vdw_rad_cosmo, get_vdw_rad_bondi&
     &, get_vdw_rad_smd
   use numsa_lebedev, only : get_angular_grid, grid_size
   use numsa_surface, only : surface_integrator, new_surface_integrator
   use numsa_version, only : get_numsa_version
   implicit none
   public

end module numsa
