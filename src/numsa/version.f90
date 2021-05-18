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

!> Version information for the library
module numsa_version
   implicit none
   private

   public :: numsa_version_string, numsa_version_compact
   public :: get_numsa_version


   !> String representation of the numsa version
   character(len=*), parameter :: numsa_version_string = "0.1.0"

   !> Numeric representation of the numsa version
   integer, parameter :: numsa_version_compact(3) = [0, 1, 0]


contains


!> Getter function to retrieve numsa version
subroutine get_numsa_version(major, minor, patch, string)
   !> Major version number of the numsa version
   integer, intent(out), optional :: major
   !> Minor version number of the numsa version
   integer, intent(out), optional :: minor
   !> Patch version number of the numsa version
   integer, intent(out), optional :: patch
   !> String representation of the numsa version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = numsa_version_compact(1)
   end if
   if (present(minor)) then
      minor = numsa_version_compact(2)
   end if
   if (present(patch)) then
      patch = numsa_version_compact(3)
   end if
   if (present(string)) then
      string = numsa_version_string
   end if

end subroutine get_numsa_version


end module numsa_version
