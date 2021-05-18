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

!> Auxilary function for searching a list of integers
module numsa_search
   implicit none
   private

   public :: list_bisection

contains

!> Integer case for bisection search
pure function list_bisection(list, val) result(pos)
   !> Array of values in monotonic order to search through
   integer, intent(in) :: list(:)
   !> Value to locate pos for
   integer, intent(in) :: val
   !> Located element such that list(pos) <= val < list(pos+1)
   integer :: pos

   integer :: n, lower, upper, current

   n = size(list)
   if (n == 0) then
      pos = 0
      return
   end if

   if (val < list(1)) then
      pos = 0
   else if (val == list(1)) then
      pos = 1
   else if (val == list(n)) then
      pos = n - 1
   else if (val > list(n)) then
      pos = n
   else
      lower = 0
      current = n+1
      do while ((current - lower) > 1)
         upper = (current + lower)/2
         if ((list(n) >= list(1)).eqv.(val >= list(upper))) then
            lower = upper
         else
            current = upper
         end if
      end do
      pos = lower
   end if

end function list_bisection

end module numsa_search
