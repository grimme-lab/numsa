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

!> Van-der-Waals radii data
module numsa_data
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use mctc_io_convert, only : aatoau
   implicit none
   private

   public :: get_vdw_rad_d3, get_vdw_rad_cosmo, get_vdw_rad_bondi


  !> In case no van-der-Waals value is provided
  real(wp), parameter :: missing = -1.0_wp

   !> Get van-der-Waals radius for a species
   interface get_vdw_rad_d3
      module procedure :: get_vdw_rad_d3_symbol
      module procedure :: get_vdw_rad_d3_number
   end interface get_vdw_rad_d3

   !> D3 pairwise van-der-Waals radii (only homoatomic pairs present here)
   real(wp), parameter :: vdw_rad_d3(94) = aatoau * [&
      & 1.09155_wp, 0.86735_wp, 1.74780_wp, 1.54910_wp, &  ! H-Be
      & 1.60800_wp, 1.45515_wp, 1.31125_wp, 1.24085_wp, &  ! B-O
      & 1.14980_wp, 1.06870_wp, 1.85410_wp, 1.74195_wp, &  ! F-Mg
      & 2.00530_wp, 1.89585_wp, 1.75085_wp, 1.65535_wp, &  ! Al-S
      & 1.55230_wp, 1.45740_wp, 2.12055_wp, 2.05175_wp, &  ! Cl-Ca
      & 1.94515_wp, 1.88210_wp, 1.86055_wp, 1.72070_wp, &  ! Sc-Cr
      & 1.77310_wp, 1.72105_wp, 1.71635_wp, 1.67310_wp, &  ! Mn-Ni
      & 1.65040_wp, 1.61545_wp, 1.97895_wp, 1.93095_wp, &  ! Cu-Ge
      & 1.83125_wp, 1.76340_wp, 1.68310_wp, 1.60480_wp, &  ! As-Kr
      & 2.30880_wp, 2.23820_wp, 2.10980_wp, 2.02985_wp, &  ! Rb-Zr
      & 1.92980_wp, 1.87715_wp, 1.78450_wp, 1.73115_wp, &  ! Nb-Ru
      & 1.69875_wp, 1.67625_wp, 1.66540_wp, 1.73100_wp, &  ! Rh-Cd
      & 2.13115_wp, 2.09370_wp, 2.00750_wp, 1.94505_wp, &  ! In-Te
      & 1.86900_wp, 1.79445_wp, 2.52835_wp, 2.59070_wp, &  ! I-Ba
      & 2.31305_wp, 2.31005_wp, 2.28510_wp, 2.26355_wp, &  ! La-Nd
      & 2.24480_wp, 2.22575_wp, 2.21170_wp, 2.06215_wp, &  ! Pm-Gd
      & 2.12135_wp, 2.07705_wp, 2.13970_wp, 2.12250_wp, &  ! Tb-Er
      & 2.11040_wp, 2.09930_wp, 2.00650_wp, 2.12250_wp, &  ! Tm-Hf
      & 2.04900_wp, 1.99275_wp, 1.94775_wp, 1.87450_wp, &  ! Ta-Os
      & 1.72280_wp, 1.67625_wp, 1.62820_wp, 1.67995_wp, &  ! Ir-Hg
      & 2.15635_wp, 2.13820_wp, 2.05875_wp, 2.00270_wp, &  ! Tl-Po
      & 1.93220_wp, 1.86080_wp, 2.53980_wp, 2.46470_wp, &  ! At-Ra
      & 2.35215_wp, 2.21260_wp, 2.22970_wp, 2.19785_wp, &  ! Ac-U
      & 2.17695_wp, 2.21705_wp]                            ! Np-Pu

  !> Get van-der-Waals radius for a species
  interface get_vdw_rad_cosmo
    module procedure :: get_vdw_rad_cosmo_symbol
    module procedure :: get_vdw_rad_cosmo_number
  end interface get_vdw_rad_cosmo


   !> Default value for unoptimized van-der-Waals radii
   real(wp), parameter :: cosmostub = 2.223_wp

   !> COSMO optimized van-der-Waals radii
   real(wp), parameter :: vdw_rad_cosmo(94) = aatoau * [ &
       & 1.3000_wp, 1.6380_wp, 1.5700_wp, 1.0530_wp, &   ! h-be
       & 2.0480_wp, 2.0000_wp, 1.8300_wp, 1.7200_wp, &   ! B-O
       & 1.7200_wp, 1.8018_wp, 1.8000_wp, 1.6380_wp, &   ! F-Mg
       & 2.1530_wp, 2.2000_wp, 2.1060_wp, 2.1600_wp, &   ! Al-S
       & 2.0500_wp, 2.2000_wp, 2.2230_wp, cosmostub, &   ! Cl-Ca
       & cosmostub, 2.2930_wp, cosmostub, cosmostub, &   ! Sc-Cr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Mn-Ni
       & cosmostub, 1.6260_wp, cosmostub, 2.7000_wp, &   ! Cu-Ge
       & 2.3500_wp, 2.2000_wp, 2.1600_wp, 2.3630_wp, &   ! As-Kr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rb-Zr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Nb-Ru
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rh-Cd
       & 2.2580_wp, 2.5500_wp, 2.4100_wp, 2.4100_wp, &   ! In-Te
       & 2.3200_wp, 2.5270_wp, cosmostub, cosmostub, &   ! I-Ba
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! La-Nd
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Pm-Gd
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tb-Er
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tm-Hf
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ta-Os
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ir-Hg
       & cosmostub, 2.3600_wp, 2.4220_wp, 2.3050_wp, &   ! Tl-Po
       & 2.3630_wp, 2.5740_wp, cosmostub, cosmostub, &   ! At-Ra
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ac-U
       & cosmostub, cosmostub]                           ! Np-Pu


  !> Get van-der-Waals radius for a species
  interface get_vdw_rad_bondi
    module procedure :: get_vdw_rad_bondi_Symbol
    module procedure :: get_vdw_rad_bondi_Number
  end interface get_vdw_rad_bondi

  !> Van-der-Waals radii from
  !> Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J. Cramer,
  !> and Donald G. Truhlar, Consistent van der Waals Radii for the Whole Main Group,
  !> J. Phys. Chem. A 2009, 113, 19, 5806–5812. https://doi.org/10.1021/jp8111556
  real(wp), parameter :: vdw_rad_Bondi(88) = aatoau * [ &
      & 1.10_wp, 1.40_wp, 1.81_wp, 1.53_wp, 1.92_wp, 1.70_wp, 1.55_wp, 1.52_wp, &  ! H-O
      & 1.47_wp, 1.54_wp, 2.27_wp, 1.73_wp, 1.84_wp, 2.10_wp, 1.80_wp, 1.80_wp, &  ! F-S
      & 1.75_wp, 1.88_wp, 2.75_wp, 2.31_wp, missing, missing, missing, missing, &  ! Cl-Cr
      & missing, missing, missing, missing, missing, missing, 1.87_wp, 2.11_wp, &  ! Mn-Ge
      & 1.85_wp, 1.90_wp, 1.83_wp, 2.02_wp, 3.03_wp, 2.49_wp, missing, missing, &  ! As-Zr
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Nb-Cd
      & 1.93_wp, 2.17_wp, 2.06_wp, 2.06_wp, 1.98_wp, 2.16_wp, 3.43_wp, 2.68_wp, &  ! I-Ba
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! La-Gd
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Tb-Hf
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Ta-Hg
      & 1.96_wp, 2.02_wp, 2.07_wp, 1.97_wp, 2.02_wp, 2.20_wp, 3.48_wp, 2.83_wp]    ! Tl-Ra

contains

!> Get van-der-Waals radius for species with a given symbol
elemental function get_vdw_rad_d3_symbol(symbol) result(radius)
   !> Element symbol
   character(len=*), intent(in) :: symbol
   !> van-der-Waals radius
   real(wp) :: radius

   radius = get_vdw_rad_d3(to_number(symbol))

end function get_vdw_rad_d3_symbol

!> Get van-der-Waals radius for species with a given atomic number
elemental function get_vdw_rad_d3_number(number) result(radius)
   !> Atomic number
   integer, intent(in) :: number
   !> van-der-Waals radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(vdw_rad_d3, dim=1)) then
      radius = vdw_rad_d3(number)
   else
      radius = missing
   end if

end function get_vdw_rad_d3_number

!> Get van-der-Waals radius for species with a given symbol
elemental function get_vdw_rad_cosmo_symbol(symbol) result(radius)
   !> Element symbol
   character(len=*), intent(in) :: symbol
   !> van-der-Waals radius
   real(wp) :: radius

   radius = get_vdw_rad_cosmo(to_number(symbol))

end function get_vdw_rad_cosmo_symbol

!> Get van-der-Waals radius for species with a given atomic number
elemental function get_vdw_rad_cosmo_number(number) result(radius)
   !> Atomic number
   integer, intent(in) :: number
   !> van-der-Waals radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(vdw_rad_cosmo, dim=1)) then
      radius = vdw_rad_cosmo(number)
   else
      radius = missing
   end if

end function get_vdw_rad_cosmo_number

!> Get van-der-Waals radius for species with a given symbol
elemental function get_vdw_rad_bondi_symbol(symbol) result(radius)
   !> Element symbol
   character(len=*), intent(in) :: symbol
   !> van-der-Waals radius
   real(wp) :: radius

   radius = get_vdw_rad_bondi(to_number(symbol))

end function get_vdw_rad_bondi_symbol

!> Get van-der-Waals radius for species with a given atomic number
elemental function get_vdw_rad_bondi_number(number) result(radius)
   !> Atomic number
   integer, intent(in) :: number
   !> van-der-Waals radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(vdw_rad_bondi, dim=1)) then
      radius = vdw_rad_bondi(number)
   else
      radius = missing
   end if

end function get_vdw_rad_bondi_number

end module numsa_data
