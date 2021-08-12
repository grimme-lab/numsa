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

program main_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, write_structure, &
      & filetype, get_filetype, to_symbol
   use mctc_io_convert, only : aatoau
   use numsa, only : get_numsa_version, surface_integrator, new_surface_integrator, &
      & get_vdw_rad_d3, get_vdw_rad_cosmo, get_vdw_rad_bondi, get_vdw_rad_smd
   use numsa_output, only : ascii_surface_area
   use smd_init, only: smd_param, init_smd
   use smd_sigma, only: smd_surft, calc_surft
   use smd_cds, only: calc_cds
   use smd_output, only: ascii_cds
   implicit none
   character(len=*), parameter :: prog_name = "numsa"

   type :: driver_config
      character(len=:), allocatable :: input
      character(len=:), allocatable :: solvent
      integer, allocatable :: input_format
      integer, allocatable :: grid_size
      real(wp), allocatable :: probe
      real(wp), allocatable :: offset
      real(wp), allocatable :: smoothing
      character(len=:), allocatable :: rad_type
   end type driver_config

   type(driver_config) :: config
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   type(surface_integrator) :: sasa
   type(smd_param) :: param
   type(smd_surft) :: surft
   real(wp), allocatable :: rad(:), surface(:), dsdr(:, :, :), cds(:)
   real(wp) :: cds_sm
   integer :: stat, unit
   logical :: exist

   call get_arguments(config, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) config%input_format = filetype%xyz
      call read_structure(mol, input_unit, config%input_format, error)
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (.not.allocated(config%probe)) config%probe = 1.4_wp * aatoau
   if (.not.allocated(config%grid_size)) config%grid_size = 110
   if (.not.allocated(config%rad_type)) config%rad_type = "d3"

   allocate(surface(mol%nat), dsdr(3, mol%nat, mol%nat))

   select case(config%rad_type)
   case default
      write(error_unit, '(a)') "Unknown radii type '"//config%rad_type//"' selected"
      error stop
   case("d3")
      rad = get_vdw_rad_d3(mol%num)
   case("bondi")
      rad = get_vdw_rad_bondi(mol%num)
   case("smd")
      rad = get_vdw_rad_smd(mol%num)
   case("cosmo")
      rad = get_vdw_rad_cosmo(mol%num)
   end select

   if (any(rad <= 0.0_wp)) then
      write(error_unit, '(a)') "Unknown species '"//mol%sym(minloc(rad))//"' present"
      error stop
   end if

   call new_surface_integrator(sasa, mol%id, rad, config%probe, config%grid_size, &
      & config%offset, config%smoothing)

   call sasa%get_surface(mol%id, mol%xyz, surface, dsdr)

   call ascii_surface_area(output_unit, mol, surface, rad(mol%id)+config%probe)

   if (allocated(config%solvent)) then
      call init_smd(param,config%solvent)
      call calc_surft(mol%xyz,mol%id,mol%sym,param,surft)
      call calc_cds(surft,surface,cds,cds_sm)
      call ascii_cds(output_unit,mol,cds,cds_sm)
   end if

contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "Calculates the surface area for a molecular structure input.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "    --probe <real>", "Probe radius for the integration (default: 1.4Å)", &
      "    --offset <real>", "Offset for the real space cutoff (default: 2.0Å)", &
      "    --smoothing <real>", "Smoothing parameter for the integration (default: 0.3Å)", &
      "    --rad-type <str>", "Select van-der-Waals radii (d3, cosmo or bondi)", &
      "-i, --input <format>", "Hint for the format of the input file", &
      "    --version", "Print program version and exit", &
      "    --help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_numsa_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine get_argument_as_real(iarg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg
   !> Real value
   real(wp), intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


subroutine get_argument_as_int(iarg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg
   !> Real value
   integer, intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read integer value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read integer value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_int


subroutine get_arguments(config, error)

   type(driver_config), intent(out) :: config
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   real(wp) :: val
   character(len=:), allocatable :: arg

   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("--help")
         call help(output_unit)
         stop
      case("--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)
      case("-s", "--smd")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for Solvent")
            exit
         end if
         config%solvent = arg
      case("--rad-type")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         call move_alloc(arg, config%rad_type)
      case("--grid")
         iarg = iarg + 1
         allocate(config%grid_size)
         call get_argument_as_int(iarg, config%grid_size, error)
         if (allocated(error)) exit
      case("--probe")
         iarg = iarg + 1
         call get_argument_as_real(iarg, val, error)
         if (allocated(error)) exit
         config%probe = val * aatoau
      case("--offset")
         iarg = iarg + 1
         call get_argument_as_real(iarg, val, error)
         if (allocated(error)) exit
         config%offset = val * aatoau
      case("--smoothing")
         iarg = iarg + 1
         call get_argument_as_real(iarg, val, error)
         if (allocated(error)) exit
         config%smoothing = val * aatoau
      end select
   end do

   if (.not.(allocated(config%input))) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments


end program main_driver
