Surface area calculator
=======================

Numerical surface area integrator for molecular inputs.
This project is based on routines from [``xtb``](https://github.com/grimme-lab/xtb) and [``dftb+``](https://github.com/dftbplus/dftbplus).


## Installing

To compile this version of ``numsa`` the following programs are needed
(the number in parentheses specifies the tested versions).

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008 (GFortran and Intel Fortran are known to work)
- [meson](https://mesonbuild.com) version 0.53 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer

Optional dependencies are
- asciidoctor to build the manual page

Setup a build with

```sh
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile and run the projects testsuite use

```sh
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```sh
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.


## Usage

Surface area calculations can be performed with the ``numsa`` executable.
To calculate the surface area for an input run:

```sh
numsa coord
```

For an overview over all command line arguments use the ``--help`` argument or checkout the [``numsa(1)``](man/numsa.1.adoc) manpage.


### Meson project

To use the ``numsa`` library in your meson project, include it as dependency with

```meson
numsa_dep = dependency('numsa', ['numsa', 'numsa_dep'])
```

and add ``numsa.wrap`` to your subprojects directory

```ini
[wrap-git]
directory = numsa
url = https://github.com/grimme-lab/numsa.git
revision = head
```


### Fpm project

This project can be used with the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)).
Include ``numsa`` as dependency in your package manifest ``fpm.toml`` with

```toml
[dependencies]
numsa.git = "https://github.com/grimme-lab/numsa.git"
```


### Library usage

This library can be used in Fortran projects by importing the ``numsa`` module.
It provides a ``surface_integrator`` type with a ``get_surface`` method to perform the actual integration.
Input to the constructor are the van-der-Waals radii for each species, the probe radius of the solvent molecule and the number of grid points for each atom.

```f90
!> Example implementation to calculate surface area for a molecule input
subroutine get_surface_area(species, symbols, coord, probe, surface, dsdr)
   use mctc_env, only : wp
   use numsa, only : surface_integrator, new_surface_integrator, get_vdw_rad_bondi, grid_size
   !> Unique chemical species in the input structure, shape: [nat]
   integer, intent(in) :: species(:)
   !> Element symbol for each chemical species, shape: [nsp]
   character(len=*), intent(in) :: symbols(:)
   !> Cartesian coordinates in Bohr, shape: [3, nat]
   real(wp), intent(in) :: coord(:, :)
   !> Probe radius for surface area integration in Bohr
   real(wp), intent(in) :: probe
   !> Accessible surface area in Bohr², shape: [nat]
   real(wp), intent(out) :: surface(:)
   !> Derivative of surface area w.r.t. atomic displacements, shape: [3, nat, nat]
   real(wp), intent(out) :: dsdr(:, :, :)

   type(surface_integrator) :: sasa
   real(wp), allocatable :: rad(:)

   rad = get_vdw_rad_bondi(symbols)
   call new_surface_integrator(sasa, species, rad, probe, grid_size(8))

   call sasa%get_surface(species, xyz, surface, dsdr)

end subroutine get_surface_area
```

For convenience the ``numsa`` module supports access to different van-der-Waals radii, including the DFT-D3 van-der-Waals radii, Bondi radii and COSMO radii.
Also, the 32 supported angular grid sizes are available with the ``grid_size`` parameter.


## References

- Angular integration grids:
  V.I. Lebedev, and D.N. Laikov, A quadrature formula for the sphere of the 131st algebraic order of accuracy, *Doklady Mathematics*, Vol. 59, No. 3, **1999**, pp. 477–481.

- Smooth numerical integration:
  W. Im, M.S. Lee, and C.L. Brooks III, Generalized Born model with a simple smoothing function, *J. Comput. Chem.*, Vol. 24, No. 14, **2003**, pp. 1691–1702.

- DFT-D3 van-der-Waals radii:
  S. Grimme, J. Antony, S. Ehrlich, and H. Krieg, A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu, *J. Chem. Phys.*, Vol. 132, **2010**, p. 154104.

- Bondi van-der-Waals radii:
  M. Mantina, A.C. Chamberlin, R. Valero, C.J. Cramer, and D.G. Truhlar, Consistent van der Waals Radii for the Whole Main Group, *J. Phys. Chem. A*, Vol. 113, No. 19, **2009**, pp.. 5806–5812.


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
Lesser GNU General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Lesser GNU General Public license, shall be licensed as above, without any
additional terms or conditions.
