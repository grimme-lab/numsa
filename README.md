Surface area calculator
=======================

Numerical surface area integrator for molecular inputs.


## Installing

To compile this version of ``numsa`` the following programs are needed
(the number in parentheses specifies the tested versions).

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
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
