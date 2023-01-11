# This file is part of numsa.
# SPDX-Identifier: LGPL-3.0-or-later
#
# numsa is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numsa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with numsa.  If not, see <https://www.gnu.org/licenses/>.

project(
  'numsa',
  'fortran',
  version: '0.2.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
lib_deps = []
inc_dirs = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
numsa_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
numsa_inc = numsa_lib.private_dir_include()
numsa_dep = declare_dependency(
  link_with: numsa_lib,
  include_directories: numsa_inc,
  dependencies: lib_deps,
)

# Add executable targets
subdir('app')

# Package the license files
numsa_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    numsa_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    numsa_lib,
    description: 'Solvent accessible surface area calculations',
    subdirs: ['', module_id],
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man/numsa.1.adoc'),
        output: '@BASENAME@',
      )
    )
  endif
endif

# add the testsuite
subdir('test')
