### What is SMAT? ###

SMAT is a simple C program with a configuration file-driven interface
supporting matrix multiplication and QR decomposition implemented in Scalapack
for the purposes of "burning in" an arbitrary number of cluster nodes while
testing performance.

### Future Work ###

Future improvements will include more LAPACK routines and implementations
supporting heterogeneous architectures such as GPGPUs (MAGMA, clMAGMA) and
Intel MICs (MAGMA MIC).

### License ###

LAMMPS-HESSIAN Copyright (C) 2014 Anthony B. Costa, Numerical Solutions, Inc.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

### Citation ###

Publications making substantial use of SMAT or its derivatives should cite this
software package and freely available online repository. An example citation is
given as:

    Costa, A. B., "SMAT", Package Version <Insert Version Number>,
    http://bitbucket.org/numericalsolutions/smat.

### Contact ###

Anthony B. Costa, anthony.costa@numericalsolutions.org
