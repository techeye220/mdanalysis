##################
Loading your data
##################


Basic Usage
===========

MDAnalysis aims to read any and all molecular simulation data,
from whatever files you provide to it.
This is done via the `Universe` object,
which accepts paths to files as its arguments.
Usually two files are required to create a `Universe` object:
A topology file, which provides information about the names and types of atoms,
and a trajectory file, which gives the positions of atoms over time.


This generally corresponds to the file required to set up the simulation acting as the topology,
while the results file of the simulation provides the trajectory.
For example to load results from a CHARMM simulation,
we provide a path to the PSF file to act as a topology,
and a path to the DCD results to act as the trajectory
::

   import MDAnalysis as mda

   u = mda.Universe('adk.psf', 'adk_dims.dcd')


.. note:: If a file which also provides coordinates is used as a topology, no trajectory
	  information is read from this file.  Ie the first frame will come from the trajectory
	  unless the `all_coordinates` keyword is set to ``True``. 

   
Single file Universes
---------------------

Occaisionally a file may contain both topology and trajectory information,
in such cases it is sufficient to provide only a single filename to Universe
::

   import MDAnalysis as mda

   u = mda.Universe('myfile.pdb')


Concatenating trajectories
--------------------------
   
It is also possible to read multiple consecutive trajectories,
(for example if a simulation was restarted),
by passing a list of trajectory filenames.
In this example, frames from `traj1.trr` and `traj2.trr` will be concatenated when iterating through the trajectory.
::

   import MDAnalysis as mda

   u = mda.Universe('topol.tpr', ['traj1.trr', 'traj2.trr'])


Supported formats and further details
=====================================

This table lists all currently supported file formats in MDAnalysis,
whether they can act as either Topology or Trajectory files,
as well as links to the relevant documentation pages.
In general MDAnalysis will try and extract all available data from a
given file, for full details of what is extracted from each file format
consult the relevant documentation page.


Specifying format
-----------------

Generally the format of a file is automatically detected from the extension,
for example a file called `system.xyz` is recognised as an XYZ file.
This can be overriden by supplying the `topology_format` and `format` keyword
arguments to Universe.
A full list of valid values for these keywords are given in the below table.


.. note:: It is possible to pass tarballed/zipped versions of files.  The
	  format detection will work around this.


.. _Supported formats:

Table of supported formats
--------------------------

+------------------+-----------+----------+--------------+-----+
| Source           | Extension | Topology | Trajectory   | I/O |
+==================+===========+==========+==============+=====+
| Generic          | XYZ       | Yes      | Yes          | r/w |
|                  +-----------+----------+--------------+-----+
|                  | PQR       | Yes      | Yes          | r/w |
+------------------+-----------+----------+--------------+-----+
| Protein          | MMTF      | Yes      | Yes          | r   |
| Databank         +-----------+----------+--------------+-----+
|                  | PDB,      | Yes      | Yes          | r/w |
|                  | ENT,      |          |              |     |
|                  | XPDB      |          |              |     |
+------------------+-----------+----------+--------------+-----+
| Autodock         | PDBQT     | Yes      | Yes          | r   |
+------------------+-----------+----------+--------------+-----+
| Amber            | TOP,      | Yes      | Yes          | r   |
|                  | PRMTOP,   |          |              |     |
|                  | PARM7     |          |              |     |
|                  +-----------+----------+--------------+-----+
|                  | TRJ,      | No       | Yes          | r   |
|                  | MDCRD     |          |              |     |
|                  +-----------+----------+--------------+-----+
|                  | INPCRD,   | No       | Yes          | r   |
|                  | RESTRT    |          |              |     |
|                  +-----------+----------+--------------+-----+
|                  | NCDF,     | Minimal  | Yes          | r/w |
|                  | NC        |          |              |     |
+------------------+-----------+----------+--------------+-----+
| CHARMM           | PSF       | Yes      | No           | r   |
|                  +-----------+----------+--------------+-----+
|                  | DCD       | Minimal  | Yes          | r/w |
|                  +-----------+----------+--------------+-----+
|                  | CRD       | No       | Yes          | r   |
+------------------+-----------+----------+--------------+-----+
| Desmond          | DMS       | No       | Yes          | r   |
+------------------+-----------+----------+--------------+-----+
| DL Poly          | CONFIG    | Yes      | Yes          | r   |
|                  +-----------+----------+--------------+-----+
|                  | HISTORY   | Yes      | Yes          | r   |
+------------------+-----------+----------+--------------+-----+
| GAMESS           | GMS,      |          | Yes          | r   |
|                  | LOG,      |          |              |     |
|                  | OUT       |          |              |     |
+------------------+-----------+----------+--------------+-----+
| Gromacs          | GRO       | Yes      | Yes          | r/w |
|                  +-----------+----------+--------------+-----+
|                  | TPR       | Yes      | No           | r   |
|                  +-----------+----------+--------------+-----+
|                  | TRR       | Minimal  | Yes          | r/w |
|                  +-----------+----------+--------------+-----+
|                  | XTC       | Minimal  | Yes          | r/w |
+------------------+-----------+----------+--------------+-----+
| Hoomd            | GSD       | Yes      | No           | r   |
|                  +-----------+----------+--------------+-----+
|                  | XML       | Yes      | No           | r   |
+------------------+-----------+----------+--------------+-----+
| IBIsCO / YASP    | TRZ       | Minimal  | Yes          | r/w |
+------------------+-----------+----------+--------------+-----+
| LAMMPS           | DATA      | Yes      | Yes          | r/w |
|                  +-----------+----------+--------------+-----+
|                  | DCD       | Minimal  | Yes          | r/w |
+------------------+-----------+----------+--------------+-----+
| Tinker           | TXYZ      | Yes      | Yes          | r   |
+------------------+-----------+----------+--------------+-----+
| Tripos           | MOL2      | Yes      | Yes          | r/w |
+------------------+-----------+----------+--------------+-----+

