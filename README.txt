PREREQUISITES
-------------

You will need to check that the following libraries are available on your system.

- GNU C compiler (gcc/g++)
- GNU scientific library (http://www.gnu.org/software/gsl/)
- Boost C++ libraries (http://www.boost.org)

These are usually pre-installed on modern linux distributions, and if not, a simple command such as 'yum install gsl-devel' and 'yum install boost-devel' (if YUM is the package manager) will suffice.

BUILDING
--------

You will then want to assign the path of the installations in locations.mk

To compile, run:

make

CONFIGURATION
-------------

You will now want to set some optional tuning parameters in the XML file that is directly relevant to the desired analysis.  Valid analysis strings are "tumor_cnv" and "seq_cnv".  The latter analysis is currently under development and should not be used yet.

The XML file must be named <analysis string>.xml (e.g. tumor_cnv.xml for the tumor CNV method).

Tumor CNV configuration
-----------------------

tumor_cnv.xml expects in one of its field a filename that contains a list of signal intensity input files.  This manifest file has a header describing the following columns: the path of the input file, a friendly name that will be pre-pended to the results to identify the dataset, and the stromal contamination level of the sample.

Each input file has a header describing the following columns: the SNP identifier, the chromosome (all records must contain the same chromosome), the base pair position (this must be sorted in increasing numerical order), the B allele frequency, and the log2 R ratio.

EXECUTION
---------

You can now execute the program by running

./analyzer <analysis string>

An example is shown in examples directory.
