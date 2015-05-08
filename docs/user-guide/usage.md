## Overview

You will now want to set some optional tuning parameters in the XML file that is directly relevant to the desired analysis.  Valid analysis strings are "tumor_cnv" and "seq_cnv".  The latter analysis is currently under development.

The XML file must be named <analysis string>.xml (e.g. tumor_cnv.xml for the tumor CNV method).

## Tumor CNV configuration

`tumor_cnv.xml` expects in one of its field a filename that contains a list of signal intensity input files.  This manifest file has a header describing the following columns: the path of the input file, a friendly name that will be pre-pended to the results to identify the dataset, and the stromal contamination level of the sample.

Each input file has a header describing the following columns: the SNP identifier, the chromosome (all records must contain the same chromosome), the base pair position (this must be sorted in increasing numerical order), the B allele frequency, and the log2 R ratio.

## Execution

You can now execute the program by running

```
./analyzer <analysis string>
```

An example is shown in examples directory.
