`gpps`: An ILP-based approach for inferring cancer progression with mutation losses from single cell data
===================================================

ILP step
--------
We provide *gpps*,
that can be used to infer cancer progressions from single cell data.
Differently from the previous tool, *gpps* employs a maximum likelihood search
to find the best tree that explain the input, starting from single cell data.

The tool can be run with the following arguments:

```
usage: gpps_ilp.py [-h] -f FILE [-k K] -t TIME -o OUTDIR [-d MAXDEL] [-e] -b
                   FALSEPOSITIVE -a FALSENEGATIVE [--mps]

gpps- ILP

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  path of the input file.
  -k K                  k-value of the selected model. Eg: Dollo(k)
  -t TIME, --time TIME  maximum time allowed for the computation. Type 0 to
                        not impose a limit.
  -o OUTDIR, --outdir OUTDIR
                        output directory.
  -d MAXDEL, --maxdel MAXDEL
                        maximum number of deletion allowed
  -e, --exp             set -e to get experimental-format results.
  -b FALSEPOSITIVE, --falsepositive FALSEPOSITIVE
                        set -b False positive probability.
  -a FALSENEGATIVE, --falsenegative FALSENEGATIVE
                        set -a False negative probability.
  --mps                 This will output the model in MPS format instead of
                        running the solver

```

Where `-a` and `-b` are respectively the false negative and false positive rates for the
Single Cell Sequencing.


Hill Climbing step
--------------------

```
sage: hill_climbing.py [-h] -i ILPFILE -s SCSFILE -k K -o OUTDIR -b
                        FALSEPOSITIVE -a FALSENEGATIVE --ns NS --mi MI
                        [--names NAMES]

gpps- Hill Climber

optional arguments:
  -h, --help            show this help message and exit
  -i ILPFILE, --ilpfile ILPFILE
                        path of the ILP output file.
  -s SCSFILE, --scsfile SCSFILE
                        path of the SCS input file. (same input feeded to the
                        ILP)
  -k K                  k-value of the selected model. Eg: Dollo(k)
  -o OUTDIR, --outdir OUTDIR
                        output directory.
  -b FALSEPOSITIVE, --falsepositive FALSEPOSITIVE
                        set -b False positive probability.
  -a FALSENEGATIVE, --falsenegative FALSENEGATIVE
                        set -a False negative probability.
  --ns NS               Hill climbing neighbourhood size.
  --mi MI               Hill climbing maximum iterations.
  --names NAMES         Mutation names.
```