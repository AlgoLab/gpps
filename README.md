# `gpps`: An ILP-based approach for inferring cancer progression with mutation losses from single cell data

## Running the entire `gpps` pipeline



## ILP step

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
                        not impose a limit (default).
  -o OUTDIR, --outdir OUTDIR
                        output directory.
  -d MAXDEL, --maxdel MAXDEL
                        maximum number of deletion allowed
  -c CORES, --cores CORES
                        Total number of cores (default 4).
  -b FALSEPOSITIVE, --falsepositive FALSEPOSITIVE
                        set -b False positive probability.
  -a FALSENEGATIVE, --falsenegative FALSENEGATIVE
                        set -a False negative probability.
  --mps                 This will output the model in MPS format instead of
                        running the solver

```

Where `-a` and `-b` are respectively the false negative and false positive rates for the
Single Cell Sequencing.

Running the example `python3 gpps_ilp.py -f data/examples/ex_scs.txt -a 0.1 -b 10e-4 -k 1 -d 4 -t 3600 -o example`.

## Hill Climbing step

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
  --ns NS               Hill climbing neighbourhood size.(default 30)
  --mi MI               Hill climbing maximum iterations. (default 100)
  --names NAMES         Mutation names.
```

Running the example `python3 gpps_hc.py -s data/examples/ex_scs.txt -i example/ex_scs.ilp.extended.out -k 1 -a 0.1 -b 10e-4 -k 1 -o example`.