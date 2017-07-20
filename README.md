# flu-sampler project

This program can be used for estimating epidemiological parameters.
This is achieved with Gibbs sampling.

## compiling the code

Simply 'cd' to the 'flu-sampler' folder, and run 'make':
```
#!bash
$ cd flu-sampler-project/flu-sampler/
$ make
```

## running the program

From the 'flu-sampler' folder, run
```
#!bash
$ ./flu-sampler [-s <seed>] [-v <vid>] [-m <mode>] [-l <length>] [-t <thinning>] [-i <imitations>]
```

where <seed> is a seed (non-negative integer) for the random number generator;
<vid> is an identifier for the output files;
<mode> can be "mcmc" for MCMC, "wbic" for WBIC, "help" for help;
<length> is the length of the chain (pre-thinning)
<thinning> is the number of states skipped between stored states
<imitations> is the number of simulated timeseries


**TODO: new features added** 

## data files

The following files need to be present in the data/essential folder:
Data needed for MCMC:

- `m-matrix-small.txt` -- contact matrix for the 6 age classes
- `sortedIliData.txt` -- simplified and combined the ILI data
- `seasons-daynumbers.txt` -- seasons defined by an interval of day numbers (from 01/01/1970 onwards)

Data needed for infection histories

- epitope-virus-testmatrix.dat -- which epitopes are present in the viruses (TODO: this is just a test matrix)
- epitope-hla-testmatrix.dat -- which epitopes are presented by the HLAs (TODO: this is just a test matrix)

If all goes well, the following files will be created in the data/ folder 
after running the program:

- full-ili-chain-<vid>.xml -- the thinned chain (posterior)
- sili-<vid>.xml -- simulated data, parameters sampled from the posterior

where <vid> is an identifier given by the uses (defaults to "test"

**TODO: new IO files added** 

## Jupyter notebooks

Jupyter notebooks with the extension `.ipython.cln` are included in the repository.
These are *clean* versions of the notebooks (i.e. without output).
In order to make the actual notebooks, run from the folder `flu-sampler-project/flu-sampler/`
```
#!bash
make -f nbs.makefile
```
The following notebook can be used for analysing the chain, plotting data,
producing easy-to-parse data files:

- `ili-data-and-fit-age-stratified.ipynb` -- load the chain, and analyse
- `plot-simulated-data-age-stratified.ipynb` -- load simulated data (sampled from posterior) and plot
- `sort-raw-ili-data.ipynb` -- make data files (requires original ILI data)
- TODO

Before you want to commit modified notebooks to the repository, run
```
#!bash
make -f cln.makefile
```
This will write cleared copies of the notebooks `x.ipynb` to the files `x.ipynb.cln`