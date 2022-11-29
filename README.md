# Population
**New generation population dynamics model with a dynamic population structure**

This is the third-generation implementation of the dynamically-structured matrix population model ([sPop](https://doi.org/10.12688/f1000research.15824.3) and [sPop2](https://doi.org/10.1038/s41598-022-15806-2)). As in [Population.jl](https://github.com/kerguler/Population.jl), population new accepts a list of processes applied sequentially on each member-class. This version implements both age-dependent and accumulative processes in a uniform structure.

## Installation

The library has been developed and tested on Linux and MacOS systems, but not on Windows. To install, download the source file `population-x.y.z.tar.gz`, unzip and untar (`tar -xvzf population-x.y.z.tar.gz`), and execute the `configure` - `make` - `make install` sequence.

**Dependencies**

The library is dependent on the GNU Scientific Library ([GSL](https://www.gnu.org/software/gsl/)), which is not included in this distribution. Please make sure GSL is installed before configuring the `population` library.

## Using the library

Please refer to the `population` [pages](https://kerguler.github.io/Population/) for usage instructions and examples.


