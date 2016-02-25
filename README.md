# lattice-playground

Simple gauge theories on the [lattice](https://xkcd.com/505/).

## details

following original `roots' of lattice qcd: [Mike Creutz et al, PRL 1979](http://thy.phy.bnl.gov/~creutz/mypubs/pub031.pdf)

## observables

Calculating the average plaquette energy is perhaps the most straighforward _measurement_. The action is normalised so that the average plaquette energy is 1 for beta=0 and 0 for low temperature.

The Polyakov loop is a closed path on the lattice, only in the temporal direction. It features as an order parameter for confinment and from it one can extract the static quark potential.

## visualisation

data is written to ``~/out/data/``. Pretty pictures using [GLE](http://glx.sourceforge.net/), can compile
directly by running ``make plotter`` after the main executable has been run. See the
makefile for details.

