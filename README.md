# dihbsa - EIC BRANCH
dihadron beam spin asymmetry for CLAS12

* this branch is specifically for the yellow report dihadron partial wave
  asymmetry projections; it is based off analysis code for CLAS12 dihadrons,
  and therefore contains a lot of extra/unused code
* it also requires some CLAS12 dependencies
* if you wish to run this code to produce projections, ask the developer to
  make this repository more stand-alone and streamlined for EIC purposes,
  removing CLAS12 dependencies etc.; for now the code is stored as is, for
  preservation purposes

# Procedure
  * first produce ROOT files using
    [`eicsim`](https://gitlab.com/c.dilks/eicsim) plugin `dihadron`; follow
    documentation there
    * `eicsim` commit
      [`c4505314`](https://gitlab.com/c.dilks/eicsim/-/tree/c45053146d9cc14ad4b3af80ec4a4700c62e28c2)
      should be sufficient; later commits might no longer be compatible
    * you need to run an `escalate` docker container; yellow report plots
      were generated using `escalate 1.0.0`, though `1.1.0` should also work
    * see documentation in `eicsim` for setup; here's the procedure outline:
      * create docker image, cd to `dihadronAnalysis`
      * run `pythia`, or whatever generator you prefer; check
        steer files carefully to make sure the settings are correct,
        such as minimum Q2
      * fast simulation is done by `eic-smear` in the `eJANA` framework
      * ROOT files contain the TTree `tree`
* return back to this repository; the `outroot/` subdirectory should contain
  these ROOT files; it is useful to link `outroot -> /path/to/eicsim/dihadronAnalysis/dataRoot`, since that's the default output area for `eicsim`
* adjust cut settings in `src/EventTree.cxx`, such as the minimum pT cut;
  check the code carefully to see if the cuts are acceptable for you
  * generator cuts, such as minimum generated Q2, are applied in `eicsim`
* to make acceptance plots, use `acceptanceEIC.cpp`
* to make projections, follow usual spin asymmetry procedure:
  * `buildSpinroot.cpp` reads trees in `outroot/` and produces root files
    in `spinroot/`
    * if you have multiple root files, use a wrapper to process them in 
      parallel (you may have to create log output directories):
      * condor wrapper: `loopBuildSpinroot.sh`
      * slurm wrapper: `slurm.buildSpinroot.sh`
  * `catSpinroot.cpp` to "`hadd`" the root files together
  * `asymFit.cpp` to fit the desired asymmetry (see code for arguments and
    their meaning)
* postprocessing
  * use `ProjectorEIC_yellowreport.C`
    * it needs two `asym*.root` files produced from `asymFit.cpp`: one
      for pT>100 MeV cut, and another for pT>300 MeV cut
    * make sure the `error scale factor` is properly calculated for your
      needs; you need to supply the cross section for your beam energy setting,
      the number of events you generated, and the luminosity to which
      you want to project; check the print outs of all relevant values
      to make sure the numbers are correct
    * check that the proper value of polarization is considered
      (currently set to 0.7)
    * check that the depolarization ratios are correct for each bin
      * the numbers are hardcoded (sorry); they are the mean depolarization
        ratios for each bin
      * I obtained these numbers from plots produced by `diagnostics.cpp`;
        they are all named `kfVs...`; you will need to modify the code
        to put in your own bin boundaries and determine the mean depolarization
        factors
      * this could be automated
    * execute the macro
    * the result is `dihadronPWprojection.png`
  * the vector file `dihadronPWprojection.svg` links directly to
    the `dihadronPWprojection.png` file, and contains overlays
    such as the LaTeX labels; edit it with Inkscape
    * export the SVG file as a PNG, with `dpi=100`
    

# Build instructions
* the `Clas12Tool` dependency should obviously be removed...
make sure `Clas12Tool` is at `../Clas12Tool`
make sure `Clas12Tool/lib` is in `$LD_LIBRARY_PATH`
type `make` to build


# Constants.h
This header contains several variables, such as pion mass, dihadron pair names, etc.
In the source code, if you can't find where a particular variable is defined, try
looking here, since variables and functions in Constants.h are used frequently
