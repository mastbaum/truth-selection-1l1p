Truth-Based 1l1p Selection
==========================
A. Mastbaum <mastbaum@uchicago.edu>, 2017/08/13

Loop through MCTracks and MCShowers, doing PID based on dE/dx and other
observable parameters, to build a list of tracks and showers in the event,
then select based on that.

**N.B. This documentation refers to an older version.**

Tools
-----
### TSSelection: Truth-Based Selection ###
This module is the main truth-based selection code. It scans through Monte
Carlo events to extract MCTracks and MCShowers, assigns them particle IDs
based on dE/dx information (using PDFs from `TSPDFGen`), applies selection
cuts (e.g. kinetic energy), and classifies events as 1e1p and 1mu1p. The
selected events are written into an output analysis tree.

**Usage:**

    $ python scripts/selection.py OUTPUT.root DATASET "INPUT*.root"

The output file `OUTPUT.root` will contain the output analysis trees. The
`DATASET` is an integer flag which is written into the output trees, as a
convenience for bookkeeping when working with multiple MC samples. The
final argument is a set of input art ROOT files.

### TSCovariance: Covariance Matrix Generator ###
This module takes analysis trees with selected event samples and Monte Carlo
event weights, and computes a covariance matrix for a user-specified subset
of weights. The output is a weighted event histogram, a graph with systematic
uncertainty bands, and the covariance and correlation matrices.

**Usage:**

    $ python scripts/cov.py input.root

Here, the input is an analysis tree ROOT file, i.e. the output of
`TSSelection`.

The user should modify `cov.py` to select an appropriate combinations of
weights.

### TSPDFGen: PDF Generator ###
The PDF Generator produces probability distributions used by `TSSelection` for
track and particle identification. It constructs a 2D histogram in dE/dx vs.
residual range for tracks, and a 1D histogram of dE/dx for showers, based on
the Monte Carlo truth information in the events in the input dataset.

**Usage:**

    $ python scripts/make_pdfs.py "input*.root"

The input is a set of art ROOT files.

This will create a file `ana_out.root` containing histograms named like
`htrackdedx_<PDG>` for track-like PDG codes, and `hshowerdedx_<PDG>` for
shower-like PDGs.

A set of PDFs (which is used by `TSSelection` is provided in
`data/dedx_pdfs.root`. To update the PDFs for the selection, run `TSPDFGen`
and replace that file.

### TSUtil: Utilities ###
This file contains common utilities used by the truth-based selection modules,
which are in namespace `tsutil`.

### Scripts ###
The `scripts` directory contains:

* Python scripts used to run each module (see the Tools section)
* Utilities in `scripts/util`
  * File handling and MC POT counting
* Analysis scripts in `scripts/analysis`

See header blocks at the top of each script for more information, including
usage.

Workflow
--------
The workflow begins with raw art ROOT files, produced by LArSoft (uboonecode).
In practice we use two Monte Carlo samples: one fully inclusive BNB neutrino
sample with MC cosmic ray overlays, and one with only nue events (to
improve the statistics for nues, since they constitute only a small fraction
of the BNB inclusive sample). Additionally, there are MC samples produced for
various signal models.

In general, MC weighting functions that assign weights based on variations in
the systematics parameters are not run on the production MC files, so we
must run them manually using grid-based LArSoft jobs that run the EventWeight
module (and can also drop most data products, since this analysis only looks
at MC truth). In this way we obtain copies of the MC datasets with only MC
information, imbued with a set of systematics weights attached to each event.
The POT counting script in `scripts/util/pot.py` can be used to extract the
equivalent POT for these datasets, which is necessary for scaling to various
exposures.

The selection script (see TSSelection above) operates on a list of these
weighted art ROOT files, applying the truth-based selection. Typically this
would be run three times -- for the electron, inclusive, and signal samples --
and produce three output files. The electron and inclusive samples can now be
added together, using `scripts/util/add.C`. This (in contrast to `hadd`) takes
CC nue and nuebar events from the electron sample and *everything else* from
the inclusive sample. Splitting in this way simplifies the reweighting, as
we combine samples representing different POT exposures. After the addition,
we have two ROOT files -- a combined 1e1p + 1mu1p background file, and a signal
file -- containing simple analysis ROOT trees.

We can now run `scripts/cov.py` on the background ROOT file. This script will
iterate through the chosen weights (the set of which comprise a systematics
"universe") and produce a number of weighted spectra with and without
systematic error bands, plus the covariance and correlation matrices. The
script allows the user to toggle on and off different combinations of
parameters, so that one can produce a covariance for e.g. the flux
uncertainties alone, one particular GENIE uncertainty, or everything together.
The `plotcov.py` script is provided for convenience, to draw these plots into
reasonable-looking PDF files.

The `scripts/analysis` directory provides some basic analysis code which
operates on either the analysis tree files, or the histograms output from
`cov.py` The main script, `chi2.py` reads in the latter and can make a
variety of plots, and compute a chi^2 significance for observing the signal
prediction under (background-only) null hypothesis using an analytic
approximation or a Monte Carlo approach, which represents the final product
of this analysis chain.

