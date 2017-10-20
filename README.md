Truth-Based 1l1p Selection
==========================
A. Mastbaum <mastbaum@uchicago.edu>, 2017/08/13

Loop through MCTracks and MCShowers, doing PID based on dE/dx and other
observable parameters, to build a list of tracks and showers in the event,
then select based on that.

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

Scripts
-------
The `scripts` directory contains:

* Python scripts used to run each module (see the Tools section)
* Utilities in `scripts/util`
  * File handling and MC POT counting
* Analysis scripts in `scripts/analysis`

See header blocks at the top of each script for more information, including
usage.

