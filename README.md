Truth-Based 1l1p Selection
==========================
A. Mastbaum <mastbaum@uchicago.edu>, 2017/08/13

Loop through MCTracks and MCShowers, doing PID based on dE/dx and other
observable parameters, to build a list of tracks and showers in the event,
then select based on that.

For tracks:
* MCTrack must not be empty
* Energy at start of track > 60 MeV
* Track starts within 10 cm of the neutrino vertex

Then, we do a dE/dx-based PID. For this, we build a dE/dx vs. residual range
joint distribution for all tracks of a given PID, and store them in a ROOT
file. Then, for each event we compute a "KS" probability against each
distribution to find the best match, and assign that PID. Anything called a
proton with a track length shorter than 1.2 cm or longer than 8 cm is marked
as unknown. Then anything where the KS test failed to find a match is called
a proton. This turns out to be pretty good.

For showers, we just require that they come from the vertex and are >60 MeV.
Then, the shower is called an electron or gamma based on whether the dE/dx
falls below or above a cut at 3.5 MeV/cm. This lets in a small fraction of
gammas.

