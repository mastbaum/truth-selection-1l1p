/**
 * \file tsutil.h
 * \brief Utilities
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef GALLERY_FMWK_TSUTIL_H
#define GALLERY_FMWK_TSUTIL_H

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"

namespace galleryfmwk {

// Check if shower endpoint is within the fiducial volume
bool inFV(const sim::MCShower& show);

// Check if track is contained (N.B. just checking endpoint for now)
bool inFV(const sim::MCTrack& track);

// Check if shower is from vertex (within distance tolerance)
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance=50.0);

// Check if track is from vertex (within distance tolerance)
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance=50.0);

// Calculate the CCQE energy from lepton four-momentum
double eccqe(const TLorentzVector& v);

}  // namespace galleryfmwk

#endif  // GALLERY_FMWK_TSUTIL_H

