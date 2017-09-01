#include "tsutil.h"

namespace galleryfmwk {

bool inFV(const sim::MCShower& s) {
  // From J. Zennamo
  const double x_lo =    0.0, x_hi =  256.35;
  const double y_lo = -116.5, y_hi =  116.50;
  const double z_lo =    0.0, z_hi = 1036.80;
  
  return (s.End().X() > x_lo && s.End().X() < x_hi &&
          s.End().Y() > y_lo && s.End().Y() < y_hi &&
          s.End().Z() > z_lo && s.End().Z() < z_hi);
}


bool inFV(const sim::MCTrack& t) {
  const double x_lo =    0.0, x_hi =  256.35;
  const double y_lo = -116.5, y_hi =  116.50;
  const double z_lo =    0.0, z_hi = 1036.80;
  
  return (t.End().X() > x_lo && t.End().X() < x_hi &&
          t.End().Y() > y_lo && t.End().Y() < y_hi &&
          t.End().Z() > z_lo && t.End().Z() < z_hi);
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector showStart = show.Start().Position();
  return (showStart - nuVtx).Mag() < distance;
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector trkStart = track.Start().Position();
  return (trkStart - nuVtx).Mag() < distance;
}

double eccqe(const TLorentzVector& v) {
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 939.565; // MeV/c^2
  double M_p = 938.272; // MeV/c^2
  double M_e = 0.511; // MeV/c^2
  double bindingE = 30.0; // MeV

  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  double l_energy = v.E();

  double l_mom = sqrt(l_energy*l_energy - me2);

  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));

  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));

  return enu_top / enu_bot;
}

}  // namespace gallery_fmwk

