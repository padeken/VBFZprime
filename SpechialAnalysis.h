#ifndef analysis_h
#define analysis_h

struct Analyzer;
#include "src/Analyzer.h"

class SpechialAnalysis {
public:
  SpechialAnalysis(Analyzer* _a);
  void init();

  void begin_run();
  void analyze();
  void end_run();
  void fill_particle(int);
  void make_tau_tree();
  void make_TTLAna();

private:
  Analyzer* a;
  const std::string particle_names[4] = {"Ele", "Muon", "Tau", "MET"};
  const std::string particleSymbols[4]       = {"e", "#mu", "#tau", "E_{T}^{miss}"};
  const vector<CUTS> cuts={CUTS::eRElec1,CUTS::eRMuon1,CUTS::eRTau1,CUTS::eMET};
  Particle* particles[3];
  const int nstage;
  unordered_map< string,float > zTauTree;
};

#endif
