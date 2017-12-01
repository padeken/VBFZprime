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
  int getGoodLeptonCombosForQCD(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst);

private:
  Analyzer* a;
  const std::string particle_names[4] = {"Ele", "Muon", "Tau", "MET"};
  const std::string particleSymbols[4]       = {"e", "#mu", "#tau", "E_{T}^{miss}"};
  const vector<CUTS> cuts={CUTS::eRElec1,CUTS::eRMuon1,CUTS::eRTau1,CUTS::eMET};
  Particle* particles[3];
  const int nstage;
  unordered_map< string,float > zTauTree;
  unordered_map< string,vector<float>* > zTauTreeV;
  TH2D* ttl_rato_one;
  TH2D* ttl_rato_two;
  Histogramer sphisto;
};

#endif
