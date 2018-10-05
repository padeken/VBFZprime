//#ifndef: If the macro isn't defined by a #define statement, then the code following the #ifndef statement will compile.
#ifndef analysis_h
//Putting a #define statement without a replacement string defines the macro to be empty; however, the macro is still defined (i.e. an empty macro is different than one that doesn't exist)
#define analysis_h

//struct command is used to define "Analyzer" as a structure. A data structure is a list of elements grouped together under one name. According to the website that I found, structures are defined as follows:
//  struct type_name {
//  member_type1 member_name1;
//  member_type2 member_name2;
//  member_type3 member_name3;
//  .
//  .
//  } object_names;
//This means that the data elements (called members) can have different types and lengths. The declaration makes a new type, type_name, which can be used to define other objects.
//I tried to understand why this (and other) struct commands have ONLY a type_name, and no defined members. For now, I am taking this to mean that it creates the structure called Analyzer, but it is empty (it simply defines Analyzer as a data structure). However, there are other header files involved: notably, Analyzer.h is brought in right after, and Analyzer.h itself has a lot of header files in there as well. Maybe these give Analyzer structure elements?
struct Analyzer;
#include "src/Analyzer.h"

//The class statement creates a new class: in this case, SpechialAnalysis. Classes are user-defined data types; a description of a class is organized in two files: a header file and an implementation file. This is the header file; the .cc file implements the class.
class SpechialAnalysis {
  // This access specifier means that the following list is visible to users of the class.
public:
  // This is a constructor that creates objects according to the given argument
  SpechialAnalysis(Analyzer* _a);
  // Void=nothing. It's used when the function isn't supposed to return anything, I think.
  // Init(); is a contructor that increases the internal static constructor by one; If the value of the internal counter was zero, the standard iostream objects are constructed and initialized; if they have not yet been constructed and initialized.
  // I supposed that the rest of these contructors are defined in another header file?
  void init();
  void begin_run();
  void analyze();
  void end_run();
  void fill_particle(int);
  void make_tau_tree();
  void make_TTLAna();
  //Getters allow for controlled access to parts of an object.
  int getGoodLeptonCombosForQCD(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst);

  //Private: All the following members are private, i.e. they are not visible to users of the object or classes that inherit from this class
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
// I'm not sure where the #if is, but #if, #else, and #endif are used as instructions for the compiler so that it only compiles what is between the #if and the #endif.
#endif
