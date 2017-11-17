#include "SpechialAnalysis.h"
#include "HistClass.h"
#include <TFile.h>
#include <Compression.h>

#define BIG_NUM 46340

SpechialAnalysis::SpechialAnalysis(Analyzer* _a):nstage(9) {
  a=_a;
}

void SpechialAnalysis::init() {
  

  //this is just to store the histograms if they are to much to hold in memory
  // string safeFileName = "SpecialHistos.root";
  // file1               = new TFile(safeFileName.c_str(), "RECREATE");

  // number of events, saved in a histogram
  HistClass::CreateHistoUnchangedName("h_counters", 10, 0, 11, "N_{events}");

  for (unsigned int i = 0; i < 4; i++) {
    for (unsigned int j = 0; j < 6; j++) {
      HistClass::CreateHisto(nstage, "num", particle_names[i].c_str(), 40, 0, 39, TString::Format("N_{%s}", particleSymbols[i].c_str()));
      HistClass::CreateHisto(nstage, ("pt_"+to_string(j)).c_str(), particle_names[i].c_str(), 5000, 0, 5000, TString::Format("p_{T}^{%s} (GeV)", particleSymbols[i].c_str()));
      HistClass::CreateHisto(nstage, ("eta_"+to_string(j)).c_str(), particle_names[i].c_str(), 80, -4, 4, TString::Format("#eta_{%s}", particleSymbols[i].c_str()));
      HistClass::CreateHisto(nstage, ("phi_"+to_string(j)).c_str(), particle_names[i].c_str(), 40, -3.2, 3.2, TString::Format("#phi_{%s} (rad)", particleSymbols[i].c_str()));
    }

    if (not a->isData) {
      HistClass::CreateHisto(1, "num_Gen", particle_names[i].c_str(), 40, 0, 39, TString::Format("N_{%s}", particleSymbols[i].c_str()));
      HistClass::CreateHisto(1, "pt_Gen", particle_names[i].c_str(), 5000, 0, 5000, TString::Format("p_{T}^{%s} (GeV)", particleSymbols[i].c_str()));
      HistClass::CreateHisto(1, "eta_Gen", particle_names[i].c_str(), 80, -4, 4, TString::Format("#eta_{%s}", particleSymbols[i].c_str()));
      HistClass::CreateHisto(1, "phi_Gen", particle_names[i].c_str(), 40, -3.2, 3.2, TString::Format("#phi_{%s} (rad)", particleSymbols[i].c_str()));
    }
  }
  
  for (unsigned int j = 0; j < 1; j++) {
    int i=2;
    HistClass::CreateHisto(14, "num", particle_names[i].c_str(), 40, 0, 39, TString::Format("N_{%s}", particleSymbols[i].c_str()));
    HistClass::CreateHisto(14, ("pt_"+to_string(j)).c_str(), particle_names[i].c_str(), 5000, 0, 5000, TString::Format("p_{T}^{%s} (GeV)", particleSymbols[i].c_str()));
    HistClass::CreateHisto(14, ("eta_"+to_string(j)).c_str(), particle_names[i].c_str(), 80, -4, 4, TString::Format("#eta_{%s}", particleSymbols[i].c_str()));
    HistClass::CreateHisto(14, ("phi_"+to_string(j)).c_str(), particle_names[i].c_str(), 40, -3.2, 3.2, TString::Format("#phi_{%s} (rad)", particleSymbols[i].c_str()));
  }

  

  for(size_t jpart =0; jpart<6; jpart++ ){
    zTauTree["tau_"+to_string(jpart)+"_pt"]       = 0;
    zTauTree["tau_"+to_string(jpart)+"_eta"]      = 0;
    zTauTree["tau_"+to_string(jpart)+"_phi"]      = 0;
    zTauTree["tau_"+to_string(jpart)+"_charge"]   = 0;
    zTauTree["cosDphi_"+to_string(jpart)]         = 0;
    zTauTree["mt_tau_"+to_string(jpart)]          = 0;
    zTauTree["tau_"+to_string(jpart)+"_iso"]      = 0;
  }
  zTauTree["tau_n"]          =0;
  zTauTree["tauIso_n"]       =0;
  zTauTree["tauNonIso_n"]    =0;
  
  zTauTree["met"]           = 0;
  zTauTree["muo1_pt"]       = 0;
  zTauTree["muo1_eta"]      = 0;
  zTauTree["muo1_phi"]      = 0;
  zTauTree["muo1_charge"]   = 0;
  zTauTree["mt_muo1"]       = 0;
  zTauTree["ele1_pt"]       = 0;
  zTauTree["ele1_eta"]      = 0;
  zTauTree["ele1_phi"]      = 0;
  zTauTree["ele1_charge"]   = 0;
  zTauTree["mt_ele1"]       = 0;
  zTauTree["jet_n"]         = 0;
  zTauTree["jet1_pt"]       = 0;
  zTauTree["jet1_eta"]      = 10;
  zTauTree["jet1_phi"]      = 10;
  zTauTree["jet2_pt"]       = 0;
  zTauTree["jet2_eta"]      = 10;
  zTauTree["jet2_phi"]      = 10;
  zTauTree["b_jet_n"]       = 0;
  zTauTree["b_jet_pt"]      = 0;
  zTauTree["jet_mass"]      = 0;
  zTauTree["weight"]        = 0;
  
  HistClass::CreateTree(&zTauTree,"TauIsoTree");
}

void SpechialAnalysis::begin_run() {
  particles[0]=a->_Electron;
  particles[1]=a->_Muon;
  particles[2]=a->_Tau;
}

void SpechialAnalysis::analyze() {
  //just because it is convinient read in  the cuts:
  // const unordered_map<string,pair<int,int> >* cut_info = a->histo.get_cuts();
  // const vector<string>* cut_order = a->histo.get_cutorder();
  
  if(a->active_part->at(a->cut_num.at("NRecoTriggers1"))->size()==0 )
    return;
  if(a->active_part->at(a->cut_num.at("NRecoVertex"))->size()==0 )
    return;
  
  fill_particle(0);
  
  vector<string>::iterator isobool =find(a->_Tau->pstats["Tau1"].bset.begin(),a->_Tau->pstats["Tau1"].bset.end(),"DoDiscrByIsolation");
  if(isobool!=a->_Tau->pstats["Tau1"].bset.end()){
    a->_Tau->pstats["Tau1"].bset.erase(isobool);
  }
  a->active_part->at(CUTS::eRTau1)->clear();
  a->getGoodRecoLeptons(*a->_Tau, CUTS::eRTau1, CUTS::eGTau, a->_Tau->pstats["Tau1"],0);
  //restore current version
  a->_Tau->pstats["Tau1"].bset.push_back("DoDiscrByIsolation");
  
  make_TTLAna();
  make_tau_tree();
  

}

void SpechialAnalysis::make_tau_tree() {
  
  
  int j1      = -1;
  int j2      = -1;
  double mass = 0;
  for (auto it : *a->active_part->at(CUTS::eDiJet)) {
    int j1tmp = (it) / a->_Jet->size();
    int j2tmp = (it) % a->_Jet->size();
    if (a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "") > mass) {
      j1   = j1tmp;
      j2   = j2tmp;
      mass = a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "");
    }
  }
  
  int iitau=0;
  int isotau=0;
  int non_isotau=0;
  
  for(size_t jpart =0; jpart<6; jpart++ ){
    zTauTree["tau_"+to_string(jpart)+"_pt"]       = 0;
    zTauTree["tau_"+to_string(jpart)+"_eta"]      = 0;
    zTauTree["tau_"+to_string(jpart)+"_phi"]      = 0;
    zTauTree["tau_"+to_string(jpart)+"_charge"]   = 0;
    zTauTree["cosDphi_"+to_string(jpart)]         = 0;
    zTauTree["mt_tau_"+to_string(jpart)]          = 0;
    zTauTree["tau_"+to_string(jpart)+"_iso"]      = 0;
  }
  for(size_t jpart =0; jpart<a->active_part->at(CUTS::eRTau1)->size(); jpart++ ){
    if(iitau>5) continue;
    int itau = a->active_part->at(CUTS::eRTau1)->at(jpart);
    if(a->_Tau->minIso.first->at(itau)<0.5){
      continue;
    }
    zTauTree["tau_"+to_string(jpart)+"_pt"]   = a->_Tau->pt(itau);
    zTauTree["tau_"+to_string(jpart)+"_eta"]  = a->_Tau->eta(itau);
    zTauTree["tau_"+to_string(jpart)+"_phi"]  = a->_Tau->phi(itau);
    zTauTree["tau_"+to_string(jpart)+"_charge"]  = a->_Tau->charge(itau);
    zTauTree["cosDphi_"+to_string(jpart)]  = absnormPhi(a->_Tau->phi(itau) - a->_MET->phi());
    zTauTree["mt_tau_"+to_string(jpart)]   = a->calculateLeptonMetMt(a->_Tau->p4(itau));
    zTauTree["tau_"+to_string(jpart)+"_iso"]   = a->_Tau->maxIso.first->at(itau);
    iitau++;
    if(a->_Tau->maxIso.first->at(itau)>0.5)
      isotau++;
    else
      non_isotau++;
  }
  if(isotau==0){
    //we do not care about these:
    return;
  }
  zTauTree["tau_n"]       = iitau;
  zTauTree["tauIso_n"]    = isotau;
  zTauTree["tauNonIso_n"] = non_isotau;
  zTauTree["met"]       = a->_MET->pt();
  
  
  if(a->active_part->at(CUTS::eRMuon1)->size()>0){
    int imuo=a->active_part->at(CUTS::eRMuon1)->at(0);
    zTauTree["muo1_pt"]       = a->_Muon->pt(imuo);
    zTauTree["muo1_eta"]      = a->_Muon->eta(imuo);
    zTauTree["muo1_phi"]      = a->_Muon->phi(imuo);
    zTauTree["muo1_charge"]   = a->_Muon->charge(imuo);
    zTauTree["mt_muo1"]       = a->calculateLeptonMetMt(a->_Muon->p4(imuo));
    
  }else{
    zTauTree["muo1_pt"]       = 0;
    zTauTree["muo1_eta"]      = 0;
    zTauTree["muo1_phi"]      = 0;
    zTauTree["muo1_charge"]   = 0;
    zTauTree["mt_muo1"]       = 0;
  }
  if(a->active_part->at(CUTS::eRElec1)->size()>0){
    int iele=a->active_part->at(CUTS::eRElec1)->at(0);
    zTauTree["ele1_pt"]       = a->_Electron->pt(iele);
    zTauTree["ele1_eta"]      = a->_Electron->eta(iele);
    zTauTree["ele1_phi"]      = a->_Electron->phi(iele);
    zTauTree["ele1_charge"]   = a->_Electron->charge(iele);
    zTauTree["mt_ele1"]       = a->calculateLeptonMetMt(a->_Electron->p4(iele));
  }else{
    zTauTree["ele1_pt"]       = 0;
    zTauTree["ele1_eta"]      = 0;
    zTauTree["ele1_phi"]      = 0;
    zTauTree["ele1_charge"]   = 0;
    zTauTree["mt_ele1"]       = 0;
  }
  
  zTauTree["jet_n"]   = a->active_part->at(CUTS::eRJet1)->size();
  if(j1>=0){
    zTauTree["jet1_pt"]   = a->_Jet->pt(j1);
    zTauTree["jet1_eta"]  = a->_Jet->eta(j1);
    zTauTree["jet1_phi"]  = a->_Jet->phi(j1);
  }else{
    zTauTree["jet1_pt"]   = 0;
    zTauTree["jet1_eta"]  = 10;
    zTauTree["jet1_phi"]  = 10;
  }
  if(j2>=0){
    zTauTree["jet2_pt"]   = a->_Jet->pt(j2);
    zTauTree["jet2_eta"]  = a->_Jet->eta(j2);
    zTauTree["jet2_phi"]  = a->_Jet->phi(j2);
  }else{           
    zTauTree["jet2_pt"]   = 0;
    zTauTree["jet2_eta"]  = 10;
    zTauTree["jet2_phi"]  = 10;
  }
  
  zTauTree["b_jet_n"] = a->active_part->at(CUTS::eRBJet)->size();
  if(a->active_part->at(CUTS::eRBJet)->size()>0){
    int ibjet=a->active_part->at(CUTS::eRBJet)->at(0);
    zTauTree["b_jet_pt"]   = a->_Jet->pt(ibjet);
  }else{
    zTauTree["b_jet_pt"]   = 0;
  }
  zTauTree["jet_mass"]  = mass;
  zTauTree["weight"]    = a->wgt;
  
  HistClass::FillTree("TauIsoTree");
}

void SpechialAnalysis::make_TTLAna() {
  
  
  //for(string i :  a->_Tau->pstats["Tau1"].bset ){
    //cout<<i<<endl;
  //}
  
  
  
  
  
  int j1      = -1;
  int j2      = -1;
  double mass = 0;
  for (auto it : *a->active_part->at(CUTS::eDiJet)) {
    int j1tmp = (it) / a->_Jet->size();
    int j2tmp = (it) % a->_Jet->size();
    if (a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "") > mass) {
      j1   = j1tmp;
      j2   = j2tmp;
      mass = a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "");
    }
  }
  
  int iitau=0;
  int isotau=0;
  int non_isotau=0;
  vector<int> isolated_taus;
  vector<int> non_isolated_taus;
  for(size_t jpart =0; jpart<a->active_part->at(CUTS::eRTau1)->size(); jpart++ ){
    if(iitau>5) continue;
    int itau = a->active_part->at(CUTS::eRTau1)->at(jpart);
    if(a->_Tau->minIso.first->at(itau)<0.5){
      continue;
    }
    
    zTauTree["tau_"+to_string(jpart)+"_pt"]   = a->_Tau->pt(itau);
    zTauTree["tau_"+to_string(jpart)+"_eta"]  = a->_Tau->eta(itau);
    zTauTree["tau_"+to_string(jpart)+"_phi"]  = a->_Tau->phi(itau);
    zTauTree["tau_"+to_string(jpart)+"_charge"]  = a->_Tau->charge(itau);
    zTauTree["cosDphi_"+to_string(jpart)]  = absnormPhi(a->_Tau->phi(itau) - a->_MET->phi());
    zTauTree["mt_tau_"+to_string(jpart)]   = a->calculateLeptonMetMt(a->_Tau->p4(itau));
    zTauTree["tau_"+to_string(jpart)+"_iso"]   = a->_Tau->maxIso.first->at(itau);
    iitau++;
    //cout<<a->_Tau->maxIso.first->at(itau)<<endl;
    if(a->_Tau->maxIso.first->at(itau)>0.5){
      HistClass::Fill(4, (particle_names[2]+"_pt_"+ to_string(isotau)).c_str(),a->_Tau->pt(itau),a->wgt);
      HistClass::Fill(4, (particle_names[2]+"_eta_"+to_string(isotau)).c_str(),a->_Tau->eta(itau),a->wgt);
      HistClass::Fill(4, (particle_names[2]+"_phi_"+to_string(isotau)).c_str(),a->_Tau->phi(itau),a->wgt);
      isolated_taus.push_back(itau);
      isotau++;
    }else{
      HistClass::Fill(5, (particle_names[2]+"_pt_"+ to_string(non_isotau)).c_str(),a->_Tau->pt(itau),a->wgt);
      HistClass::Fill(5, (particle_names[2]+"_eta_"+to_string(non_isotau)).c_str(),a->_Tau->eta(itau),a->wgt);
      HistClass::Fill(5, (particle_names[2]+"_phi_"+to_string(non_isotau)).c_str(),a->_Tau->phi(itau),a->wgt);
      non_isolated_taus.push_back(itau);
      non_isotau++;
    }
  }
  HistClass::Fill(4, (particle_names[3]+"_pt_0").c_str(),a->_MET->pt(),a->wgt);
  HistClass::Fill(4, (particle_names[3]+"_eta_0").c_str(),a->_MET->eta(),a->wgt);
  HistClass::Fill(4, (particle_names[3]+"_phi_0").c_str(),a->_MET->phi(),a->wgt);
  
  HistClass::Fill(3, (particle_names[2]+"_num").c_str(),iitau,a->wgt);
  HistClass::Fill(4, (particle_names[2]+"_num").c_str(),isotau,a->wgt);
  HistClass::Fill(5, (particle_names[2]+"_num").c_str(),non_isotau,a->wgt);
  
  if(isolated_taus.size()==0){
    for(int &itau : non_isolated_taus){
        HistClass::Fill(7, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(7, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(7, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
  }
  if(isolated_taus.size()==1){
    for(int &itau : isolated_taus){
        HistClass::Fill(8, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(8, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(8, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
    for(int &itau : non_isolated_taus){
        HistClass::Fill(9, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(9, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(9, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
  }
  if(isolated_taus.size()==2){
    for(int &itau : isolated_taus){
        HistClass::Fill(10, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(10, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(10, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
    for(int &itau : non_isolated_taus){
        HistClass::Fill(11, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(11, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(11, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
  }
  if(isolated_taus.size()>2){
    for(int &itau : isolated_taus){
        HistClass::Fill(12, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(12, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(12, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
    for(int &itau : non_isolated_taus){
        HistClass::Fill(13, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(13, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(13, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
    }
  }
 
  
  //zTauTree["tau_n"]       = iitau;
  //zTauTree["tauIso_n"]    = isotau;
  //zTauTree["tauNonIso_n"] = non_isotau;
  //zTauTree["met"]       = a->_MET->pt();
  
  
  //zTauTree["b_jet_n"] = a->active_part->at(CUTS::eRBJet)->size();
  //if(a->active_part->at(CUTS::eRBJet)->size()>0){
    //int ibjet=a->active_part->at(CUTS::eRBJet)->at(0);
    //zTauTree["b_jet_pt"]   = a->_Jet->pt(ibjet);
  //}else{
    //zTauTree["b_jet_pt"]   = 0;
  //}
  //zTauTree["jet_mass"]  = mass;
  //zTauTree["weight"]    = a->wgt;
  
  //HistClass::FillTree("TauIsoTree");
}

void SpechialAnalysis::fill_particle(int stage) {
  
  
  for (unsigned int ipart = 0; ipart < 3; ipart++) {
    HistClass::Fill(stage, (particle_names[ipart]+"_num").c_str(),a->active_part->at(cuts[ipart])->size(),a->wgt);
    for( size_t jpart =0; jpart<a->active_part->at(cuts[ipart])->size(); jpart++ ){
      if(jpart>5) continue;
      int iipart=a->active_part->at(cuts[ipart])->at(jpart);
      HistClass::Fill(stage, (particle_names[ipart]+"_pt_"+to_string(jpart)).c_str(),particles[ipart]->pt(iipart),a->wgt);
      HistClass::Fill(stage, (particle_names[ipart]+"_eta_"+to_string(jpart)).c_str(),particles[ipart]->eta(iipart),a->wgt);
      HistClass::Fill(stage, (particle_names[ipart]+"_phi_"+to_string(jpart)).c_str(),particles[ipart]->phi(iipart),a->wgt);
    }
  }
  HistClass::Fill(stage, (particle_names[3]+"_pt_0").c_str(),a->_MET->pt(),a->wgt);
  HistClass::Fill(stage, (particle_names[3]+"_eta_0").c_str(),a->_MET->eta(),a->wgt);
  HistClass::Fill(stage, (particle_names[3]+"_phi_0").c_str(),a->_MET->phi(),a->wgt);
  
}
void SpechialAnalysis::end_run() {
  
  TFile* outfile = new TFile(a->histo.outname.c_str(), "UPDATE", a->histo.outname.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));
  
  outfile->cd();
  outfile->mkdir("Spechial");
  for(int i=0; i<13; i++){
    outfile->cd();
    outfile->mkdir(("Spechial/Stage_"+to_string(i)).c_str());
    outfile->cd(("Spechial/Stage_"+to_string(i)).c_str());
    HistClass::WriteAll(("h1_"+to_string(i)+"_").c_str());
  }
  outfile->cd();
  outfile->cd("Spechial");
  HistClass::WriteAllTrees();
  outfile->Close();
}

