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
  
  HistClass::CreateHisto(14,"Tau_pt_pt", 1000,  0, 1000, 1000, 0, 1000, "tau 1 pt", "tau 2 pt");

  

  zTauTreeV["tau_pt"]       = new vector<float>;
  zTauTreeV["tau_lpt"]       = new vector<float>;
  zTauTreeV["tau_eta"]      = new vector<float>;
  zTauTreeV["tau_phi"]      = new vector<float>;
  zTauTreeV["tau_charge"]   = new vector<float>;
  zTauTreeV["tau_cosDphi"]      = new vector<float>;
  zTauTreeV["tau_mt"]       = new vector<float>;
  zTauTreeV["tau_iso"]      = new vector<float>;
  zTauTreeV["tau_dz"]      = new vector<float>;
  zTauTreeV["tau_mass"]      = new vector<float>;
  
  //zTauTreeV["tau_noiso_pt"]       = new vector<float>;
  //zTauTreeV["tau_noiso_eta"]      = new vector<float>;
  //zTauTreeV["tau_noiso_phi"]      = new vector<float>;
  //zTauTreeV["tau_noiso_charge"]   = new vector<float>;
  //zTauTreeV["tau_noiso_cosDphi"]      = new vector<float>;
  //zTauTreeV["tau_noiso_mt"]       = new vector<float>;
  //zTauTreeV["tau_noiso_iso"]      = new vector<float>;
  
  
  //zTauTree["tau_n"]          =0;
  //zTauTree["tauIso_n"]       =0;
  //zTauTree["tauNonIso_n"]    =0;
  
  zTauTree["met"]           = 0;
  zTauTree["met_phi"]           = 0;
  //zTauTree["muo1_pt"]       = 0;
  //zTauTree["muo1_eta"]      = 0;
  //zTauTree["muo1_phi"]      = 0;
  //zTauTree["muo1_charge"]   = 0;
  //zTauTree["mt_muo1"]       = 0;
  //zTauTree["ele1_pt"]       = 0;
  //zTauTree["ele1_eta"]      = 0;
  //zTauTree["ele1_phi"]      = 0;
  //zTauTree["ele1_charge"]   = 0;
  //zTauTree["mt_ele1"]       = 0;
  //zTauTree["jet_n"]         = 0;
  //zTauTree["jet1_pt"]       = 0;
  //zTauTree["jet1_eta"]      = 10;
  //zTauTree["jet1_phi"]      = 10;
  //zTauTree["jet2_pt"]       = 0;
  //zTauTree["jet2_eta"]      = 10;
  //zTauTree["jet2_phi"]      = 10;
  //zTauTree["b_jet_n"]       = 0;
  //zTauTree["b_jet_pt"]      = 0;
  //zTauTree["jet_mass"]      = 0;
  zTauTree["weight"]        = 0;
  
  HistClass::CreateTree(&zTauTree,&zTauTreeV,"TauIsoTree");
  
  TFile ttl_file("Pileup/ttl_rate.root","READ");
  ttl_rato_one = (TH2D*) ttl_file.Get("one_fake");
  ttl_rato_two = (TH2D*) ttl_file.Get("two_fake");
  vector<string> cr_variables;
  sphisto = Histogramer(1, a->filespace+"Hist_entries.in", a->filespace+"Cuts.in", a->histo.outname, a->isData, cr_variables);
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
  
  //for(int &itau : *(a->active_part->at(CUTS::eRTau1))){
    //if(!a->_Tau->decayModeFinding->at(itau)){
        //cout<<"strange event"<<endl;
        //cout<<a->_Tau->decayMode->at(itau)<<"  "<<a->_Tau->nProngs->at(itau)<<"   "<<   a->_Tau->decayModeFindingNewDMs->at(itau) <<endl;
    //}
  //}
  
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

  double backup_wgt=a->wgt;
  const vector<string>* groups = sphisto.get_groups();
  
  vector<int> jet1_backup;
  vector<int> jet2_backup;
  for(int ijet : *a->active_part->at(CUTS::eRJet1)){
    jet1_backup.push_back(ijet);
  }
  for(int ijet : *a->active_part->at(CUTS::eRJet2)){
    jet2_backup.push_back(ijet);
  }
  
  isobool =find(a->_Tau->pstats["Tau1"].bset.begin(),a->_Tau->pstats["Tau1"].bset.end(),"DoDiscrByIsolation");
  if(isobool!=a->_Tau->pstats["Tau1"].bset.end()){
    a->_Tau->pstats["Tau1"].bset.erase(isobool);
  }
  a->active_part->at(CUTS::eRTau1)->clear();
  a->getGoodRecoLeptons(*a->_Tau, CUTS::eRTau1, CUTS::eGTau, a->_Tau->pstats["Tau1"],0);
  //restore current version
  a->_Tau->pstats["Tau1"].bset.push_back("DoDiscrByIsolation");
  
  isobool =find(a->_Tau->pstats["Tau2"].bset.begin(),a->_Tau->pstats["Tau2"].bset.end(),"DoDiscrByIsolation");
  if(isobool!=a->_Tau->pstats["Tau2"].bset.end()){
    a->_Tau->pstats["Tau2"].bset.erase(isobool);
  }
  a->active_part->at(CUTS::eRTau2)->clear();
  a->getGoodRecoLeptons(*a->_Tau, CUTS::eRTau2, CUTS::eGTau, a->_Tau->pstats["Tau2"],0);
  //restore current version
  a->_Tau->pstats["Tau2"].bset.push_back("DoDiscrByIsolation");

  a->active_part->at(CUTS::eDiTau)->clear();
  int tau_combo_type=getGoodLeptonCombosForQCD(*a->_Tau, *a->_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, a->distats["DiTau"],0);
  if(tau_combo_type==0 or a->active_part->at(CUTS::eDiTau)->size()!=1){
    return;
  }
  int combindex=a->active_part->at(CUTS::eDiTau)->at(0);
  int p1= (combindex) / BIG_NUM;
  int p2= (combindex) % BIG_NUM;
  double ttl_weight=1.;
  double pt1=a->_Tau->pt(p1);
  double pt2=a->_Tau->pt(p2);
  if(pt2>pt1){
    double tmp=pt1;
    pt1=pt2;
    pt2=tmp;
  }
  if(tau_combo_type==1){
    ttl_weight=ttl_rato_one->GetBinContent(ttl_rato_one->FindBin(pt1,pt2));
  }
  if(tau_combo_type==2){
    ttl_weight=ttl_rato_two->GetBinContent(ttl_rato_two->FindBin(pt1,pt2));
  }
  a->wgt*= ttl_weight;
  
  a->active_part->at(CUTS::eRJet1)->clear();
  for(int ijet : jet1_backup){
    if(a->_Tau->p4(p1).DeltaR(a->_Jet->p4(ijet))>0.3){
      a->active_part->at(CUTS::eRJet1)->push_back(ijet);
    }
  }
  a->active_part->at(CUTS::eRJet2)->clear();
  for(int ijet : jet2_backup){
    if(a->_Tau->p4(p2).DeltaR(a->_Jet->p4(ijet))>0.3){
      a->active_part->at(CUTS::eRJet2)->push_back(ijet);
    }
  }
  
  a->fillCuts(false);
  for(auto it: *groups) {
    a->fill_Folder(it, a->maxCut, sphisto, false);
  }
  
  a->wgt=backup_wgt;
  for(Particle* ipart: a->allParticles) ipart->setCurrentP(0);
  a->_MET->setCurrentP(0);
  a->active_part = &(a->goodParts);
  //make_tau_tree();
  
  a->active_part->at(CUTS::eRJet1)->clear();
  for(int ijet : jet1_backup){
    a->active_part->at(CUTS::eRJet1)->push_back(ijet);
  }
  a->active_part->at(CUTS::eRJet2)->clear();
  for(int ijet : jet2_backup){
    a->active_part->at(CUTS::eRJet2)->push_back(ijet);
  }
  

}

void SpechialAnalysis::make_tau_tree() {
  
  
  //int j1      = -1;
  //int j2      = -1;
  //double mass = 0;
  //for (auto it : *a->active_part->at(CUTS::eDiJet)) {
    //int j1tmp = (it) / a->_Jet->size();
    //int j2tmp = (it) % a->_Jet->size();
    //if (a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "") > mass) {
      //j1   = j1tmp;
      //j2   = j2tmp;
      //mass = a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "");
    //}
  //}
  
  for(auto &myv : zTauTreeV){
    myv.second->clear();
  }
  vector<int> selTaus;
  for(size_t itau =0; itau<a->_Tau->size(); itau++ ){
    if(a->_Tau->leadChargedCandPt->at(itau)<70 or a->_Tau->chargedIsoPtSum->at(itau)/a->_Tau->leadChargedCandPt->at(itau)>0.03 or a->_Tau->leadChargedCandDz_pv->at(itau)>0.01 ){
      continue;
    }
    if(a->isOverlaping(a->_Tau->p4(itau), *a->_Muon, CUTS::eRMuon1, 0.3)){
      continue;
    }
    if(a->isOverlaping(a->_Tau->p4(itau), *a->_Electron, CUTS::eRElec1, 0.3)){
      continue;
    }
    zTauTreeV["tau_pt"]->push_back(a->_Tau->pt(itau));
    zTauTreeV["tau_lpt"]->push_back(a->_Tau->leadChargedCandPt->at(itau));
    zTauTreeV["tau_eta"]->push_back(a->_Tau->eta(itau));
    zTauTreeV["tau_phi"]->push_back(a->_Tau->phi(itau));
    zTauTreeV["tau_charge"]->push_back(a->_Tau->charge(itau));
    zTauTreeV["tau_cosDphi"]->push_back(absnormPhi(a->_Tau->phi(itau) - a->_MET->phi()));
    zTauTreeV["tau_mt"]->push_back(a->calculateLeptonMetMt(a->_Tau->p4(itau)));
    zTauTreeV["tau_chiso"]->push_back(a->_Tau->chargedIsoPtSum->at(itau));
    zTauTreeV["tau_neutiso"]->push_back(a->_Tau->neutralIsoPtSum->at(itau));
    zTauTreeV["tau_puiso"]->push_back(a->_Tau->puCorrPtSum->at(itau));
    
    
    zTauTreeV["tau_dz"]->push_back(a->_Tau->leadChargedCandDz_pv->at(itau));
    selTaus.push_back(itau);
    
  }
  //if( (zTauTreeV["tau_pt"]->size() +zTauTreeV["tau_noiso_pt"]->size()) ==2){
    //////we do not care about these:
    //return;
  //}
  if( (zTauTreeV["tau_pt"]->size() <2) or (zTauTreeV["tau_pt"]->size() >4) ){
    return;
  }
  for(int itau : selTaus){
    for(int jtau : selTaus){
      if(itau>=jtau){
        continue;
      }
      zTauTreeV["tau_mass"]->push_back((a->_Tau->p4(itau)+a->_Tau->p4(jtau)).M());
    }
  }
  //zTauTree["tau_n"]       = iitau;
  //zTauTree["tauIso_n"]    = isotau;
  //zTauTree["tauNonIso_n"] = non_isotau;
  zTauTree["met"]       = a->_MET->pt();
  zTauTree["met_phi"]       = a->_MET->phi();
  
  
  //if(a->active_part->at(CUTS::eRMuon1)->size()>0){
    //int imuo=a->active_part->at(CUTS::eRMuon1)->at(0);
    //zTauTree["muo1_pt"]       = a->_Muon->pt(imuo);
    //zTauTree["muo1_eta"]      = a->_Muon->eta(imuo);
    //zTauTree["muo1_phi"]      = a->_Muon->phi(imuo);
    //zTauTree["muo1_charge"]   = a->_Muon->charge(imuo);
    //zTauTree["mt_muo1"]       = a->calculateLeptonMetMt(a->_Muon->p4(imuo));
    
  //}else{
    //zTauTree["muo1_pt"]       = 0;
    //zTauTree["muo1_eta"]      = 0;
    //zTauTree["muo1_phi"]      = 0;
    //zTauTree["muo1_charge"]   = 0;
    //zTauTree["mt_muo1"]       = 0;
  //}
  //if(a->active_part->at(CUTS::eRElec1)->size()>0){
    //int iele=a->active_part->at(CUTS::eRElec1)->at(0);
    //zTauTree["ele1_pt"]       = a->_Electron->pt(iele);
    //zTauTree["ele1_eta"]      = a->_Electron->eta(iele);
    //zTauTree["ele1_phi"]      = a->_Electron->phi(iele);
    //zTauTree["ele1_charge"]   = a->_Electron->charge(iele);
    //zTauTree["mt_ele1"]       = a->calculateLeptonMetMt(a->_Electron->p4(iele));
  //}else{
    //zTauTree["ele1_pt"]       = 0;
    //zTauTree["ele1_eta"]      = 0;
    //zTauTree["ele1_phi"]      = 0;
    //zTauTree["ele1_charge"]   = 0;
    //zTauTree["mt_ele1"]       = 0;
  //}
  
  //zTauTree["jet_n"]   = a->active_part->at(CUTS::eRJet1)->size();
  //if(j1>=0){
    //zTauTree["jet1_pt"]   = a->_Jet->pt(j1);
    //zTauTree["jet1_eta"]  = a->_Jet->eta(j1);
    //zTauTree["jet1_phi"]  = a->_Jet->phi(j1);
  //}else{
    //zTauTree["jet1_pt"]   = 0;
    //zTauTree["jet1_eta"]  = 10;
    //zTauTree["jet1_phi"]  = 10;
  //}
  //if(j2>=0){
    //zTauTree["jet2_pt"]   = a->_Jet->pt(j2);
    //zTauTree["jet2_eta"]  = a->_Jet->eta(j2);
    //zTauTree["jet2_phi"]  = a->_Jet->phi(j2);
  //}else{           
    //zTauTree["jet2_pt"]   = 0;
    //zTauTree["jet2_eta"]  = 10;
    //zTauTree["jet2_phi"]  = 10;
  //}
  
  //zTauTree["b_jet_n"] = a->active_part->at(CUTS::eRBJet)->size();
  //if(a->active_part->at(CUTS::eRBJet)->size()>0){
    //int ibjet=a->active_part->at(CUTS::eRBJet)->at(0);
    //zTauTree["b_jet_pt"]   = a->_Jet->pt(ibjet);
  //}else{
    //zTauTree["b_jet_pt"]   = 0;
  //}
  //zTauTree["jet_mass"]  = mass;
  zTauTree["weight"]    = a->wgt;
  
  HistClass::FillTree("TauIsoTree");
}

void SpechialAnalysis::make_TTLAna() {
  
  
  //for(string i :  a->_Tau->pstats["Tau1"].bset ){
    //cout<<i<<endl;
  //}
  
  
  
  
  
  //int j1      = -1;
  //int j2      = -1;
  //double mass = 0;
  //for (auto it : *a->active_part->at(CUTS::eDiJet)) {
    //int j1tmp = (it) / a->_Jet->size();
    //int j2tmp = (it) % a->_Jet->size();
    //if (a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "") > mass) {
      //j1   = j1tmp;
      //j2   = j2tmp;
      //mass = a->diParticleMass(a->_Jet->p4(j1tmp), a->_Jet->p4(j2tmp), "");
    //}
  //}
  //if(not a->_MET->pt()<80){
  if(a->passCutRange(a->_MET->pt(), a->distats["Run"].pmap.at("MetCut")))
    return;
  
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
      
      if(a->_Tau->decayMode->at(itau)==1 and !a->_Tau->decayModeFinding->at(itau)){
        
        cout<<"strange event"<<endl;
      }
      
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
  
  vector<int> used_tau;
  if(isolated_taus.size()==0 && non_isolated_taus.size()==2){
    for(int &itau : non_isolated_taus){
        HistClass::Fill(7, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(7, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(7, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    if(a->_Tau->pt(used_tau[0])>a->_Tau->pt(used_tau[1])){
      HistClass::Fill(7, "Tau_pt_pt",a->_Tau->pt(used_tau[0]),a->_Tau->pt(used_tau[1]),a->wgt);
    }else{
      HistClass::Fill(7, "Tau_pt_pt",a->_Tau->pt(used_tau[1]),a->_Tau->pt(used_tau[0]),a->wgt);
    }
  }
  
  
  if(isolated_taus.size()==1 && non_isolated_taus.size()==1){
    for(int &itau : isolated_taus){
        HistClass::Fill(8, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(8, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(8, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    for(int &itau : non_isolated_taus){
        HistClass::Fill(9, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(9, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(9, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    if(a->_Tau->pt(used_tau[0])>a->_Tau->pt(used_tau[1])){
      HistClass::Fill(8, "Tau_pt_pt",a->_Tau->pt(used_tau[0]),a->_Tau->pt(used_tau[1]),a->wgt);
    }else{
      HistClass::Fill(8, "Tau_pt_pt",a->_Tau->pt(used_tau[1]),a->_Tau->pt(used_tau[0]),a->wgt);
    }
  }
  if(isolated_taus.size()==2 && non_isolated_taus.size()==0){
    for(int &itau : isolated_taus){
        HistClass::Fill(10, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(10, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(10, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    for(int &itau : non_isolated_taus){
        HistClass::Fill(11, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(11, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(11, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    if(a->_Tau->pt(used_tau[0])>a->_Tau->pt(used_tau[1])){
      HistClass::Fill(10, "Tau_pt_pt",a->_Tau->pt(used_tau[0]),a->_Tau->pt(used_tau[1]),a->wgt);
    }else{
      HistClass::Fill(10, "Tau_pt_pt",a->_Tau->pt(used_tau[1]),a->_Tau->pt(used_tau[0]),a->wgt);
    }
  }
  if(non_isolated_taus.size()+isolated_taus.size()!=2){
    for(int &itau : isolated_taus){
        HistClass::Fill(12, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(12, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(12, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    for(int &itau : non_isolated_taus){
        HistClass::Fill(13, (particle_names[2]+"_pt_0").c_str(),a->_Tau->pt(itau),a->wgt);
        HistClass::Fill(13, (particle_names[2]+"_eta_0").c_str(),a->_Tau->eta(itau),a->wgt);
        HistClass::Fill(13, (particle_names[2]+"_phi_0").c_str(),a->_Tau->phi(itau),a->wgt);
        used_tau.push_back(itau);
    }
    if(used_tau.size()>1){
      HistClass::Fill(12, "Tau_pt_pt",a->_Tau->pt(used_tau[0]),a->_Tau->pt(used_tau[1]),a->wgt);
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
  sphisto.fill_histogram("QCD");
  
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
  outfile->cd();
  outfile->cd("Spechial");
  outfile->mkdir("Spechial/2D");
  outfile->cd("Spechial/2D");
  HistClass::WriteAll2();
  outfile->Close();
  
}



/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
///like in the analyzer, but with non isolated leptons so we will build the combinations here
///it will return the type of pair that was found
int SpechialAnalysis::getGoodLeptonCombosForQCD(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst) {
  if(! a->neededCuts.isPresent(ePosFin)) return 0;
  string systname = a->syst_names.at(syst);

  if(!lep1.needSyst(syst) && !lep2.needSyst(syst)) {
    a->active_part->at(ePosFin)=a->goodParts[ePosFin];
    return 0;
  }

  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2;
  
  bool two_non_isolated_taus=false;
  bool one_non_isolated_tau=false;

  for(auto i1 : *a->active_part->at(ePos1)) {
    part1 = lep1.p4(i1);
    for(auto i2 : *a->active_part->at(ePos2)) {
      if(sameParticle && i2 <= i1) continue;
      // check for isolation here
      if(a->_Tau->minIso.first->at(i1)<0.5 or a->_Tau->minIso.first->at(i2)<0.5){
        continue;
      }
      if( (a->_Tau->maxIso.first->at(i1)>0.5 and a->_Tau->maxIso.first->at(i2)<0.5) or 
          (a->_Tau->maxIso.first->at(i2)>0.5 and a->_Tau->maxIso.first->at(i1)<0.5)){
        one_non_isolated_tau=true;
      }
      if( (a->_Tau->maxIso.first->at(i1)<0.5 and a->_Tau->maxIso.first->at(i2)<0.5) ){
        two_non_isolated_taus=true;
      }
      if( !( one_non_isolated_tau or two_non_isolated_taus ) ){
        continue;
      }
      
      
      part2 = lep2.p4(i2);
      bool passCuts = true;
      for(auto cut : stats.bset) {
        if(!passCuts) break;
        else if (cut == "DiscrByDeltaR") passCuts = passCuts && (part1.DeltaR(part2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && a->passCutRange(cos(absnormPhi(part1.Phi() - part2.Phi())), stats.pmap.at("CosDphiCut"));
        else if(cut == "DiscrByDeltaPt") passCuts = passCuts && a->passCutRange(part1.Pt() - part2.Pt(), stats.pmap.at("DeltaPtCutValue"));
        else if(cut == "DiscrByCDFzeta2D") {
          double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * a->getPZeta(part1, part2).first
            + stats.dmap.at("PZetaVisCutCoefficient") * a->getPZeta(part1, part2).second;
          passCuts = passCuts && a->passCutRange(CDFzeta, stats.pmap.at("CDFzeta2DCutValue"));
        }
        else if(cut == "DiscrByDeltaPtDivSumPt") {
          double ptDiv = (part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt());
          passCuts = passCuts && a->passCutRange(ptDiv, stats.pmap.at("DeltaPtDivSumPtCutValue"));
        }
        else if (cut == "DiscrByMassReco") {
          double diMass = a->diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
          passCuts = passCuts && a->passCutRange(diMass, stats.pmap.at("MassCut"));
        }
        else if(cut == "DiscrByCosDphiPtAndMet"){
          double CosDPhi1 = cos(absnormPhi(part1.Phi() - a->_MET->phi()));
          passCuts = passCuts && a->passCutRange(CosDPhi1, stats.pmap.at("CosDphiPtAndMetCut"));
        }


        else cout << "cut: " << cut << " not listed" << endl;
      }
      if (stats.bfind("DiscrByOSLSType")){
        //   if it is 1 or 0 it will end up in the bool map!!
        if(stats.bfind("DiscrByOSLSType") && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
      }else if (stats.dmap.find("DiscrByOSLSType") != stats.dmap.end() ){
        if(lep1.charge(i1) * lep2.charge(i2) > 0) continue;
      }else if (stats.smap.find("DiscrByOSLSType") != stats.smap.end() ){
        if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
        else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge(i1) * lep2.charge(i2) >= 0)) continue;
      }

      ///Particles that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2
      if(passCuts)
        a->active_part->at(ePosFin)->push_back(i1*BIG_NUM + i2);
    }
  }
  
  if(one_non_isolated_tau)
    return 1;
  if(two_non_isolated_taus)
    return 2;
  return 0;
}
