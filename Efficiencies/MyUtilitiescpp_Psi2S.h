#ifndef MYUTILITIESCPP_PSI2S_H
#define MYUTILITIESCPP_PSI2S_H

#include <RooPlot.h>
#include "RooFit.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"

using namespace RooFit;
using namespace std;

//[1.1, 2.2, 2.98, 3.21, 3.59, 3.782, 4.18, 4.8]
//[1.2100000000000002, 4.840000000000001, 8.8804, 10.3041, 12.8881, 14.303524, 17.472399999999997, 23.04]
std::vector<double> Psi2S_Effy_Pre(UInt_t year=2022, UInt_t sample=7)
 {
  //**************************************************
  // prefilter Effy
  //**************************************************
  //MC with NO filter cuts
  TString MCGenModel = "PHP";
  TChain tree("ntuple");
  if (year==2022 && sample==6){
    tree.Add("../MC_Psi2S_2022/*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==7){
    tree.Add("/home/oscar/Documents/Thesis/Efficencies/MC/Psi2S/MC_Psi2SK_OnlyGen.root");
    MCGenModel = "SVS";
  }
  
  //Datframe with NO filter
  ROOT::RDataFrame DFPrei(tree);
  //DFPreI.Describe().Print();
  auto nentriesPreI = DFPrei.Count();
  cout << " Total entries in MC with NO filter " << *nentriesPreI <<endl;

  //23, 3.59, 3.782
  auto DFPreI = DFPrei.Define("massJ","gen_jpsi_p4.M()");
  auto hRes_pre1 = DFPreI.Histo1D( {"massJ", "", 23, 3.59, 3.782},"massJ");
  //return;
  
  //Datframe with filter cuts by hand
  auto DFPreF = DFPreI.Filter("gen_muon1_p4.Pt()>3.5 && gen_muon2_p4.Pt()>3.5 && fabs(gen_muon1_p4.Eta())<=2.5 && fabs(gen_muon2_p4.Eta())<=2.5 && gen_pion3_p4.Pt()>0.5 && fabs(gen_pion3_p4.Eta())<=2.5", "selfilters");
  auto nentriesPreF = DFPreF.Count();
  cout << " Total entries in MC with filters by hand " << *nentriesPreF <<endl;

  //auto hRes_pre2 = DFPreF.Histo1D( {"massJ", "", 76, 1.1, 4.8},"massJ");
  auto hRes_pre2 = DFPreF.Histo1D( {"massJ", "", 23, 3.59, 3.782},"massJ");

  //Divide
  auto hRes_preD = new TH1F("hRes_preD","hRes_preD",23, 3.59, 3.782);
  hRes_preD->Divide( &(*hRes_pre2), &(*hRes_pre1),1,1);

  TCanvas* canv_histopre = CreateCanvasth1("canv_histopre", hRes_preD, "Pre-filter efficiency", "m_{#mu#mu}"); 
  canv_histopre->SaveAs("plotsEffy/Resonant_prefilter_Psi2S.png");

  //double preEffyVal = (double)(*nentriesPreF) / (double)(*nentriesPreI);
  //double preEffyVal = static_cast< double >(*nentriesPreF) / (*nentriesPreI);
  double EFi = (double)(*nentriesPreI);
  double EFie = sqrt(EFi);
  double EFf = (double)(*nentriesPreF);
  double EFfe = sqrt(EFf);
  double preEffyVal = EFf/EFi;
  double preEffyVale = (EFf/EFi)*sqrt( (EFie/EFi)*(EFie/EFi) + (EFfe/EFf)*(EFfe/EFf) );
  //return preEffyVal;

  std::vector<double> myeffv;
  myeffv.push_back( preEffyVal );
  myeffv.push_back( preEffyVale );
  //std::cout << "Effy Prefilter = " << preEffyVal << ", +/- " << preEffyVale << std::endl;
  //std::cout << " " << std::endl;
  return myeffv;

}

std::vector<double> Psi2S_Effy_Reco(UInt_t year=2022, UInt_t sample=2)
 {
  //**************************************************
  // Reco Effy
  //**************************************************
  //MC with filter cuts
  TString MCGenModel = "PHP";
  TChain tree("ntuple");
  if (year==2022 && sample==6){
    tree.Add("../MC_Psi2S_2022/*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==7){
    tree.Add("/home/oscar/Documents/Thesis/Efficencies/MC/Psi2S/MC_Psi2SK_RecoGen.root");
    MCGenModel = "SVS";
  }

  //Datframe with filter cuts
  ROOT::RDataFrame DFRecii(tree);
  auto nentriesRecII = DFRecii.Count();
  cout << " Total entries in MC with filters " << *nentriesRecII << endl;

  //Datframe with filter and cuts by hand
  auto DFReci = DFRecii.Filter("gen_muon1_p4.Pt()>3.5 && gen_muon2_p4.Pt()>3.5 && fabs(gen_muon1_p4.Eta())<=2.5 && fabs(gen_muon2_p4.Eta())<=2.5 && gen_pion3_p4.Pt()>0.5 && fabs(gen_pion3_p4.Eta())<=2.5", "selfilters");
  auto nentriesRecI = DFReci.Count();
  cout << " Total entries in MC with filters by hand " << *nentriesRecI <<endl;

  auto DFRecI = DFReci.Define("massJ","gen_jpsi_p4.M()");
  auto hRes_rec1 = DFRecI.Histo1D( {"massJ", "", 23, 3.59, 3.782},"massJ");
  //return;


  //https://root.cern/doc/v626/namespaceROOT_1_1RDF.html
  //auto DFRecFi = ROOT::RDF::MakeCsvDataFrame(Form("../MyReweiting/results_v5/weighted_MCntuple_Jpsiv0_FromSplot_Year%1i_Sample%1i.csv",year,sample));
  //auto DFRecF = DFRecFi.Filter("massB > 5.0 && massB < 5.700 && massJ>3.59 && massJ<3.782 && massv0>0.4876 && massv0<0.5076", "selfilters");
  //auto DFRecF = ROOT::RDF::MakeCsvDataFrame(Form("../MyReweiting/MCntuple_Psi2Sv0_FromSplot_Year%1i_Sample%1i.csv",year,sample));

  TChain tree2("treeBu");
  //TChain tree2("ntuple");


  if (year==2022 && sample==3){
    tree2.Add("../MC_NOres2022/Rootuple_Bdtomumupipi_MCNoRes2022_MCRecoGen_bacht*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==7){
    tree2.Add("/home/oscar/Documents/Thesis/Efficencies/MC/Psi2S/MC_Psi2SK.root");
    MCGenModel = "SVS";
  }


  ROOT::RDataFrame DFRecF(tree2);  
  
  //auto DFRecF = ROOT::RDF::MakeCsvDataFrame(Form("/cms/home/hcrotte/Reweighting/Results_Dec23_JPsi/MCntuple_Jpsiv0_FromSplot_Year%1i_Sample%1i_BDT.csv",year,sample));
  auto nentriesRecF = DFRecF.Count();
  cout << " Total entries in MC with filters " << *nentriesRecF <<endl;
  
  auto DFRecF2 = DFRecF.Filter("XGB_auc > 0.9983-0.003", "selfilters2");
  auto nentriesRecF2 = DFRecF2.Count();
  cout << " Total entries in MC with Filters and XGB cut " << *nentriesRecF2 <<endl;

  //auto hRes_rec2 = DFRecF.Histo1D( {"massJ", "", 76, 1.1, 4.8},"massJ");
  //auto hRes_rec2 = DFRecF.Histo1D( {"massJ", "", 23, 3.59, 3.782},"massJ");
  // ************************************************
  // *** be careful in case weights become needed ***
  // ************************************************
  // *** https://root.cern.ch/doc/master/df005__fillAnyObject_8C.html ***
  // *** https://root-forum.cern.ch/t/rdataframe-with-weight/43992  ***
  // *** https://root-forum.cern.ch/t/rdataframe-histo1d-multiple-weight-branches/31288  ***
  auto hRes_rec2 = DFRecF2.Histo1D( {"massJ", "", 23,  3.59, 3.782},"massJ");//Fill with weights
  //auto hRes_rec2 = DFRecF.Define("w","splotW*5.0").Histo1D( {"massJ", "", 23,  3.59, 3.782},"massJ","w");//Fill with weights

  //Divide
  auto hRes_recD = new TH1F("hRes_recD","hRes_recD",23, 3.59, 3.782);
  hRes_recD->Divide( &(*hRes_rec2), &(*hRes_rec1),1,1);

  TCanvas* canv_historec = CreateCanvasth1("canv_historec", hRes_recD, "Reco efficiency", "m_{#mu#mu}"); 
  canv_historec->SaveAs("plotsEffy/Resonant_recor_Psi2S.png");
  //double recEffyVal = static_cast< double >(*nentriesRecF) / (*nentriesRecI);
  //return recEffyVal;

  double EFi = (double)(*nentriesRecI);
  double EFie = sqrt(EFi);
  //double EFf = (double)(*nentriesRecF);
  cout << " Total entries in MC with filters and Full selection " << *nentriesRecF <<endl;
  cout << " Total entries in MC with filters and Full selection TH1 " << hRes_rec2->Integral() <<endl;
  double EFf = (double)(hRes_rec2->Integral());
  double EFfe = sqrt(EFf);
  double recEffyVal = EFf/EFi;
  double recEffyVale = (EFf/EFi)*sqrt( (EFie/EFi)*(EFie/EFi) + (EFfe/EFf)*(EFfe/EFf) );

  std::vector<double> myeffv;
  myeffv.push_back( recEffyVal );
  myeffv.push_back( recEffyVale );
  //std::cout << "Effy Recfilter = " << recEffyVal << ", +/- " << recEffyVale << std::endl;
  //std::cout << " " << std::endl;
  return myeffv;
  
}

#endif
