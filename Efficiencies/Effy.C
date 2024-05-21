#include <iostream>
#include <string>
#include <TMath.h>
#include <math.h>
#include <Math/Vector4D.h>
#include "Math/GenVector/Boost.h"

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooBernstein.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGenericPdf.h"

#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"

//#include <TError.h>
#include <ROOT/RDataFrame.hxx>
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"

//#include "Utilities.h"
#include "MyUtilitiescpp.h"
#include "MyUtilitiescpp_Psi2S.h"

using namespace std;
using namespace ROOT;
using namespace RooFit;

const Int_t Nbin = 9;  
//Double_t bindiv[Nbin+1] = {1.1, 2.2, 2.98, 3.210, 3.590, 3.782, 4.18, 4.8};
Double_t bindiv[Nbin+1] = {1.1, 2.0, 4.3, 8.68, 10.09, 12.86, 14.18, 16, 18, 22};

TH1* NoRes_Effy_Pre(UInt_t year=2022, UInt_t sample=5)
 {
  //**************************************************
  // prefilter Effy
  //**************************************************

  //MC with NO filter cuts
  TString MCGenModel = "PHP";
  TChain tree("ntuple");
  if (year==2022 && sample==3){
    tree.Add("../MC_NOres2022/Rootuple_Bdtomumupipi_MCNoRes2022_OnlyGen_bacht1.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==5){
    tree.Add("/home/oscar/Documents/Thesis/Efficencies/MC/MC_uuK_OnlyGen.root");
    MCGenModel = "SVS";
  }
  
  //Datframe with NO filter
  ROOT::RDataFrame DFPrei(tree);
  //DFPreI.Describe().Print();
  auto nentriesPreI = DFPrei.Count();
  cout << " Total entries in MC with NO filter " << *nentriesPreI <<endl;

  auto DFPreI = DFPrei.Define("massJ","gen_jpsi_p4.M()").Define("massJ2","gen_jpsi_p4.M()*gen_jpsi_p4.M()");
  //auto hNoRes_pre1 = DFPreI.Histo1D( {"massJ", "", Nbin, bindiv},"massJ");
  auto hNoRes_pre1 = DFPreI.Histo1D( {"massJ2", "", Nbin, bindiv},"massJ2");
  //return;
  
  //Datframe with filter cuts by hand
  auto DFPreF = DFPreI.Filter("gen_muon1_p4.Pt()>3.5 && gen_muon2_p4.Pt()>3.5 && fabs(gen_muon1_p4.Eta())<=2.5 && fabs(gen_muon2_p4.Eta())<=2.5 && gen_pion3_p4.Pt()>0.5 && fabs(gen_pion3_p4.Eta())<=2.5", "selfilters");
  auto nentriesPreF = DFPreF.Count();
  cout << " Total entries in MC with filters by hand " << *nentriesPreF <<endl;

  //auto hNoRes_pre2 = DFPreF.Histo1D( {"massJ", "", Nbin, bindiv},"massJ");
  auto hNoRes_pre2 = DFPreF.Histo1D( {"massJ2", "", Nbin, bindiv},"massJ2");

  //Divide
  auto hNoRes_preD = new TH1F("hNoRes_preD","hNoRes_preD",Nbin, bindiv);
  hNoRes_preD->Sumw2();
  hNoRes_preD->Divide( &(*hNoRes_pre2), &(*hNoRes_pre1), 1, 1);

  TCanvas* canv_histopre = CreateCanvasth1("canv_histopre", hNoRes_preD, "Pre-filter efficiency", "q^{2}"); 
  canv_histopre->SaveAs("plotsEffy/NoResonant_prefilter.png");

   /*
  for(int ii=0; ii<Nbin; ii++){
    cout << hNoRes_pre1->GetBinContent(ii+1) << " , " << hNoRes_pre1->GetBinError(ii+1)  << "," <<  endl;
    cout << hNoRes_pre2->GetBinContent(ii+1) << " , " << hNoRes_pre2->GetBinError(ii+1)  << "," <<  endl;
    cout << hNoRes_preD->GetBinContent(ii+1) << " , " << hNoRes_preD->GetBinError(ii+1)  << ", \n" <<  endl;
  }
  */
  
  return hNoRes_preD;
}

TH1* NoRes_Effy_Reco(UInt_t year=2022, UInt_t sample=5)
 {
  //**************************************************
  // Reco Effy
  //**************************************************
  //MC with filter cuts
  TString MCGenModel = "PHP";
  TChain tree("ntuple");
  if (year==2022 && sample==3){
    tree.Add("../MC_NOres2022/Rootuple_Bdtomumupipi_MCNoRes2022_MCRecoGen_bacht*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==5){
    tree.Add("/home/oscar/Documents/Thesis/Efficencies/MC/MC_RecoGen_NonRes.root");
    MCGenModel = "SVS";
  }
  
  //Datframe with filter cuts
  ROOT::RDataFrame DFRecii(tree);
  auto nentriesRecII = DFRecii.Count();
  cout << " Total entries in MC with filters " << *nentriesRecII <<endl;

  auto DFReci = DFRecii.Filter("gen_muon1_p4.Pt()>3.5 && gen_muon2_p4.Pt()>3.5 && fabs(gen_muon1_p4.Eta())<=2.5 && fabs(gen_muon2_p4.Eta())<=2.5 && gen_pion3_p4.Pt()>0.5 && fabs(gen_pion3_p4.Eta())<=2.5", "selfilters");
  auto nentriesRecI = DFReci.Count();
  cout << " Total entries in MC with filters by hand (3p5 both muons) " << *nentriesRecI <<endl;
  
  auto DFRecI = DFReci.Define("massJ","gen_jpsi_p4.M()").Define("massJ2","gen_jpsi_p4.M()*gen_jpsi_p4.M()");
  auto hNoRes_rec1 = DFRecI.Histo1D( {"massJ2", "", Nbin, bindiv},"massJ2");
  //return;

  //Http://root.cern/doc/v626/namespaceROOT_1_1RDF.html
  /*
  auto BMcut = "massB>(4.9) && massB<(5.7) ";
  auto Low_Cut = "massJ>1.1";
  auto JP_cut = "massJ<2.980 || massJ>3.210 ";
  auto PP_cut = "massJ<3.590 || massJ>3.782";
  auto V0M_cut = "massv0>0.4876 && massv0<0.5076";
  auto antiRad_JP1 = " ( fabs( (massB-5.27934)- (massJ-3.0969) ) > 0.137 ) || (massJ>3.0969) ";
  auto antiRad_JP2 = " ( fabs( (massB-5.27934)- (massJ-3.0969) ) > 0.134 ) || (massJ>3.43) || (massJ<3.0969)";

  auto antiRad_PP1 = " ( fabs( (massB-5.27934)- (massJ-3.6861) ) > 0.097 ) || (massJ>3.6861) ";
  auto antiRad_PP2 = " ( fabs( (massB-5.27934)- (massJ-3.6861) ) > 0.044 ) || (massJ>3.92) || (massJ<3.6861)";
  auto xgbcut = " (XGB_aucpr > 0.9865) ";

  //auto DFRecFi = ROOT::RDF::MakeCsvDataFrame(Form("../MyXGboostTry/CSV_files/MC_sllBall_Jhovanny_sample%1i_V5.csv",sample));
  //auto DFRecFi2 = DFRecFi.Filter(BMcut,"BMcut").Filter(Low_Cut,"Low_Cut").Filter(JP_cut,"JP_cut").Filter(PP_cut,"PP_cut").Filter(V0M_cut,"V0M_cut")
  //  .Filter(antiRad_JP1,"antiRad_JP1").Filter(antiRad_JP2,"antiRad_JP2").Filter(antiRad_PP1,"antiRad_PP1").Filter(antiRad_PP2,"antiRad_PP2").Filter(xgbcut,"xgbcut");
  */

  //auto DFRecFi2 = ROOT::RDF::MakeCsvDataFrame(Form("../MyReweiting/MCntuple_NoRes_FromSplot_Year%1i_Sample%1i.csv",year,sample));
  TChain tree2("treeBu");
  //TChain tree2("ntuple");


  if (year==2022 && sample==3){
    tree2.Add("../MC_NOres2022/Rootuple_Bdtomumupipi_MCNoRes2022_MCRecoGen_bacht*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==5){
    tree2.Add("/home/oscar/Documents/Thesis/Efficencies/MC/ntuple_BuJpsiK_ForXGboost_Year2022_Sample3_v0.root");
    MCGenModel = "SVS";
  }


  ROOT::RDataFrame DFRecFi2(tree2);  

  //auto DFRecFi2 = ROOT::RDF::MakeCsvDataFrame(Form("/cms/home/hcrotte/Reweighting/Results_Dec23_JPsi/MCntuple_NoRes_FromSplot_Year%1i_Sample%1i_BDT.csv",year,sample));
  auto DFRecF = DFRecFi2.Define("massJ2","massJ*massJ");
  auto nentriesRecF = DFRecF.Count();
  cout << " Total entries in MC with filters and Full selection " << *nentriesRecF <<endl;

  //auto hNoRes_rec2 = DFRecF.Histo1D( {"massJ", "", Nbin, bindiv},"massJ");
  //auto hNoRes_rec2 = DFRecF.Histo1D( {"massJ2", "", Nbin, bindiv},"massJ2");
  // ************************************************
  // *** be careful in case weights become needed ***
  // ************************************************
  // *** https://root.cern.ch/doc/master/df005__fillAnyObject_8C.html ***
  // *** https://root-forum.cern.ch/t/rdataframe-with-weight/43992  ***
  // *** https://root-forum.cern.ch/t/rdataframe-histo1d-multiple-weight-branches/31288  ***
  //auto hNoRes_rec2 = DFRecF.Histo1D( {"massJ2", "",  Nbin, bindiv},"massJ2","splotW");//Fill with weights
  auto hNoRes_rec2 = DFRecF.Histo1D( {"massJ2", "",  Nbin, bindiv},"massJ2");//Fill with weights
  //auto hNoRes_rec2 = DFRecF.Define("w","splotW*5.0").Histo1D( {"massJ", "", Nbin, bindiv},"massJ","w");//Fill with weights
  cout << " Total entries in MC with filters and Full selection TH1 " << hNoRes_rec2->Integral() <<endl;

  //Divide
  auto hNoRes_recD = new TH1F("hNoRes_recD","hNoRes_recD", Nbin, bindiv);
  hNoRes_recD->Sumw2();
  hNoRes_recD->Divide( &(*hNoRes_rec2), &(*hNoRes_rec1),1,1);

  TCanvas* canv_historec = CreateCanvasth1("canv_historec", hNoRes_recD, "Reco efficiency", "q^{2}"); 
  canv_historec->SaveAs("plotsEffy/NoResonant_reco.png");

  for(int ii=0; ii<Nbin; ii++){
    cout << hNoRes_rec1->GetBinContent(ii+1) << " , " << hNoRes_rec1->GetBinError(ii+1)  << "," <<  endl;
    cout << hNoRes_rec2->GetBinContent(ii+1) << " , " << hNoRes_rec2->GetBinError(ii+1)  << "," <<  endl;
    cout << hNoRes_recD->GetBinContent(ii+1) << " , " << hNoRes_recD->GetBinError(ii+1)  << ", \n" <<  endl;
  }
  

  return hNoRes_recD;
}

//****************************************************************************************************************************
// Notabene: Very Imporntant
// "year" is 2022 or 20223 
// "sample" is 1 if is Data, 2 if is MCResonante(J/psi) PHSP, 3 if is MC NoResonante(mumu) PHSP,  
//             4 if is MCResonante(J/psi) SVS, 5 if is MC NoResonante(mumu) sllBall.
//****************************************************************************************************************************
void Effy(UInt_t Year=2022, UInt_t Sample=2)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  
  //gErrorIgnoreLevel = 2001;
  //gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");// Ignore Warnings
  ROOT::EnableImplicitMT(8);// Tell ROOT you want to go parallel
  
  // **************************************************
  // Jpsi prefilter Effy
  // **************************************************
  
  std::vector<double> Jpsi_PreEffyVal = Jpsi_Effy_Pre(Year, 4);
  cout<< "J/psi prefilter Effy = " << Jpsi_PreEffyVal[0] << " +/- " << Jpsi_PreEffyVal[1] << " \n"  << endl;

  // **************************************************
  // Jpsi Reco Effy
  // **************************************************
  
  std::vector<double> Jpsi_RecoEffyVal = Jpsi_Effy_Reco(Year, 4);
  cout<< "J/psi Reco Effy = " << Jpsi_RecoEffyVal[0] << " +/- " << Jpsi_RecoEffyVal[1] << " \n" << endl;
  
  // **************************************************
  // Jpsi Total Effy
  // **************************************************
  
  double Jpsi_TotalEffyVal = Jpsi_PreEffyVal[0]*Jpsi_RecoEffyVal[0];
  double Jpsi_TotalEffyVale = (Jpsi_TotalEffyVal)*sqrt( (Jpsi_RecoEffyVal[1]/Jpsi_RecoEffyVal[0])*(Jpsi_RecoEffyVal[1]/Jpsi_RecoEffyVal[0])
							+ (Jpsi_PreEffyVal[1]/Jpsi_PreEffyVal[0])*(Jpsi_PreEffyVal[1]/Jpsi_PreEffyVal[0]) );

  cout<< "J/psi Total Effy = " << Jpsi_TotalEffyVal << " +/- " << Jpsi_TotalEffyVale << " \n"  << endl;
  
  // **************************************************
  // Psi2S prefilter Effy
  // **************************************************
  
  std::vector<double> Psi2S_PreEffyVal = Psi2S_Effy_Pre(Year, 7);
  cout<< "Psi2S prefilter Effy = " << Psi2S_PreEffyVal[0] << " +/- " << Psi2S_PreEffyVal[1] << " \n"  << endl;
  
  // **************************************************
  // Psi2S Reco Effy
  // **************************************************
  
  std::vector<double> Psi2S_RecoEffyVal = Psi2S_Effy_Reco(Year, 7);
  cout<< "Psi2S Reco Effy = " << Psi2S_RecoEffyVal[0] << " +/- " << Psi2S_RecoEffyVal[1] << " \n" << endl;
  
  // **************************************************
  // Psi2S Total Effy
  // **************************************************
  
  double Psi2S_TotalEffyVal = Psi2S_PreEffyVal[0]*Psi2S_RecoEffyVal[0];
  double Psi2S_TotalEffyVale = (Psi2S_TotalEffyVal)*sqrt( (Psi2S_RecoEffyVal[1]/Psi2S_RecoEffyVal[0])*(Psi2S_RecoEffyVal[1]/Psi2S_RecoEffyVal[0])
							+ (Psi2S_PreEffyVal[1]/Psi2S_PreEffyVal[0])*(Psi2S_PreEffyVal[1]/Psi2S_PreEffyVal[0]) );

  cout<< "Psi2S Total Effy = " << Psi2S_TotalEffyVal << " +/- " << Psi2S_TotalEffyVale << " \n"  << endl;
  
  // **************************************************
  // NoRes prefilter Effy
  // **************************************************
  
  //TH1* NoRes_Preh = NoRes_Effy_Pre(Year, 5);
  
  // **************************************************
  // NoRes Reco Effy
  // **************************************************
  
  //TH1* NoRes_Recoh = NoRes_Effy_Reco(Year, 5);
  
  // **************************************************
  // Nores Total Effy
  // **************************************************
  /*
  TH1* NoRes_Totalh = new TH1F("NoRes_Totalh","NoRes_Totalh", Nbin, bindiv);
  NoRes_Totalh->Sumw2();    
  NoRes_Totalh->Multiply(NoRes_Preh,NoRes_Recoh,1,1);

  TCanvas* canv_NoRes_Totalh = CreateCanvasth1("canv_NoRes_Totalh", NoRes_Totalh, "Total efficiency", "q^{2}"); 
  canv_NoRes_Totalh->SaveAs("plotsEffy/NoResonant_Total.png");
  */

  ofstream salida ("plotsEffy/Totaleffy_q2.csv");

  //cout << "NoRes_Totalh->GetXaxis()->GetNbins() = " << NoRes_Totalh->GetXaxis()->GetNbins() << endl; //7
  salida.is_open();
  for(int ii=0; ii<Nbin; ii++){
    /*NoRes_Preh->GetBinContent(ii+1) << "," << NoRes_Preh->GetBinError(ii+1)  << "," << 
      NoRes_Recoh->GetBinContent(ii+1) << "," << NoRes_Recoh->GetBinError(ii+1)  << "," <<
      NoRes_Totalh->GetBinContent(ii+1) << "," << NoRes_Totalh->GetBinError(ii+1) << "," << */
    
    salida <<  Jpsi_PreEffyVal[0] << "," << Jpsi_PreEffyVal[1] << "," << Jpsi_RecoEffyVal[0] << "," << Jpsi_RecoEffyVal[1] << "," << Jpsi_TotalEffyVal << "," << Jpsi_TotalEffyVale << "," << endl ;
      //Psi2S_PreEffyVal[0] << "," << Psi2S_PreEffyVal[1] << "," << Psi2S_RecoEffyVal[0] << "," << Psi2S_RecoEffyVal[1] << "," << Psi2S_TotalEffyVal << "," << Psi2S_TotalEffyVale << endl;

  }
  salida.close(); 
  
  
}// End of main funcion

