#ifndef MYUTILITIESCPP_H
#define MYUTILITIESCPP_H

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

TCanvas* CreateCanvasth1(TString cname, TH1 *h1, TString yhname="yname", TString xhname="xname") 
{
 int H = 600;
 int W = 800;
 TCanvas* canv = new TCanvas(cname,cname,50,50,W,H);
 canv->cd();
 canv->SetLeftMargin(0.15);
 canv->SetRightMargin(0.06);
 canv->SetTopMargin(0.09);
 canv->SetBottomMargin(0.14);
 gPad->SetLogy();
 //h1->Scale(1.0/h1->Integral());

 h1->Draw("");
 h1->SetMarkerStyle(20);
 //h1->SetMarkerSize(1.5);
 h1->SetMarkerColor(1);
 h1->SetLineColor(1);
 h1->SetLineWidth(2);
 h1->GetYaxis()->CenterTitle(true);
 //h1->SetYTitle("Pt(J/#psi K^{+}) [GeV]");
 h1->SetYTitle(yhname);
 //h1->SetYTitle("Distribution normalized."+yhname);
 h1->GetXaxis()->CenterTitle(true); 
 //h1->SetXTitle("|#eta(J/#psi )|"); 
 h1->SetXTitle(xhname); 
 h1->SetTitleSize(35,"XY"); 
 h1->SetLabelSize(30,"XY");
 h1->SetTitleOffset(1.1,"Y");
 h1->SetTitleOffset(1.0,"X");
 h1->SetLabelFont(43,"XY");  
 h1->SetTitleFont(43,"XY");
 h1->SetMinimum(0.0001);
 //h1->SetMinimum(0.0);
 h1->SetMaximum(1.0);
 //h1->SetMaximum(1.0);
 //h1->SetMaximum(0.7);

 TLatex *   tex1 = new TLatex(0.94,0.926,"MC simulation");
 tex1->SetNDC();
 tex1->SetTextAlign(31);
 tex1->SetTextFont(42);
 tex1->SetTextSize(0.04); 
 tex1->SetLineWidth(2);
 
 TLatex *tex2 = new TLatex(0.15,0.926,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.04); 
 tex2->SetLineWidth(2);

 TLatex *tex3 = new TLatex(0.23,0.926,"Preliminary");
 tex3->SetNDC();
 tex3->SetTextFont(52);
 tex3->SetTextSize(0.04); 
 tex3->SetLineWidth(2); 

 tex1->Draw();
 tex2->Draw();
 tex3->Draw();

 return canv;
}

std::vector<double> Jpsi_Effy_Pre(UInt_t year=2022, UInt_t sample=4)
 {
  //**************************************************
  // prefilter Effy
  //**************************************************
  //MC with NO filter cuts
  TString MCGenModel = "PHP";
  TChain tree("ntuple");
  if (year==2022 && sample==2){
    tree.Add("../MC_Jpsi2022/YESpions/Rootuple_Bd_toJpipi_2022_MiniAOD_MConlyGen_bacht*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==4){
    tree.Add("/home/oscar/Documents/Thesis/Efficencies/MC/JPsi/MC_JPsiK_OnlyGen.root");
    MCGenModel = "SVS";
  }
  
  //Datframe with NO filter
  ROOT::RDataFrame DFPrei(tree);
  //DFPreI.Describe().Print();
  auto nentriesPreI = DFPrei.Count();
  cout << " Total entries in MC with NO filter " << *nentriesPreI <<endl;

  //40, 2.9, 3.3
  //23, 2.980, 3.210
  auto DFPreI = DFPrei.Define("massJ","gen_jpsi_p4.M()");
  auto hRes_pre1 = DFPreI.Histo1D( {"massJ", "", 23, 2.980, 3.210},"massJ");
  //return;
  
  //Datframe with filter cuts by hand
  auto DFPreF = DFPreI.Filter("gen_muon1_p4.Pt()>3.5 && gen_muon2_p4.Pt()>3.5 && fabs(gen_muon1_p4.Eta())<=2.5 && fabs(gen_muon2_p4.Eta())<=2.5 && gen_pion3_p4.Pt()>0.5 && fabs(gen_pion3_p4.Eta())<=2.5", "selfilters");
  auto nentriesPreF = DFPreF.Count();
  cout << " Total entries in MC with filters by hand " << *nentriesPreF <<endl;

  //auto hRes_pre2 = DFPreF.Histo1D( {"massJ", "", 76, 1.1, 4.8},"massJ");
  auto hRes_pre2 = DFPreF.Histo1D( {"massJ", "", 23, 2.980, 3.210},"massJ");

  //Divide
  auto hRes_preD = new TH1F("hRes_preD","hRes_preD",23, 2.980, 3.210);
  hRes_preD->Divide( &(*hRes_pre2), &(*hRes_pre1),1,1);

  TCanvas* canv_histopre = CreateCanvasth1("canv_histopre", hRes_preD, "Pre-filter efficiency", "m_{#mu#mu}"); 
  canv_histopre->SaveAs("plotsEffy/Resonant_prefilter.png");

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

std::vector<double> Jpsi_Effy_Reco(UInt_t year=2022, UInt_t sample=4)
 {
  //**************************************************
  // Reco Effy
  //**************************************************
  //MC with filter cuts
  TString MCGenModel = "PHP";
  TChain tree("ntuple");
  if (year==2022 && sample==2){
    tree.Add("../MC_Jpsi2022/YESpions/Rootuple_Bd_toJpipi_2022_MiniAOD_MCRecoGen_bacht*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==4){
    tree.Add("/home/oscar/Documents/Thesis/Efficencies/MC/JPsi/MC_JPsiK_RecoGen.root");
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
  auto hRes_rec1 = DFRecI.Histo1D( {"massJ", "", 23, 2.980, 3.210},"massJ");
  //return;

  //https://root.cern/doc/v626/namespaceROOT_1_1RDF.html
  //auto DFRecFi = ROOT::RDF::MakeCsvDataFrame(Form("../MyReweiting/results_v5/weighted_MCntuple_Jpsiv0_FromSplot_Year%1i_Sample%1i.csv",year,sample));
  //auto DFRecF = DFRecFi.Filter("massB > 5.0 && massB < 5.700 && massJ>2.980 && massJ<3.210 && massv0>0.4876 && massv0<0.5076", "selfilters");
  //auto DFRecF = ROOT::RDF::MakeCsvDataFrame(Form("../MyReweiting/MCntuple_Jpsiv0_FromSplot_Year%1i_Sample%1i.csv",year,sample));
  
  TChain tree2("treeBu");
  //TChain tree2("ntuple");


  if (year==2022 && sample==3){
    tree2.Add("../MC_NOres2022/Rootuple_Bdtomumupipi_MCNoRes2022_MCRecoGen_bacht*.root/rootuple/ntuple");
  }
  else if (year==2022 && sample==4){
    tree2.Add("/home/oscar/Documents/Thesis/Efficencies/MC/JPsi/MC_JPsiK_hyperOpt_v3_Resonant.root");
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
  //auto hRes_rec2 = DFRecF.Histo1D( {"massJ", "", 23, 2.980, 3.210},"massJ");
  // ************************************************
  // *** be careful in case weights become needed ***
  // ************************************************
  // *** https://root.cern.ch/doc/master/df005__fillAnyObject_8C.html ***
  // *** https://root-forum.cern.ch/t/rdataframe-with-weight/43992  ***
  // *** https://root-forum.cern.ch/t/rdataframe-histo1d-multiple-weight-branches/31288  ***
  auto hRes_rec2 = DFRecF2.Histo1D( {"massJ", "", 23, 2.980, 3.210},"massJ");//Fill with weights
  //auto hRes_rec2 = DFRecF.Define("w","splotW*5.0").Histo1D( {"massJ", "", 23, 2.980, 3.210},"massJ","w");//Fill with weights

  //Divide
  auto hRes_recD = new TH1F("hRes_recD","hRes_recD",23, 2.980, 3.210);
  hRes_recD->Divide( &(*hRes_rec2), &(*hRes_rec1),1,1);

  TCanvas* canv_historec = CreateCanvasth1("canv_historec", hRes_recD, "Reco efficiency", "m_{#mu#mu}"); 
  canv_historec->SaveAs("plotsEffy/Resonant_reco.png");
  //double recEffyVal = static_cast< double >(*nentriesRecF) / (*nentriesRecI);
  //return recEffyVal;

  double EFi = (double)(*nentriesRecI);
  double EFie = sqrt(EFi);
  //double EFf = (double)(*nentriesRecF);
  cout << " Total entries in MC with filters and Full selection " << *nentriesRecF2 <<endl;
  cout << " Total entries in MC with filters and Full selection TH1 " << hRes_rec2->Integral() <<endl;
  //double EFf = (double)(hRes_rec2->Integral());
  double EFf = (double)(*nentriesRecF2);
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
