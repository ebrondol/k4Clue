#include "tdrstyle.C"
#include "library.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TProfile.h>

float enMin = 0.03;
float enMax = 5;

std::vector< TProfile*> returnEnergyProfiles(TString fileName, float ymax ){

  std::cout << "Read file: " << fileName << std::endl;
  std::vector<TProfile*> product{};
  auto f = TFile::Open(fileName);
  if (!f || f->IsZombie()) {
     return product;
  }

  TTree* treeHits = (TTree*)f->Get("CLUEHits");
  std::vector<int> *status = 0;
  std::vector<float> *mcenergy = 0;
  std::vector<int> *region = 0;
  treeHits->SetBranchAddress("status",&status);
  treeHits->SetBranchAddress("MCEnergy",&mcenergy);
  treeHits->SetBranchAddress("region",&region);
  Int_t nentries = (Int_t)treeHits->GetEntries();
  std::cout << nentries << " events to read " << std::endl;

  float xmin = float(enMin);
  float xmax = float(enMax);
  int bins   = 50;

  TString profileNameSeeds = "p_layer_seeds";
  TProfile *p_var_Seeds = new TProfile("p_layer_seeds", "p_layer_seeds", bins, xmin, xmax, 0, ymax);
  TProfile *p_var_Followers = new TProfile("p_layer_followers", "p_layer_followers", bins, xmin, xmax, 0, ymax);
  TProfile *p_var_Outliers = new TProfile("p_layer_outliers", "p_layer_outliers", bins, xmin, xmax, 0, ymax);

  for (Int_t i = 0; i < nentries; i++){
    //if(i > 10) break;
    std::vector<int> seedsPerLayer(bins, 0);
    std::vector<int> outliersPerLayer(bins, 0);
    std::vector<int> followersPerLayer(bins, 0);
    std::cout << "event #" << i << std::endl; 
    treeHits->GetEntry(i);

    int countS = std::count(status->begin(), status->end(), 2);
    int countO = std::count(status->begin(), status->end(), 0);
    int countF = std::count(status->begin(), status->end(), 1);
    p_var_Seeds->Fill(mcenergy->at(0), countS, 1);
    p_var_Outliers->Fill(mcenergy->at(0), countO, 1);
    p_var_Followers->Fill(mcenergy->at(0), countF, 1);

/*
    for (auto ihit = 0; ihit < int(layer->size()); ihit++){
        int hitLayer = layer->at(ihit);
        if( hitLayer > 40 ) {
          hitLayer = hitLayer%40;
        }
        if(int(status->at(ihit))==2){
          seedsPerLayer[hitLayer] += 1;
        } else if(int(status->at(ihit))==1) {
          followersPerLayer[hitLayer] += 1;
        } else if(int(status->at(ihit))==0){
          outliersPerLayer[hitLayer] += 1;
        }
    }
    //std::for_each(seedPerLayer.begin(), seedPerLayer.end(),
    //             [&] (int elem) {
    //               std::cout << elem << std::endl;
    //             });
    for (auto ilay = 0; ilay < int(seedsPerLayer.size()); ilay++){
      p_var_Seeds->Fill(ilay, seedsPerLayer[ilay], 1);
      p_var_Followers->Fill(ilay, followersPerLayer[ilay], 1);
      p_var_Outliers->Fill(ilay, outliersPerLayer[ilay], 1);
    }
*/
  }

  product.push_back(p_var_Outliers);
  product.push_back(p_var_Followers);
  product.push_back(p_var_Seeds);

  return product;

}

void profile(TString detector, TString variable, TString energy, 
             TString dc, TString rhoc, TString of,
             TString var_short, float ymax){

  TString settings = energy+"GeV_"+dc+"dc_"+rhoc+"rhoc_"+of+"of";

  //gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  int totLayers = nLayersFromDet(detector);
  TString detectorLabel = detlabelFromDet(detector);
  TString folder = folderFromDet(detector);
  TString fileName = folder+"/"+filenameFromDet(detector, energy, dc, rhoc, of);

  float xmin = 0;
  float xmax = 10;
  if(variable=="layer"){
    xmax = float(totLayers);
  } else if(variable=="energy"){
    xmin = float(enMin);
    xmax = float(enMax);
  }

  std::vector<TProfile*> profiles = {};
  std::cout << variable << std::endl;
  if(variable=="layer"){
    profiles = returnLayerProfiles(fileName, ymax, totLayers);
  } else if(variable=="energy"){
    profiles = returnEnergyProfiles(fileName, ymax);
  }
  if(profiles.empty())
    return;

  profiles[0]->SetLineColor(kRed+1);
  profiles[0]->SetLineWidth(3);
  profiles[0]->SetLineStyle(3);
  profiles[0]->SetMarkerColor(kRed+1);
  profiles[0]->SetMarkerSize(0.9);
  profiles[0]->SetMarkerStyle(25);
  profiles[0]->GetYaxis()->SetTitle("<Hits>");
  profiles[0]->GetYaxis()->SetTitleOffset(1.4);
  profiles[0]->SetMaximum(ymax);
  if(variable=="layer"){
    profiles[0]->GetXaxis()->SetTitle("Layer Number");
  }

  profiles[1]->SetLineColor(kBlue+1);
  profiles[1]->SetLineWidth(2);
  profiles[1]->SetLineStyle(2);
  profiles[1]->SetMarkerColor(kBlue+1);
  profiles[1]->SetMarkerSize(1.2);
  profiles[1]->SetMarkerStyle(33);

  profiles[2]->SetLineColor(kGreen+3);
  profiles[2]->SetLineWidth(1);
  profiles[2]->SetLineStyle(1);
  profiles[2]->SetMarkerColor(kGreen+3);

  TLegend* t1 = new TLegend(0.506689,0.68,0.8996656,0.906087);
  TString header = "#splitline{"+detectorLabel+", "+energy+" GeV gamma}";
  header += "{o_{f} = "+of+", #rho_{c} = "+rhoc+"}";
  t1->SetHeader(header);
  t1->AddEntry(profiles[2], "Seeds", "pl");
  t1->AddEntry(profiles[1], "Followers", "pl");
  t1->AddEntry(profiles[0], "Outliers", "pl");

  TLine *l1 = new TLine(xmin, 1, xmax, 1);
  l1->SetLineColor(kOrange+2);

  TString canvasName = "cProf_"+var_short+"_"+settings;
  TCanvas *c_var = new TCanvas(canvasName,canvasName, 50, 50, 600, 600);//,200,10,700,780);
  c_var->SetRightMargin(0.06);
  profiles[0]->Draw();
  profiles[1]->Draw("same");
  profiles[2]->Draw("same");
  t1->Draw("same");
  l1->Draw("same");
  c_var->Draw(); 

  TString plotName = folder + "/" + variable+"Profile_" + settings;
  c_var->SaveAs(plotName+".png", "png");
  c_var->SaveAs(plotName+".eps", "eps");
  c_var->SaveAs(plotName+".pdf", "pdf");

  return;
}

void profileComparison(TString detector, TString hitStatus, TString variable, TString var_short, TString energy,
                       TString dc, TString rhoc1, TString of1, float ymax,
                       TString rhoc2 = "", TString of2 = "", TString rhoc3 = "", TString of3 = ""){

  //TString settings = energy+"GeV_"+hitStatus+"Comparison_"+rhoc1+"rhoc_"+of1+"of";
  TString settings = energy+"GeV_"+rhoc1+"rhoc_"+of1+"of";
  if(rhoc2 != "" && of2 != "") {
    settings += "_"+rhoc2+"rhoc_"+of2+"of";
  }
  if(rhoc3 != "" && of3 != "") {
    settings += "_"+rhoc3+"rhoc_"+of3+"of";
  }

  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  int totLayers = nLayersFromDet(detector);
  TString detectorLabel = detlabelFromDet(detector);
  TString folder = folderFromDet(detector);
  TString fileName1 = folder+"/"+filenameFromDet(detector, energy, dc, rhoc1, of1);
  TString fileName2 = folder+"/"+filenameFromDet(detector, energy, dc, rhoc2, of2);
  TString fileName3 = folder+"/"+filenameFromDet(detector, energy, dc, rhoc3, of3);

  float xmin = 0;
  float xmax = float(totLayers);

  std::vector<TProfile*> profiles1 = returnLayerProfiles(fileName1, ymax, totLayers);
  if(profiles1.empty())
    return;

  std::vector<TProfile*> profiles2;
  if(rhoc2 != "" && of2 != "") {
    profiles2 = returnLayerProfiles(fileName2, ymax, totLayers);
    if(profiles2.empty())
      return;
  }
  std::vector<TProfile*> profiles3;
  if(rhoc3 != "" && of3 != "") {
    profiles3 = returnLayerProfiles(fileName3, ymax, totLayers);
    if(profiles3.empty())
      return;
  }

  int status = 0;
  if(hitStatus == "outliers"){
    status = 0;
    profiles1[status]->GetYaxis()->SetTitle("<Outliers>");
    profiles1[status]->SetLineColor(kRed-4);
    profiles1[status]->SetMarkerColor(kRed-4);
    if(rhoc2 != "" && of2 != "") {
      profiles2[status]->SetLineColor(kRed+1);
      profiles2[status]->SetMarkerColor(kRed+1);
      profiles2[status]->SetMarkerSize(0.9);
    }
    if(rhoc3 != "" && of3 != "") {
      profiles3[status]->SetLineColor(kRed+3);
      profiles3[status]->SetMarkerColor(kRed+3);
      profiles3[status]->SetMarkerSize(1.2);
      profiles3[status]->SetMarkerStyle(33);
    }
  } else if(hitStatus == "followers"){
    status = 1;
    profiles1[status]->GetYaxis()->SetTitle("<Followers>");
    profiles1[status]->SetLineColor(kBlue-4);
    profiles1[status]->SetMarkerColor(kBlue-4);
    if(rhoc2 != "" && of2 != "") {
      profiles2[status]->SetLineColor(kBlue+1);
      profiles2[status]->SetMarkerColor(kBlue+1);
      profiles2[status]->SetMarkerSize(0.9);
    }
    if(rhoc3 != "" && of3 != "") {
      profiles3[status]->SetLineColor(kBlue+3);
      profiles3[status]->SetMarkerColor(kBlue+3);
      profiles3[status]->SetMarkerSize(1.2);
      profiles3[status]->SetMarkerStyle(33);
    }
  } else if(hitStatus == "seeds"){
    status = 2;
    profiles1[status]->GetYaxis()->SetTitle("<Seeds>");
    profiles1[status]->SetLineColor(kGreen-2);
    profiles1[status]->SetMarkerColor(kGreen-2);
    if(rhoc2 != "" && of2 != "") {
      profiles2[status]->SetLineColor(kGreen+1);
      profiles2[status]->SetMarkerColor(kGreen+1);
      profiles2[status]->SetMarkerSize(0.9);
    }
    if(rhoc3 != "" && of3 != "") {
      profiles3[status]->SetLineColor(kGreen+3);
      profiles3[status]->SetMarkerColor(kGreen+3);
      profiles3[status]->SetMarkerSize(1.2);
      profiles3[status]->SetMarkerStyle(33);
    }
  } else {
    std::cout << "Status does not exist. Check spelling" << std::endl; 
    return;
  }

  profiles1[status]->SetMarkerSize(0.9);
  profiles1[status]->SetMarkerStyle(25);
  profiles1[status]->SetMaximum(ymax);
  profiles1[status]->GetYaxis()->SetTitleOffset(1.4);
  profiles1[status]->SetMaximum(ymax);
  if(variable=="layer"){
    profiles1[status]->GetXaxis()->SetTitle("Layer Number");
  }

  TLegend* t1 = new TLegend(0.506689,0.68,0.8996656,0.906087);
  TString header = detectorLabel+", "+energy+" GeV gamma";
  t1->SetHeader(header);
  t1->AddEntry(profiles1[status], "o_{f} = "+of1+", #rho_{c} = "+rhoc1[0]+"."+rhoc1[2]+rhoc1[3], "pl");

  TString canvasName = "cProf_"+var_short+"_"+settings;
  TCanvas *c_var = new TCanvas(canvasName,canvasName, 50, 50, 600, 600);//,200,10,700,780);
  c_var->SetRightMargin(0.06);
  profiles1[status]->Draw();
  if(rhoc2 != "" && of2 != "") {
    profiles2[status]->Draw("same");
    t1->AddEntry(profiles2[status], "o_{f} = "+of2+", #rho_{c} = "+rhoc2[0]+"."+rhoc2[2]+rhoc2[3], "pl");
  }
  if(rhoc3 != "" && of3 != "") {
    profiles3[status]->Draw("same");
    t1->AddEntry(profiles3[status], "o_{f} = "+of3+", #rho_{c} = "+rhoc3[0]+"."+rhoc3[2]+rhoc3[3], "pl");
  }

  t1->Draw("same");

  if(hitStatus == "seeds"){
    TLine *l1 = new TLine(xmin, 1, xmax, 1);
    l1->SetLineColor(kOrange+2);
    l1->Draw("same");
  }
  c_var->Draw();

  TString plotName = folder + "/" + hitStatus + "_" + settings;
  c_var->SaveAs(plotName+".png", "png");
  c_var->SaveAs(plotName+".eps", "eps");
  c_var->SaveAs(plotName+".pdf", "pdf");

  return;

}

void distribution(TString detector, TString variable, TString energy, TString dc, TString rhoc, TString of,
                  TString var_short, bool endcap, int bins   = 40, float xmin = 0, float xmax = 40, float ymax = 1200){

  TString settings = energy+"GeV_"+rhoc+"rhoc_"+of+"of";

  TString regionCut = "region == 0";
  if(endcap){
    regionCut = "region == 1";
    if(variable=="layer"){
      xmin = 0;
      xmax = 80;
      bins = 80;
    }
  }

  //gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString detectorLabel = detlabelFromDet(detector);
  TString folder = folderFromDet(detector);
  TString fileName = folder+"/"+filenameFromDet(detector, energy, dc, rhoc, of);

  // Open input file
  auto f = TFile::Open(fileName);
  if (!f || f->IsZombie()) {
     return;
  }

  TNtuple *ntupleHits = (TNtuple*)f->Get("CLUEHits");
  TString histoNameFollowers = "h_"+var_short+"_Followers_"+settings;
  TH1F* h_var_Followers = new TH1F(histoNameFollowers, histoNameFollowers, bins, xmin, xmax);
  float ratio = (xmax-xmin)/bins;
  h_var_Followers->GetXaxis()->SetTitle(variable);
  if(variable=="layer"){
    h_var_Followers->GetXaxis()->SetTitle("Layer Number");
  } else if(variable=="energy"){
    h_var_Followers->GetXaxis()->SetTitle("Hit Energy [GeV]");
  }
  h_var_Followers->GetYaxis()->SetTitle("<Hits>");
  h_var_Followers->GetYaxis()->SetTitleOffset(1.4);
  h_var_Followers->SetMaximum(ymax);
  h_var_Followers->SetLineColor(kBlue+1);
  h_var_Followers->SetLineWidth(2);
  h_var_Followers->SetLineStyle(3);

  TString histoNameOutliers = "h_"+var_short+"_Outliers_"+settings;
  TH1F* h_var_Outliers = new TH1F(histoNameOutliers, histoNameOutliers, bins, xmin, xmax);
  h_var_Outliers->SetLineColor(kRed+1);
  h_var_Outliers->SetLineWidth(2);
  h_var_Outliers->SetLineStyle(2);

  TString histoNameSeeds = "h_"+var_short+"_Seeds_"+settings;
  TH1F* h_var_Seeds = new TH1F(histoNameSeeds, histoNameSeeds, bins, xmin, xmax);
  h_var_Seeds->SetLineColor(kGreen+2);
  h_var_Seeds->SetLineWidth(2);
  h_var_Seeds->SetLineStyle(1);

  TCanvas *cTrash = new TCanvas();//,200,10,700,780);
  ntupleHits->Draw(variable + " >> "+histoNameFollowers, "status==1 && "+regionCut);
  ntupleHits->Draw(variable + " >> "+histoNameOutliers, "status==0  && "+regionCut);
  ntupleHits->Draw(variable + " >> "+histoNameSeeds, "status==2     && "+regionCut);

  //TLegend* t1 = new TLegend(0.60,0.23,0.90,0.38);
  TLegend* t1 = new TLegend(0.506689,0.68,0.8996656,0.906087);
  TString header = "#splitline{"+detectorLabel+", "+energy+" GeV gamma";
  if(endcap){
    header += ", endcap}";
  } else {
    header += ", barrel}";
  }
  header += "{o_{f} = "+of+", #rho_{c} = "+rhoc[0]+"."+rhoc[2]+rhoc[3]+"}";
  t1->SetHeader(header);
  t1->AddEntry(h_var_Seeds, "Seeds", "pl");
  t1->AddEntry(h_var_Followers, "Followers", "pl");
  t1->AddEntry(h_var_Outliers, "Outliers", "pl");

  TLine *t500 = new TLine(xmin, 500, xmax, 500);
  t500->SetLineColor(kOrange+2);

  TString canvasName = "cDistr_"+var_short+"_"+settings;
  TCanvas *c_var = new TCanvas(canvasName,canvasName, 50, 50, 600, 600);//,200,10,700,780);
  gStyle->SetOptStat(0);
  if(variable=="energy"){
    c_var->SetLogy();
  }
  c_var->SetRightMargin(0.06);
  c_var->cd();
  h_var_Followers->Draw();
  h_var_Outliers->Draw("same");
  h_var_Seeds->Draw("same");
  t1->Draw("same");
  if(variable=="layer")
    t500->Draw("same");
  c_var->Draw();
  //TString plotName = "../plots/" + variable+"Distr_" + settings;
  TString plotName = folder + "/" + variable+"Distr_" + settings;
  if(endcap){
    plotName += "_endcap";
  } else {
    plotName += "_barrel";
  }
  plotName += ".png";
// +"_log.png"
  c_var->SaveAs(plotName, "png");

  return;

}

void hitsPlots(){

//  profile("energy", "5", "40.0", "0.03", "3.0", "en", 100);

//  distribution("cld", "energy", "10", "15.0", "0.02", "2.0", "en", false, 100, 0.0, 1.0, 100000);
//  distribution("energy", "10", "002", "2", "en", true, 100, 0.0, 1.0, 5500);
//  distribution("energy", "10", "002", "3", "en", false, 100, 0.0, 1.0, 100000);
//  distribution("energy", "10", "003", "2", "en", false, 100, 0.0, 1.0, 100000);




//  profileComparison("cld", "outliers", "layer", "layOut",  "10", "15.0", "0.02", "1.0",  4, "0.02", "2.0", "0.02", "3.0");
//  profileComparison("cld", "seeds", "layer", "layOut",     "10", "15.0", "0.02", "1.0",  4, "0.02", "2.0", "0.02", "3.0");
//  profileComparison("cld", "followers", "layer", "layOut", "10", "15.0", "0.02", "1.0", 20, "0.02", "2.0", "0.02", "3.0");

//  profileComparison("cld", "outliers", "layer", "layOut",  "10", "15.0", "0.01", "3.0",  4, "0.02", "3.0", "0.03", "3.0");
//  profileComparison("cld", "seeds", "layer", "layOut",     "10", "15.0", "0.01", "3.0",  4, "0.02", "3.0", "0.03", "3.0");
//  profileComparison("cld", "followers", "layer", "layOut", "10", "15.0", "0.01", "3.0", 20, "0.02", "3.0", "0.03", "3.0");

  profile("cld", "layer",     "10",  "15.0", "0.02", "3.0", "lay", 20);
  profile("cld", "layer",     "100",  "15.0", "0.02", "3.0", "lay", 100);
//  profile("clicdet", "layer", "10",  "15.0", "0.02", "3.0", "lay", 20);
//  profile("clicdet", "layer", "100", "15.0", "0.02", "3.0", "lay", 100);





//  profileComparison("followers:outliers", "layer", "layOut", "10", "40.0", "0.03", "3.0", 70, "0.03", "2.0", "0.02", "3.0");
 // profileComparison("seeds", "layer", "layOut", "10", "40.0", "0.03", "3.0", 4, "0.03", "2.0", "0.02", "3.0");
  //profileComparison("followers", "layer", "layOut", "10", "40.0", "0.03", "3.0", 70, "0.03", "2.0", "0.02", "3.0");

  //profileComparison("seeds", "layer", "layOut", "10", "40.0", "0.03", "3.0", 10, "0.03", "3.0");

  //profile("layer", "10", "40.0", "0.03", "3.0", "lay", 70);

  return;
}
