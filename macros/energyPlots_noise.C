#include "tdrstyle.C"
#include "library.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

void distributionComparison(TString detector, TString variable, TString var_short, TString energy,
                            TString dc, TString rhoc, TString of,
                            float xmin, float xmax, float ymax, int bins,
                            bool filterAndNoiseTopo = true,
                            TString dc2 = "", TString rhoc2 = "", TString of2 = "", 
                            TString dc3 = "", TString rhoc3 = "", TString of3 = ""){

  TString settings = energy+"GeV_"+rhoc+"rhoc_"+of+"of";
  if(dc2 != "" && rhoc2 != "" && of2 != "" && dc2 != dc && rhoc2 != rhoc && of2 != of) {
    settings += "_"+rhoc2+"rhoc_"+of2+"of";
  }
  if(dc3 != "" && rhoc3 != "" && of3 != "" && dc3 != dc && rhoc3 != rhoc && of3 != of) {
    settings += "_"+rhoc3+"rhoc_"+of3+"of";
  }

  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString detectorLabel = detlabelFromDet(detector);
  TString folder = folderFromDet(detector);
  //TString fileName3 = folder+"/"+filenameFromDet(detector, energy, dc3, rhoc3, of3);
  TString fileName1= folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc+"dc_"+rhoc+"rhoc_"+of+"of.root";
  TString fileName2= folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc2+"dc_"+rhoc2+"rhoc_"+of2+"of_clue_filter2sigma.root";
  TString fileName3= folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc3+"dc_"+rhoc3+"rhoc_"+of3+"of_noise.root";
  settings += "_noise";

  if(filterAndNoiseTopo){
    fileName1= folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc+"dc_"+rhoc+"rhoc_"+of+"of_clue_filter2sigma.root";
    fileName2= folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc2+"dc_"+rhoc2+"rhoc_"+of2+"of_noise.root";
    fileName3= folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc3+"dc_"+rhoc3+"rhoc_"+of3+"of_noise.root";
  } else {
    settings += "_only_clue";
  }

  TString ntupleName = "CLUEClusters";
  if( detector == "lar"){
    ntupleName = "CaloClueClusters";
  }

  TH1F* histo1 = histosFromNtuples(fileName1, variable, ntupleName, xmin, xmax, bins);
  float ratio = (xmax-xmin)/bins;
  histo1->GetXaxis()->SetTitle(variable);
  if(variable=="totEnergyHits"){
    //histo1->GetXaxis()->SetTitle(TString::Format("Total energy [GeV/%.2f]", ratio));
    histo1->GetXaxis()->SetTitle("Total energy [GeV]");
  } else if(variable=="totEnergyHits/MCEnergy"){
    histo1->GetXaxis()->SetTitle("E_{reco}/E_{MC}");
  }
  histo1->GetYaxis()->SetTitle("a.u.");
  histo1->GetYaxis()->SetTitleOffset(1.4);
  if(filterAndNoiseTopo){
    histo1->SetLineColor(kRed+2);
    histo1->SetMarkerColor(kRed+2);
    histo1->SetMarkerStyle(25);
  } else {
    histo1->SetLineColor(kRed);
    histo1->SetMarkerColor(kRed);
    histo1->SetMarkerStyle(21);
  }
  histo1->SetLineWidth(3);
  histo1->SetLineStyle(1);
  histo1->SetStats(0);
  histo1->Scale(1./histo1->Integral());

  TH1F* histo2;
  if(dc2 != "" && rhoc2 != "" && of2 != "") {
    histo2 = histosFromNtuples(fileName2, variable, ntupleName, xmin, xmax, bins);
    if(filterAndNoiseTopo){
      histo2->SetLineColor(kBlue);
      histo2->SetMarkerColor(kBlue);
      histo2->SetMarkerStyle(33);
    } else { 
      histo2->SetLineColor(kRed+2);
      histo2->SetMarkerColor(kRed+2);
      histo2->SetMarkerStyle(25);
    }
    histo2->SetLineWidth(2);
    histo2->SetLineStyle(1);
    histo2->Scale(1./histo2->Integral());
  }

  TH1F* histo3;
  if(dc3 != "" && rhoc3 != "" && of3 != "") {
    if( filterAndNoiseTopo ){
      ntupleName = "CaloTopoClusters";
    }
    histo3 = histosFromNtuples(fileName3, variable, ntupleName, xmin, xmax, bins);
    if(filterAndNoiseTopo){
      histo3->SetLineColor(kCyan+4);
      histo3->SetMarkerColor(kCyan+4);
      histo3->SetMarkerStyle(22);
    } else {
      histo3->SetLineColor(kBlue);
      histo3->SetMarkerColor(kBlue);
      histo3->SetMarkerStyle(33);
    }
    histo3->SetMarkerSize(1.3);
    histo3->SetLineWidth(2);
    histo3->SetLineStyle(1);
    histo3->Scale(1./histo3->Integral());
  }

  histo1->SetMaximum(ymax);

  TString canvasName = "cProf_"+var_short+"_"+settings;
  TCanvas *c_var = new TCanvas(canvasName,canvasName, 50, 50, 600, 600);
  gStyle->SetOptStat(0);
  c_var->SetRightMargin(0.06);
  c_var->cd();
  histo1->Draw("PLE");

  histo1->Fit("gaus", "", "", xmin, xmax);
  TF1 *fit1 = histo1->GetFunction("gaus");
  fit1->SetLineColor(kRed); 
  if(filterAndNoiseTopo)  fit1->SetLineColor(kRed+2); 
  fit1->SetLineStyle(3);
  fit1->SetLineWidth(3);
  fit1->Draw("same");

  //TLegend* t1 = new TLegend(0.187291,0.7495652,0.8996656,0.9165217);
  TLegend* t1 = new TLegend(0.5919732,0.4730435,0.9130435,0.9286957);
  //t1->SetNColumns(3);
  TString header = detectorLabel+", "+energy+" GeV gamma";

  TString title1 = "No noise ";
  if(filterAndNoiseTopo)  title1 = "No noise, CLUE, >2#sigma_{noise} ";
  t1->SetHeader(header);
  t1->AddEntry(histo1, title1, "pl");
  t1->AddEntry((TObject*)0, TString::Format("#mu = %.3f #pm %.3f", fit1->GetParameter(1), fit1->GetParError(1)), "");
  t1->AddEntry((TObject*)0, TString::Format("#sigma = %.3f #pm %.3f", fit1->GetParameter(2), fit1->GetParError(2)), "");

  if(dc2 != "" && rhoc2 != "" && of2 != "") {
    histo2->Draw("samePLE");
    histo2->Fit("gaus", "", "", xmin, xmax);
    TF1 *fit2 = histo2->GetFunction("gaus");
    fit2->SetLineColor(kRed+2); 
    if(filterAndNoiseTopo)  fit2->SetLineColor(kBlue); 
    fit2->SetLineStyle(2);
    fit2->SetLineWidth(2);
    fit2->Draw("same");
    TString title2 = "No noise, >2#sigma_{noise} ";
    //title2 = "With noise, Sliding Window";
    if(dc2 != dc || (dc2 != dc3 && dc3 != "") ){
      title2 += ", d_{c} = "+dc2;
    }
    t1->AddEntry(histo2, title2, "pl");
    t1->AddEntry((TObject*)0, TString::Format("#mu = %.3f #pm %.3f", fit2->GetParameter(1), fit2->GetParError(1)), "");
    t1->AddEntry((TObject*)0, TString::Format("#sigma = %.3f #pm %.3f", fit2->GetParameter(2), fit2->GetParError(2)), "");
  }
  if(dc3 != "" && rhoc3 != "" && of3 != "") {
    histo3->Draw("samePLE");
    histo3->Fit("gaus", "", "", xmin, xmax);
    TF1 *fit3 = histo3->GetFunction("gaus");
    fit3->SetLineColor(kBlue);
    if(filterAndNoiseTopo)  fit3->SetLineColor(kCyan+4);
    fit3->SetLineStyle(4);
    fit3->SetLineWidth(1);
    fit3->Draw("same");
    TString title3 = " o_{f} = "+of3+", #rho_{c} = "+rhoc;
    title3 = "With noise, >2#sigma_{noise} ";
    if(filterAndNoiseTopo)  title3 = "With noise, topological";
    if(dc3 != dc || dc3 != dc2 ){
      title3 += ", d_{c} = "+dc3;
    }
    t1->AddEntry(histo3, title3, "pl");
    t1->AddEntry((TObject*)0, TString::Format("#mu = %.3f #pm %.3f", fit3->GetParameter(1), fit3->GetParError(1)), "");
    t1->AddEntry((TObject*)0, TString::Format("#sigma = %.3f #pm %.3f", fit3->GetParameter(2), fit3->GetParError(2)), "");
  }

  t1->Draw("same");
  c_var->Draw();

  TString plotName = folder + "/" + variable+"_" + settings;
  c_var->SaveAs(plotName+".pdf", "pdf");
  c_var->SaveAs(plotName+".eps", "eps");
  c_var->SaveAs(plotName+".png", "png");

  return;

}


void comparisonWithOtherClusteringAlgorithms(TString detector, TString variable, TString var_short, 
                                             TString energy, TString dc, TString rhoc, TString of,
                                             float xmin, float xmax, int bins, float ymax = 120, 
                                             TString comparedTo = "PandoraClusters", TString comparedToSecond = ""){

  TString settings = energy+"GeV_"+rhoc+"rhoc_"+of+"of_clue_filter2sigma";

  setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString detectorLabel = detlabelFromDet(detector);
  TString folder = folderFromDet(detector);
  TString fileName = folder+"/"+filenameFromDet(detector, energy, dc, rhoc, of);
  //TString folder = "/eos/home-e/ebrondol/SWAN_projects/CLUE/data/"+detector+"/70_110_theta/";
  //TString folder = "/eos/home-e/ebrondol/SWAN_projects/CLUE/data/"+detector+"/gpsProduction/";
  //TString fileName = folder + "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc+"dc_"+rhoc+"rhoc_"+of+"of_clue_filter2sigma.root";

  // Open input file
  auto f = TFile::Open(fileName);
  if (!f || f->IsZombie()) {
     return;
  }

  TString ntupleName = "CLUEClusters";
  if( detector == "lar"){
    ntupleName = "CaloClueClusters";
  }

  TNtuple *ntupleCLUE = (TNtuple*)f->Get(ntupleName);
  TNtuple *ntupleOther = (TNtuple*)f->Get(comparedTo);
  TNtuple *ntupleOther2 = (TNtuple*)f->Get(comparedToSecond);

  TString histoNameCLUE = "h_"+var_short+"_CLUE_"+energy+"GeV";
  TH1F* h_var_CLUE = new TH1F(histoNameCLUE, histoNameCLUE, bins, xmin, xmax);
  float ratio = (xmax-xmin)/bins;
  h_var_CLUE->GetXaxis()->SetTitle(variable);
  if(variable=="totEnergyHits"){
    h_var_CLUE->GetXaxis()->SetTitle("Total energy [GeV]");
    //h_var_CLUE->GetXaxis()->SetTitle(TString::Format("Total energy [GeV/%.2f]", ratio));
  } else if(variable=="totEnergyHits/MCEnergy"){
    h_var_CLUE->GetXaxis()->SetTitle("E_{reco}/E_{MC}");
  }
  h_var_CLUE->GetYaxis()->SetTitle("a.u.");
  h_var_CLUE->GetYaxis()->SetTitleOffset(1.4);
  h_var_CLUE->SetLineColor(kRed+1);
  h_var_CLUE->SetMarkerColor(kRed+1);
  h_var_CLUE->SetLineWidth(3);
  h_var_CLUE->SetLineStyle(1);
  h_var_CLUE->SetMarkerStyle(21);
  h_var_CLUE->SetStats(0);

  TString histoNameOther = "h_"+var_short+"_Other_"+energy+"GeV";
  TH1F* h_var_Other = new TH1F(histoNameOther, histoNameOther, bins, xmin, xmax);
  h_var_Other->SetLineColor(kCyan+2);
  h_var_Other->SetMarkerColor(kCyan+2);
  h_var_Other->SetLineWidth(2);
  h_var_Other->SetLineStyle(1);
  h_var_Other->SetMarkerStyle(25);
  h_var_Other->SetStats(0);

  TString histoNameOther2 = "h_"+var_short+"_Other2_"+energy+"GeV";
  TH1F* h_var_Other2 = new TH1F(histoNameOther2, histoNameOther2, bins, xmin, xmax);
  h_var_Other2->SetLineColor(kCyan+4);
  h_var_Other2->SetMarkerColor(kCyan+4);
  h_var_Other2->SetLineWidth(2);
  h_var_Other2->SetLineStyle(1);
  h_var_Other2->SetMarkerStyle(22);
  h_var_Other2->SetStats(0);

  TCanvas *cTrash = new TCanvas();//,200,10,700,780);
  if(comparedToSecond != "")
    ntupleOther2->Draw(variable + " >> "+histoNameOther2);
  ntupleOther->Draw(variable + " >> "+histoNameOther);
  ntupleCLUE->Draw(variable + " >> "+histoNameCLUE);

  h_var_CLUE->Scale(1./h_var_CLUE->Integral());
  h_var_Other->Scale(1./h_var_Other->Integral());
  h_var_Other2->Scale(1./h_var_Other2->Integral());
  h_var_CLUE->SetMaximum(ymax);

  h_var_Other->Fit("gaus", "", "", xmin, xmax);
  TF1 *fit_var_Other = h_var_Other->GetFunction("gaus");
  fit_var_Other->SetLineColor(kCyan+2); 
  fit_var_Other->SetLineStyle(3);
  fit_var_Other->SetLineWidth(2);
  fit_var_Other->Draw("same");

  h_var_CLUE->Fit("gaus", "", "", xmin, xmax);
  TF1 *fit_var_CLUE = h_var_CLUE->GetFunction("gaus");
  fit_var_CLUE->SetLineColor(kRed+2); 
  fit_var_CLUE->SetLineStyle(2);
  fit_var_CLUE->SetLineWidth(3);
  fit_var_CLUE->Draw("same");

  TLegend* t1 = new TLegend(0.632107,0.5486957,0.9130435,0.8973913);
  if(comparedToSecond != ""){
    t1 = new TLegend(0.5518395,0.4591304,0.9130435,0.8973913);
  }
  TString header = "#splitline{"+detectorLabel+", "+energy+" GeV gamma}";
  header += "{o_{f} = "+of+", #rho_{c} = "+rhoc+"}";
  t1->SetHeader(header);
  t1->AddEntry(h_var_CLUE, "CLUE", "pl");
  TString fitResults = "#splitline{"+TString::Format("#mu = %.3f #pm %.3f", fit_var_CLUE->GetParameter(1), fit_var_CLUE->GetParError(1))+"}";
  fitResults = fitResults + "{"+TString::Format("#sigma = %.3f #pm %.3f", fit_var_CLUE->GetParameter(2), fit_var_CLUE->GetParError(2))+"}";
  t1->AddEntry((TObject*)0, fitResults, "");
  if(comparedTo == "CaloClusters"){
    t1->AddEntry(h_var_Other, "Sliding Window", "pl");
  } else {
    t1->AddEntry(h_var_Other, comparedTo, "pl");
  }
  TString fitResultsOther = "#splitline{"+TString::Format("#mu = %.3f #pm %.3f", fit_var_Other->GetParameter(1), fit_var_Other->GetParError(1))+"}";
  fitResultsOther = fitResultsOther + "{"+TString::Format("#sigma = %.3f #pm %.3f", fit_var_Other->GetParameter(2), fit_var_Other->GetParError(2))+"}";
  t1->AddEntry((TObject*)0, fitResultsOther, "");
  if(comparedToSecond != ""){
    h_var_Other2->Fit("gaus", "", "", xmin, xmax);
    TF1 *fit_var_Other2 = h_var_Other2->GetFunction("gaus");
    fit_var_Other2->SetLineColor(kBlue+2); 
    fit_var_Other2->SetLineStyle(4);
    fit_var_Other2->SetLineWidth(1);
    fit_var_Other2->Draw("same");

    if(comparedToSecond == "CaloTopoClusters"){
      t1->AddEntry(h_var_Other2, "Topological", "pl");
    } else {
      t1->AddEntry(h_var_Other2, comparedToSecond, "pl");
    }
    TString fitResultsOther2 = "#splitline{"+TString::Format("#mu = %.3f #pm %.3f", fit_var_Other2->GetParameter(1), fit_var_Other2->GetParError(1))+"}";
    fitResultsOther2 = fitResultsOther2 + "{"+TString::Format("#sigma = %.3f #pm %.3f", fit_var_Other2->GetParameter(2), fit_var_Other2->GetParError(2))+"}";
    t1->AddEntry((TObject*)0, fitResultsOther2, "");
  }

  TString canvasName = "cProf_"+var_short+"_"+settings;
  TCanvas *c_var = new TCanvas(canvasName,canvasName, 50, 50, 600, 600);
  gStyle->SetOptStat(0);
  c_var->SetRightMargin(0.06);
  c_var->cd();
  h_var_CLUE->Draw("");
  h_var_Other->Draw("same");
  if(comparedToSecond != ""){
    h_var_Other2->Draw("same");
  }
  t1->Draw("same");
  c_var->Draw();

  TString plotName = folder + "/" + variable;
  if(comparedToSecond != "")
    plotName += "_"+comparedToSecond+"_";
  plotName += "_"+comparedTo+"Comparison_" + settings;
  c_var->SaveAs(plotName+".png", "png");
  c_var->SaveAs(plotName+".pdf", "pdf");
  c_var->SaveAs(plotName+".eps", "eps");

  return;

}

void energyPlots_noise(){
  distributionComparison("lar", "totEnergyHits", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 0.4, 25, false, "40.0", "0.03", "3.0", "40.0", "0.03", "3.0");
  distributionComparison("lar", "totEnergyHits", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 0.4, 25, true, "40.0", "0.03", "3.0", "40.0", "0.03", "3.0");
  distributionComparison("lar", "totEnergyHits", "ene", "2",  "40.0", "0.03", "3.0", 1.2,  2.5, 0.35, 25, false, "40.0", "0.03", "3.0", "40.0", "0.03", "3.0");
  distributionComparison("lar", "totEnergyHits", "ene", "2",  "40.0", "0.03", "3.0", 1.2,  2.5, 0.35, 25, true, "40.0", "0.03", "3.0", "40.0", "0.03", "3.0");

  //distributionComparison("totEnergyHits", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 0.4, 25, "0.03", "3.0", "0.03", "3.0");
  //distributionComparison("totEnergyHits", "ene", "2", "40.0", "0.03", "3.0", 1.2, 2.5, 0.4, 25, "0.03", "3.0", "0.03", "3.0");
  //comparisonWithOtherClusteringAlgorithms("totEnergyHits", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 25, 0.4, "CaloClusters", "CaloTopoClusters");
  //comparisonWithOtherClusteringAlgorithms("totEnergyHits", "ene", "2", "40.0", "0.03", "3.0", 1.2, 2.5, 25, 0.4, "CaloClusters", "CaloTopoClusters");
  //comparisonWithOtherClusteringAlgorithms("totEnergy", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 25, 200, "CaloClusters");
  //comparisonWithOtherClusteringAlgorithms("totEnergy", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 25, 200, "CorrectedCaloClusters");
  //comparisonWithOtherClusteringAlgorithms("totEnergyHits", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 25, 200, "CaloTopoClusters");
  //comparisonWithOtherClusteringAlgorithms("totEnergy", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 25, 200, "CaloTopoClusters");
  //comparisonWithOtherClusteringAlgorithms("totEnergy", "ene", "10", "40.0", "0.03", "3.0", 7.5, 12.5, 25, 200, "CorrectedCaloTopoClusters");

  return;
}
