int nLayersFromDet(TString detector){
  if(detector == "cld" || detector == "clicdet")
    return 40;
  else if(detector == "lar")
    return 12;
  else
    return 1;
}

TString folderFromDet(TString detector){
  TString currentFold = "/eos/home-e/ebrondol/SWAN_projects/CLUE/data/"+detector+"/";
  if(detector == "cld")
    //return currentFold+"/gpsProduction/tuningPars/";
    return currentFold+"/gpsProduction/";
  else if(detector == "clicdet")
    return currentFold+"/gpsProduction/";
  else if(detector == "lar")
    //return currentFold+"/tuningPars/";
    //return currentFold+"/timingStudies/";
    return currentFold+"/50_130_theta_I_think/";
  else
    return currentFold;
}

TString filenameFromDet(TString detector, TString energy, TString dc, TString rhoc, TString of){
  if(detector == "cld" || detector == "clicdet")
    return "k4clue_"+detector+"_output_gamma_"+energy+"GeV_uniform_500events_"+dc+"dc_"+rhoc+"rhoc_"+of+"of.root";
  else if(detector == "lar")
    return "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc+"dc_"+rhoc+"rhoc_"+of+"of.root";
    //return "k4clue_"+detector+"_output_gamma_"+energy+"GeV_barrel_500events_"+dc+"dc_"+rhoc+"rhoc_"+of+"of_timingStudies16.root";
  else
    return "k4clue_output.root";
}

TString detlabelFromDet(TString detector){
  TString detLabel = "Detector";
  if(detector == "cld")
    return "CLD";
  else if(detector == "clicdet")
    return "CLICdet";
  else if(detector == "lar")
    return "LAr";
  else
    return detLabel;
}


std::vector< TProfile*> returnLayerProfiles(TString fileName, float ymax, int nLayers){

  std::cout << "Read file: " << fileName << std::endl;
  std::vector<TProfile*> product{};
  auto f = TFile::Open(fileName);
  if (!f || f->IsZombie()) {
     return product;
  }

  TTree* treeHits = (TTree*)f->Get("CLUEHits");
  std::vector<int> *status = 0;
  std::vector<int> *layer = 0;
  std::vector<int> *region = 0;
  treeHits->SetBranchAddress("status",&status);
  treeHits->SetBranchAddress("layer",&layer);
  treeHits->SetBranchAddress("region",&region);
  Int_t nentries = (Int_t)treeHits->GetEntries();
//  std::cout << nentries << " events to read " << std::endl;

  float xmin = 0;
  float xmax = float(nLayers);
  int bins   = nLayers;

  TString profileNameSeeds = "p_layer_seeds";
  TProfile *p_var_Seeds = new TProfile("p_layer_seeds", "p_layer_seeds", bins, xmin, xmax, 0, ymax);
  TProfile *p_var_Followers = new TProfile("p_layer_followers", "p_layer_followers", bins, xmin, xmax, 0, ymax);
  TProfile *p_var_Outliers = new TProfile("p_layer_outliers", "p_layer_outliers", bins, xmin, xmax, 0, ymax);

  for (Int_t i = 0; i < nentries; i++){
    //if(i > 10) break;
    std::vector<int> seedsPerLayer(40, 0);
    std::vector<int> outliersPerLayer(40, 0);
    std::vector<int> followersPerLayer(40, 0);
//    std::cout << "event #" << i << std::endl; 
    treeHits->GetEntry(i);

    int countS = std::count(status->begin(), status->end(), 2);
    int countO = std::count(status->begin(), status->end(), 0);
    int countF = std::count(status->begin(), status->end(), 1);
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
    for (auto ilay = 0; ilay < int(seedsPerLayer.size()); ilay++){
      p_var_Seeds->Fill(ilay, seedsPerLayer[ilay], 1);
      p_var_Followers->Fill(ilay, followersPerLayer[ilay], 1);
      p_var_Outliers->Fill(ilay, outliersPerLayer[ilay], 1);
    }
  }

  product.push_back(p_var_Outliers);
  product.push_back(p_var_Followers);
  product.push_back(p_var_Seeds);

  return product;
}

TH1F* histosFromNtuples(TString fileName, TString variable, TString ntupleName,
                        float xmin, float xmax, int bins){

  // Open input file
  std::cout << "Read file: " << fileName << std::endl;
  auto f = TFile::Open(fileName);

  std::cout << "Ntuple name: " << ntupleName << std::endl;
  TNtuple *ntupleCLUE = (TNtuple*)f->Get(ntupleName);
  TString histoNameCLUE = "h_CLUE";
  TH1F* h_var_CLUE = new TH1F(histoNameCLUE, histoNameCLUE, bins, xmin, xmax);

  TCanvas *cTrash = new TCanvas();//,200,10,700,780);
  ntupleCLUE->Draw(variable + " >> "+histoNameCLUE);
  //ntupleCLUE->Draw(" energy >> "+histoNameCLUE, "energy > 1");

  return h_var_CLUE;
}

