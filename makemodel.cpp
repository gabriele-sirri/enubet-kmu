// usage:
//   $ root model.root
//   root[ ] makemodel(w)
//

int makemodel(string filepath)
{
  // file containing templates
  TFile* inputFile = new TFile(filepath.c_str());

  // get neutrino spectrum
  TH1F* hnu_e;
  inputFile->GetObject("hNuE", hnu_e);

  int NBinsNu = hnu_e->GetNbinsX();

  double tot_evts = hnu_e->Integral();

  // get templates from file

  vector<TH2F*> hmu_sgn_templ;

  for (auto t = 0; t < NBinsNu; t++) {
    TH2F* h_temp;
    inputFile->GetObject(Form("hMu_sgn_obs_nuBin%d", t + 1), h_temp);

    if (h_temp->Integral() != 0) hmu_sgn_templ.push_back(h_temp);
  }

  // total number of parameters: number of templates minus one
  // since we are not considering the extended fit
  int n_params = hmu_sgn_templ.size() - 1;

  // check template validity and determine actual parameters for model
  double integr(0);
  vector<double> f_val;
  for (auto t = 0; t < n_params; t++) {
    integr = hmu_sgn_templ[t]->Integral();

    f_val.push_back(integr / tot_evts);
  }

  // save number of bins in variable
  RooInt evis_bins(hmu_sgn_templ[0]->GetNbinsX());
  evis_bins.SetName("evis_bins");

  RooInt z_bins(hmu_sgn_templ[0]->GetNbinsY());
  z_bins.SetName("z_bins");

  // define observables
  RooRealVar evis{"evis", "visible energy", 0, 400, "MeV"};
  RooRealVar z{"z", "Impact Position", 0, 4400, "cm"};

  // define parameters
  RooRealVar* f[n_params];

  for (auto p = 0; p < n_params; p++)
    f[p] = new RooRealVar{Form("f%d", p), Form("f%d", p), f_val[p], 0., 1.};

  // template RooDataHist: needed to build pdf
  RooDataHist* dhist[n_params + 1];

  for (auto p = 0; p < n_params + 1; p++)
    dhist[p] = new RooDataHist{Form("dh%d", p), Form("dh%d", p),
                               RooArgSet{evis, z}, hmu_sgn_templ[p]};

  // pdf of evis,z from histograms
  RooHistPdf* pdf[n_params + 1];

  for (auto p = 0; p < n_params + 1; p++)
    pdf[p] = new RooHistPdf{Form("pdf%d", p), Form("pdf%d", p),
                            RooArgSet{evis, z}, *dhist[p], 2};

  RooArgSet pdfset;
  for (auto p = 0; p < n_params + 1; p++) pdfset.add(*pdf[p]);

  RooArgSet parset;
  for (auto p = 0; p < n_params; p++) parset.add(*f[p]);

  // define the model (not extended)
  RooAddPdf model{"model", "model", pdfset, parset};

  //
  RooWorkspace w{"w", "w"};
  w.import(evis_bins);
  w.import(z_bins);
  w.import(model);
  w.writeToFile("model.root");

  for (auto p = 0; p < n_params + 1; p++) {
    if (p < n_params) delete f[p];

    delete dhist[p];
    delete pdf[p];
  }

  return 0;
}
