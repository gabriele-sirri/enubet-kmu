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
  //
  // signal
  vector<TH2F*> hmu_sgn_templ;

  for (auto t = 0; t < NBinsNu; t++) {
    TH2F* h_temp;
    inputFile->GetObject(Form("hMu_sgn_obs_nuBin%d", t + 1), h_temp);

    if (h_temp->Integral() != 0) hmu_sgn_templ.push_back(h_temp);
  }

  // background
  TH2F* hmu_bkg_templ;
  inputFile->GetObject("hMu_bkg_obs", hmu_bkg_templ);

  // must add bkg events to tot number of events
  tot_evts += hmu_bkg_templ->Integral();

  // total number of parameters: number of templates minus one
  // since we are not considering the extended fit
  int n_params = hmu_sgn_templ.size();

  // check template validity and determine actual parameters for model
  double integr(0);
  vector<double> f_val;
  for (auto t = 0; t < n_params; t++) {
    integr = hmu_sgn_templ[t]->Integral();

    f_val.push_back(integr / tot_evts);
  }

  integr = hmu_bkg_templ->Integral();
  f_val.push_back(integr / tot_evts);

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
  RooDataHist* dhist[n_params];

  for (auto p = 0; p < n_params; p++) {
    hmu_sgn_templ[p]->Scale(1. / hmu_sgn_templ[p]->Integral());

    dhist[p] = new RooDataHist{Form("dh%d", p), Form("dh%d", p),
                               RooArgSet{evis, z}, hmu_sgn_templ[p]};
  }

  // signal pdfs from RooDataHist
  RooHistPdf* pdf[n_params];

  // for (auto p = 0; p < n_params + 1; p++)
  for (auto p = 0; p < n_params; p++)
    pdf[p] = new RooHistPdf{Form("pdf%d", p), Form("pdf%d", p),
                            RooArgSet{evis, z}, *dhist[p]};

  // background pdf: model shape systematic
  //
  // nominal value
  TH2F* hmu_bkg_nom = (TH2F*)hmu_bkg_templ->Clone("hmu_bkg_nom");
  hmu_bkg_nom->Scale(1. / hmu_bkg_nom->Integral());
  RooDataHist dh_mu_bkg_nom{"dh_mu_bkg_nom", "dh_mu_bkg_nom",
                            RooArgSet{evis, z}, hmu_bkg_nom};
  RooHistPdf bkg_0("bkg_0", "bkg_0", RooArgSet{evis, z}, dh_mu_bkg_nom);
  //
  // bkg +1 sigma
  TH2F* hmu_bkg_p;
  inputFile->GetObject("hMu_bkg_p_templ", hmu_bkg_p);
  hmu_bkg_p->Scale(1. / hmu_bkg_p->Integral());
  RooDataHist dh_mu_bkg_p{"dh_mu_bkg_p", "dh_mu_bkg_p", RooArgSet{evis, z},
                          hmu_bkg_p};
  RooHistPdf bkg_p("bkg_p", "bkg_p", RooArgSet{evis, z}, dh_mu_bkg_p);
  //
  // bkg -1 sigma
  TH2F* hmu_bkg_m;
  inputFile->GetObject("hMu_bkg_m_templ", hmu_bkg_m);
  hmu_bkg_m->Scale(1. / hmu_bkg_m->Integral());
  RooDataHist dh_mu_bkg_m{"dh_mu_bkg_m", "dh_mu_bkg_m", RooArgSet{evis, z},
                          hmu_bkg_m};
  RooHistPdf bkg_m("bkg_m", "bkg_m", RooArgSet{evis, z}, dh_mu_bkg_m);
  //
  // build interpolation between -1 -> +1 sigma
  RooRealVar alpha{"alpha", "alpha bkg sys", 0, -5, 5};
  PiecewiseInterpolation bkg_sys("bkg_sys", "bkg_sys", bkg_0, bkg_m, bkg_p,
                                 alpha);

  // constrain alpha
  RooGaussian gaus_alpha("gaus_alpha", "gaus_alpha", alpha, RooFit::RooConst(0),
                         RooFit::RooConst(1));

  // assemble all components to fnal model
  RooArgList pdflist;
  for (auto p = 0; p < n_params; p++) pdflist.add(*pdf[p]);

  pdflist.add(bkg_sys);

  RooArgList parlist;
  for (auto p = 0; p < n_params; p++) parlist.add(*f[p]);

  // define the model (not extended)
  RooRealSumPdf sb{"sb", "sb", pdflist, parlist};
  RooProdPdf model{"model", "model", sb, gaus_alpha};

  //
  RooWorkspace w{"w", "w"};
  w.import(evis_bins);
  w.import(z_bins);
  w.import(model);
  w.writeToFile("model.root");

  for (auto p = 0; p < n_params; p++) {
    delete f[p];
    delete pdf[p];
    delete dhist[p];
  }

  return 0;
}
