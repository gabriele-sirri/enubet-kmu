// usage:
//   $ root model.root
//   root[ ] testmodel(w)
//

void setTrueFrac(int n_par, double sigma, RooRealVar* f[], TMatrixDSym& cov_f);

int testmodel(RooWorkspace& w)
{
  gROOT->SetBatch(kTRUE);

  // extract elements from workspace (as pointers)
  auto evis = w.var("evis");
  auto z = w.var("z");

  // nuisance parameter
  auto alpha = w.var("alpha");

  RooInt* evis_bins = (RooInt*)w.genobj("evis_bins");
  RooInt* z_bins = (RooInt*)w.genobj("z_bins");

  evis->setBins(Int_t(*evis_bins));
  z->setBins(Int_t(*z_bins));

  auto model = (RooProdPdf*)w.pdf("model");

  // list of template pdfs: to determine number of parameters
  auto sb = (RooRealSumPdf*)w.pdf("sb");
  auto pdf_list = sb->funcList();

  auto n_pdf = pdf_list.getSize();
  auto n_params = n_pdf - 1;

  // get model parameters
  RooRealVar* f[n_params];
  for (auto p = 0; p < n_params; p++) f[p] = w.var(Form("f%d", p));

  // copy initial parameters values: this will be used for initialization,
  // they are the actual MC fractions
  vector<double> vf_in;
  for (auto p = 0; p < n_params; p++) vf_in.push_back(f[p]->getValV());

  // simulate b coefficients (weights) covariance matrix (and corresponding f
  // covariance matrix) and from this the actual b coefficients true value
  double sigma_perc = 0.2;
  TMatrixDSym Vf(n_params);

  setTrueFrac(n_params, sigma_perc, f, Vf);

  // true values of parameters for toy data generation
  //
  // f true values
  cout << "\n========== f fractions true values =========\n" << endl;
  for (auto p = 0; p < n_params; p++)
    cout << "f" << p << " = " << f[p]->getValV() << endl;

  // alpha true value
  TRandom r(0);
  double alpha_true = r.Gaus(0, 1);

  alpha->setVal(alpha_true);
  // alpha->setConstant(kTRUE);

  cout << "\n========== nuisance true values =========\n" << endl;
  cout << "alpha = " << alpha_true << endl;

  // Simulate a data-taking with a given number of POT
  double N_data_pot = 4.5e10;     // total POT in 1 y data-taking @ SPS
  double N_TL_pot = 7.97468e+07;  // total POT in TL production (current number
                                  // from TL_r5_v5_clean)
  double N_TL_K =
      16650;  // total K in TL production (current number from TL_r5_v5_clean)
  double N_prod_K = 16650;  // total K in production for systematics study (1x
                            // production from TL_r5_v5_clean)
  double N_TL_evt = 4095;   // total events in TL production (integral of
                            // distributions for neutrino and background)

  double N_data = N_TL_evt * (N_data_pot / N_TL_pot) *
                  (N_TL_K / N_prod_K);  // normalization factor to scale
                                        // variables to actual data-taking POT

  // generate toy data
  cout << "\nWill generate a total of " << N_data
       << " events, corresponding to " << N_data_pot << " POT\n"
       << endl;

  RooDataHist* data =
      model->generateBinned(RooArgSet{*evis, *z}, N_data, RooFit::Verbose(1));

  // set initial parameters values to MC fraction: used to perform the fit
  //
  // initialize f
  for (auto p = 0; p < n_params; p++) f[p]->setVal(vf_in[p]);

  // initialize alpha
  alpha->setVal(0.);

  // plot data and model overlaid before fitting
  //
  // evis var
  auto evis_frame_pref = evis->frame(RooFit::Title("Initialized model"));

  // draw data
  data->plotOn(evis_frame_pref);

  // draw total model
  model->plotOn(evis_frame_pref, RooFit::MoveToBack());

  // draw model components
  //
  // draw backgrouund first
  model->plotOn(evis_frame_pref, RooFit::Components("bkg_sys"),
                RooFit::FillColor(1), RooFit::FillStyle(3005), RooFit::VLines(),
                RooFit::DrawOption("F"), RooFit::MoveToBack());
  model->plotOn(evis_frame_pref, RooFit::Components("bkg_sys"),
                RooFit::FillColor(10), RooFit::Name("filled"),
                RooFit::DrawOption("F"), RooFit::MoveToBack());

  // draw signal
  string comp{"bkg_sys"};
  for (auto p = 0; p < n_pdf - 1; p++) {
    comp.append(Form(",pdf%d", p));

    model->plotOn(evis_frame_pref, RooFit::Components(comp.c_str()),
                  RooFit::FillColor(p + 40), RooFit::DrawOption("F"),
                  RooFit::MoveToBack());
  }

  // z var
  auto z_frame_pref = z->frame(RooFit::Title("Initialized model"));

  data->plotOn(z_frame_pref);

  model->plotOn(z_frame_pref, RooFit::MoveToBack());

  model->plotOn(z_frame_pref, RooFit::Components("bkg_sys"),
                RooFit::FillColor(1), RooFit::FillStyle(3005),
                RooFit::DrawOption("F"), RooFit::MoveToBack());
  model->plotOn(z_frame_pref, RooFit::Components("bkg_sys"),
                RooFit::FillColor(10), RooFit::Name("filled"),
                RooFit::DrawOption("F"), RooFit::MoveToBack());

  comp.clear();
  comp = "bkg_sys";
  for (auto p = 0; p < n_pdf - 1; p++) {
    comp.append(Form(",pdf%d", p));

    model->plotOn(z_frame_pref, RooFit::Components(comp.c_str()),
                  RooFit::FillColor(p + 40), RooFit::DrawOption("F"),
                  RooFit::MoveToBack());
  }

  // Prepare for fit
  //
  // parameters mean values: needed for constraint
  RooConstVar* mu_f[n_params];
  for (auto p = 0; p < n_params; p++)
    mu_f[p] = new RooConstVar{Form("mu_f%d", p), Form("mu_f%d", p), vf_in[p]};

  RooArgSet fset, mu_fset;
  for (auto p = 0; p < n_params; p++) {
    fset.add(*f[p]);
    mu_fset.add(*mu_f[p]);
  }

  // constraint on model parameters
  RooMultiVarGaussian c_Vf("c_Vf", "c_Vf", fset, mu_fset, Vf);

  // fit model to data
  auto fitresult = model->fitTo(*data, RooFit::ExternalConstraints(c_Vf),
                                RooFit::Constrain(*alpha), RooFit::Save(),
                                RooFit::PrintLevel(3));

  // display results
  fitresult->Print();

  // plot data and model overlaid after fitting
  //
  // evis var
  auto evis_frame_postf = evis->frame(RooFit::Title("Fitted model"));

  data->plotOn(evis_frame_postf);

  model->plotOn(evis_frame_postf, RooFit::MoveToBack());

  model->plotOn(evis_frame_postf, RooFit::Components("bkg_sys"),
                RooFit::FillColor(1), RooFit::FillStyle(3005),
                RooFit::DrawOption("F"), RooFit::MoveToBack());
  model->plotOn(evis_frame_postf, RooFit::Components("bkg_sys"),
                RooFit::FillColor(10), RooFit::Name("filled"),
                RooFit::DrawOption("F"), RooFit::MoveToBack());

  comp.clear();
  comp = "bkg_sys";
  for (auto p = 0; p < n_pdf - 1; p++) {
    comp.append(Form(",pdf%d", p));

    model->plotOn(evis_frame_postf, RooFit::Components(comp.c_str()),
                  RooFit::FillColor(p + 40), RooFit::DrawOption("F"),
                  RooFit::MoveToBack());
  }

  // z var
  auto z_frame_postf = z->frame(RooFit::Title("Fitted model"));

  data->plotOn(z_frame_postf);

  model->plotOn(z_frame_postf, RooFit::MoveToBack());

  model->plotOn(z_frame_postf, RooFit::Components("bkg_sys"),
                RooFit::FillColor(1), RooFit::FillStyle(3005),
                RooFit::DrawOption("F"), RooFit::MoveToBack());
  model->plotOn(z_frame_postf, RooFit::Components("bkg_sys"),
                RooFit::FillColor(10), RooFit::Name("filled"),
                RooFit::DrawOption("F"), RooFit::MoveToBack());

  comp.clear();
  comp = "bkg_sys";
  for (auto p = 0; p < n_pdf - 1; p++) {
    comp.append(Form(",pdf%d", p));

    model->plotOn(z_frame_postf, RooFit::Components(comp.c_str()),
                  RooFit::FillColor(p + 40), RooFit::DrawOption("F"),
                  RooFit::MoveToBack());
  }

  // Draw frame on canvas
  //
  // evis var
  TCanvas* c_evis = new TCanvas("c_evis", "Data vs model", 800, 400);
  c_evis->Divide(2);

  c_evis->cd(1);

  gPad->SetLeftMargin(0.15);
  evis_frame_pref->GetYaxis()->SetTitleOffset(1.4);
  evis_frame_pref->Draw();

  c_evis->cd(2);

  gPad->SetLeftMargin(0.15);
  evis_frame_postf->GetYaxis()->SetTitleOffset(1.4);
  evis_frame_postf->Draw();

  c_evis->SaveAs("plots/evis_observable.pdf");

  // z var
  TCanvas* c_z = new TCanvas("c_z", "Data vs model", 800, 400);
  c_z->Divide(2);

  c_z->cd(1);

  gPad->SetLeftMargin(0.15);
  z_frame_pref->GetYaxis()->SetTitleOffset(1.4);
  z_frame_pref->Draw();

  c_z->cd(2);

  gPad->SetLeftMargin(0.15);
  z_frame_postf->GetYaxis()->SetTitleOffset(1.4);
  z_frame_postf->Draw();

  c_z->SaveAs("plots/z_observable.pdf");

  // Plot likelihood profiles
  //

  // Construct likelihood
  //
  // -> it will be binned likelihood, since data are binned (see:
  // https://root.cern/doc/master/classRooAbsPdf.html#a44f3ac4e956a5ba4abb0b74841f22704)
  RooAbsReal* nll = model->createNLL(*data, RooFit::NumCPU(20));

  // Construct profile likelihood for the parameters
  for (auto p = 0; p < n_params; p++) {
    // set appropriate range to scan
    RooRealVar* par =
        (RooRealVar*)(fitresult->floatParsFinal()).find(Form("f%d", p));

    double val = par->getValV();
    double er = par->getError();

    f[p]->setRange(val - 4 * er, val + 4 * er);

    RooAbsReal* pll = nll->createProfile(*f[p]);

    RooPlot* frame = f[p]->frame(RooFit::Title("Profile likelihood"));

    pll->plotOn(frame, RooFit::ShiftToZero(), RooFit::LineColor(kRed));

    frame->SetMinimum(0);
    frame->SetMaximum(8);

    TCanvas* c_pll = new TCanvas("c_pll", "Profile likelihood", 800, 400);
    c_pll->cd();

    frame->Draw();

    c_pll->SaveAs(Form("plots/pll_par_f%d.pdf", p));

    delete c_pll;
  }

  return 1;
}

void setTrueFrac(int n_par, double sigma_perc, RooRealVar* f[],
                 TMatrixDSym& cov_f)
{
  // simulate values for cov matrix
  for (auto i = 0; i < n_par; i++)
    for (auto j = 0; j < n_par; j++) {
      double sigma_i = sigma_perc * f[i]->getValV();
      double sigma_j = sigma_perc * f[j]->getValV();

      if (i == j)
        cov_f(i, j) = 1.0 * sigma_i * sigma_i;
      else if (j == i + 1 || j == i - 1)
        cov_f(i, j) = 0.2 * sigma_i * sigma_j;
      else if (j == i + 2 || j == i - 2)
        cov_f(i, j) = 0.05 * sigma_i * sigma_j;
      else
        cov_f(i, j) = 0.;
    }

  // We want to get variables y with covariance matrix M
  //
  // BK decomposition: matrix M = U*V*U^T: where
  //
  // - M is the covariance  matriix of y variables;
  // - U is the transformation matrx to apply to the the uncorrelated x
  // variables;
  // - V is the diagonal matrix containing the standard deviatons of the
  // uncorrelated
  //   variables x;
  //
  // To obtain the y, we generate the x varriables with standard deviaton given
  // by V, then transform these using matrix U
  //

  // Decompose the M matrix
  TDecompSVD* svd = new TDecompSVD(cov_f);
  bool ret = svd->Decompose();

  // get the U matrix
  TMatrixD* U = new TMatrixD(n_par, n_par);
  *U = svd->GetU();

  // get the S matrix: in this case is equal to the U matrix, since M is
  // simmetric
  TVectorD* S = new TVectorD(n_par);
  *S = svd->GetSig();

  double w[n_par];

  TRandom* rnd = new TRandom(0);

  // generate y variables using the procedure described
  for (int p = 0; p < n_par; p++)
    w[p] = rnd->Gaus(0, TMath::Sqrt((double)(*S)[p]));

  TVectorD wDec(0, n_par - 1, w);
  // The vector of correlated coefficient to be used to simulate data
  TVectorD fcor = *U * wDec;

  for (int p = 0; p < n_par; p++) fcor[p] += f[p]->getValV();

  // set true values for fractions

  for (int p = 0; p < n_par; p++) f[p]->setVal(fcor[p]);

  // FIXME: for debuging purposes
  // f[0]->setVal(0.0107923);
  // f[1]->setVal(0.286132);
  // f[2]->setVal(0.0977188);

  delete svd;
  delete U;
  delete S;
  delete rnd;
}
