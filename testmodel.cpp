// usage:
//   $ root model.root
//   root[ ] testmodel(w)
//

TVectorD getTrueWeights(int n_par, double sigma, RooRealVar* f[],
                        TMatrixDSym& cov_f);

int testmodel(RooWorkspace& w)
{
  gROOT->SetBatch(kTRUE);

  // extract elements from workspace (as pointers)

  auto evis = w.var("evis");
  auto z = w.var("z");

  RooInt* evis_bins = (RooInt*) w.genobj("evis_bins");
  RooInt* z_bins = (RooInt*) w.genobj("z_bins");

  evis->setBins(Int_t(*evis_bins));
  z->setBins(Int_t(*z_bins));

  auto model = (RooAddPdf*) w.pdf("model");

  // list of model pdfs: useful for drawing
  auto pdf_list = model->pdfList();

  auto n_pdf = pdf_list.getSize();
  auto n_params = n_pdf - 1;

  RooArgSet pdfset;
  for (auto p = 0; p < pdf_list.getSize(); p++)
    pdfset.add(*((RooHistPdf*) pdf_list.find(Form("pdf%d", p))));

  // get model parameters
  RooRealVar* f[n_params];
  for (auto p = 0; p < n_params; p++) f[p] = w.var(Form("f%d", p));

  // copy initial parameters values: this will be used for initialization,
  // they are the actual MC fractions
  vector<double> vf_in;
  for (auto p = 0; p < n_params; p++) vf_in.push_back(f[p]->getValV());

  // simulate b coefficients (weights) covariance matrix (and corresponding f
  // covariance matrix) and from this the actual b coefficients true value
  double sigma_b = 0.2;
  TMatrixDSym Vf(n_params);
  TVectorD b_coef(n_params);

  // need to ensure that sum of fractions is less than 1
  int n_trial(1);
  double sum_f(0);
  bool condition(false);
  while (!condition) {
    cout << "Get b-coefficients: trial = " << n_trial << endl;

    b_coef = getTrueWeights(n_params, sigma_b, f, Vf);

    for (auto p = 0; p < n_params; p++) sum_f += b_coef[p] * vf_in[p];

    if (sum_f < 1.) condition = true;

    sum_f = 0;
    n_trial++;
  }

  // print b coefficients true values
  cout << " ========== b-coefficients true values =========" << endl;
  for (auto p = 0; p < n_params; p++)
    cout << "b" << p << " = " << b_coef[p] << endl;

  // set parameters value for toy data generation: weight the MC fractions with
  // b_coef true values
  cout << " ========== f fractions true values =========" << endl;
  for (auto p = 0; p < n_params; p++) {
    f[p]->setVal(b_coef[p] * vf_in[p]);
    cout << "f" << p << " = " << f[p]->getValV() << endl;
  }

  // generate toy data
  auto data =
      model->generateBinned(RooArgSet{*evis, *z}, 160000, RooFit::Verbose(1));

  // set initial values for parameters to MC fraction: used to perform the fit
  for (auto p = 0; p < n_params; p++) f[p]->setVal(vf_in[p]);

  // plot data and model overlaid before fitting
  //
  // evis var
  auto evis_frame_pref = evis->frame(RooFit::Title("Pre-fit model"));

  data->plotOn(evis_frame_pref);

  model->plotOn(evis_frame_pref);

  RooArgSet pdf_set;
  for (auto p = 0; p < n_pdf; p++) {
    pdf_set.add(*((RooHistPdf*) pdf_list.find(
        Form("pdf%d", p))));  // in order to draw stacked
    model->plotOn(evis_frame_pref, RooFit::Components(pdf_set),
                  RooFit::LineStyle(kDotted), RooFit::FillColor(p + 40),
                  RooFit::DrawOption("F"), RooFit::MoveToBack());
  }

  // z var
  auto z_frame_pref = z->frame(RooFit::Title("Pre-fit model"));

  data->plotOn(z_frame_pref);

  model->plotOn(z_frame_pref);

  pdf_set.removeAll();
  for (auto p = 0; p < n_pdf; p++) {
    pdf_set.add(*((RooHistPdf*) pdf_list.find(
        Form("pdf%d", p))));  // in order to draw stacked
    model->plotOn(z_frame_pref, RooFit::Components(pdf_set),
                  RooFit::LineStyle(kDotted), RooFit::FillColor(p + 40),
                  RooFit::DrawOption("F"), RooFit::MoveToBack());
  }

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
                                RooFit::Save(), RooFit::PrintLevel(3));

  // display results
  fitresult->Print();

  // plot data and model overlaid after fitting
  //
  // evis var
  auto evis_frame_postf = evis->frame(RooFit::Title("Post-fit model"));

  data->plotOn(evis_frame_postf);

  model->plotOn(evis_frame_postf);

  pdf_set.removeAll();
  for (auto p = 0; p < n_pdf; p++) {
    pdf_set.add(*((RooHistPdf*) pdf_list.find(
        Form("pdf%d", p))));  // in order to draw stacked
    model->plotOn(evis_frame_postf, RooFit::Components(pdf_set),
                  RooFit::LineStyle(kDotted), RooFit::FillColor(p + 40),
                  RooFit::DrawOption("F"), RooFit::MoveToBack());
  }

  // z var
  auto z_frame_postf = z->frame(RooFit::Title("Post-fit model"));

  data->plotOn(z_frame_postf);

  model->plotOn(z_frame_postf);

  pdf_set.removeAll();
  for (auto p = 0; p < n_pdf; p++) {
    pdf_set.add(*((RooHistPdf*) pdf_list.find(
        Form("pdf%d", p))));  // in order to draw stacked
    model->plotOn(z_frame_postf, RooFit::Components(pdf_set),
                  RooFit::LineStyle(kDotted), RooFit::FillColor(p + 40),
                  RooFit::DrawOption("F"), RooFit::MoveToBack());
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
  RooAbsReal* nll = model->createNLL(*data, RooFit::NumCPU(2));

  // Construct profile likelihood for the parameters
  for (auto p = 0; p < n_params; p++) {
    // set appropriate range to scan
    RooRealVar* par =
        (RooRealVar*) (fitresult->floatParsFinal()).find(Form("f%d", p));

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

TVectorD getTrueWeights(int n_par, double sigma, RooRealVar* f[],
                        TMatrixDSym& cov_f)
{
  // simulate values for cov matrix
  TMatrixDSym cov_b(n_par);
  for (auto i = 0; i < n_par; i++)
    for (auto j = 0; j < n_par; j++) {
      if (i == j)
        cov_b(i, j) = 1.0;
      else if (j == i + 1 || j == i - 1)
        cov_b(i, j) = 0.2;
      else if (j == i + 2 || j == i - 2)
        cov_b(i, j) = 0.05;
      else
        cov_b(i, j) = 0.;
    }

  for (auto i = 0; i < n_par; i++)
    for (auto j = 0; j < n_par; j++) {
      double fij = (f[i]->getValV()) * (f[j]->getValV());
      if (i == j)
        cov_f(i, j) = 1.0 * fij;
      else if (j == i + 1 || j == i - 1)
        cov_f(i, j) = 0.2 * fij;
      else if (j == i + 2 || j == i - 2)
        cov_f(i, j) = 0.05 * fij;
      else
        cov_f(i, j) = 0. * fij;
    }

  cov_f *= pow(sigma, 2);

  // Simulate observed data: to do this, simulate a given vaalue for the vector
  // of normalization coefficients b taking into account the flux correlation
  // matrix. Use these values to normalize the templates and sum them to obtain
  // total spetra for the observables:
  //
  // let's assume we have the above correlation matrix: starting from
  // uncorrelated values of the normalization coefficients b of the bins in the
  // neutrino energy spectrum, we can obtain correlated ones by applying the
  // change of coordinates matrix obtained from the Cholesky decomposition.
  //
  // Cholesky decomposition of covariance flux matrix (NOTE: actually this is
  // the matrix divided by sigma^2)
  TDecompChol* ch = new TDecompChol(cov_b);
  bool ret = ch->Decompose();

  TMatrixD* U = new TMatrixD(n_par, n_par);
  *U = ch->GetU();

  // weights
  double w[n_par];

  TRandom* rnd = new TRandom(0);
  for (int bn = 0; bn < n_par; bn++) w[bn] = rnd->Gaus(0, sigma);

  TVectorD wDec(0, n_par - 1, w);
  // The vector of correlated coefficient to be used to simulate data
  TVectorD wCor = *U * wDec;

  wCor += 1.;

  cov_b *= pow(sigma, 2);

  cov_b.Print();

  delete ch;
  delete U;
  delete rnd;

  return wCor;
}
