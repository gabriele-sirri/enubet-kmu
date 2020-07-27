// usage:
//   $ root model.root
//   root[ ] testmodel(w)
//

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

  RooArgSet pdfset;
  for (auto p = 0; p < pdf_list.getSize(); p++)
    pdfset.add(*((RooHistPdf*) pdf_list.find(Form("pdf%d", p))));

  // generate toy data
  auto data =
      model->generateBinned(RooArgSet{*evis, *z}, 160000, RooFit::Verbose(1));

  // plot data and model overlaid before fitting
  //
  // evis var
  auto evis_frame_pref = evis->frame(RooFit::Title("Pre-fit model"));

  data->plotOn(evis_frame_pref);

  model->plotOn(evis_frame_pref);

  RooArgSet pdf_set;
  for (auto p = 0; p < pdf_list.getSize(); p++) {
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
  for (auto p = 0; p < pdf_list.getSize(); p++) {
    pdf_set.add(*((RooHistPdf*) pdf_list.find(
        Form("pdf%d", p))));  // in order to draw stacked
    model->plotOn(z_frame_pref, RooFit::Components(pdf_set),
                  RooFit::LineStyle(kDotted), RooFit::FillColor(p + 40),
                  RooFit::DrawOption("F"), RooFit::MoveToBack());
  }

  // fit model to data
  auto fitresult = model->fitTo(*data, RooFit::Save(), RooFit::PrintLevel(3));

  // display results
  fitresult->Print();

  // plot data and model overlaid after fitting
  //
  // evis var
  auto evis_frame_postf = evis->frame(RooFit::Title("Post-fit model"));

  data->plotOn(evis_frame_postf);

  model->plotOn(evis_frame_postf);

  pdf_set.removeAll();
  for (auto p = 0; p < pdf_list.getSize(); p++) {
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
  for (auto p = 0; p < pdf_list.getSize(); p++) {
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
  // list of normalization coefficients
  RooRealVar* f[pdf_list.getSize() - 1];
  for (auto p = 0; p < pdf_list.getSize() - 1; p++)
    f[p] = w.var(Form("f%d", p));

  // Construct likelihood
  //
  // -> it will be binned likelihood, since data are binned (see:
  // https://root.cern/doc/master/classRooAbsPdf.html#a44f3ac4e956a5ba4abb0b74841f22704)
  RooAbsReal* nll = model->createNLL(*data, RooFit::NumCPU(2));

  // Construct profile likelihood for the parameters
  for (auto p = 0; p < pdf_list.getSize() - 1; p++) {
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
  }

  return 1;
}
