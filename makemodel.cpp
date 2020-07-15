// usage: 
//   $ root model.root
//   root[ ] makemodel(w)
//

int makemodel() 
{
  // define observables
  RooRealVar evis{"evis", "visible energy", 0,	400, "MeV"};
  RooRealVar z{"z", "Impact Position", 0, 4400, "cm"};
  
  // define parameters 
  RooRealVar f0{"f0", "f0", 0.1, 0, 1};
  RooRealVar f1{"f1", "f1", 0.1, 0, 1};
  RooRealVar f2{"f2", "f2", 0.1, 0, 1};
  RooRealVar f3{"f3", "f3", 0.1, 0, 1};
  RooRealVar f4{"f4", "f4", 0.1, 0, 1};
  RooRealVar f5{"f5", "f5", 0.1, 0, 1};
  RooRealVar f6{"f6", "f6", 0.1, 0, 1};
  RooRealVar f7{"f7", "f7", 0.1, 0, 1};
  RooRealVar f8{"f8", "f8", 0.1, 0, 1};
  RooRealVar f9{"f9", "f9", 0.1, 0, 1};
  
  // pdf of evis,z from histograms
  RooHistPdf pfd00{"pfd00", "pfd00", ... };  
  RooHistPdf pfd01{"pfd01", "pfd01", ... };
  RooHistPdf pfd02{"pfd02", "pfd02", ... };  
  RooHistPdf pfd03{"pfd03", "pfd03", ... };
  RooHistPdf pfd04{"pfd04", "pfd04", ... };  
  RooHistPdf pfd05{"pfd05", "pfd05", ... };
  RooHistPdf pfd06{"pfd06", "pfd06", ... };  
  RooHistPdf pfd07{"pfd07", "pfd07", ... };
  RooHistPdf pfd08{"pfd08", "pfd08", ... };  
  RooHistPdf pfd09{"pfd09", "pfd09", ... };  
  RooHistPdf pfd10{"pfd10", "pfd10", ... };
  
  RooArgSet pdfset{pfd00, pfd01, pdf02, pdf03, pdf04, pdf05, pdf06, pdff07, pdf08,  pdf09};
  pdfset.add(pdf10);
  
  RooArgSet parset{f0, f1, f2, f3, f4, f5, f6, f7, f8, f9};
   
  // define the model (not extended)
  RooAddPdf model{"model", "model", pdfset, parset};
  
  // 
  RooWorkspace w{"w", "w"};
  w.import(model);
  w.writeToFile("model.root");
  
  return 0;
}

  
  