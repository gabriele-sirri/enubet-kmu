// usage: 
//   $ root model.root
//   root[ ] makemodel(w)
//

int makemodel(string filepath) 
{
	// file containing templates
	TFile* inputFile = new TFile(filepath.c_str());

	// get templates from file
	int NBinsNu = 10;
	TH2F* hmu_sgn_templ[NBinsNu];

	for(auto t=0; t<NBinsNu; t++)
		inputFile -> GetObject(Form("hMu_sgn_obs_nuBin%d", t+1), hmu_sgn_templ[t]);

	// define observables
	RooRealVar evis{"evis", "visible energy", 0,	400, "MeV"};
	RooRealVar z{"z", "Impact Position", 0, 4400, "cm"};

	// define parameters 
	RooRealVar f0{"f0", "f0", 	  0,     0,  0.1};
	RooRealVar f1{"f1", "f1", 	  0,    10,  300};
	RooRealVar f2{"f2", "f2",  1000,  1500,  2500};
	RooRealVar f3{"f3", "f3",  8000, 12000, 20000};
	RooRealVar f4{"f4", "f4", 15000, 28000, 40000};
	RooRealVar f5{"f5", "f5", 20000, 30000, 50000};
	RooRealVar f6{"f6", "f6", 20000, 30000, 50000};
	RooRealVar f7{"f7", "f7", 15000, 20000, 40000};
	RooRealVar f8{"f8", "f8",   500,  1400,  4000};

	// template RooDataHist: needed to build pdf
	RooDataHist dhist00{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[0]};
	RooDataHist dhist01{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[1]};
	RooDataHist dhist02{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[2]};
	RooDataHist dhist03{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[3]};
	RooDataHist dhist04{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[4]};
	RooDataHist dhist05{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[5]};
	RooDataHist dhist06{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[6]};
	RooDataHist dhist07{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[7]};
	RooDataHist dhist08{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[8]};
	RooDataHist dhist09{"dh00", "dh00", RooArgSet{evis, z}, hmu_sgn_templ[9]};

	// pdf of evis,z from histograms
	RooHistPdf pdf00{"pdf00", "pdf00", RooArgSet{evis, z}, dhist00, 2};  
	RooHistPdf pdf01{"pdf01", "pdf01", RooArgSet{evis, z}, dhist01, 2};
	RooHistPdf pdf02{"pdf02", "pdf02", RooArgSet{evis, z}, dhist02, 2};  
	RooHistPdf pdf03{"pdf03", "pdf03", RooArgSet{evis, z}, dhist03, 2};
	RooHistPdf pdf04{"pdf04", "pdf04", RooArgSet{evis, z}, dhist04, 2};  
	RooHistPdf pdf05{"pdf05", "pdf05", RooArgSet{evis, z}, dhist05, 2};
	RooHistPdf pdf06{"pdf06", "pdf06", RooArgSet{evis, z}, dhist06, 2};  
	RooHistPdf pdf07{"pdf07", "pdf07", RooArgSet{evis, z}, dhist07, 2};
	RooHistPdf pdf08{"pdf08", "pdf08", RooArgSet{evis, z}, dhist08, 2};  

	RooArgSet pdfset{pdf00, pdf01, pdf02, pdf03, pdf04, pdf05, pdf06, pdf07, pdf08};

	RooArgSet parset{f0, f1, f2, f3, f4, f5, f6, f7, f8};

	// define the model (not extended)
	RooAddPdf model{"model", "model", pdfset, parset};

	// 
	RooWorkspace w{"w", "w"};
	w.import(model);
	w.writeToFile("model.root");

	return 0;
}



