// usage: 
//   $ root model.root
//   root[ ] testmodel(w)
//

int testmodel(RooWorkspace& w)
{
  // extract elements from workspace (as pointers)
  
  auto evis = w.var("evis");
  auto z = w.var("z");
  
  auto model = w.pdf("model");
  
  // generate toy data
  auto data = model->generate(RooArgSet{*evis, *z}, 160000);
  
  // fit model to data
  auto fitresult = model->fitTo(*data, Save());
  
  // display results
  fitresult->Print();
  
  return 1;
}
