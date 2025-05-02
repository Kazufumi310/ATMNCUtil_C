#include <myheader/myRoot.h>
#include <myheader/myStandard.h>
#include <FileHandler.h>

int main(int argc, char** argv){
  if(argc<3){
    std::cout<<"arg error"<<std::endl;
    std::cout<<"usage: sumRoot input1.root input2.root ... inputN.root -o output.root"<<std::endl;
    std::cout<<"usage: sumRoot -i inputList.txt -o output.root"<<std::endl;
    return 0;
  }

  std::vector<std::string> inputRoots;
  std::string outputRoot="out.root";

  for(int iarg=1; iarg < argc; iarg++){
    std::string arg = argv[iarg];
    if(arg=="-o"){
      iarg++;
      outputRoot = argv[iarg];
    }else if(arg=="-i"){
      iarg++;
      std::string inputList = argv[iarg];
      std::ifstream ifs(inputList.c_str());
      std::string inputRoot;
      while(ifs>>inputRoot){
	if(ifs.eof()) break;
	inputRoots.push_back(inputRoot);
      }
      
    }else{
      inputRoots.push_back(argv[iarg]);
    }
  }

  FileHandler filer;
  filer.sum(inputRoots,outputRoot);
}



/*
void recordEachGraph(Histogrammer &hist,TFile &fout){
  int const nazim = hist.getNAzim();
  int const ncosz = hist.getNCosz();
  int const ncosz_HAKKM = 20;
  int const dicosz = ncosz/ncosz_HAKKM-1;
  int const isite = 0;
  int const nkind = 4;
  int const color[nkind]={4,3,2,6};
  fout.cd();
  for(int ik=0;ik<nkind;ik++){
    SpectrumHandler spec;
    std::cout<<"ik:"<<ik<<std::endl;
    for(int ic=0;ic<ncosz;ic+=dicosz+1){
      for(int ia=0;ia<nazim;ia++){
	hist.averageFlux(ik,isite,ic,ic+dicosz,ia,ia,spec);
	int nn = spec.getN();
	double *E = spec.getEneArray();
	double rate[1024],rateErr[1024];
	double flux[1024],fluxErr[1024];
	for(int i=0;i<nn;i++){
	  rate[i] = spec.getRate(E[i]);
	  rateErr[i] = spec.getRateErr(E[i]);
	  flux[i] = spec.getFlux(E[i]);
	  fluxErr[i] = spec.getFluxErr(E[i]);
	  //	  std::cout<<"i,r,f:"<<i<<" "<<rate[i]<<" "<<flux[i]<<std::endl;
	}
	TGraph *grate = new TGraphErrors(nn,E,rate,0,rateErr);
	grate->SetMarkerStyle(20);
	grate->SetMarkerSize(0.7);
	grate->SetName(Form("gmyrate_k%d_c%d_a%d",ik,ic/(dicosz+1),ia));
	grate->Write();

	TGraph *gflux = new TGraphErrors(nn,E,flux,0,fluxErr);
	gflux->SetMarkerStyle(20);
	gflux->SetMarkerSize(0.7);
	gflux->SetName(Form("gmyflux_k%d_c%d_a%d",ik,ic/(dicosz+1),ia));
	gflux->Write();
	delete grate;
	delete gflux;
      }
    }
  }
}

void recordAverageGraph(Histogrammer &hist,TFile &fout){
  int const nazim = hist.getNAzim();
  int const ncosz = hist.getNCosz();
  int const ncosz_HAKKM = 20;
  int const dicosz = ncosz/ncosz_HAKKM-1;
  int const isite = 0;
  int const nkind = 4;
  int const color[nkind]={4,3,2,6};
  fout.cd();
  for(int ik=0;ik<nkind;ik++){
    SpectrumHandler spec;
    std::cout<<"ik:"<<ik<<std::endl;
    for(int ic=0;ic<ncosz;ic+=dicosz+1){
      hist.averageFlux(ik,isite,ic,ic+dicosz,0,nazim-1,spec);
      int nn = spec.getN();
      double *E = spec.getEneArray();
      double rate[1024],rateErr[1024];
      double flux[1024],fluxErr[1024];
      for(int i=0;i<nn;i++){
	rate[i] = spec.getRate(E[i]);
	rateErr[i] = spec.getRateErr(E[i]);
	flux[i] = spec.getFlux(E[i]);
	fluxErr[i] = spec.getFluxErr(E[i]);
      }
      TGraph *grate = new TGraphErrors(nn,E,rate,0,rateErr);
      grate->SetMarkerStyle(20);
      grate->SetMarkerSize(0.7);
      grate->SetName(Form("gmyrateAve_k%d_c%d",ik,ic/(dicosz+1)));
      grate->Write();
      
      TGraph *gflux = new TGraphErrors(nn,E,flux,0,fluxErr);
      gflux->SetMarkerStyle(20);
      gflux->SetMarkerSize(0.7);
      gflux->SetName(Form("gmyfluxAve_k%d_c%d",ik,ic/(dicosz+1)));
      gflux->Write();
      delete grate;
      delete gflux;
    }
  }

  for(int ik=0;ik<nkind;ik++){
    SpectrumHandler spec;
    hist.averageFlux(ik,isite,0,ncosz-1,0,nazim-1,spec);
    int nn = spec.getN();
    double *E = spec.getEneArray();
    double rate[1024],rateErr[1024];
    double flux[1024],fluxErr[1024];
    for(int i=0;i<nn;i++){
      rate[i] = spec.getRate(E[i]);
      rateErr[i] = spec.getRateErr(E[i]);
      flux[i] = spec.getFlux(E[i]);
      fluxErr[i] = spec.getFluxErr(E[i]);
    }
    TGraph *grate = new TGraphErrors(nn,E,rate,0,rateErr);
    grate->SetMarkerStyle(20);
    grate->SetMarkerSize(0.7);
    grate->SetName(Form("gmyrateTot_k%d",ik));
    grate->Write();
    
    TGraph *gflux = new TGraphErrors(nn,E,flux,0,fluxErr);
    gflux->SetMarkerStyle(20);
    gflux->SetMarkerSize(0.7);
    gflux->SetName(Form("gmyfluxTot_k%d",ik));
    gflux->Write();
    delete grate;
    delete gflux;
  }
}
*/
