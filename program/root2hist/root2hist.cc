#include <myheader/myRoot.h>
#include <myheader/myStandard.h>
#include "atmnc-constants.h"
#include "atmnc-particle-code.h"
#include "ContainerClass.h"
#include "Histogrammer.h"


class ArgHandler{
public:
  ArgHandler();
  std::string ifname(){return m_ifname;}
  std::string wfname(){return m_wfname;}
  std::string ofname(){return m_ofname;}

  double getEMin(){return m_emin;}
  int getNCos(){return m_ncos;}
  double getCMin(){return m_cmin;}

  bool getArg(int &argc,char** argv);
  bool isUserDarea(){ return m_defineUserDarea;}
  bool needToCountTequiv(){ return m_needTequiv;}
  
private:
  bool m_defineUserDarea;
  double m_emin;
  int m_ncos;
  double m_cmin;
  std::string m_ifname;
  std::string m_wfname;
  std::string m_ofname;
  bool m_needTequiv;
};


class EventWeighter{
public:
  EventWeighter(std::string ifname);
  double getWeight(double ek1ry);
private:
  std::vector<double> m_threE;
  std::vector<double> m_weight;
};

////////////////////////
/// begin of main program
/////////////////////////
int main(int argc, char** argv){
  std::cout<<"////////////////////////////"<<std::endl;
  std::cout<<"program to make flux histgrams from root files with full-hit information."<<std::endl;
  std::cout<<"usage: bin/root2hist -i input_file -o output_file (-w weight_file)"<<std::endl;
  std::cout<<"////////////////////////////"<<std::endl;

  ArgHandler gArg;  
  gArg.getArg(argc,argv);

  std::string ifname = gArg.ifname();
  std::string wfname = gArg.wfname();
  std::string ofname = gArg.ofname();

  std::cout<<"ifname:"<<ifname<<std::endl;
  std::cout<<"ofname:"<<ofname<<std::endl;
  std::cout<<"wfname:"<<wfname<<std::endl;


  //prepare output file & trees
  TFile fout(ofname.c_str(),"recreate");
  std::cout<<"create "<<ofname<<std::endl;
  TTree *trout = new TTree("trhist","");
  TTree *troutInfo = new TTree("trinfo","");

  //prepare input trees
  TChain chain("tr");
  TChain chinfo("trinfo");
  chain.Add(ifname.c_str());
  chinfo.Add(ifname.c_str());

  int const nEntry = chain.GetEntries();
  if(nEntry==0){
    std::cout<<"can't find "<<ifname<<std::endl;
    return 0;
  }else{
    std::cout<<"read "<<ifname<<" # events:"<<nEntry<<std::endl;
  }

  ATMNCInputInfo input;
  input.setBranchAddress(&chinfo);
  ATMNCStatInfo statInfo;
  statInfo.setBranchAddress(&chinfo);
  
  HitInfo hit;
  hit.setBranchAddress(&chain);
  int ievt;
  chain.SetBranchAddress("ievt",&ievt);


  //prepare and initialize Histogrammer
  chinfo.GetEntry(0);
  if(gArg.getEMin()!=DBL_MAX){
    input.etmi = gArg.getEMin();
  }
  if(gArg.getNCos()>0){
    input.nczo = gArg.getNCos();
    input.comi = gArg.getCMin();
  }
  Histogrammer hister(input);
  double tequiv = (gArg.needToCountTequiv())?statInfo.tequiv:0;
  hister.setTequiv(tequiv);

  if(gArg.isUserDarea()){
    std::cout<<"!!!!!!!!!!!!!  WARNING  !!!!!!!!!!!!!!!!"<<std::endl;
    std::cout<<"special criteria defined by user will be used to judge the detector hit."<<std::endl;
    std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  }



  // record header information
  input.branch(troutInfo);
  troutInfo->Branch("tequiv",&tequiv,"tequiv/D");
  troutInfo->Fill();
  troutInfo->Write();

  

  // event loop to fill histograms with applying event weight
  EventWeighter weighter(wfname);
  for(int ientry=0;ientry<nEntry;ientry++){
    chain.GetEntry(ientry);
    if(hit.nhit==0) continue;
    int pid1ry = hit.pid[0];
    double p1ry = sqrt(pow(hit.bornP[0][0],2)+pow(hit.bornP[0][1],2)+pow(hit.bornP[0][2],2));

    double m1ry = ParticleCode::mass[pid1ry];
    double ek1ry =  sqrt(p1ry*p1ry+m1ry*m1ry) - m1ry;
    if(pid1ry==ParticleCode::kalpha) ek1ry /= 4;
    double totWeight = weighter.getWeight(ek1ry);
    
    hister.addEvent(hit,totWeight);
  }

  
  // record histograms into TTree
  ATMNCHistInfo output;
  output.branch(trout);
  
  bool isFirstTime = true;
  int nebin;
  double ebins[4096];
  
  for(output.isite=0;output.isite<input.nsite;output.isite++){
    for(output.icosz=0;output.icosz<input.nczo;output.icosz++){
      for(output.iazim=0;output.iazim<input.nazo;output.iazim++){
	output.nkind = input.nkind;
	for(int ikind=0;ikind<input.nkind;ikind++){
	  
	  hister.getHistogram(output.kind[ikind], output.isite, 
			      output.icosz, output.iazim,
			      nebin, ebins, output.nobs[ikind], output.nobsErr[ikind]);
	  if(isFirstTime){
	    isFirstTime = false;
	    if(nebin!=output.nbin){
	      std::cout<<"Unkonwn Error. #bin of histgram = "<<output.nbin<<std::endl;
	      return 0;
	    }
	    for(int i=0;i<nebin;i++){
	      output.ebin[i] = ebins[i];
	    }
	  }else{
	    for(int i=0;i<nebin;i++){
	      if(output.ebin[i] != ebins[i]){
		std::cout<<"Error. #bin of histgram is not consisten."<<std::endl;
		return 0;
	      }
	    }
	  }
	  SpectrumHandler spec;
	  hister.getSpectrum(output.kind[ikind], output.isite, 
			     output.icosz, output.iazim,spec);
	  for(int i=0;i<output.nbin;i++){
	    output.flux[ikind][i] = spec.getFlux(output.ebin[i]);
	  }
	  trout->Fill();
	}
      }
    }
  }
  trout->Write();
  fout.Close();
  std::cout<<"finish"<<std::endl;
}
////////////////////////
/// end of main program
/////////////////////////




ArgHandler::ArgHandler(){
  m_defineUserDarea = false;
  m_emin = DBL_MAX;
  m_ncos = -1;
  m_cmin = -1;
  m_needTequiv = true;
}

bool ArgHandler::getArg(int &argc, char** argv){
  std::cout<<"command: ";
  for(int i=0;i<argc;i++){
    std::cout<<argv[i]<<" ";
  }
  std::cout<<std::endl;
  
  for(int iarg=1;iarg<argc;iarg++){
    std::string arg = argv[iarg];
    if(arg[0]=='-'){ //option
      if(arg.size()!=2){
	std::cout<<"argument error! unknown option "<<arg<<": option should be 1 charactor. "<<std::endl;
	return false;
      }
      char opt = arg[1];
      /*      if(opt=='d'){ // use user's special shape detector
	char* optarg = argv[++iarg];
	/////// should edit user_hitTestFunc.icc ///////////
	m_defineUserDarea = atoi(optarg);
	}else */ if(opt=='e'){ 
	char* optarg = argv[++iarg];
	m_emin = atof(optarg);
      }else if(opt=='t'){  // cost range
	char* optarg = argv[++iarg];
	m_ncos = atoi(optarg);
	optarg = argv[++iarg];
	m_cmin = atof(optarg);
	std::cout<<"cost range specified : ncosz = "<<m_ncos<<" cmin = "<<m_cmin<<std::endl;
      }else if(opt=='t'){ 
	char* optarg = argv[++iarg];
	m_needTequiv = atoi(optarg);
      }else if(opt=='i'){  // cost range
	char* optarg = argv[++iarg];
	m_ifname = optarg;
      }else if(opt=='o'){  // cost range
	char* optarg = argv[++iarg];
	m_ofname = optarg;
      }else if(opt=='w'){  // cost range
	char* optarg = argv[++iarg];
	m_wfname = optarg;
      }else{
	std::cout<<"argument error! unknown option:"<<opt<<std::endl;
	return false;
      }
    }
    if(m_ifname.empty() || m_ofname.empty()){
      std::cout<<"Arg error. you need to assign input file and output file"<<std::endl;
      std::cout<<"usage: bin/root2hist -i input_file -o output_file (-w weight_file)"<<std::endl;
    }
  }
  return true;
}




EventWeighter::EventWeighter(std::string ifname){
  m_weight.clear();
  m_threE.clear();
  std::ifstream ifs(ifname.c_str());
  if(!ifs){ 
    std::cout<<"EventWeighter : setting weight data: can't find "<<ifname<<std::endl;
    exit(1);
  }else{
    std::cout<<"read "<<ifname<<std::endl;
  }
  
  std::string line;
  while(std::getline(ifs,line)){
    std::stringstream ss(line);
    double epmi,weight;
    ss>>epmi>>weight;
    m_threE.push_back(epmi);
    m_weight.push_back(weight);
  }
  
  int const nBnd = m_threE.size();
  for(int ibnd=0;ibnd<nBnd;ibnd++){
    std::cout<<"event weight="<<m_weight[ibnd]<<" for Ekine > "<<m_threE[ibnd]<<std::endl;
  }
}

double EventWeighter::getWeight(double ekin){
  int const nBnd = m_threE.size();
  if(nBnd==0) return 1.;
  

  for(int ibnd = nBnd-1;ibnd>=0;ibnd--){
    if(m_threE[ibnd]<ekin) return m_weight[ibnd];
  }
  
  static bool alreadyWarn = false;
  if(!alreadyWarn){
    alreadyWarn = true;
    std::cout<<"Warning: EventWeight::getWeight: no data for primary particle with ek ="<<ekin<<std::endl;
  }
  return 1.;
};
