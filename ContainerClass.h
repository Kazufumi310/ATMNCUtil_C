int const MAX_ARRAY_SIZE = 32768;

struct ATMNCStatInfo{
  ATMNCStatInfo();
  void branch(TTree* tree);
  void setBranchAddress(TTree* tree);
  double crsamp;
  double crsim;
  double fintsum;
  double fintEach[4]; // element size is dummy. acturally it is [n1ry]
  double sarea;
  double omega;
  double tequiv;
  double fact1ry;
};


struct ATMNCInputInfo{
  ATMNCInputInfo();
  void init();
  void dump();
  void read(std::string inputFile);
  void branch(TTree *tree);
  void setBranchAddress(TTree *tree);

  // for the simulation
  //  static int const fnameLength = 100;
  //  char pref[fnameLength];
  std::string* m_inputFile;
  std::string* pref;
  int n1ry;
  int kind1ry[16];          // [n1ry] kind1ry
  //  char f1ryFile[16][fnameLength]; // [n1ry] data table file
  std::string* f1ryFile[16];
  //  char cfil[fnameLength]; // cutoff table
  //  char gfil[fnameLength]; // geomagnetic field table (IGRF)
  //  char afil[fnameLength]; //  Air density table
  std::string* cfil; // cutoff table
  std::string* gfil; // geomagnetic field table (IGRF)
  std::string* afil; //  Air density table
  int mnth;   // start_month*100 + month_end
  double dmyr; // Year from IGRF date
  double epmi; // Min energy for primary CRs
  double lat0; // Center Lattitide of CR injection
  double lon0; // Center Longitude of CR injection
  double rsmp; // Sampling radius in degree
  double czmi; // Sampling minimum cos(zenith)
  double hinj; // CR injection Height
  double ehmi; // Minimum hadron int. energy
  double eemi; // Minimum elemag int. energy
  double eprp; // Minimum Energy for all particle propagation
  double tlim; // Maximum time for p,n,e,g propagation
  double rbnd; // Maximum distance from earth center
  double epht; // below this use jam int model
  bool noem;   // swith off elemag int.
  bool auti;  // Auto generation of random seed
  int nwup;  // Call rndc before simulation
  bool dbg1;  // Record random seed for each event
  bool dbg2;  // Record random seed for each particle
  int seed[2];  // Seeds in case for auti
  //char vfil[fnameLength]; //
  //  char vpfl[fnameLength]; //
  std::string* vfil;
  std::string* vpfl;
  int neve; // CR// to sample
  double fcrs; // Factor for cross section
  double rpft; // Production ratio for Jam
  double jobt; // Maximum Jobtime
  
  // for the observation
  //  char rfil[fnameLength];
  //  char ofil[fnameLength];
  std::string* rfil;
  std::string* ofil;
  bool is1d;
  int nsite;
  double siteParam[16][4]; //[nsite][nparam(=4)]
  double comi;  // Minimum cos(zenith) to obs
  int nczo;  // Observation cos(zenith) devision #
  int nazo;  // Observation azimuth devision #
  int nkind;
  int obsp[16]; //[nkind]
  double eomi; // Minimum energy to obs
  double etmi; // Minimum energy in table
  bool alpr; // Print when observed with 1ry
  int nstp[2]; // print every [0] events and stop at [1] prints
};

class ATMNCFileHandler{
 public:
  ATMNCFileHandler();
  void ascii2root(std::string inputAscii, std::string outputRoot);
  void sum(std::string Root1, std::string Root2, std::string newRoot="");
  
};



ATMNCStatInfo::ATMNCStatInfo(){
  crsamp = crsim = 0;
  sarea =  omega = fintsum = -1;
  tequiv = 0;
  fact1ry = -1;
  for(int i=0;i<4;i++) fintEach[i] = 0;
}

void ATMNCStatInfo::branch(TTree* tree){
  tree->Branch("crsim",&crsim,"crsim/D");
  tree->Branch("crsamp",&crsamp,"crsamp/D");
  tree->Branch("fintsum",&fintsum,"fintsum/D");
  tree->Branch("fintEach",fintEach,"fintEach[4]/D");
  tree->Branch("sarea",&sarea,"sarea/D");
  tree->Branch("omega",&omega,"omega/D");
  tree->Branch("tequiv",&tequiv,"tequiv/D");
  tree->Branch("fact1ry",&fact1ry,"fact1ry/D");
}

void ATMNCStatInfo::setBranchAddress(TTree* tree){
  tree->SetBranchAddress("crsim",&crsim);
  tree->SetBranchAddress("crsamp",&crsamp);
  tree->SetBranchAddress("fintsum",&fintsum);
  tree->SetBranchAddress("fintEach",fintEach);
  tree->SetBranchAddress("sarea",&sarea);
  tree->SetBranchAddress("omega",&omega);
  tree->SetBranchAddress("tequiv",&tequiv);
  tree->SetBranchAddress("fact1ry",&fact1ry);
}

ATMNCInputInfo::ATMNCInputInfo(){
  init();
}

void ATMNCInputInfo::dump(){
  std::cout<<"dump"<<std::endl;
  std::cout<<"pref:"<<*pref<<std::endl;
  std::cout<<"n1ry:"<<n1ry<<std::endl;
  for(int i=0;i<n1ry;i++){
    std::cout<<kind1ry[i]<<" "<<*f1ryFile[i]<<std::endl;    
  }
  std::cout<<"cfil:"<<*cfil<<std::endl;
  std::cout<<"gfil:"<<*gfil<<std::endl;
  std::cout<<"afil:"<<*afil<<std::endl;
  std::cout<<"mnth:"<<mnth<<std::endl; 
  std::cout<<"dmyr:"<<dmyr<<std::endl; 
  std::cout<<"epmi:"<<epmi<<std::endl; 
  std::cout<<"lat0:"<<lat0<<std::endl; 
  std::cout<<"lon0:"<<lon0<<std::endl; 
  std::cout<<"rsmp:"<<rsmp<<std::endl; 
  std::cout<<"czmi:"<<czmi<<std::endl; 
  std::cout<<"hinf:"<<hinj<<std::endl; 
  std::cout<<"ehmi:"<<ehmi<<std::endl; 
  std::cout<<"eemi:"<<eemi<<std::endl; 
  std::cout<<"eprp:"<<eprp<<std::endl; 
  std::cout<<"tlim:"<<tlim<<std::endl; 
  std::cout<<"rbnd:"<<rbnd<<std::endl; 
  std::cout<<"epht:"<<epht<<std::endl; 
  std::cout<<"noem:"<<noem<<std::endl;  
  std::cout<<"auti:"<<auti<<std::endl;  
  std::cout<<"nwup:"<<nwup<<std::endl;  
  std::cout<<"dbg1:"<<dbg1<<std::endl;  
  std::cout<<"dbg2:"<<dbg2<<std::endl;  
  std::cout<<"seed:"<<seed[0]<<" "<<seed[1]<<std::endl; 
  std::cout<<"vfil:"<<*vfil<<std::endl; 
  std::cout<<"vpfl:"<<*vpfl<<std::endl; 
  std::cout<<"neve:"<<neve<<std::endl; 
  std::cout<<"fcrs:"<<fcrs<<std::endl; 
  std::cout<<"rpft:"<<rpft<<std::endl; 
  std::cout<<"jobt:"<<jobt<<std::endl; 
 
  // for the observation
  std::cout<<"rfil:"<<*rfil<<std::endl;
  std::cout<<"ofil:"<<*ofil<<std::endl;
  std::cout<<"nsite:"<<nsite<<std::endl;
  std::cout<<"siteParam:"<<std::endl;
  for(int i=0;i<nsite;i++){
    std::cout<<i<<" "<<siteParam[i][0]<<" "<<siteParam[i][1]<<" "<<siteParam[i][2]<<" "<<siteParam[i][3]<<std::endl;
  }
  std::cout<<"comi:"<<comi<<std::endl; 
  std::cout<<"nczo:"<<nczo<<std::endl; 
  std::cout<<"nazo:"<<nazo<<std::endl; 
  std::cout<<"nkind:"<<nkind<<std::endl;
  std::cout<<"obsp:";
  for(int i=0;i<nkind;i++){std::cout<<obsp[i]<<" ";}
  std::cout<<std::endl;
  std::cout<<"eomi:"<<eomi<<std::endl; 
  std::cout<<"etmi:"<<etmi<<std::endl; 
  std::cout<<"alpr:"<<alpr<<std::endl; 
  std::cout<<"nstp:"<<nstp[0]<<" "<<nstp[1]<<std::endl;
}

void ATMNCInputInfo::init(){
  //  sprintf(pref,"Null");
  //  sprintf(vfil,"Null");
  //  sprintf(vpfl,"Null");
  //  sprintf(cfil,"Null");
  //  sprintf(rfil,"Null");
  //  sprintf(gfil,"irgf05.d");
  //  sprintf(afil,"airdensity-std.d");
  m_inputFile = new std::string;
  pref = new std::string;
  vfil = new std::string;
  vpfl = new std::string;
  cfil = new std::string;
  rfil = new std::string;
  ofil = new std::string;
  gfil = new std::string;
  afil = new std::string;
  for(int i=0;i<16;i++) f1ryFile[i] = new std::string;  
  *m_inputFile = "dummy";
  *pref = "Null";
  *vfil = "Null";
  *vpfl = "Null";
  *cfil = "Null";
  *rfil = "Null";
  *gfil = "irgf05.d";
  *afil = "airdensity-std.d";
  dmyr   =    0.0;
  epmi =    1.0;
  n1ry  =    0;
  nkind = 0;
  nsite = 0;
  is1d=false;

  lat0 =   36.4258333;
  lon0 =  137.3102777;
  rsmp =   1.0;
  czmi =   0.0;

  epht = 31.6227766;
  rpft = 0.95; // key = rpft
  neve = 20000000.;
  hinj = 1.e5;
  ehmi = 0.2;
  eemi = 0.1;
  eprp = 0.001;
  tlim = 100.0;
  rbnd = 6378.140e4;
  auti =false;
  noem = true;
  dbg1 = false;
  dbg2 = false;
  seed[0] = seed[1] = -1;
  nwup = 0;
  mnth = 0;
  fcrs = 1.0;

  eomi = 1e-3;
  etmi = 1e-1;
  nczo = 1;
  nazo = 1;
  comi = -1;
  nstp[0] = 200000;
  nstp[1] = 10;
  alpr = false;
}

void ATMNCInputInfo::read(std::string inputFile){
  std::string const parName[ParticleCode::nParType]={
    "", "E+", "E-","GAMM","K+",
    "K-","K0L","K0S","N+","N-",
    "P+","P-","PI+","PI-","PI0",
    "MU+","MU-","NUE","NUEB","NUM",
    "NUMB","ALPH"};


  *m_inputFile = inputFile;
  std::ifstream ifs(inputFile.c_str());
  std::string line;
  while(std::getline(ifs,line)){
    // replace "," to " "  
    while(1){
      std::string::size_type pos=line.find(",");
      if(pos==std::string::npos){
	break;
      }
      line.erase(pos,1);
      line.insert(pos," ");
    }

    std::stringstream ss(line);
    std::string first;
    ss>>first;
    std::string arg=first.substr(0,4);

    if(arg=="pref" ){
      ss>>*pref;
    }
    else if(arg=="cfil" ){
      ss>>*cfil;
    }
    else if(arg =="gfil" ){
      ss>>*gfil;
    }
    else if(arg=="afil" ){
      ss>>*afil;
    }
    else if(arg=="vfil" ){ 
      ss>>*vfil;
    }
    else if(arg=="vpfl" ){
      ss>>*vpfl;
    }
    else if(arg=="rfil" ){
      ss>>*rfil;
    }
    else if(arg=="ofil"){
      ss>>*ofil;
    }
    /*  
    if(arg=="pref" || arg=="cfil" || arg =="gfil" || arg=="afil" || arg=="vfil" || arg=="vpfl" || arg=="rfil" || arg=="ofil"){
      std::string fname;
      ss >> fname;
      
      if(arg=="pref" ){
	if(fname.length()<fnameLength-1) sprintf(pref,fname.c_str());
	else sprintf(pref,fname.substr(0,fnameLength-1).c_str());
      } 
      if(arg=="cfil" ){
	if(fname.length()<fnameLength-1) sprintf(cfil,fname.c_str());
	else sprintf(cfil,fname.substr(0,fnameLength-1).c_str());
      }
      if(arg =="gfil" ){
	if(fname.length()<fnameLength-1) sprintf(gfil,fname.c_str());
	else sprintf(gfil,fname.substr(0,fnameLength-1).c_str());
      }
      if(arg=="afil" ){
	if(fname.length()<fnameLength-1) sprintf(afil,fname.c_str());
	else sprintf(afil,fname.substr(0,fnameLength-1).c_str());
      }
      if(arg=="vfil" ){ 
	if(fname.length()<fnameLength-1) sprintf(vfil,fname.c_str());
	else sprintf(vfil,fname.substr(0,fnameLength-1).c_str());
      }
      if(arg=="vpfl" ){
	if(fname.length()<fnameLength-1) sprintf(vpfl,fname.c_str());
	else sprintf(vpfl,fname.substr(0,fnameLength-1).c_str());
      }
      if(arg=="rfil" ){
	if(fname.length()<fnameLength-1) sprintf(rfil,fname.c_str());
	else sprintf(rfil,fname.substr(0,fnameLength-1).c_str());
      }
      if(arg=="ofil"){
	if(fname.length()<fnameLength-1) sprintf(ofil,fname.c_str());
	else sprintf(ofil,fname.substr(0,fnameLength-1).c_str());
      }
    }
    */
    else if(arg=="f1ry"){
      std::string kind,fname;
      ss >> kind >> fname;
      //      std::cout<<"fname:"<<fname<<" "<<fname.length()<<std::endl;
      //      sprintf(f1ryFile[n1ry],fname.c_str());   
      //      kind1ry[n1ry] = 0;
      *f1ryFile[n1ry] = fname;
      for(int i=0;i<nParType;i++){
	if(kind==parName[i]) kind1ry[n1ry] = i;
      }
      n1ry++;
    }
    else if(arg=="mnth") ss >> mnth; 
    else if(arg=="dmyr") ss >> dmyr; 
    else if(arg=="epmi") ss >> epmi; 
    else if(arg=="lat0") ss >> lat0; 
    else if(arg=="lon0") ss >> lon0; 
    else if(arg=="rsmp") ss >> rsmp; 
    else if(arg=="czmi") ss >> czmi; 
    else if(arg=="hinj") ss >> hinj; 
    else if(arg=="ehmi") ss >> ehmi; 
    else if(arg=="eemi") ss >> eemi; 
    else if(arg=="eprp") ss >> eprp; 
    else if(arg=="tlim") ss >> tlim; 
    else if(arg=="rbnd") ss >> rbnd; 
    else if(arg=="epht" || arg=="ejam") ss >> epht; 
    else if(arg=="noem" || arg=="auti" || arg=="dbg1" || arg=="dbg2" || arg=="alpr" || arg=="is1d"){
      std::string tf;
      ss >> tf;       
      bool isTrue;
      if(tf=="T" || tf=="t" || tf=="Yes" || tf=="yes" || tf=="Y" || tf=="y")
	isTrue=true;
      else isTrue = false;
      if(arg=="noem") noem = isTrue;
      if(arg=="auti") auti = isTrue;
      if(arg=="dbg1") dbg1 = isTrue;
      if(arg=="dbg2") dbg2 = isTrue;
      if(arg=="alpr") alpr = isTrue;
      if(arg=="is1d") is1d = isTrue;
    }
    else if(arg=="nwup") ss >> nwup; 
    else if(arg=="seed"){
      ss>>seed[0]>>seed[1];
    }
    else if(arg=="neve") ss >> neve;
    else if(arg=="fcrs") ss >> fcrs;
    else if(arg=="rpft") ss >> rpft;
    else if(arg=="jobt") ss >> jobt;
  
  // for the observation
    //    else if(arg=="rfil") ss >> rfil;
    //    else if(arg=="ofil") ss >> ofil;
    else if(arg=="site"){
      //      std::cout<<line<<std::endl;
      for(int i=0;i<4;i++){
	ss>>siteParam[nsite][i];
	//	std::cout<<" "<<siteParam[nsite][i]<<std::endl;
      }
      nsite++;
    }
    else if(arg=="comi") ss >> comi; 
    else if(arg=="nczo") ss >> nczo; 
    else if(arg=="nazo") ss >> nazo; 
    else if(arg=="obsp"){
      std::string kind;
      while(ss>>kind){
	if(kind.at(0)=='#'){
	  //	  std::cout<<kind<<std::endl;
	  break;
	}
	obsp[nkind] = 0;
	for(int i=0;i<nParType;i++){
	  if(kind==parName[i]) obsp[nkind] = i;
	}
	nkind++;
      }
    }
    else if(arg=="eomi") ss >> eomi;       
    else if(arg=="etmi") ss >> etmi;
    else if(arg=="nstp"){
      ss >> nstp[0]>>nstp[1];       
    }
  }
}
  
void ATMNCInputInfo::branch(TTree *tree){
  tree->Branch("inputFile",&m_inputFile);
  tree->Branch("pref",&pref);
  tree->Branch("n1ry",&n1ry,"n1ry/I");
  tree->Branch("kind1ry",&kind1ry,"kind1ry[n1ry]/I");          
  for(int i=0;i<n1ry;i++){
    tree->Branch(Form("f1ryFile%d",i),&f1ryFile[i]); 
  }
  tree->Branch("cfil",&cfil);
  tree->Branch("gfil",&gfil);
  tree->Branch("afil",&afil);
  tree->Branch("mnth",&mnth,"mnth/I");
  tree->Branch("dmyr",&dmyr,"dmyr/D");
  tree->Branch("epmi",&epmi,"epmi/D");
  tree->Branch("lat0",&lat0,"lat0/D");
  tree->Branch("lon0",&lon0,"lon0/D");
  tree->Branch("rsmp",&rsmp,"rsmp/D");
  tree->Branch("czmi",&czmi,"czmi/D");
  tree->Branch("hinj",&hinj,"hinj/D");
  tree->Branch("ehmi",&ehmi,"ehmi/D");
  tree->Branch("eemi",&eemi,"eemi/D");
  tree->Branch("eprp",&eprp,"eprp/D");
  tree->Branch("tlim",&tlim,"tlim/D");
  tree->Branch("rbnd",&rbnd,"rbnd/D");
  tree->Branch("epht",&epht,"epht/D");
  tree->Branch("noem",&noem,"noem/O");
  tree->Branch("auti",&auti,"auti/O");
  tree->Branch("nwup",&nwup,"nwup/I");
  tree->Branch("dbg1",&dbg1,"dbg1/O");
  tree->Branch("dbg2",&dbg2,"dbg2/O");
  tree->Branch("seed",seed,"seed[2]/I");
  //  leafName = Form("vfil[%d]/C",fnameLength);
  //  tree->Branch("vfil",vfil,leafName.c_str());
  tree->Branch("vfil",&vfil);
  //  leafName = Form("vpfl[%d]/C",fnameLength);
  //  tree->Branch("vpfl",vpfl,leafName.c_str());
  tree->Branch("vpfl",&vpfl);
  tree->Branch("neve",&neve,"neve/I");
  tree->Branch("fcrs",&fcrs,"fcrs/D");
  tree->Branch("rpft",&rpft,"rpft/D");
  tree->Branch("jobt",&jobt,"jobt/D");
  
  // for the observation
  //  leafName = Form("rfil[%d]/C",fnameLength);
  //  tree->Branch("rfil",rfil,leafName.c_str());
  tree->Branch("rfil",&rfil);
  //  leafName = Form("ofil[%d]/C",fnameLength);
  //  tree->Branch("ofil",ofil,leafName.c_str());
  tree->Branch("ofil",&ofil);
  tree->Branch("is1d",&is1d,"is1d/O");
  tree->Branch("nsite",&nsite,"nsite/I");
  tree->Branch("siteParam",siteParam,"siteParam[nsite][4]/D");
  tree->Branch("comi",&comi,"comi/D");
  tree->Branch("nczo",&nczo,"nczo/I");
  tree->Branch("nazo",&nazo,"nazo/I");
  tree->Branch("nkind",&nkind,"nkind/I");
  tree->Branch("obsp",&obsp,"obsp[nkind]/I");
  tree->Branch("eomi",&eomi,"eomi/D");
  tree->Branch("etmi",&etmi,"etmi/D");
  tree->Branch("alpr",&alpr,"alpr/O");
  tree->Branch("nstp",&nstp,"nstp[2]/I");
}

void ATMNCInputInfo::setBranchAddress(TTree *tree){
  tree->SetBranchAddress("inputFile",&m_inputFile);         
  tree->SetBranchAddress("pref",&pref);         
  tree->SetBranchAddress("n1ry",&n1ry);         
  tree->SetBranchAddress("kind1ry",&kind1ry);   
  for(int i=0;i<16;i++){
    std::string bname = Form("f1ryFile%d",i);
    if(tree->GetBranch(bname.c_str()))
      tree->SetBranchAddress(bname.c_str(),&f1ryFile[i]);
  }

  tree->SetBranchAddress("cfil",&cfil);
  tree->SetBranchAddress("gfil",&gfil);
  tree->SetBranchAddress("afil",&afil);
  tree->SetBranchAddress("mnth",&mnth); 
  tree->SetBranchAddress("dmyr",&dmyr); 
  tree->SetBranchAddress("epmi",&epmi); 
  tree->SetBranchAddress("lat0",&lat0); 
  tree->SetBranchAddress("lon0",&lon0); 
  tree->SetBranchAddress("rsmp",&rsmp); 
  tree->SetBranchAddress("czmi",&czmi); 
  tree->SetBranchAddress("hinj",&hinj); 
  tree->SetBranchAddress("ehmi",&ehmi); 
  tree->SetBranchAddress("eemi",&eemi); 
  tree->SetBranchAddress("eprp",&eprp); 
  tree->SetBranchAddress("tlim",&tlim); 
  tree->SetBranchAddress("rbnd",&rbnd); 
  tree->SetBranchAddress("epht",&epht); 
  tree->SetBranchAddress("noem",&noem); 
  tree->SetBranchAddress("auti",&auti); 
  tree->SetBranchAddress("nwup",&nwup); 
  tree->SetBranchAddress("dbg1",&dbg1); 
  tree->SetBranchAddress("dbg2",&dbg2); 
  tree->SetBranchAddress("seed",seed);  
  tree->SetBranchAddress("vfil",&vfil);
  tree->SetBranchAddress("vpfl",&vpfl);
  tree->SetBranchAddress("neve",&neve);
  tree->SetBranchAddress("fcrs",&fcrs);
  tree->SetBranchAddress("rpft",&rpft);
  tree->SetBranchAddress("jobt",&jobt);
  
  // for the observation
  tree->SetBranchAddress("rfil",&rfil);
  tree->SetBranchAddress("ofil",&ofil);
  tree->SetBranchAddress("nsite",&nsite);
  tree->SetBranchAddress("siteParam",siteParam);
  tree->SetBranchAddress("comi",&comi); 
  tree->SetBranchAddress("nczo",&nczo); 
  tree->SetBranchAddress("nazo",&nazo); 
  tree->SetBranchAddress("nkind",&nkind);     
  tree->SetBranchAddress("obsp",&obsp);       
  tree->SetBranchAddress("eomi",&eomi);       
  tree->SetBranchAddress("etmi",&etmi);       
  tree->SetBranchAddress("alpr",&alpr);       
  tree->SetBranchAddress("nstp",&nstp);       
}


