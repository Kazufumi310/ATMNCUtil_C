#include <ContainerClass.h>

ATMNCHistInfo::ATMNCHistInfo(){
  nkind = 0;
}

void ATMNCHistInfo::branch(TTree* tree){
  tree->Branch("isite",&isite,"isite/I");
  tree->Branch("icosz",&icosz,"icosz/I");
  tree->Branch("iazim",&iazim,"iazim/I");
  tree->Branch("Darea",&Darea,"Darea/D");
  tree->Branch("Domega",&Domega,"Domega/D");
  tree->Branch("coszFrom",&coszFrom,"coszFrom/D");
  tree->Branch("coszTo",&coszTo,"coszTo/D");
  tree->Branch("azimFrom",&azimFrom,"azimFrom/D");
  tree->Branch("azimTo",&azimTo,"azimTo/D");

  tree->Branch("nkind",&nkind,"nkind/I");
  tree->Branch("kind",kind,"kind[nkind]/I");
  tree->Branch("ebin",ebin,Form("ebin[%d]/D",nbin));
  tree->Branch("nobs",nobs,Form("nobs[nkind][%d]/D",nbin));
  tree->Branch("nobsErr",nobsErr,Form("nobsErr[nkind][%d]/D",nbin));
  tree->Branch("flux",flux,Form("flux[nkind][%d]/D",nbin));  
}

void ATMNCHistInfo::setBranchAddress(TTree* tree){
  tree->SetBranchAddress("isite",&isite);
  tree->SetBranchAddress("icosz",&icosz);
  tree->SetBranchAddress("iazim",&iazim);
  tree->SetBranchAddress("Darea",&Darea);
  tree->SetBranchAddress("Domega",&Domega);
  tree->SetBranchAddress("coszFrom",&coszFrom);
  tree->SetBranchAddress("coszTo",&coszTo);
  tree->SetBranchAddress("azimFrom",&azimFrom);
  tree->SetBranchAddress("azimTo",&azimTo);

  tree->SetBranchAddress("nkind",&nkind);
  tree->SetBranchAddress("kind",kind);
  tree->SetBranchAddress("ebin",ebin);
  tree->SetBranchAddress("nobs",nobs);
  tree->SetBranchAddress("nobsErr",nobsErr);
  tree->SetBranchAddress("flux",flux);
}


ATMNCStatInfo::ATMNCStatInfo(){
  crsamp = crsim = 0;
  sarea =  omega = fintsum = -1;
  tequiv = 0;
  fact1ry = -1;
  for(int i=0;i<2;i++){
    for(int k=0;k<mtbl;k++){
      ebin1ry[i][k] = 0;
    }
    for(int j=0;j<4;j++){
      for(int k=0;k<mtbl;k++){
	nobs1ry[i][j][k] = dflux1ry[i][j][k] = 0;
      }
    }
  }
}


void ATMNCStatInfo::branch(TTree* tree){
  tree->Branch("crsim",&crsim,"crsim/D");
  tree->Branch("crsamp",&crsamp,"crsamp/D");
  tree->Branch("fintsum",&fintsum,"fintsum/D");
  tree->Branch("sarea",&sarea,"sarea/D");
  tree->Branch("omega",&omega,"omega/D");
  tree->Branch("tequiv",&tequiv,"tequiv/D");
  tree->Branch("fact1ry",&fact1ry,"fact1ry/D");


  tree->Branch("ebin1ry",ebin1ry,Form("ebin1ry[2][%d]/D",mtbl));
  tree->Branch("nobs1ry",nobs1ry,Form("nobs1ry[2][4][%d]/D",mtbl));
  tree->Branch("dflux1ry",dflux1ry,Form("dflux1ry[2][4][%d]/D",mtbl));
}

void ATMNCStatInfo::setBranchAddress(TTree* tree){
  tree->SetBranchAddress("crsim",&crsim);
  tree->SetBranchAddress("crsamp",&crsamp);
  tree->SetBranchAddress("fintsum",&fintsum);
  tree->SetBranchAddress("sarea",&sarea);
  tree->SetBranchAddress("omega",&omega);
  tree->SetBranchAddress("tequiv",&tequiv);
  tree->SetBranchAddress("fact1ry",&fact1ry);

  tree->SetBranchAddress("ebin1ry",ebin1ry);
  tree->SetBranchAddress("nobs1ry",nobs1ry);
  tree->SetBranchAddress("dflux1ry",dflux1ry);
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



