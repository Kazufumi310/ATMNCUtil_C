#include <FileHandler.h>

FileHandler::FileHandler(){
  ;
};

void FileHandler::sum(std::vector<std::string> inputFileList, std::string outputFile){
  
  TFile fout(outputFile.c_str(),"recreate");
  
  TChain chhist("tr");
  TChain chinfo("trinfo");
  int nFile = 0;
  for(auto input : inputFileList){
    chinfo.Add(input.c_str());
    chhist.Add(input.c_str());
    std::cout<<"read "<<input<<std::endl;
    nFile++;
    if(chinfo.GetEntries()!=nFile){
      std::cout<<"warning: "<<input<<" seems empty "<<std::endl;
    }
  }

  ATMNCInputInfo inputInfo;
  inputInfo.setBranchAddress(&chinfo);
  ATMNCStatInfo statInfo;
  statInfo.setBranchAddress(&chinfo);
  ATMNCHistInfo histInfo;
  histInfo.setBranchAddress(&chhist);

  chinfo.GetEntry(0);
  Histogrammer hister(inputInfo);

  //////////// 1ry particle
  double tequivTot = 0;
  ATMNCStatInfo statOut;
  for(int ientry=0;ientry<chinfo.GetEntries();ientry++){
    chinfo.GetEntry(ientry);
    tequivTot+=statInfo.tequiv;
    for(int i=0;i<ATMNCStatInfo::mtbl;i++){
      for(int iab=0;iab<2;iab++){
	for(int i1ry=0;i1ry<inputInfo.n1ry;i1ry++){
	  statOut.nobs1ry[iab][i1ry][i]+=statInfo.nobs1ry[iab][i1ry][i];
	  statOut.dflux1ry[iab][i1ry][i]+=statInfo.dflux1ry[iab][i1ry][i];
	}
      }
    }
  }
  std::cout<<"tequivTot:"<<tequivTot<<std::endl;
  {  //
    SpectrumHandler spec;
    for(int iab=0;iab<2;iab++){
      for(int i1ry=0;i1ry<inputInfo.n1ry;i1ry++){
	double rate[ATMNCStatInfo::mtbl];
	statInfo.fact1ry = 1/(tequivTot*statInfo.sarea*statInfo.omega);
	for(int i=0;i<ATMNCStatInfo::mtbl;i++){
	  rate[i] = statOut.nobs1ry[iab][i1ry][i]*statInfo.fact1ry;
	}
	spec.setRate(ATMNCStatInfo::mtbl,statInfo.ebin1ry[iab],rate);
	double *flux = spec.getFluxArray();

	for(int i=0;i<ATMNCStatInfo::mtbl;i++){
	  statInfo.nobs1ry[iab][i1ry][i] = statOut.nobs1ry[iab][i1ry][i];
	  statInfo.dflux1ry[iab][i1ry][i] = flux[i];
	}
      }
    }
  }
  
  hister.setTequiv(tequivTot);
  std::cout<<"normalization factor:"<<hister.getNormFactor(0,0)<<std::endl;

  int const nEntry = chhist.GetEntries();
  std::cout<<"nEntry:"<<nEntry<<std::endl;
  for(int ientry=0;ientry<nEntry;ientry++){
    if(ientry % int(nEntry/100.) == 0) std::cout<<ientry / int(nEntry/100.)<<"%"<<std::endl;
    chhist.GetEntry(ientry);
    for(int ik=0;ik<histInfo.nkind;ik++){
      hister.addHistogram(ik, histInfo.isite, histInfo.icosz, histInfo.iazim,
			  histInfo.nbin, histInfo.ebin, histInfo.nobs[ik], histInfo.nobsErr[ik]);
    }
  }

  TTree* trinfo = new TTree("trinfo","");
  TTree* trout = new TTree("tr","");
  statInfo.tequiv = tequivTot;
  inputInfo.branch(trinfo);
  statInfo.branch(trinfo);
  trinfo->Fill();
  
  histInfo.branch(trout);

  int const nazim = hister.getNAzim();
  int const ncosz = hister.getNCosz();
  int const nsite = hister.getNSite();
  int const nkind = hister.getNKind();
  histInfo.nkind = nkind;
  for(int is=0;is<nsite;is++){
    histInfo.isite = is;
    for(int ic=0;ic<ncosz;ic++){
      histInfo.icosz =ic;
      histInfo.coszFrom = hister.getCosz(ic);
      histInfo.coszTo = hister.getCosz(ic+1);
      for(int ia=0;ia<nazim;ia++){
	histInfo.iazim = ia;
	histInfo.azimFrom = hister.getAzim(ia);
	histInfo.azimTo = hister.getAzim(ia+1);
	
	for(int ik=0;ik<nkind;ik++){
	  histInfo.kind[ik] = hister.getPID(ik);
	  int num;
	  hister.getHistogram(ik, is, ic, ia,
			      num, histInfo.ebin, histInfo.nobs[ik], histInfo.nobsErr[ik]);
	  if(num!=histInfo.nbin){
	    std::cout<<"Error : FileHandler::sum: unexpected bin number:"<<num<<std::endl;
	    exit(1);
	  }
	  SpectrumHandler spec;
	  hister.getSpectrum(ik, is, ic, ia, spec);
	  double *flux = spec.getFluxArray();
	  for(int i=0;i<num;i++){
	    histInfo.flux[ik][i] = flux[i];
	  }
	}
	trout->Fill();
      }
    }
  }
  trinfo->Write();
  trout->Write();
  fout.Close();
}


void FileHandler::ascii2root(std::string inputAscii, std::string outputRoot){
  std::ifstream ifs(inputAscii);
  if(!ifs){
    std::cout<<"Error: can't open "<<inputAscii<<std::endl;
    exit(1);
  }else{
    std::cout<<"read "<<inputAscii<<std::endl;
  }
  TFile fout(outputRoot.c_str(),"recreate");
  TTree* trinfo = new TTree("trinfo","simulation setting");
  TTree* tr = new TTree("tr","observed spectrum");

  ATMNCInputInfo inputInfo;
  inputInfo.branch(trinfo);
  ATMNCStatInfo statInfo;
  statInfo.branch(trinfo);
  ATMNCHistInfo histInfo;
  histInfo.branch(tr);

  readInputInfo(inputInfo,ifs);
  read1ryInfo(statInfo,ifs);
  trinfo->Fill();
  while( read1Histogram(histInfo,ifs) ){
    tr->Fill();
  }
  
  trinfo->Write();
  tr->Write();
  fout.Close();
}

bool FileHandler::read1Histogram( ATMNCHistInfo& memo, std::ifstream &ifs){
  std::string line;
  { // header 1
    std::getline(ifs,line); //
    if(ifs.eof()) return false;
    std::stringstream ss(line);
    ss>>memo.isite;
    ss>>memo.icosz;
    ss>>memo.iazim;
    memo.isite -= 1; // Fortran -> C style
    memo.icosz -= 1; // Fortran -> C style
    memo.iazim -= 1; // Fortran -> C style
    ss>>memo.Darea;
    ss>>memo.Domega;
  }
  { // header 2
    std::getline(ifs,line); //
    if(ifs.eof()){
      std::cout<<"warning: in FileHandler::read1Histogram: Unexpected EOF"<<std::endl;
      return false;
    }
    std::stringstream ss(line);
    std::string trash;
    for(int i=0;i<6;i++) ss>>trash;
    ss>>trash>>memo.coszFrom>>trash>>memo.coszTo;
    ss>>trash>>memo.azimFrom>>trash>>memo.azimTo;
  }
  { // histogram title
    std::getline(ifs,line); //
    std::stringstream ss(line);
    std::string trash;
    ss>>trash>>trash>>trash; // 1st column title
    std::string kind;
    memo.nkind = 0;
    while(ss>>kind>>trash){
      if(!ss) break;
      for(int i=0;i<ParticleCode::nParType;i++){
	if(kind==ParticleCode::parName[i]) memo.kind[memo.nkind] = i;
	}
      memo.nkind++;
    }
  }
  
  for(int i=0;i<memo.nbin;i++){
    std::getline(ifs,line);
    std::stringstream ss(line);
    ss>>memo.ebin[i];
    for(int ik=0;ik<memo.nkind;ik++){
      ss>>memo.nobs[ik][i] >> memo.flux[ik][i];
      memo.nobsErr[ik][i] = sqrt(memo.nobs[ik][i]);
      if(!ss){
	std::cout<<"Error in FileHandler : reading line :"<<line<<std::endl;
      }
    }
  }
  return true;
}


void FileHandler::read1ryInfo( ATMNCStatInfo& memo, std::ifstream &ifs){
  std::string line;
  {
    std::getline(ifs,line); //
    std::stringstream ss(line);
    ss>>memo.crsim;
    ss>>memo.crsamp;
    ss>>memo.fintsum;
    ss>>memo.sarea;
    ss>>memo.omega;
    ss>>memo.tequiv;
    memo.fact1ry = 1./(memo.sarea * memo.omega * memo.tequiv);
  }

  for(int iaft=0;iaft<2;iaft++){ // before and after cut-off,
    std::getline(ifs,line); // header of 1ry histogram
    int n1ry=0;
    line = line.substr(line.find("cutoff)")+7);
    {
      std::stringstream ss(line);
      while(1){
	std::string trash;
	ss>>trash;
	if(!ss) break;
	n1ry++;
      }
    }
    for(int i=0;i<memo.mtbl;i++){
      std::getline(ifs,line);
      std::stringstream ss(line);
      ss>>memo.ebin1ry[iaft][i];
      for(int ik=0;ik<n1ry;ik++){
	ss>>memo.nobs1ry[iaft][ik][i] >> memo.dflux1ry[iaft][ik][i];
      }
    }
  }
}

void FileHandler::readInputInfo( ATMNCInputInfo& memo, std::ifstream &ifs){
  
  std::string line;
  while(std::getline(ifs,line)){
    if(line.empty()) continue;
    if(line.substr(0,10)=="#---------"){ // end of header
      break; 
    }else if(line[0]=='#') continue; // comment
    
    std::stringstream ss(line);
    std::string first;
    ss>>first;
    std::string arg=first.substr(0,4);

    if(arg=="pref" ){
      ss>>*memo.pref;
    }
    else if(arg=="cfil" ){
      ss>>*memo.cfil;
    }
    else if(arg =="gfil" ){
      ss>>*memo.gfil;
    }
    else if(arg=="afil" ){
      ss>>*memo.afil;
    }
    else if(arg=="vfil" ){ 
      ss>>*memo.vfil;
    }
    else if(arg=="vpfl" ){
      ss>>*memo.vpfl;
    }
    else if(arg=="rfil" ){
      ss>>*memo.rfil;
    }
    else if(arg=="ofil"){
      ss>>*memo.ofil;
    }
    else if(arg=="f1ry"){
      std::string kind,fname;
      ss >> kind >> fname;
      //      std::cout<<"fname:"<<fname<<" "<<fname.length()<<std::endl;
      //      sprintf(f1ryFile[n1ry],fname.c_str());   
      //      kind1ry[n1ry] = 0;
      *memo.f1ryFile[memo.n1ry] = fname;
      for(int i=0;i<ParticleCode::nParType;i++){
	if(kind==ParticleCode::parName[i]) memo.kind1ry[memo.n1ry] = i;
      }
      memo.n1ry++;
    }
    else if(arg=="mnth") ss >> memo.mnth; 
    else if(arg=="dmyr") ss >> memo.dmyr; 
    else if(arg=="epmi") ss >> memo.epmi; 
    else if(arg=="lat0") ss >> memo.lat0; 
    else if(arg=="lon0") ss >> memo.lon0; 
    else if(arg=="rsmp") ss >> memo.rsmp; 
    else if(arg=="czmi") ss >> memo.czmi; 
    else if(arg=="hinj") ss >> memo.hinj; 
    else if(arg=="ehmi") ss >> memo.ehmi; 
    else if(arg=="eemi") ss >> memo.eemi; 
    else if(arg=="eprp") ss >> memo.eprp; 
    else if(arg=="tlim") ss >> memo.tlim; 
    else if(arg=="rbnd") ss >> memo.rbnd; 
    else if(arg=="epht" || arg=="ejam") ss >> memo.epht; 
    else if(arg=="noem" || arg=="auti" || arg=="dbg1" || arg=="dbg2" || arg=="alpr" || arg=="is1d"){
      std::string tf;
      ss >> tf;       
      bool isTrue;
      if(tf=="T" || tf=="t" || tf=="Yes" || tf=="yes" || tf=="Y" || tf=="y")
	isTrue=true;
      else isTrue = false;
      if(arg=="noem") memo.noem = isTrue;
      if(arg=="auti") memo.auti = isTrue;
      if(arg=="dbg1") memo.dbg1 = isTrue;
      if(arg=="dbg2") memo.dbg2 = isTrue;
      if(arg=="alpr") memo.alpr = isTrue;
      if(arg=="is1d") memo.is1d = isTrue;
    }
    else if(arg=="nwup") ss >> memo.nwup; 
    else if(arg=="seed"){
      ss>>memo.seed[0]>>memo.seed[1];
    }
    else if(arg=="neve") ss >> memo.neve;
    else if(arg=="fcrs") ss >> memo.fcrs;
    else if(arg=="rpft") ss >> memo.rpft;
    else if(arg=="jobt") ss >> memo.jobt;
  
    else if(arg=="site"){
      for(int i=0;i<4;i++){
	ss>>memo.siteParam[memo.nsite][i];
      }
      memo.nsite++;
    }
    else if(arg=="comi") ss >> memo.comi; 
    else if(arg=="nczo") ss >> memo.nczo; 
    else if(arg=="nazo") ss >> memo.nazo; 
    else if(arg=="obsp"){
      std::string kind;
      while(ss>>kind){
	if(kind.at(0)=='#'){
	  break;
	}
	memo.obsp[memo.nkind] = 0;
	for(int i=0;i<ParticleCode::nParType;i++){
	  if(kind==ParticleCode::parName[i]) memo.obsp[memo.nkind] = i;
	}
	memo.nkind++;
      }
    }
    else if(arg=="eomi") ss >> memo.eomi;       
    else if(arg=="etmi") ss >> memo.etmi;
    else if(arg=="nstp"){
      ss >> memo.nstp[0]>>memo.nstp[1];
    }
  }
}
