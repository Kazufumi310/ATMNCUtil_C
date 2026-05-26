#include <myheader/myRoot.h>
#include <myheader/myStandard.h>
#include "atmnc-constants.h"
#include "atmnc-particle-code.h"
#include "ContainerClass.h"

int main(int argc, char** argv){
  if(argc!=3) return 0;
  std::string inputList = argv[1];
  std::string outputFile = argv[2];

  std::ifstream ifs(inputList.c_str());
  if(!ifs){
    std::cout<<"can't find inputList:"<<inputList<<std::endl;
    return 0;
  }

  std::map<double,double> teqSum;
  std::string fname;
  while(std::getline(ifs,fname)){
    TChain ch("trinfo");
    ch.Add(fname.c_str());
    if(ch.GetEntries()==1){
      std::cout<<"read "<<fname<<std::endl;
    }else{
      if(ch.GetEntries()==0){
	std::cout<<"can't find "<<fname<<std::endl;
      }else{
	std::cout<<"strange. The number of entries != 1 in "<<fname<<std::endl;
      }
      return 0;
    }
    ATMNCInputInfo input;
    input.setBranchAddress(&ch);
    ATMNCStatInfo statInfo;
    statInfo.setBranchAddress(&ch);
    ch.GetEntry(0);
    double tequiv = statInfo.tequiv;
    double epmi = input.epmi;
    teqSum[epmi]+=tequiv;
  }

  int nbnd=0;
  double epmi[16],tequiv[16];
  for(std::map<double,double>::iterator it = teqSum.begin();
      it!=teqSum.end();it++){
    epmi[nbnd] = it->first;
    tequiv[nbnd] = it->second;
    nbnd++;
  }
  for(int ibnd=1;ibnd<nbnd;ibnd++){
    tequiv[ibnd]+=tequiv[ibnd-1];
  }

  std::ofstream ofs(outputFile.c_str());
  for(int ibnd=0;ibnd<nbnd;ibnd++){
    ofs<<epmi[ibnd]<<" "<<tequiv[0]/tequiv[ibnd]<<std::endl;
  }
  ofs.close();
}

