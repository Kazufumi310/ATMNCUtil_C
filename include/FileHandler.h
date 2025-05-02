#ifndef FileHandler_h__
#define FileHandler_h__

#include <myheader/myRoot.h>
#include <myheader/myStandard.h>
#include <ContainerClass.h>
#include <Histogrammer.h>
#include <atmnc-constants.h>
#include <atmnc-particle-code.h>

class FileHandler{
 public:
  FileHandler();
  void ascii2root(std::string inputAscii, std::string outputRoot);
  void sum(std::vector<std::string> inputFileList, std::string outputFile);

private:
  // ascii2root
  std::ifstream m_ifs;
  void readInputInfo( ATMNCInputInfo& memo, std::ifstream &ifs );
  void read1ryInfo( ATMNCStatInfo& memo, std::ifstream &ifs );
  bool read1Histogram( ATMNCHistInfo &memo, std::ifstream &ifs );
};

#endif // FileHandler_h__
