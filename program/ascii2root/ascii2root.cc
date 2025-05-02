#include <FileHandler.h>

int main(int argc, char** argv){
  if(argc!=3){
    std::cout<<"arg error"<<std::endl;
    std::cout<<"usage: ascii2root input.d output.root"<<std::endl;
    return 0;
  }
  std::string input = argv[1];
  std::string output = argv[2];

  FileHandler filer;
  filer.ascii2root(input,output);
}
