#include <iostream>
#include <fstream>
#include "stdlib.h"

int makeFillBunches(int bunchNumOffset)
{
  std::ofstream bunchFile;
  bunchFile.open("bunchFile.csv");

  int bunchNum = bunchNumOffset;
  int sumNum[2] = {5, 10};
  
  for(int iter = 0; iter < 18; iter++){
    for(int iter2 = 0; iter2 < 12; iter2++){
      for(int iter3 = 0; iter3 < 2; iter3++){
	bunchFile << bunchNum << ",\n";
	bunchNum += sumNum[iter3];
      }
    }
    bunchNum += 19;
  }
  bunchFile.close();

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc != 2){
    std::cout << "Usage: makeFillBunches <bunchNumOffset>" << std::endl;

    std::cout << "argc: " << argc << std::endl;
    for(int iter = 0; iter < argc; iter++){
      std::cout << "argv[" << iter << "]: " << argv[iter] << std::endl;
    }

    return 1;
  }

  int rStatus = -1;
  rStatus = makeFillBunches(atoi(argv[1]));
  return rStatus;
}
