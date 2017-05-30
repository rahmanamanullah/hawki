#include <iostream>
#include <string>
#include <vector>
#include <CCfits/CCfits>

using namespace std;
using namespace CCfits;

int main(int argc, char **argv) {
  
  string frFileName = string(argv[1]);
  FITS* frFile = new FITS(frFileName,Read,true);
  PHDU& header = frFile->pHDU();
  header.readAllKeys();

  for (int n=2; n<argc; n++) {
    string toFileName = string(argv[n]);
    FITS* toFile = new FITS(toFileName,Write,true);
    
    PHDU& theader = toFile->pHDU();
    

  }

  return 0;
}

// EOF

