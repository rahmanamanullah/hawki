#include <string>
#include "hawki_image.h"

using namespace std;

int main(int argc, char **argv) {

  for (int i=1; i<argc; i++) {
    string dbimdir = string(argv[i]);

    HAWKIimage dbim(dbimdir);
    dbim.flatFieldImage();
  }

  return 0;
}
