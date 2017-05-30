/* hawki_expandmask.cc
 * 
 * USAGE
 *  
 *   hawki_expandmask <imask.fits> <omask.fits> [<regions.reg>]
 *
 * ABOUT
 *
 *   The recipe from the HAWKI reduction pipeline for object detection
 *   does a poor job in determining the masking area which results in
 *   negative sky subtraction in the vicinity of bright objects.
 *
 *   This program will expand the regions of the detected objects, or
 *   if a ds9 regions file is given, it will mask out all objects
 *   specified in this.
 *
 */
#include <string>
#include <vector>
#include <CCfits/CCfits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace CCfits;

float string2float(string str) {
  istringstream foo(str);
  float bar;
  foo >> bar;
  return bar;
}

float string2int(string str) {
  istringstream foo(str);
  int bar;
  foo >> bar;
  return bar;
}

struct PixCoord {
  double x;
  double y;

  PixCoord () : x(0.0), y(0.0) {};
  PixCoord (double tx, double ty) : x(tx), y(ty) {};
  PixCoord (const PixCoord &pc) {
    x = pc.x;
    y = pc.y;
  }

  PixCoord operator=(const PixCoord &c) {
    x = c.x;
    y = c.y;
    return *this;
  }

  PixCoord operator+(const PixCoord &c) {
    x += c.x;
    y += c.y;
    return *this;
  }

  PixCoord operator-(const PixCoord &c) {
    x -= c.x;
    y -= c.y;
    return *this;
  }

  ~PixCoord() {}
} Q1(0,0), Q2(2048+153,0+3), Q3(2048+157,2048+144), Q4(0+5,2048+142);


// chipOffsets.push_back(PixCoord());
// chipOffsets.resize(4);
// chipOffsets[0] = PixCoord();



void Tokenize(const string& str, vector<string>& tokens,
	      const string& delimiters);

const int MAX_LINE = 256;
long hw = 10;

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr,
	    "Usage: hawki_expandmask <imask.fits> <omask.fits> [<regions>]\n");
    exit(0);
  }

  string inFileName(argv[1]);
  string outFileName(argv[2]);

  // HAWKI chip offsets
  //
  vector<PixCoord> chipOffsets(4);
  chipOffsets[0] = Q1;
  chipOffsets[1] = Q2;
  chipOffsets[2] = Q3;
  chipOffsets[3] = Q4;

  // Use a regions file where the objects are defined
  //
  bool useRegions = false;
  if (argc == 4) useRegions = true;
  vector<float> x, y, dx, dy;

  if (useRegions) {
    char *regionsFile = argv[3];

    ifstream regions;
    regions.open(regionsFile);
    while(!regions.eof()) {
      char foo[MAX_LINE];
      regions.getline(foo,MAX_LINE);
      string object(foo);
      if (object.find("box",0) == 0) {  // Object starts with box
	string bar = object.substr(4);
	vector<string> coords;
	Tokenize(object.substr(4), coords, ",");

	x.push_back(string2float(coords[0]));
	y.push_back(string2float(coords[1]));
	dx.push_back(0.5*string2float(coords[2]));
	dy.push_back(0.5*string2float(coords[3]));
      }
    }
    regions.close();
  }

  remove(outFileName.c_str());          // Delete old outfile if it exists

  FITS inFile(inFileName,Read,true);
  FITS outFile(outFileName,inFile);

  for (int c=1; c<=4; c++) {
    ExtHDU& imageExt = inFile.extension(c);
    imageExt.readAllKeys();
    
    valarray<float> mask, newmask;
    imageExt.read(mask);
    newmask.resize(mask.size());
    
    long nx(imageExt.axis(0));
    long ny(imageExt.axis(1));

    if (useRegions) {
      newmask = mask;     // Keep the original mask, just add the regions
      
      int chip = imageExt.name().at(4) - '1'; // Chip ID

      // Loop over objects
      //
      for (unsigned int n=0; n<x.size(); n++) {
	int cx = x[n] - chipOffsets[chip].x;
	int cy = y[n] - chipOffsets[chip].y;

	// Object is on mask
	//
	if ( cx >= 0 && cy >= 0 && cx < nx && cy < ny ) {

	  // Loop over all pixels in the vicinity of the object
	  // centre
	  //
	  for (int yy=cy-dy[n]; yy<=cy+dy[n]; yy++) {
	    for (int xx=cx-dx[n]; xx<=cx+dx[n]; xx++) {
	      
	      // Angle and distance between the object centre
	      // and the current pixel coordinate
	      //
	      float th = atan2(yy-cy,xx-cx);
	      float r = sqrt((xx-cx)*(xx-cx) + (yy-cy)*(yy-cy));

	      // Mask pixel if it is on the image and within
	      // the ellipse defined by the region x and y widths.
	      //
	      //	      cout << xx << ' ' << yy <<  ' ' << r << endl;
	      if (xx >= 0 && yy >= 0 && xx < nx && yy < ny &&
		  r < sqrt( (dx[n]*cos(th))*(dx[n]*cos(th)) +
			    (dy[n]*sin(th))*(dy[n]*sin(th)) ) )
		newmask[yy*nx+xx] = 1;
	    }
	  } 
	}	
      }
    } else { // Enlarge masked areas
      // Start with a zero mask
      //
      for (unsigned int n=0; n < newmask.size(); n++) newmask[n] = 0;

      // Loop over pixels in mask
      //
      for (long j=0; j<ny; j++) {
	for (long i=0; i<nx; i++) {
	  
	  // Loop over surrounding pixels
	  //
	  if (mask[j*nx+i] == 1) {
	    for (long jj=j-hw; jj<=j+hw; jj++) {
	      for (long ii=i-hw; ii<=i+hw; ii++) {
		if (ii >= 0 && jj >= 0 && ii < nx && jj < ny)
		  newmask[jj*nx+ii] = 1;
	      }
	    }
	  }
	  
	}
      }
    }
    
    vector<long> extAx(2);
    extAx[0] = nx;
    extAx[1] = ny;
    long fpixel(1);

    ExtHDU* newext = outFile.addImage(imageExt.name(),FLOAT_IMG,extAx);
    newext->copyAllKeys(&imageExt);
    newext->write(fpixel,newmask.size(),newmask);
    //    newext->write(fpixel,mask.size(),mask);
    // cout << outFile.pHDU() << endl;


    //    free(&newmask);
    //    free(&mask);
    //    free(newext);
    //    free(&imageExt);
   }

  inFile.destroy();
  outFile.destroy();

  return 0;
}


void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}


