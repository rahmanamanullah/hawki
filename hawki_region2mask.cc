/* hawki_region2mask.cc
 * 
 * USAGE
 *  
 *   hawki_region2mask <dbimg1> <dbimg2> ...
 *
 * ABOUT
 *
 *   The recipe from the HAWKI reduction pipeline for object detection
 *   does a poor job in determining the masking area which results in
 *   negative sky subtraction in the vicinity of bright objects.
 *
 *   This program will convert the ds9 regions file called 'objects.reg'
 *   in each dbimage to an 'objectmask.fits' file.
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
  if (argc < 2) {
    fprintf(stderr,
	    "Usage: hawki_region2mask <dbimag1> <dbimag2> ...\n");
    exit(0);
  }

  // HAWKI chip offsets
  //
  vector<PixCoord> chipOffsets(4);
  chipOffsets[0] = Q1;
  chipOffsets[1] = Q2;
  chipOffsets[2] = Q3;
  chipOffsets[3] = Q4;

  for (int imn=1; imn < argc; imn++) {
    string dbimage(argv[imn]);
    
    string regnFileName(dbimage+"/objects.reg");
    string maskFileName(dbimage+"/objectmask.fits");
    
    vector<float> x, y, a, b, theta;
    
    // Read regions file
    //
    const char *regionsFile = regnFileName.c_str();
    ifstream regions;
    regions.open(regionsFile);
    while(!regions.eof()) {
      char foo[MAX_LINE];
      regions.getline(foo,MAX_LINE);
      string object(foo);
      if (object.find("ellipse",0) == 0) {  // Object starts with 'ellipse'
	vector<string> coords;
	int i1 = object.find_first_of("(",0);
	int i2 = object.find_first_of(")",i1);

	Tokenize(object.substr(i1+1,i2-i1-1), coords, ",");

	x.push_back(string2float(coords[0]));
	y.push_back(string2float(coords[1]));
	a.push_back(1.5*string2float(coords[2]));
	b.push_back(1.5*string2float(coords[3]));
	theta.push_back(3.14159265*string2float(coords[4])/180.0);
      }
    }
    regions.close();

    remove(maskFileName.c_str());          // Delete old outfile if it exists

    FITS bpmFile(dbimage+"/bpm.fits");
    FITS outFile(maskFileName,bpmFile);
    
    for (int c=1; c<=4; c++) {
      ExtHDU& imageExt = bpmFile.extension(c);
      imageExt.readAllKeys();
      
      valarray<float> mask, newmask;
      imageExt.read(mask);
      newmask.resize(mask.size());
      
      long nx(imageExt.axis(0));
      long ny(imageExt.axis(1));
      
      newmask = 0.0;
      int chip = imageExt.name().at(4) - '1'; // Chip ID
      
      // Loop over objects
      //
      for (unsigned int n=0; n<x.size(); n++) {
	int cx = x[n] - chipOffsets[chip].x;
	int cy = y[n] - chipOffsets[chip].y;
	
	// Object is on mask or not (use the semi-major axis and the object
	// centre to determine this).
	//
	if ( cx+a[n] >= 0 && cy+a[n] >= 0 && cx-a[n] < nx && cy-a[n] < ny ) {
	  
	  // Loop over all pixels in the vicinity of the object
	  // centre (use the semi-major axis).
	  //
	  for (int yy=cy-a[n]; yy<=cy+a[n]; yy++) {
	    for (int xx=cx-a[n]; xx<=cx+a[n]; xx++) {
	      
	      // Angle and distance between the object centre
	      // and the current pixel coordinate
	      //
	      float th = atan2(yy-cy,xx-cx);
	      float r2 = (xx-cx)*(xx-cx) + (yy-cy)*(yy-cy);
	      
	      // Mask pixel if it is on the image and within
	      // the ellipse defined by the region x and y widths.
	      //
	      //	      cout << xx << ' ' << yy <<  ' ' << r << endl;
	      float foo = cos(th-theta[n])/a[n];
	      float bar = sin(th-theta[n])/b[n];
	      float re2 = 1.0/( foo*foo + bar*bar );
	      
	      //	      cout << n << ' ' << r2 << ' ' << re2 << endl;

	      if (xx >= 0 && yy >= 0 && xx < nx && yy < ny && r2 < re2 ) {
		newmask[yy*nx+xx] = 1;
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
    bpmFile.destroy();
    outFile.destroy();
  }
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


