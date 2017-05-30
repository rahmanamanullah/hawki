CC = gcc
CPP = gcc -E
CPPFLAGS = -g -Wall
CXX = g++
CXXCPP = g++ -E
#CXXFLAGS = -g -O2 -Wall -pedantic -I/usr/include/CCfits
#CXXLIBS = -L/usr/lib -lCCfits
CXXFLAGS = -g -O2 -Wall -I$(SNOVASOFTARCH)/include -I$(SNOVASOFTARCH)/include/CCfits
CXXLIBS = -L$(SNOVASOFTARCH)/lib -lCCfits

all : hawki_buildskymap hawki_buildskymap2 hawki_smoothskymap \
	hawki_skysubtract hawki_flatfield hawki_expandmask \
	hawki_region2mask hawki_splitfits

hawki_splitfits : hawki_splitfits.cc Makefile dbimage.h dbimage.cc
	$(CXX) $(CXXFLAGS) -o $@ $@.cc $(CXXLIBS)

hawki_region2mask : hawki_region2mask.cc Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc $(CXXLIBS)

hawki_expandmask : hawki_expandmask.cc Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc $(CXXLIBS)

hawki_buildskymap : hawki_buildskymap.cc \
	hawki_image.o dbimage.o stats.h Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc hawki_image.o dbimage.o \
	$(CXXLIBS)

hawki_buildskymap2 : hawki_buildskymap2.cc \
	hawki_image.o dbimage.o stats.h Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc hawki_image.o dbimage.o \
	$(CXXLIBS)

hawki_smoothskymap : hawki_smoothskymap.cc \
	hawki_image.o dbimage.o stats.h Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc hawki_image.o dbimage.o \
	$(CXXLIBS)

hawki_skysubtract : hawki_skysubtract.cc \
	hawki_image.o dbimage.o Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc hawki_image.o dbimage.o $(CXXLIBS)

hawki_flatfield : hawki_flatfield.cc \
	hawki_image.o dbimage.o Makefile
	$(CXX) $(CXXFLAGS) -o $@ $@.cc hawki_image.o dbimage.o $(CXXLIBS)

hawki_image.o : hawki_image.cc hawki_image.h dbimage.o
	$(CXX) $(CXXFLAGS) -c $< -o $@

dbimage.o : dbimage.cc dbimage.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

runningmean.o : runningmean.c
	$(CC) -c $< -o $@

clean :
	rm -rf *.o hawki_flatfield hawki_skysubtract hawki_expandmask hawki_splitfits hawki_buildskymap hawki_region2mask
