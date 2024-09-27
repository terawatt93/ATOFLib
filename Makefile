CXX = `root-config --cxx`
CXXFLAGS = `root-config --cflags` -lSpectrum -fPIC -g -Wall 
ROOTLIBS = `root-config --glibs` -lSpectrum
SHARED = -shared
SRCS = ATOFLibDict.cxx ATOFLib.cpp
HDRS = ATOFLibLinkDef.h ATOFLib.hh
PROGRAM = ATOFLib.so

all: $(PROGRAM)

ATOFLibDict.cxx: $(HDRS) ATOFLibLinkDef.h
	@echo "Generating dictionary ..."
	#@rootcint -f $@ -c -p $^
	@rootcling -f ATOFLibDict.cxx -rml ATOFLib.so -rmf ATOFLib.rootmap ATOFLib.hh ATOFLibLinkDef.h

$(PROGRAM): $(SRCS)
	@echo "Building $(PROGRAM) ..."
	@rm -f $(PROGRAM)
	@$(CXX) $(CXXFLAGS) $(SHARED) -o $@ $^ $(ROOTLIBS)
	@echo "done"
#options:
clean:; @rm -rf core *.so *.rootmap *.cxx *.pcm
