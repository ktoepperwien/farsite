# Builds TestFARSITE on Linux and macOS with no external libraries.
#
# Deliberately free of platform detection: every source file compiles with the
# same flags on both platforms. Use CXX=clang++ (or any other compiler) on the
# command line to override.

CXX ?= g++
INC =
# -MMD -MP (added at the compile rule) emit .d files so that editing a header
# rebuilds its dependents.
#
# -O2 is ~2.5x faster than the unoptimised build this Makefile used to produce.
# It gives identical results on reproducible cases, but FARSITE's spread front is
# threshold-driven, so on floating-point-sensitive cases the optimisation level
# can change output. To reproduce the numbers an unoptimised build produces:
#
#     make clean && make CXXFLAGS="-std=c++11 -g -Wall -DUNIX -Wno-deprecated"
#
CXXFLAGS ?= -std=c++11 -O2 -g -Wall -DUNIX -Wno-deprecated
CPPFLAGS ?=

SOURCES=barriers.cpp burnupw.cpp cdtlib.cpp \
	crossthread.cpp dbfopen.cpp Far_BSG.cpp Far_Cond.cpp Far_FPC.cpp \
	FARSITE.cpp fsairatk.cpp fsxpfront.cpp fsxwaccl.cpp \
	fsxwatm.cpp fsxwattk.cpp fsxwbar.cpp fsxwburn4.cpp fsxwcrwn.cpp \
	fsxwenvt2.cpp fsxwfms2.cpp fsxwfotp.cpp fsxwignt.cpp fsxwmech.cpp \
	fsxwrast.cpp fsxwshap.cpp fsxwspot.cpp fsxwutil.cpp fsxwvect.cpp \
	gridthem.cpp newclip.cpp newfms.cpp nn_cls.cpp PerimeterData.cpp \
	polygon.cpp shpopen.cpp vec.cpp icf_input.cpp icf_chk.cpp \
	deadfuelmoisture.cpp FMC_CFMC.cpp Farsite5.cpp\
	FMC_CI.cpp FMC_CI2.cpp FMC_CI3.cpp FMC_CI4.cpp FMC_FE2.cpp \
	FMC_FE22.cpp FMC_FE3.cpp semtime.cpp rand3.cpp TestFARSITE.cpp

OBJECTS=$(SOURCES:.cpp=.o)
DEPS=$(OBJECTS:.o=.d)
EXECUTABLE=TestFARSITE

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@

%.o: %.cpp
	$(CXX) $(INC) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f *.o *.d $(EXECUTABLE)

-include $(DEPS)
