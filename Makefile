include Makefile.common


all:  verticalsegment segment


AtmTestGrISU: AtmTestGrISU.o atmosphere.o kasatmGrISU.o
	$(CXX) $(LDFLAGS) $^ -o AtmTestGrISU 

test2: test2.o atmosphere.o kasatmGrISU.o
	$(CXX) $(LDFLAGS)  $(ROOTLIBS) $^ -o test2

test2.o: test2.cpp
	$(CXX) $(CXXFLAGS) -c test2.cpp -o test2.o

test3: test3.o atmosphere.o kasatmGrISU.o
	$(CXX) $(LDFLAGS) $(ROOTLIBS) $^ -o test3

segment: segment.o atmosphere.o kasatmGrISU.o
	$(CXX) $(LDFLAGS) $(ROOTLIBS) $^ -o segment

verticalsegment: verticalsegment.o atmosphere.o kasatmGrISU.o
	$(CXX) $(LDFLAGS) $(LIBS) $^ -o verticalsegment

kasatmGrISU.o: kasatmGrISU.cpp kasatmGrISU.h 
	$(CXX) $(CXXFLAGS) -c kasatmGrISU.cpp -o kasatmGrISU.o


atmosphere.o: atmosphere.h atmosphere.cpp
	$(CXX) $(CXXFLAGS) -c atmosphere.cpp -o atmosphere.o

AtmTestGrISU.o: AtmTestGrISU.cpp
	$(CXX) $(CXXFLAGS) -c AtmTestGrISU.cpp -o AtmTestGrISU.o

clean:
	rm *.o
	rm AtmTestGrISU test2 test3 verticalsegment segment

.PHONY: flags
flags:
	@echo "CXXFLAGS:  $(CXXFLAGS)"
	@echo "LDFLAGS:   $(LDFLAGS)"
