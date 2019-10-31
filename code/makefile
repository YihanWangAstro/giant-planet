all: single-single single-binary binary-single binary-binary

VPATH =${MY_REPO_PATH}/SpaceX:${MY_REPO_PATH}/SpaceHub/src/
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall  
OBJS = dynamicSystem.o InitFileCreator.o \
	 chain.o chainSystem.o tools.o\
	 IO.o symplectic.o regularizator.o BSIterator.o\
	particleSystem.o integrator.o seqIterator.o Wrapper.o
	
INCL = dynamicSystem.h macros.h chain.h\
	regularizator.h particleSystem.h interaction.h chainSystem.h \
	integrator.h BSIterator.h errhand.h InitFileCreator.h\
	seqIterator.h odeIterator.h vector3d.h tools.h

single-single : $(OBJS) single-single.o
	$(CXX) -std=c++17 -O3 -o single-single.ex $(OBJS) single-single.o -pthread

binary-single : $(OBJS) binary-single.o
	$(CXX) -std=c++17 -O3 -o binary-single.ex $(OBJS) binary-single.o -pthread

single-binary : $(OBJS) single-binary.o
	$(CXX) -std=c++17 -O3 -o single-binary.ex $(OBJS) single-binary.o -pthread

binary-binary : $(OBJS) binary-binary.o
	$(CXX) -std=c++17 -O3 -o binary-binary.ex $(OBJS) binary-binary.o -pthread

$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) *.o *.ex
