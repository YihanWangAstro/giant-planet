all: binary-binary binary-single

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

binary-single : $(OBJS) single-cross.o
	$(CXX) -std=c++17 -O3 -o binary-single $(OBJS) single-cross.o -pthread
	
binary-binary : $(OBJS) binary-cross.o
	$(CXX) -std=c++17 -O3 -o binary-binary $(OBJS) binary-cross.o -pthread
	
$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) *.o 
