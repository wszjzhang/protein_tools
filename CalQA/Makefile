                                                                                                                      
CXX = g++
CXXFLAGS = -O2 -Wall

DIFF = ./sdiff
PRE = ./
MAJOR = 1
MINOR = 0

%.o:            %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp 

everything:            calGdtRmsdScore

calGdtRmsdScore_obj = calGdtRmsdScore.o TM_score.o 

calGdtRmsdScore:       $(calGdtRmsdScore_obj) 
	$(CXX) -o $@ $(calGdtRmsdScore_obj) -L. -lm

calGdtRmsdScore.o:     calGdtRmsdScore.cpp
TM_score.o:            TM_score.cpp

calGdtRmsdScore.txx:     calGdtRmsdScore
	$(PRE) calGdtRmsdScore > calGdtRmsdScore.txt
	$(DIFF) calGdtRmsdScore.txt calGdtRmsdScore.txx

