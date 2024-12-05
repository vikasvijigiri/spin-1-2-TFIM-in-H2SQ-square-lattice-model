CXX = mpic++
CXXFLAGS = -O3 -std=c++17 -Wall  -pedantic-errors -g -msse2 -march=native -funroll-loops -ffast-math #-fomit-frame-pointer -fstrict-aliasing
#-O3 -g -Wall -msse2 -march=native -funroll-loops -ffast-math -fomit-frame-pointer -fstrict-aliasing
SRCS =  *.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = *.hpp

MAIN = main

all: ${MAIN}
	@echo   Compilation Successful.................

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c ${SRCS}

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
