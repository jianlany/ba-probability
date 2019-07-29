# default build mode
target ?= ubuntu

ifeq (${target},agave)
CPP=pgc++
CFLAGS=-fast -acc -ta=nvidia,cc70 -Minfo=acc -std=c++17
endif

ifeq (${target},ubuntu)
CPP=/opt/pgi/linux86-64-llvm/19.4/bin/pgc++
CFLAGS=-fast -acc -ta=nvidia,cc35 -Minfo=acc
endif

all:
	${CPP} main.cpp ${CFLAGS} -o ba-probability

clean:
	@rm ba-probability
	
