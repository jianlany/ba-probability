CPP=/opt/pgi/linux86-64-llvm/19.4/bin/pgc++
CFLAGS=-fast -acc -ta=nvidia,cc35 -Minfo=acc

all:
	${CPP} main.cpp ${CFLAGS} -o ba-probability

clean:
	@rm ba-probability
	
