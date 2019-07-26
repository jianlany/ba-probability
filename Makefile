
CPP=/opt/pgi/linux86-64-llvm/19.4/bin/pgc++

all:
	${CPP} main.cpp -fast -o ba-probability

clean:
	@rm ba-probability
	
