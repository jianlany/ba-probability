# default build mode

CFLAGS=-fast -acc -Minfo=acc -std=c++17
CPP=pgc++

# compute capability of 70, works for V100
70:
	${CPP} main.cpp ${CFLAGS} -ta=nvidia,cc70 -o ba-probability

# compute capability of 6x, works for 1080, 1080Ti
60:
	${CPP} main.cpp ${CFLAGS} -ta=nvidia,cc60 -o ba-probability-60

# compute capability of 35, works for GTX670
35:
	${CPP} main.cpp ${CFLAGS} -ta=nvidia,cc35 -o ba-probability-35

clean:
	@rm ba-probability*
