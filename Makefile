# default build mode

CFLAGS=-fast -acc -Minfo=acc -std=c++17
ifneq (,$(findstring jjo5,${HOSTNAME}))
CPP=/opt/pgi/linux86-64-llvm/19.4/bin/pgc++
else
CPP=pgc++
endif

# compute capability of 70,75,60, works for all the gpus on agave
agave:
	${CPP} main.cpp ${CFLAGS} -ta=nvidia,cc35,cc70,cc75,cc60 -o ba-probability


# compute capability of 35, works for GTX670
35:
	${CPP} main.cpp ${CFLAGS} -ta=nvidia,cc35 -o ba-probability-35

clean:
	@rm ba-probability*
