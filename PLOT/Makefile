LIB        = -L.  
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = Plot.x
CXX        = g++

${EXEC}: PLOT.c 
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} PLOT.c -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

