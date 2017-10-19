LIB        = -L.  
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = Scan.x
CXX        = g++

${EXEC}: SCAN.c 
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} SCAN.c -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

