DEBUG=
all: build query
	
build: h2h.cpp Readfile.o Heap.o Graph.h
	g++ $(DEBUG) -o build -O3 h2h.cpp lca.h Readfile.o Graph.h Heap.o

buildtime: h2htime.cpp Readfile.o Heap.o Graph.h
	g++ $(DEBUG) -o buildtime -O3 h2htime.cpp lca.h Readfile.o Graph.h Heap.o
query: query_orignal.cpp Graph.h
	g++ $(DEBUG) -o query query_orignal.cpp
queryO3: query_orignal.cpp Graph.h
	g++ $(DEBUG) -o queryO3 -O3 query_orignal.cpp
collect: collectBranches.cpp Graph.h
	g++ $(DEBUG) -o collect -O3 collectBranches.cpp
queryO3Nf: query_orignal_Nf.cpp Graph.h
	g++ $(DEBUG) -o queryO3Nf -O3 query_orignal_Nf.cpp
queryvd: queryvd.cpp Graph.h bitvector.o
	g++ $(DEBUG) -o queryvd queryvd.cpp bitvector.o -std=c++11 -mavx2 -m64
queryvdO3: queryvd.cpp Graph.h bitvector.o
	g++ $(DEBUG) -o queryvdO3 -O3 queryvd.cpp bitvector.o -std=c++11 -mavx2 -m64
querydd: querydd.cpp Graph.h
	g++ $(DEBUG) -o querydd querydd.cpp
queryddO3: querydd.cpp Graph.h
	g++ $(DEBUG) -o queryddO3 -O3 querydd.cpp
queryddO3Nf: queryddNf.cpp Graph.h
	g++ $(DEBUG) -o queryddO3Nf -O3 queryddNf.cpp

Readfile.o: Readfile.cpp Readfile.h Heap.o
	g++ $(DEBUG) -c Readfile.cpp 

Heap.o: Heap.cpp Heap.h Graph.h
	g++ $(DEBUG) -c Heap.cpp 

bitvector.o: bitvector.cpp bitvector.h macros.h
	g++ $(DEBUG) -c bitvector.cpp -mavx2 -m64 -std=c++11
qbit: querybit.cpp bitvector.o
	g++ $(DEBUG) -o qbit -O3 querybit.cpp bitvector.o -mavx2 -m64 -std=c++11 -lpthread
qbyte: query_byte.cpp bitvector.o
	g++ $(DEBUG) -o qbyte query_byte.cpp bitvector.o -mavx2 -m64 -std=c++11 -lpthread -O3
qbitp: query_bit.cpp bitvector.o
	g++ $(DEBUG) -o qbitp  query_bit.cpp bitvector.o -mavx2 -m64 -std=c++11 -lpthread -O3
qbool: querybool.cpp bitvector.o
	g++ $(DEBUG) -o qbool -O3 querybool.cpp bitvector.o -mavx2 -m64 -std=c++11 
qmul: querymulthread.cpp Graph.h
	g++ $(DEBUG) -o qmul -O3 querymulthread.cpp -lpthread -std=c++11
qnew: query_new.cpp bitvector.o
	g++ $(DEBUG) -o qnew -O3 query_new.cpp bitvector.o -std=c++11 -mavx2 -m64
qlast: query_last.cpp
	g++ $(DEBUG) -o qlast -O3 query_last.cpp
