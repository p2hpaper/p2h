DEBUG=
all: build query queryO3 queryF queryFO3
	
build: h2h.cpp Readfile.o Heap.o Graph.h
	g++ $(DEBUG) -o build2 -O3 h2h.cpp lca.h Readfile.o Graph.h Heap.o

Readfile.o: Readfile.cpp Readfile.h Heap.o
	g++ $(DEBUG) -c Readfile.cpp 

Heap.o: Heap.cpp Heap.h Graph.h
	g++ $(DEBUG) -c Heap.cpp 
