all: test.o O3Tetrahedron.o O3Triangulation.o
	g++ $$(regina-engine-config --libs) test.o O3Tetrahedron.o O3Triangulation.o -o enumerate

test.o: test.cpp
	g++ -g -c $$(regina-engine-config --cflags) test.cpp -o test.o
O3Tetrahedron.o: O3Tetrahedron.cpp O3Triangulation.h
	g++ -g -c $$(regina-engine-config --cflags) O3Tetrahedron.cpp -o O3Tetrahedron.o
O3Triangulation.o: O3Triangulation.cpp O3Tetrahedron.h
	g++ -g -c $$(regina-engine-config --cflags) O3Triangulation.cpp -o O3Triangulation.o

