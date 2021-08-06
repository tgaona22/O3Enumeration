enumerate: enumerate.o O3Tetrahedron.o O3Triangulation.o
	g++ $$(regina-engine-config --libs) enumerate.o O3Tetrahedron.o O3Triangulation.o -o enumerate

analyze: analyze.o O3Tetrahedron.o O3Triangulation.o
	g++ $$(regina-engine-config --libs) analyze.o O3Tetrahedron.o O3Triangulation.o -o analyze

analyze.o: analyze.cpp
	g++ -g -c $$(regina-engine-config --cflags) analyze.cpp -o analyze.o

enumerate.o: enumerate.cpp
	g++ -g -c $$(regina-engine-config --cflags) enumerate.cpp -o enumerate.o
O3Tetrahedron.o: O3Tetrahedron.cpp O3Triangulation.h
	g++ -g -c $$(regina-engine-config --cflags) O3Tetrahedron.cpp -o O3Tetrahedron.o
O3Triangulation.o: O3Triangulation.cpp O3Tetrahedron.h
	g++ -g -c $$(regina-engine-config --cflags) O3Triangulation.cpp -o O3Triangulation.o

