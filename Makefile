FLAST = `pkg-config --cflags --libs opencv`


all: main.o GA.o
	g++ -ggdb main.o GA.o -o main $(FLAST)
GA.o: GA.cpp GA.h
	g++ -ggdb GA.cpp -c $(FLAST)
main.o: main.cpp GA.h
	g++ main.cpp -c $(FLAST)
