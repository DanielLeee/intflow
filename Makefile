
.PHONY = all clean

CC = g++
CXX = g++
CXXFLAGS = -Wall

all: intFlow labelTransfer
	
intFlow: 3D_SIFT_FLOW.o BPFlow.o Image3D.o
	$(CC) -o $@ 3D_SIFT_FLOW.o BPFlow.o Image3D.o

labelTransfer: LabelTransfer.o Image3D.o
	$(CC) -o $@ LabelTransfer.o Image3D.o

clean:
	rm -fv *.o intFlow labelTransfer warp


