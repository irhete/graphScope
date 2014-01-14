CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -fopenmp -std=c++0x

OBJS =		GraphScope.o

LIBS = -lgomp

TARGET =	GraphScope

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
