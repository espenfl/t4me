TARGET = skw
SPGLIB_INCLUDE = $(HOME)/include
SPGLIB_LIB = $(HOME)/libs
LIBS = -L$(SPGLIB_LIB) -lsymspg $(MKLROOT)/lib/intel64/libfftw3xc_intel.a  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
INCLUDE = -I$(SPGLIB_INCLUDE) -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
CC = g++
CFLAGS = -fPIC -std=c++11
LD = g++
#DEBUG = -g -DINFO -DDEBUG_FFT
#DEBUG = -DINFO 

default: $(TARGET)

lib: $(TARGET).o
	$(LD) -shared -o lib$(TARGET).so $(TARGET).o
$(TARGET): $(TARGET).o
	$(CC) $(TARGET).o $(LIBS) -o $(TARGET)
$(TARGET).o: $(TARGET).cpp
	$(CC) $(DEBUG) $(CFLAGS) $(INCLUDE) -c $(TARGET).cpp -o $(TARGET).o


clean:
	rm -f skw.o
	rm -f skw
