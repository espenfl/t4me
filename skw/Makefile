TARGET = skw
MKLROOT = /opt/intel/mkl
SPGLIB = /home/flagre/local/spglib-1.7.4
FFTW = /home/flagre/local/fftw-3.3.4
LIBS = $(SPGLIB)/lib/libsymspg.a $(FFTW)/lib/libfftw3.a  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
#LIBS = $(SPGLIB)/lib/libsymspg.a $(FFTW)/lib/libfftw3.a  -Wl,--no-as-needed -lpthread -lm -ldl
INCLUDE = -I/home/flagre/local/spglib-1.7.4/include/spglib -I$(MKLROOT)/include
CC = g++
CFLAGS = -fPIC -std=c++11
LD = ar rvs
#DEBUG = -g -DINFO -DDEBUG_FFT
DEBUG = -DINFO 

default: $(TARGET)

lib: $(TARGET).o
	$(LD) lib$(TARGET).a $(TARGET).o
$(TARGET): $(TARGET).o
	$(CC) $(TARGET).o $(LIBS) -o $(TARGET)
$(TARGET).o: $(TARGET).cpp
	$(CC) $(DEBUG) $(CFLAGS) $(INCLUDE) -c $(TARGET).cpp -o $(TARGET).o


clean:
	rm -f skw.o
	rm -f skw
