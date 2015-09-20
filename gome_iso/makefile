

BASE_DIR=/home/home/pmills
LIB_DIR=$(BASE_DIR)/lib
INCLUDE_DIR=$(BASE_DIR)/include

CPP=g++
CFLAGS=-I$(INCLUDE_DIR) -g

LDFLAGS=-L$(LIB_DIR) -lpetey

all: convert_gome_pmd print_gome_coords convert_gome print_gome_interpolated

convert_gome_pmd: read_amsu.o convert_gome_pmd.cc
	$(CPP) $(CFLAGS) -o$@ convert_gome_pmd.cc read_amsu.o $(LDFLAGS)

convert_gome: read_gome.o convert_gome.cc rawread_gome.o
	$(CPP) $(CFLAGS) -o$@ convert_gome.cc read_gome.o rawread_gome.o $(LDFLAGS)

print_gome_coords: read_gome.o print_gome_coords.cc rawread_gome.o
	$(CPP) $(CFLAGS) -o$@ print_gome_coords.cc read_gome.o rawread_gome.o $(LDFLAGS)

print_gome_interpolated: read_gome.o print_gome_interpolated.cc rawread_gome.o
	$(CPP) $(CFLAGS) -o$@ print_gome_interpolated.cc read_gome.o rawread_gome.o $(LDFLAGS)

read_amsu.o: read_amsu.cc read_amsu.h
	$(CPP) $(CFLAGS) -c read_amsu.cc

read_gome.o: read_gome.cc read_gome.h
	$(CPP) $(CFLAGS) -c read_gome.cc

rawread_gome.o: rawread_gome.cc read_gome.h
	$(CPP) $(CFLAGS) -c rawread_gome.cc

