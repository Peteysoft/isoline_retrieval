OPT_VER = -g
CC = g++

BASE_PATH = /home/home/pmills
LIB_PATH = $(BASE_PATH)/lib
BIN_PATH = $(BASE_PATH)/bin
INCLUDE_PATH = $(BASE_PATH)/include
TRAJ_PATH = $(BASE_PATH)/my_software/ctraj

CFLAGS = $(OPT_VER) -I$(INCLUDE_PATH) -I$(TRAJ_PATH)
CFLAGS2 = $(OPT_VER) -I$(INCLUDE_PATH) 
LDFLAGS = -L$(LIB_PATH) -lagf$(OPT_VER) -ldataset$(OPT_VER) -lpetey$(OPT_VER) -lgsl -lgslcblas

EXECUTABLES=print_amsu_1c_file \
	interpolate_amsu \
	index_amsu \
	index_gome \
	retrieve_isoline_g \
	convert_gome \
	gome_colloc2 \
	print_gome_raw \
	extract_gome_colloc2 \
	o3_class \
	retrieve_gome_classes

all: $(EXECUTABLES)

install:
	cp $(EXECUTABLES) $(BIN_PATH)

clean:
	rm -f $(EXECUTABLES)
	rm -f *.o 

# retrieval:


gome_colloc: gome_colloc.cc read_amsu.o amsu_interp_obj.o amsu_interp_base.o \
		amsu_util.o global_metric.o amsu_tinterp_obj.o
	$(CC) $(CFLAGS) -o gome_colloc gome_colloc.cc \
		amsu_interp_obj.o read_amsu.o amsu_util.o global_metric.o \
			amsu_tinterp_obj.o amsu_interp_base.o \
			$(TRAJ_PATH)/traj_int_obj.o $(TRAJ_PATH)/tcoord_defs.o \
			$(LDFLAGS)

gome_colloc2: gome_colloc2.cc read_gome.o colloc_obj.o colloc_base.o \
		global_metric.o tcolloc_obj.o
	$(CC) $(CFLAGS) -o gome_colloc2 gome_colloc2.cc \
		colloc_obj.o read_gome.o global_metric.o \
			tcolloc_obj.o colloc_base.o \
			$(TRAJ_PATH)/traj_int_obj.o $(TRAJ_PATH)/tcoord_defs.o \
			$(LDFLAGS)

retrieve_gome_classes: retrieve_gome_classes.cc read_gome.o read_spec_file.cc
	$(CC) $(CFLAGS) -I$(INCLUDE_PATH) -o retrieve_gome_classes retrieve_gome_classes.cc \
		read_gome.o $(LDFLAGS)

retrieve_isoline_g: retrieve_isoline_g.cc read_amsu.o amsu_util.o global_metric.o
	$(CC) $(CFLAGS) -I$(INCLUDE_PATH) -o retrieve_isoline_g retrieve_isoline_g.cc \
		read_amsu.o amsu_util.o global_metric.o $(LDFLAGS)

index_amsu: index_amsu.cc read_amsu.o amsu_util.o global_metric.o $(INCLUDE_PATH)/linked.cc
	$(CC) $(CFLAGS2) -o index_amsu index_amsu.cc read_amsu.o amsu_util.o \
		global_metric.o $(LDFLAGS)

convert_gome: read_gome.o convert_gome.cc
	$(CC) $(CFLAGS) -o $@ convert_gome.cc read_gome.o $(LDFLAGS)

print_gome_raw: read_gome.o print_gome_raw.cc
	$(CC) $(CFLAGS) -o $@ print_gome_raw.cc read_gome.o $(LDFLAGS)

extract_gome_colloc2: read_gome.o extract_gome_colloc2.cc read_spec_file.cc
	$(CC) $(CFLAGS) -o $@ extract_gome_colloc2.cc read_gome.o $(LDFLAGS)

convert_gome_colloc: read_gome.o convert_gome_colloc.cc
	$(CC) $(CFLAGS) -o $@ convert_gome_colloc.cc read_gome.o $(LDFLAGS)

o3_class: o3_class.cc
	$(CC) $(CFLAGS) -o $@ o3_class.cc $(LDFLAGS)

index_gome: index_gome.cc read_gome.o global_metric.o $(INCLUDE_PATH)/linked.cc
	$(CC) $(CFLAGS2) -o index_gome index_gome.cc read_gome.o global_metric.o $(LDFLAGS)

interpolate_amsu: interpolate_amsu.cc read_amsu.o amsu_interp_obj.o amsu_interp_base.o \
		amsu_util.o global_metric.o amsu_tinterp_obj.o
	$(CC) $(CFLAGS) -o interpolate_amsu interpolate_amsu.cc \
		amsu_interp_obj.o read_amsu.o amsu_util.o global_metric.o \
			amsu_tinterp_obj.o amsu_interp_base.o \
			$(TRAJ_PATH)/traj_int_obj.o $(TRAJ_PATH)/tcoord_defs.o \
			$(LDFLAGS)

print_amsu_1c_file:  print_amsu_1c_file.cc read_amsu.o 
	$(CC) $(CFLAGS) -o print_amsu_1c_file print_amsu_1c_file.cc read_amsu.o $(LDFLAGS)

fswap_endian:  fswap_endian.cc 
	$(CC) $(CFLAGS) -o fswap_endian fswap_endian.cc

read_gome.o: read_gome.cc read_gome.h
	$(CC) $(CFLAGS) -c read_gome.cc

read_amsu.o: read_amsu.cc read_amsu.h
	$(CC) $(CFLAGS) -c read_amsu.cc

amsu_util.o: amsu_util.cc read_amsu.h
	$(CC) $(CFLAGS) -c amsu_util.cc

amsu_tinterp_obj.o: amsu_tinterp_obj.cc amsu_tinterp_obj.h 
	$(CC) $(CFLAGS) -c amsu_tinterp_obj.cc 

amsu_interp_obj.o: amsu_interp_obj.cc amsu_interp_obj.h amsu_interp_base.h \
			read_amsu.h global_metric.h
	$(CC) $(CFLAGS) -c amsu_interp_obj.cc

amsu_interp_base.o: amsu_interp_base.cc amsu_interp_base.h read_amsu.h
	$(CC) $(CFLAGS) -c amsu_interp_base.cc 

tcolloc_obj.o: tcolloc_obj.cc tcolloc_obj.h 
	$(CC) $(CFLAGS) -c tcolloc_obj.cc 

colloc_obj.o: colloc_obj.cc colloc_obj.h colloc_base.h \
			read_gome.h global_metric.h
	$(CC) $(CFLAGS) -c colloc_obj.cc

colloc_base.o: colloc_base.cc colloc_base.h read_gome.h
	$(CC) $(CFLAGS) -c colloc_base.cc 

global_metric.o: global_metric.cc global_metric.h 
	$(CC) $(CFLAGS) -c global_metric.cc

$(TRAJ_PATH)/traj_int_obj.o: $(TRAJ_PATH)/traj_int_obj.cc $(TRAJ_PATH)/traj_int_obj.h
	make -C $(TRAJ_PATH)

$(TRAJ_PATH)/tcoord_defs.o: $(TRAJ_PATH)/tcoord_defs.cc $(TRAJ_PATH)/tcoord_defs.h
	make -C $(TRAJ_PATH)

