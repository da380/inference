MOD=../mod/
LOBJ=../obj/
BIN=../../bin/
MV=mv -f
F95 = gfortran

# linked files from util directory
LUOBJ=../../util/obj
LUMOD=../../util/mod
UOBJ =  $(LUOBJ)/nrtype.o \
        $(LUOBJ)/nrutil.o \
        $(LUOBJ)/module_util.o \
        $(LUOBJ)/module_GLL.o  \
	$(LUOBJ)/module_ODE.o  \
	$(LUOBJ)/module_fourier.o \
	$(LUOBJ)/ran_state.o \
	$(LUOBJ)/module_maps.o \
        $(LUOBJ)/module_interp.o \
	$(LUOBJ)/module_gsht.o \
	$(LUOBJ)/module_root.o   \
	$(LUOBJ)/module_function.o \
	$(LUOBJ)/rotmx2.o          \
	$(LUOBJ)/prott.o           \
	$(LUOBJ)/module_opt.o      \
	$(LUOBJ)/module_crust.o




FFLAGS = -I../mod -I$(LUMOD) -J../mod  -O5 


SRCS =  module_mesh.f90 \
        module_eig.f90  \
	module_S2_point.f90


OBJS =  $(LOBJ)module_mesh.o 	\
	$(LOBJ)module_eig.o     \
	$(LOBJ)module_S2_point.o


MODS =  $(MOD)module_mesh.mod  \
	$(MOD)module_eig.mod   \
	$(MOD)module_S2_point.mod

all:	sobolev_disk   \
	delta_test     \
	S2_point_bound \
	S2_point_data    \
	S2_point_mnne

clean:
	rm ../obj/*
	rm ../mod/*
	rm $(BIN)sobolev_disk
	rm $(BIN)delta_test
	rm $(BIN)S2_point_bound
	rm $(BIN)S2_point_data
	rm $(BIN)S2_point_mnne



sobolev_disk: $(BIN)sobolev_disk $(OBJS) $(UOBJ) $(RAYOBJ)
$(BIN)sobolev_disk: $(OBJS) $(LOBJ)sobolev_disk.o
	$(F95) $(FFLAGS) $(LOBJ)sobolev_disk.o $(OBJS) \
	$(UOBJ) -l lapack  -o $(BIN)sobolev_disk

$(LOBJ)sobolev_disk.o: sobolev_disk.f90 
	$(F95) $(FFLAGS) sobolev_disk.f90 -c
	$(MV) sobolev_disk.o $(LOBJ)


delta_test: $(BIN)delta_test $(OBJS) $(UOBJ) $(RAYOBJ)
$(BIN)delta_test: $(OBJS) $(LOBJ)delta_test.o
	$(F95) $(FFLAGS) $(LOBJ)delta_test.o $(OBJS) \
	$(UOBJ) -l lapack  -o $(BIN)delta_test

$(LOBJ)delta_test.o: delta_test.f90 
	$(F95) $(FFLAGS) delta_test.f90 -c
	$(MV) delta_test.o $(LOBJ) 


S2_point_bound: $(BIN)S2_point_bound $(OBJS) $(UOBJ) $(RAYOBJ)
$(BIN)S2_point_bound: $(OBJS) $(LOBJ)S2_point_bound.o
	$(F95) $(FFLAGS) $(LOBJ)S2_point_bound.o $(OBJS) \
	$(UOBJ) -l lapack  -o $(BIN)S2_point_bound

$(LOBJ)S2_point_bound.o: S2_point_bound.f90 
	$(F95) $(FFLAGS) S2_point_bound.f90 -c
	$(MV) S2_point_bound.o $(LOBJ)


S2_point_data: $(BIN)S2_point_data $(OBJS) $(UOBJ) $(RAYOBJ)
$(BIN)S2_point_data: $(OBJS) $(LOBJ)S2_point_data.o
	$(F95) $(FFLAGS) $(LOBJ)S2_point_data.o $(OBJS) \
	$(UOBJ) -l lapack  -o $(BIN)S2_point_data

$(LOBJ)S2_point_data.o: S2_point_data.f90 
	$(F95) $(FFLAGS) S2_point_data.f90 -c
	$(MV) S2_point_data.o $(LOBJ)


S2_point_mnne: $(BIN)S2_point_mnne $(OBJS) $(UOBJ) $(RAYOBJ)
$(BIN)S2_point_mnne: $(OBJS) $(LOBJ)S2_point_mnne.o
	$(F95) $(FFLAGS) $(LOBJ)S2_point_mnne.o $(OBJS) \
	$(UOBJ) -l lapack  -o $(BIN)S2_point_mnne

$(LOBJ)S2_point_mnne.o: S2_point_mnne.f90 
	$(F95) $(FFLAGS) S2_point_mnne.f90 -c
	$(MV) S2_point_mnne.o $(LOBJ)


$(OBJS): $(SRCS) $(UOBJ)	
	$(F95)   $(FFLAGS) -c $(SRCS)
	$(MV)  *.o $(LOBJ) 


