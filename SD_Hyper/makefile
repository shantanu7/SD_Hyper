FCOMP = ifort
OPTS  = -c -fpe0 -traceback -CB -O3 -g
LINKOPTS = -O3 -o
OBJS =  Modules.o main.o read_input.o connectivity.o free_memory.o \
		disp_messages.o initsetup.o mapder.o update_mesh_metrics.o \
		calculate_jacobian.o getXXSolu.o write_tec_data.o \
		compresid.o compflux.o interface_flux.o BC_flux.o \
		iterations.o

main: $(OBJS) 
	  $(FCOMP) $(LINKOPTS) sd_hyper $(OBJS)

clean:
	rm *.o *.mod
	rm sd_hyper

.SUFFIXES: .o .F90 .f90

.f90.o:
	$(FCOMP) $(OPTS) $(DEFINES) $<

.F90.o:
	$(FCOMP) $(OPTS) $(DEFINES) $<
