FC      = gfortran
FFLAGS  = -Wall -Wextra -march=native -O3 -fbounds-check
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LDFLAGS = $(shell pkg-config --libs plplotd-f95)
LIBS    =

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = plot.o forces.o md_sim.o

md_sim.exe: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o:%.f95
	$(COMPILE) -c $<
