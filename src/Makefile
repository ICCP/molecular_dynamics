FC = gfortran
FFLAGS = -Wall -Wextra -O3 -march=native -ffast-math #-std=f95
LDFLAGS =
#FFLAGS += $(shell pkg-config --cflags plplotd-f95)
#LDFLAGS = $(shell pkg-config --libs plplotd-f95)
#LIBS = -llapack

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS =
#OBJS += md_plot.o
OBJS += initialization.o
OBJS += force.o
OBJS += energy.o
OBJS += const_temp.o
OBJS += myprog.o


all: myprog
LIBS=-L$(HOME)/XPS/lib -lxps -lX11
myprog: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

#md_plot.o: md_plot.f95
#	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) myprog $(OBJS) *.mod
