FC      = gfortran
FFLAGS  = -Wall -Wextra -march=native -O3 -fimplicit-none -fdefault-real-8
#FFLAGS += -pedantic -fbounds-check -fmax-errors=1 -g
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LDFLAGS = $(shell pkg-config --libs plplotd-f95)
LIBS    =

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

TARGET = md.exe			# Name of final executable to produce
OBJS = plot.o md.o	# List of object dependencies

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o:%.f95
	$(COMPILE) -c $<
