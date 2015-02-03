FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3
LDFLAGS =
#LIBS = -llapack

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS =
OBJS += myprog.o

all: myprog

myprog: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) myprog $(OBJS) *.mod
