# files to compile and where to put them
srcDir=./src
objDir=./obj
modDir=./mod
files=$(srcDir)/md.f90
objFiles=
objects=$(addprefix $(objDir)/,$(objFiles))
libs=

comp90=gfortran
comp77=gfortran
debugfl=-g -fbacktrace -fbounds-check
opt=-O3 -Wall -cpp
modd=-J $(modDir)
f90flags=$(debugfl) $(opt) -c $(modd)
f77flags=$(debugfl) $(opt) -c $(modd)
lflags=$(debugfl) $(opt)

md.exe: $(files) $(objects)
	$(comp90) $(lflags) -o md.exe $(modd) $(files) $(objects) $(libs)
$(objDir)/%.o: $(srcDir)/%.f90
	$(comp90) $(f90flags) $< -o $@ $(libs)
$(objDir)/%.o: $(srcDir)/%.f
	$(comp77) $(f77flags) $< -o $@ $(libs)

clean:
	rm -f $(objDir)/*.o
	rm -f $(modDir)/*.mod
	rm -f fort.*
	rm -f  *.exe
	rm -f *.out
