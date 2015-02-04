program test
implicit none

real::r(3), d = 1.0

r = (/ .1, .5, .7 /)
print *, r - nint(r/d)*d

end program
