
! tests floating point accuracy

real(8) a,b,c,d

real(8) e

a=5000.d0
read(*,*) e
b=a*(1-e)

write(*,*) a-b

end
