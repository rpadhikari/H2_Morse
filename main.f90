program wavefunction
  implicit none
  integer(4) i
  real(8) De, a, xe, x, x_min, x_max, mass, hg
  integer(4) nSteps, n
  real(4) lambda, alpha, z, h, Nn, Cn
  real(8), external :: fact
  real(8), parameter :: hbar=1.054571800130d-34 ! J-S
  real(8), parameter :: Ese=1.602176620898d-19 ! electron charge C
  real(8), parameter :: Lse=1.0d-09 ! Nanometer to meter factor
  open(1,file='in.dat', action='read') 
  open(2,file='psi.dat', action='write')
  read(1,*) De  ! in the unit of eV
  read(1,*) a   ! in the unit of per nm
  read(1,*) xe  ! in the unit of nm
  read(1,*) x_min  ! in the unit of nm
  read(1,*) x_max  ! in the unit of nm
  read(1,*) nSteps ! segemnts between x_min and x_max
  read(1,*) mass   ! mass in kg
  read(1,*) n      ! nth eigenfunction, starts from 0 (zero)
  close(1)
  lambda=dsqrt(2.0d0*mass*De*Ese)*a*Lse/hbar
  write(*,*) lambda
  alpha=2.0d0*(lambda-dble(n))-1.0d0
  write(*,*) alpha
  h=abs(x_max-x_min)/nSteps
  write(*,*) h
  Nn=dsqrt(fact(n)*alpha/gamma(2.0d0*lambda-dble(n)))
  write(*,*) Nn 
  Cn=gamma(alpha+dble(n+1))/(gamma(alpha+1.0d0)*gamma(dble(n+1)))
  write(*,*) Cn
  do i=0, nSteps
    x=x_min+dble(i)*h
    z=2.0d0*lambda*exp(-(x-xe))
    call chgm(-dble(n), alpha+1.0d0, z, hg)
    write(2,10) x, Nn*z**(0.50d0*alpha)*exp(-0.5d0*z)*Cn*hg
  end do
  close(2)
  10 format(f10.5, e24.15)
end program wavefunction
 
  
  
