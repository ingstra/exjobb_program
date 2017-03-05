
program main

  use module1, only: dp, print_matrix, trace, dagger, superH, identity_matrix, superG, kronecker, expm, liouvillian
 
  implicit none

  complex(dp), dimension(2,2) :: sigma_z, sigma_plus, sigma_minus, tm, rhozero, H, identity, rho
  complex(dp) :: i = complex(0,1), rhovec(4),rhovec_zero(4)
  real(dp), allocatable :: store(:)
  real(dp), parameter :: w = 1, delta_t = 1e-3
  real(dp) :: random, time, gamma, Omega
  integer, allocatable :: seed(:)
  integer, parameter :: file_unit = 1
  integer :: values(8), k, j, nruns = 10000

  character(len=15) :: file_name = 'testfil'
  
  allocate(store(nruns))

  sigma_plus = reshape ( (/ double precision :: 0,0,0.5,0/),(/2,2/) )
  sigma_minus = reshape ( (/double precision :: 0,0.5,0,0/),(/2,2/) )
  sigma_z = reshape ( (/double precision :: 0.5,0,0,-0.5/),(/2,2/) )
  rhozero = reshape ( (/0,0,0,1/),(/2,2/) )
  identity = reshape ( (/1,0,0,1/),(/2,2/) )
  rhovec_zero = (/0,0,0,1/)

   gamma = 1
   Omega = 5

  H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 + Omega*(sigma_plus + sigma_minus)/2
  time = 0

!  call print_matrix(kronecker(identity,sigma_plus+sigma_minus))


call date_and_time(VALUES=values)

call random_seed(size=k)
allocate(seed(k))
seed(:) = values
call random_seed(put=seed)

!do j=1,nruns
!   call random_number(random)
!   if (random < dt ) then
      
!end do

open(unit=file_unit, file=file_name, action="write")

do j=1,nruns
rhovec = matmul( expm(j*delta_t,liouvillian(sigma_minus,sigma_plus,H)), rhovec_zero)
rho(1,1) = rhovec(1)
rho(2,1) = rhovec(2)
rho(1,2) = rhovec(3)
rho(2,2) = rhovec(4)
write(file_unit,'(E22.7,A1,E22.7)') j*delta_t, char(9), real(trace(matmul(sigma_plus,sigma_minus)*rho ))
!print *, j*delta_t, char(9), real(trace(rho*sigma_z))
end do

!open(unit=file_unit, file=file_name, action="write")
!write(file_unit,'(E22.15)') store
close(file_unit)


   !  write (*,*) 'Hello, world!'   ! This is an inline comment
end program main

