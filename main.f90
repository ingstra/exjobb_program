
program main

  use module1, only: dp, print_matrix, trace, dagger, superH, identity_matrix, superG, kronecker, expm, liouvillian, delta_rho
 
  implicit none

  complex(dp), dimension(2,2) :: sigma_z, sigma_plus, sigma_minus, tm, rhozero, H, identity, rho, H_traj
  complex(dp) :: i = complex(0,1), rhovec(4),rhovec_zero(4)
  real(dp), allocatable :: store(:), time(:)
  real(dp), parameter ::  dt = 1e-3
  real(dp) :: random, gamma, Omega, w
  integer, allocatable :: seed(:)
  integer, parameter :: file_unit1 = 1, file_unit2 = 2
  integer :: values(8), k, j, nruns = 10000, dN, ntrajs = 10
  
  allocate(store(nruns),time(nruns))

  sigma_plus = reshape ( (/ 0,0,1,0/),(/2,2/) )
  sigma_minus = reshape ( (/ 0,1,0,0/),(/2,2/) )
  sigma_z = reshape ( (/ 1,0,0,-1/),(/2,2/) )
  rhozero = reshape ( (/1,0,0,0/),(/2,2/) )
  identity = reshape ( (/1,0,0,1/),(/2,2/) )
  rhovec_zero = (/1,0,0,0/)

   gamma = 1
   Omega = 5
   w = 1

   ! resonant fluorescence
   !H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 + Omega*(sigma_plus + sigma_minus)/2
   
   H = w*sigma_z/2  - i*gamma*matmul(sigma_plus,sigma_minus)/2

   H_traj = w*sigma_z/2 

  time = 0


call date_and_time(VALUES=values)
call random_seed(size=k)
allocate(seed(k))
seed(:) = values
call random_seed(put=seed)

store=0

!open(unit=3, file='trace.dat', action="write")

do k=1,ntrajs
   rho = rhozero
   do j=1,nruns

!      write(3,'(E22.7,A1,E22.7)') j*dt, char(9), trace(rho)

      store(j) = store(j)+ trace(matmul(sigma_plus,sigma_minus)*rho )

      call random_number(random)

      if (random < dt ) then
         dN = 1
      else
         dN = 0
      end if

      rho = rho + delta_rho(rho,H_traj,sigma_minus,sigma_plus,dt,dN)
     ! call print_matrix(delta_rho(rho,H_traj,sigma_minus,sigma_plus,dt,dN))

   end do
end do
!close(3)

open(unit=1, file='traj.dat', action="write")
 do j=1,nruns
    write(1,'(E22.7,A1,E22.7)') j*dt, char(9), store(j)/ntrajs
 end do

close(1)


open(unit=2, file='exact.dat', action="write")

do j=1,nruns
rhovec = matmul( expm(j*dt,liouvillian(sigma_minus,sigma_plus,H)), rhovec_zero)
rho(1,1) = rhovec(1)
rho(2,1) = rhovec(2)
rho(1,2) = rhovec(3)
rho(2,2) = rhovec(4)
write(2,'(E22.7,A1,E22.7)') j*dt, char(9), real(trace(matmul(sigma_plus,sigma_minus)*rho ))
end do

close(2)




end program main

