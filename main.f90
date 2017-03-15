
program main

  use module1, only: dp, print_matrix, trace, dagger, superH, superG,kronecker, expm,&
& liouvillian, delta_rho, r8_normal_01, delta_rho_homodyne, superD, homodyne_detection,&
& integrate_photocurrent, exact_solution
 
  implicit none

  complex(dp), dimension(2,2) :: sigma_z, sigma_plus, sigma_minus, rhozero, H, H_traj
  complex(dp) :: i = complex(0,1),rhovec_zero(4)
  real(dp), allocatable :: store(:), time(:)
  real(dp), parameter ::  dt = 1e-3
  real(dp) ::  gamma, Omega, w, start, finish
 
  integer ::  nruns = 10000, ntrajs = 4000

  allocate(store(nruns),time(nruns))

  sigma_plus = reshape ( (/ 0,0,1,0/),(/2,2/) )
  sigma_minus = reshape ( (/ 0,1,0,0/),(/2,2/) )
  sigma_z = reshape ( (/ 1,0,0,-1/),(/2,2/) )
  rhozero = reshape ( (/0,0,0,1/),(/2,2/) )
  rhovec_zero = (/0,0,0,1/)

   gamma = 1
   Omega = 5
   w = 1

   ! resonant fluorescence
   !H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 + Omega*(sigma_plus + sigma_minus)/2
   
   H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 + Omega*(sigma_plus + sigma_minus)/2

  ! H_traj =  w*sigma_z/2 

 H_traj =  Omega*(sigma_plus + sigma_minus)/2


call cpu_time(start)

call homodyne_detection(nruns,ntrajs,dt,rhozero,sigma_minus,sigma_plus,'traj.dat',H_traj)

call cpu_time(finish)
print '("Time = ",f10.3," seconds for trajectory solution.")',finish-start

!call integrate_photocurrent(nruns,1.0_dp,dt,rhozero,sigma_minus,sigma_plus,'current.dat',H_traj)
!call integrate_photocurrent(nruns,1.0_dp,dt,rhozero,sigma_minus,sigma_plus,'current.dat',H_traj)

call cpu_time(start)
call exact_solution(nruns,sigma_minus,sigma_plus,dt,'exact.dat',rhovec_zero, H)

call cpu_time(finish)
    print '("Time = ",f6.3," seconds for exact solution.")',finish-start


end program main

