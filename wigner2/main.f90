
program main
  
  use omp_lib
  use module1, only: dp, print_matrix, trace, dagger, superH, superG,kronecker, expm,&
& liouvillian,exact_solution, factorial,correlations

 
 
  implicit none

  complex(dp), dimension(2,2) :: sigma_z, sigma_plus, sigma_minus, rhozero, H, H_traj
  complex(dp) :: i = complex(0,1),rhovec_zero(4)
 
  real(dp), parameter ::  dt = 1e-3
  real(dp) ::  gamma, Omega, start, finish, t_start = 10
  integer ::  channels, nruns


  sigma_plus = reshape ( (/ 0,0,1,0/),(/2,2/) )
  sigma_minus = reshape ( (/ 0,1,0,0/),(/2,2/) )
  sigma_z = reshape ( (/ 1,0,0,-1/),(/2,2/) )
  rhozero = reshape ( (/1,0,0,0/),(/2,2/) )
  rhovec_zero = (/1,0,0,0/)
  
  nruns = 20000

   gamma = 1
   Omega = 0.5


   

   ! resonant fluorescence
   H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 -i*sqrt(gamma)*Omega*(sigma_plus-sigma_minus)
  ! H_traj =  w*sigma_z/2 

! H_traj = -i*sqrt(const)*Omega*(sigma_plus-sigma_minus) !Omega*(sigma_plus + sigma_minus)/2


call cpu_time(start)
call correlations(nruns,rhozero,sigma_minus,sigma_plus,dt,'g2.dat',t_start,rhovec_zero, H)
call exact_solution(nruns,sigma_minus,sigma_plus,dt,'exact.dat',rhovec_zero, H)

call cpu_time(finish)
 print '("Time = ",f6.3," seconds for exact solution.")',finish-start

end program main

