
program main
  
  use omp_lib
  use module1, only: dp, print_matrix, trace, dagger, superH, superG,kronecker, expm,&
& liouvillian, delta_rho, r8_normal_01, delta_rho_homodyne, superD, homodyne_detection,&
& integrate_current, exact_solution, factorial,reconstruct_state

 
 
  implicit none

  complex(dp), dimension(2,2) :: sigma_z, sigma_plus, sigma_minus, rhozero, H, H_traj
  complex(dp) :: i = complex(0,1),rhovec_zero(4)
 
  real(dp), parameter ::  dt = 1e-3
  real(dp) ::  gamma, Omega, w, start, finish, ompstart, ompend, const
 
  integer ::  channels, nruns = 13000, ntrajs = 500
!logical :: milstein = .true.
logical :: milstein = .false.

  sigma_plus = reshape ( (/ 0,0,1,0/),(/2,2/) )
  sigma_minus = reshape ( (/ 0,1,0,0/),(/2,2/) )
  sigma_z = reshape ( (/ 1,0,0,-1/),(/2,2/) )
  rhozero = reshape ( (/0,0,0,1/),(/2,2/) )
  rhovec_zero = (/0,0,0,1/)

   gamma = 1
   Omega = 0.5
   w = 1
   channels = 1

   if (channels == 1) then
       const = gamma
    else
       const = 0.5_dp*gamma
    end if

   

   ! resonant fluorescence
   H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 + Omega*(sigma_plus + sigma_minus)/2
   
  ! H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 -i*sqrt(0.5_dp*gamma)*Omega*(sigma_plus-sigma_minus) ! + Omega*(sigma_plus + sigma_minus)/2

  ! H_traj =  w*sigma_z/2 

 H_traj = -i*sqrt(const)*Omega*(sigma_plus-sigma_minus) !Omega*(sigma_plus + sigma_minus)/2

call cpu_time(start)

!call homodyne_detection(nruns,ntrajs,dt,rhozero,sigma_minus,sigma_plus,gamma,'traj1.dat',H_traj,milstein)

call cpu_time(finish)

print '("Time = ",f10.3," seconds for trajectory solution.")',finish-start

ompstart = omp_get_wtime()
 
call reconstruct_state(nruns,ntrajs,dt,rhozero,sigma_minus,sigma_plus,H_traj,gamma,milstein,channels)

ompend = omp_get_wtime()
print *, char(10)
print '("Time = ",f10.3," seconds for reconstruction.")',ompend-ompstart


call cpu_time(start)
!call exact_solution(nruns,sigma_minus,sigma_plus,dt,'exact.dat',rhovec_zero, H)

call cpu_time(finish)
 print '("Time = ",f6.3," seconds for exact solution.")',finish-start

end program main

