
program main
  
  use omp_lib
  use module1, only: dp, print_matrix, trace, dagger, superH, superG,kronecker, expm,&
& liouvillian, delta_rho, r8_normal_01, delta_rho_homodyne, superD, homodyne_detection,&
& integrate_current, exact_solution, factorial,reconstruct_state,correlations

 
 
  implicit none

  complex(dp), dimension(2,2) :: sigma_z, sigma_plus, sigma_minus, rhozero, H, H_traj
  complex(dp) :: i = complex(0,1),rhovec_zero(4)
 
  real(dp), parameter ::  dt = 1e-3
  real(dp) ::  gamma, Omega, w, start, finish, ompstart, ompend, const, t_start = 10
 logical :: shift_current
  integer ::  channels, nruns, ntrajs,nangles
logical :: milstein = .true.
!logical :: milstein = .false.


  sigma_plus = reshape ( (/ 0,0,1,0/),(/2,2/) )
  sigma_minus = reshape ( (/ 0,1,0,0/),(/2,2/) )
  sigma_z = reshape ( (/ 1,0,0,-1/),(/2,2/) )
  rhozero = reshape ( (/1,0,0,0/),(/2,2/) )
  rhovec_zero = (/1,0,0,0/)
  
  nruns = 20000
  ntrajs = 100
  nangles= 1

   gamma = 1
   Omega = 0
   w = 1
   channels = 1
   shift_current = .false.

   if (channels == 1) then
       const = gamma
    else
       const = 0.5_dp*gamma
    end if

   

   ! resonant fluorescence
   H =  - i*gamma*matmul(sigma_plus,sigma_minus)/2 -i*sqrt(const)*Omega*(sigma_plus-sigma_minus)
  ! H_traj =  w*sigma_z/2 

 H_traj = -i*sqrt(const)*Omega*(sigma_plus-sigma_minus) !Omega*(sigma_plus + sigma_minus)/2

call cpu_time(start)

call homodyne_detection(nruns,ntrajs,dt,rhozero,sigma_minus,sigma_plus,gamma,'traj1.dat',H_traj,milstein,channels)

call cpu_time(finish)

print '("Time = ",f10.3," seconds for trajectory solution.")',finish-start

ompstart = omp_get_wtime()
 
!call reconstruct_state(Omega,nruns,ntrajs,dt,rhozero,sigma_minus,sigma_plus,H_traj,gamma,&
!&t_start,milstein,channels, 'rho_im.dat', 'rho_re_.dat',nangles,shift_current)

ompend = omp_get_wtime()
print *, char(10)
print '("Time = ",f10.3," seconds for reconstruction.")',ompend-ompstart


call cpu_time(start)
!call correlations(nruns,rhozero,sigma_minus,sigma_plus,dt,'g2_Omega0_5.dat',t_start,rhovec_zero, H)
call exact_solution(nruns,sigma_minus,sigma_plus,dt,'exact.dat',t_start,rhovec_zero, H)

call cpu_time(finish)
 print '("Time = ",f6.3," seconds for exact solution.")',finish-start

end program main

