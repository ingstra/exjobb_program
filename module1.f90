module module1

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: pi = 3.141592653589793

contains

pure function factorial(n)
integer, intent(in) :: n

integer :: i
real(dp) :: factorial

factorial = PRODUCT((/(i, i=1,n)/))

end function factorial

  subroutine print_matrix(matrix,long_flag)
    implicit none

    !input
    !>matrix which you want to print
    complex(dp),intent(in) :: matrix(:,:)

    !>true value make subroutine print with longer format
    logical, intent(in), optional :: long_flag

    !locals
    !row index of matrix
    integer :: i

    !format string for output
    character(20) :: format

    if (present(long_flag).and.long_flag) then
       format='(20G20.12)'
    else
       format='(20G12.4)'
    end if

    do i=1,size(matrix,1) 
       write(*,format) real(matrix(i,:))
    end do

  end subroutine print_matrix

subroutine direct_detection(nruns,ntrajs,dt,rhozero,c,cdagger,filename,H)
  integer, intent(in) :: nruns, ntrajs
  real(dp), intent(in) :: dt
  complex(dp), intent(in) :: rhozero(:,:),c(:,:),cdagger(:,:),H(:,:)
  character(len=*), intent(in) :: filename

  real(dp), allocatable :: store(:)
  integer :: k, j, dN, seed
  complex(dp) :: rho(2,2)
  allocate(store(nruns))

    seed = 123456789
    call srand(seed)
!open(unit=3, file='trace.dat', action="write")
store=0
do k=1,ntrajs
   rho = rhozero
   do j=1,nruns
!      write(3,'(E22.7,A1,E22.7)') j*dt, char(9), trace(rho)
      store(j) = store(j)+ trace(matmul(cdagger,c)*rho )
      call random_number(random)
      if (random < dt ) then
         dN = 1
      else
         dN = 0
      end if
      rho = rho + delta_rho(rho,H,c,cdagger,dt,dN)
   end do
end do
!close(3)

open(unit=1, file=filename, action="write")
 do j=1,nruns
    write(1,'(E22.7,A1,E22.7)') j*dt, char(9), store(j)/ntrajs
 end do
close(1)

end subroutine direct_detection


subroutine homodyne_detection(nruns,ntrajs,dt,rhozero,c,cdagger,gamma,filename,H)

use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
  integer, intent(in) :: nruns, ntrajs
  real(dp), intent(in) :: dt, gamma
  complex(dp), intent(in) :: rhozero(:,:), c(:,:),cdagger(:,:), H(:,:)
  character(len=*), intent(in) :: filename
  real(dp), allocatable :: store(:)

  integer :: k,j,seed
  real(dp) :: dW
  complex(dp) :: rho(2,2)
  character(len=*), parameter :: carriage_return =  char(13)
  allocate(store(nruns))
  
  seed = 123456789

  !open(unit=3, file='trace.dat', action="write")
 
  do k=1,ntrajs
     rho = rhozero
     do j=1,nruns
       dW = r8_normal_01(seed)*sqrt(dt)
        !  write(3,'(E22.7,A1,E22.7)') j*dt, char(9), trace(rho)
             ! write(3, '(E22.7)') dW
        store(j) = store(j)+ trace(matmul(cdagger,c)*rho )
        rho = rho + delta_rho_homodyne(rho,H,c,dt,dW,gamma)
     end do
     write(stdout,"(2a,i10,$)") carriage_return,"Calculating trajectory: ", k
  end do
  write(stdout,*) linefeed
 ! close(3)

  open(unit=1, file=filename, action="write")
  do j=1,nruns
     write(1,'(E22.7,A1,E22.7)') j*dt, char(9), store(j)/ntrajs
  end do
  close(1)

end subroutine homodyne_detection

subroutine integrate_photocurrent(nruns,ntrajs,dt,rhozero,c,cdagger,filename,H,gamma)

use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
use hermite, only: evalHermitePoly

  integer, intent(in) :: nruns,ntrajs
  real(dp), intent(in) :: dt, gamma
  complex(dp), intent(in) :: rhozero(:,:),c(:,:),cdagger(:,:),H(:,:)
  character(len=*), intent(in) :: filename

  character(len=*), parameter :: carriage_return =  char(13)
  integer ::  seed, n=1
  real(dp) :: dW, current, t, t_tot, const,env_test
  real(dp), allocatable :: h_poly(:),current_store(:)
  
  complex(dp) :: rho(2,2)

  allocate(h_poly(ntrajs),current_store(ntrajs))
  const = 1._dp/(sqrt(2._dp**n * factorial(n)) * (2._dp*pi)**0.25_dp )
  seed = 123456789 

  write(stdout,*) 'Time to integrate: ', nruns*dt
  open(unit=4, file=filename, action="write")
  t_tot = nruns*dt ! total time

do k=1,ntrajs
   rho = rhozero
   current = 0
   t = 0
env_test = 0
   do j=1,nruns
      dW = r8_normal_01(seed)*sqrt(dt)
      current = current + (sqrt(gamma)*trace(matmul(cdagger+c,rho))*dt + dW)*envelope(gamma,t_tot,t)!/sqrt(t_tot)
      rho = rho + delta_rho_homodyne(rho,H,c,dt,dW,gamma)
      t = t + dt
      env_test = env_test + dt*envelope(gamma,t_tot,t)**2._dp
   end do
   write(stdout,"(2a,i10,$)") carriage_return,"Integrating current. Trajectory: ", k
   write(4,'(E22.7)') current
   current_store(k) = current
!print *, env_test
end do
close(4)

h_poly = evalHermitePoly(current_store/sqrt(2d0),n)
open(unit=5, file='hermite.dat', action="write")

do j=1,ntrajs
   write(5,'(E22.7,A1,E22.7)') current_store(j), char(9), &
&const*exp(-current_store(j)**2._dp/4._dp)*h_poly(j)
end do

close(5)

end subroutine integrate_photocurrent

 function envelope(gamma,t_tot,t)
  real(dp), intent(in) :: gamma, t_tot, t
  real(dp) :: envelope, N

  ! Normalization factor
  N = sqrt(gamma/(1._dp-exp(-gamma*t_tot)))

  envelope = N*exp(-0.5_dp*gamma*t)
end function envelope

subroutine exact_solution(nruns, c, cdagger, dt ,filename, rhovec_zero, H)
 integer, intent(in) :: nruns
  real(dp), intent(in) :: dt
  complex(dp), intent(in) :: rhovec_zero(:),c(:,:),cdagger(:,:),H(:,:)
  character(len=*), intent(in) :: filename

  complex(dp) :: rho(2,2), rhovec(4)

open(unit=2, file=filename, action="write")
do j=1,nruns
rhovec = matmul( expm(j*dt,liouvillian(c,H)), rhovec_zero)
rho(1,1) = rhovec(1)
rho(2,1) = rhovec(2)
rho(1,2) = rhovec(3)
rho(2,2) = rhovec(4)
write(2,'(E22.7,A1,E22.7)') j*dt, char(9), real(trace(matmul(cdagger,c)*rho ))
end do
close(2)

end subroutine exact_solution


  pure function trace(matrix)
    complex(dp), intent(in) :: matrix(:,:)
    real(dp) :: trace
    integer :: n 
    n = size(matrix,1)
    trace = 0

    do i=1,n
       trace = trace + real(matrix(i,i))
    enddo
  end function trace

  function superH(r,rho)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::superH(:,:) 

    integer :: n
    n=size(r,1)

    allocate(superH(n,n))

    superH = matmul(r,rho) + matmul(rho,dagger(r)) - trace( matmul(r,rho) + matmul(rho,dagger(r)) )*rho

  end function superH

  pure function superG(r,rho)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::superG(:,:) 
    real(dp) :: tracevalue

    integer :: n
    n=size(r,1)

    allocate(superG(n,n))

    tracevalue = real( trace( matmul(matmul(r,rho),dagger(r)) ) )

    if (tracevalue /= 0) then
       superG = matmul(matmul(r,rho),dagger(r))/tracevalue - rho
    else
       superG = 0!- rho
    end if

  end function superG

  pure function superD(c,rho)

    complex(dp), intent(in) :: c(:,:), rho(:,:)
    complex(dp), allocatable ::superD(:,:) 

    complex(dp), allocatable :: cdagger(:,:)
    integer :: n
    n=size(c,1)

    allocate(superD(n,n),cdagger(n,n))

    cdagger = conjg(transpose(c))

    superD = matmul(matmul(c,rho),cdagger) - matmul(matmul(cdagger,c),rho)/2 - matmul(matmul(rho,cdagger),c)/2

  end function superD


  pure function dagger(matrix)
    complex(dp), intent(in) :: matrix(:,:)
    complex(dp), allocatable :: dagger(:,:)
    integer :: n 
    n = size(matrix,1)
    allocate(dagger(n,n))

    dagger = conjg(transpose(matrix))
  end function dagger

  pure function identity_matrix(n)
    integer, intent(in) :: n
    real(dp), dimension(n,n) :: tmp
    complex(dp), allocatable :: identity_matrix(:,:)
    allocate(identity_matrix(n,n))

    tmp = 0
    do i=1,n
       tmp(i,i) = 1 
    end do

    identity_matrix = tmp
  end function identity_matrix

  function delta_rho(rho,H,c,cdagger,delta_t,dN)
    complex(dp), intent(in) :: rho(:,:), H(:,:), c(:,:), cdagger(:,:)
    real(dp), intent(in) :: delta_t
    integer, intent(in) :: dN
    complex(dp), allocatable ::delta_rho(:,:) 

    complex(dp) :: i 
    integer :: n
    i = complex(0,1)
    n = size(rho,1)
    allocate(delta_rho(n,n))

    delta_rho = dN*superG(c,rho)-delta_t*superH(i*H + matmul(cdagger,c)/2,rho)

  end function delta_rho

  function delta_rho_homodyne(rho,H,c,dt,dW,gamma)
    complex(dp), intent(in) :: rho(:,:), H(:,:), c(:,:)
    real(dp), intent(in) :: dt, dW, gamma
    complex(dp), allocatable ::delta_rho_homodyne(:,:) 

    complex(dp) :: i 
    integer :: n
    i = complex(0,1)
    n = size(rho,1)
    allocate(delta_rho_homodyne(n,n))

    delta_rho_homodyne = -i*(matmul(H,rho)-matmul(rho,H))*dt + gamma*superD(c,rho)*dt + sqrt(gamma)*superH(c,rho)*dW

  end function delta_rho_homodyne

  function liouvillian(c, H_eff) 
    complex(dp), intent(in) :: H_eff(:,:), c(:,:)
    complex(dp) :: i = complex(0,1)
    complex(dp), dimension(2,2) :: identity
    complex(dp), allocatable :: liouvillian(:,:)
    integer :: n 
    n = size(c,1)
    allocate(liouvillian(n**2,n**2))

    identity = reshape ( (/1,0,0,1/),(/2,2/) )

    liouvillian = -i*kronecker(identity,H_eff) + i*kronecker(conjg(H_eff),identity) +&
         &kronecker(conjg(c),c)

  end function liouvillian

  function kronecker(A,B) 

    complex(dp), dimension (:,:), intent(in)  :: A, B
    complex(dp), dimension (:,:), allocatable :: kronecker
    integer :: i = 0, j = 0, m = 0, n = 0, p = 0, q = 0

    allocate(kronecker(size(A,1)*size(B,1),size(A,2)*size(B,2)))
    kronecker = 0

    do i = 1,size(A,1)
       do j = 1,size(A,2)
          n=(i-1)*size(B,1) + 1
          m=n+size(B,1)
          p=(j-1)*size(B,2) + 1
          q=p+size(B,2) 
          kronecker(n:m,p:q) = A(i,j)*B
       enddo
    enddo

  end function kronecker


  function expm(t,H) result(expH)
    real(dp), intent(in) :: t
    complex(dp), dimension(:,:), intent(in) :: H
    complex(dp), dimension(size(H,1),size(H,2)) :: expH

    ! Expokit variables
    external :: DGPADM
    integer, parameter :: ideg = 6
    complex(dp), dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
    integer, dimension(size(H,1))  :: iwsp
    integer :: iexp, ns, iflag, n

    if (size(H,1) /= size(H,2)) then
       stop 'expm: matrix must be square'
    end if

    n = size(H,1)
    call ZGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))
  end function expm


  subroutine ZGPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

    implicit none
    real(dp) :: t
    integer ::         ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
    complex(dp) ::       H(ldh,m), wsp(lwsp)

    !-----Purpose----------------------------------------------------------|
    !
    !     Computes exp(t*H), the matrix exponential of a general complex 
    !     matrix in full, using the irreducible rational Pade approximation
    !     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
    !     combined with scaling-and-squaring.
    !
    !-----Arguments--------------------------------------------------------|
    !
    !     ideg      : (input) the degre of the diagonal Pade to be used.
    !                 a value of 6 is generally satisfactory.
    !
    !     m         : (input) order of H.
    !
    !     H(ldh,m)  : (input) argument matrix.
    !
    !     t         : (input) time-scale (can be < 0).
    !                  
    !     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
    !
    !     ipiv(m)   : (workspace)
    !
    !>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
    !                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
    !                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !                 NOTE: if the routine was called with wsp(iptr), 
    !                       then exp(tH) will start at wsp(iptr+iexph-1).
    !
    !     ns        : (output) number of scaling-squaring used.
    !
    !     iflag     : (output) exit flag.
    !                       0 - no problem
    !                      <0 - problem
    !
    !----------------------------------------------------------------------|
    !     Roger B. Sidje (rbs@maths.uq.edu.au)
    !     EXPOKIT: Software Package for Computing Matrix Exponentials.
    !     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
    !----------------------------------------------------------------------|
    !
    integer :: i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
    real(dp) :: hnorm
    complex(dp) :: cp, cq, scale, scale2, ZERO, ONE

    parameter ( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX

    !---  check restrictions on input parameters ...
    mm = m*m
    iflag = 0
    if ( ldh.lt.m ) iflag = -1
    if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
    if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGPADM)'
    !
    !---  initialise pointers ...
    !
    icoef = 1
    ih2 = icoef + (ideg+1)
    ip  = ih2 + mm
    iq  = ip + mm
    ifree = iq + mm
    !
    !---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
    !     and set scale = t/2^ns ...
    !
    do i = 1,m
       wsp(i) = ZERO
    enddo
    do j = 1,m
       do i = 1,m
          wsp(i) = wsp(i) + ABS( H(i,j) )
       enddo
    enddo
    hnorm = 0.0d0
    do i = 1,m
       hnorm = MAX( hnorm,DBLE(wsp(i)) )
    enddo
    hnorm = ABS( t*hnorm )
    if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of ZGPADM.'
    ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
    scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
    scale2 = scale*scale
    !
    !---  compute Pade coefficients ...
    !
    i = ideg+1
    j = 2*ideg+1
    wsp(icoef) = ONE
    do k = 1,ideg
       wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
    enddo
    !
    !---  H2 = scale2*H*H ...
    !
    call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
    !
    !---  initialise p (numerator) and q (denominator) ...
    !
    cp = wsp(icoef+ideg-1)
    cq = wsp(icoef+ideg)
    do j = 1,m
       do i = 1,m
          wsp(ip + (j-1)*m + i-1) = ZERO
          wsp(iq + (j-1)*m + i-1) = ZERO
       enddo
       wsp(ip + (j-1)*(m+1)) = cp
       wsp(iq + (j-1)*(m+1)) = cq
    enddo
    !
    !---  Apply Horner rule ...
    !
    iodd = 1
    k = ideg - 1
100 continue
    iused = iodd*iq + (1-iodd)*ip
    call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,wsp(ih2),m, ZERO,wsp(ifree),m )
    do j = 1,m
       wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
    enddo
    ip = (1-iodd)*ifree + iodd*ip
    iq = iodd*ifree + (1-iodd)*iq
    ifree = iused
    iodd = 1-iodd
    k = k-1
    if ( k.gt.0 )  goto 100
    !
    !---  Obtain (+/-)(I + 2*(p\q)) ...
    !
    if ( iodd.ne.0 ) then
       call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m, H,ldh, ZERO,wsp(ifree),m )
       iq = ifree
    else
       call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,H,ldh, ZERO,wsp(ifree),m )
       ip = ifree
    endif
    call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
    call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
    if ( iflag.ne.0 ) stop 'Problem in ZGESV (within ZGPADM)'
    call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
    do j = 1,m
       wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
    enddo
    iput = ip
    if ( ns.eq.0 .and. iodd.ne.0 ) then
       call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
       goto 200
    endif
    !
    !--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
    !
    iodd = 1
    do k = 1,ns
       iget = iodd*ip + (1-iodd)*iq
       iput = (1-iodd)*ip + iodd*iq
       call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,ZERO,wsp(iput),m )
       iodd = 1-iodd
    enddo
200 continue
    iexph = iput

  END subroutine ZGPADM

 function r8_uniform_01 ( seed )

    !*****************************************************************************80
    !
    !! R8_UNIFORM_01 returns a unit pseudorandom R8.
    !
    !  Discussion:
    !
    !    This routine implements the recursion
    !
    !      seed = 16807 * seed mod ( 2^31 - 1 )
    !      r8_uniform_01 = seed / ( 2^31 - 1 )
    !
    !    The integer arithmetic never requires more than 32 bits,
    !    including a sign bit.
    !
    !    If the initial seed is 12345, then the first three computations are
    !
    !      Input     Output      R8_UNIFORM_01
    !      SEED      SEED
    !
    !         12345   207482415  0.096616
    !     207482415  1790989824  0.833995
    !    1790989824  2035175616  0.947702
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 May 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Second Edition,
    !    Springer, 1987,
    !    ISBN: 0387964673,
    !    LC: QA76.9.C65.B73.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, December 1986, pages 362-376.
    !
    !    Pierre L'Ecuyer,
    !    Random Number Generation,
    !    in Handbook of Simulation,
    !    edited by Jerry Banks,
    !    Wiley, 1998,
    !    ISBN: 0471134031,
    !    LC: T57.62.H37.
    !
    !    Peter Lewis, Allen Goodman, James Miller,
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, 1969, pages 136-143.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
    !    should NOT be 0.
    !    On output, SEED has been updated.
    !
    !    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
    !    strictly between 0 and 1.
    !
    implicit none

    integer, intent(inout) :: seed
    integer ::  k
    real(dp) :: r8_uniform_01
    

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
       seed = seed + 2147483647
    end if
    !
    !  Although SEED can be represented exactly as a 32 bit integer,
    !  it generally cannot be represented exactly as a 32 bit real number!
    !
    r8_uniform_01 = seed * 4.656612875D-10


  end function r8_uniform_01

  function r8_normal_01 ( seed )

    !*****************************************************************************80
    !
    !! R8_NORMAL_01 returns a unit pseudonormal R8.
    !
    !  Discussion:
    !
    !    The standard normal probability distribution function (PDF) has
    !    mean 0 and standard deviation 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 August 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, integer  SEED, a seed for the random
    !    number generator.
    !
    !    Output, real  R8_NORMAL_01, a normally distributed
    !    random value.
    !
    implicit none

    real(dp):: r1, r2, r8_normal_01
    real(dp), parameter :: r8_pi = 3.141592653589793D+00

    integer, intent(inout) :: seed

    r1 = r8_uniform_01(seed)
    r2 = r8_uniform_01(seed)
    r8_normal_01 = sqrt( - 2.0_dp * log ( r1 ) ) * cos ( 2.0_dp * r8_pi * r2 )

    !print *, r8_normal_01 , 'bajs', seed

  end function r8_normal_01

  ! From numerical recipes
  SUBROUTINE gasdev_s(harvest)

    REAL(dp), INTENT(OUT) :: harvest
   ! Returns in harvest a normally distributed deviate with zero mean and unit variance, using the intrinsic random_number function as the source of uniform deviates.
    REAL(dp) :: rsq,v1,v2
    REAL(dp), SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.

    integer:: n, i, clock,fixseed
    integer, allocatable :: seed(:)
   CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    call srand(fixseed)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /) !+ int(rand()*10)
    CALL RANDOM_SEED(PUT = seed)



    if (gaus_stored) then
     !  We have an extra deviate handy so return it, and unset the flag.
       harvest=g
       gaus_stored=.false.     
    else
       ! We don’t have an extra deviate handy, so pick two uniform numbers in the square extending from -1 to +1 in each direction,
       do
          call random_number(v1)
          call random_number(v2)
          v1 = rand()
          v2 = rand()
         ! print *, v1, v2
          v1=2.0_dp*v1-1.0_dp
          v2=2.0_dp*v2-1.0_dp
          rsq=v1**2+v2**2
          !see if they are in the unit circle,
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
      ! otherwise try again.
       rsq=sqrt(-2.0_dp*log(rsq)/rsq)
      ! Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
       harvest=v1*rsq
       g=v2*rsq
       gaus_stored=.true.
       !Set flag.
    end if
  END SUBROUTINE gasdev_s


end module module1


