module module1

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: pi = 3.141592653589793

contains

pure function wavefunction(n,x,arraysize) result(y)

    use hermite, only: evalHermitePoly

    integer, intent(in) :: n, arraysize
    real(dp), intent(in) ::  x(:)

    real(dp) :: const, y(arraysize), h_poly(arraysize)
 
    const = 1._dp/(sqrt(2._dp**n * factorial(n)) * (pi)**0.25_dp )
    !print *, const
    h_poly = evalHermitePoly(x,n)

    y = const*exp(-x**2._dp/2._dp)*h_poly
    return
  end function wavefunction

  pure function factorial(n) result(y)
    integer, intent(in) :: n

    integer :: i
    real(dp) :: y

   y = PRODUCT((/(i, i=1,n)/))

  end function factorial

 subroutine print_matrix_real(matrix,long_flag)
    implicit none

    !input
    !>matrix which you want to print
    real(dp),intent(in) :: matrix(:,:)

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

  end subroutine print_matrix_real

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
       write(*,format) matrix(i,:)
    end do

  end subroutine print_matrix

  subroutine test_integrate()
    integer, parameter :: n = 1000
    real(dp), dimension(n) :: x,y,test
    real(dp),  dimension(n,n) :: f


    integer :: i,j
    real(dp) :: lim1x,lim2x,lim1y,lim2y,test2
    lim1x=0
    lim2x=2
    lim1y=0
    lim2y=3
    
    call linspace(lim1x,lim2x,n,x)
    call linspace(lim1y,lim2y,n,y)
    
    do i=1,n
      do j=1,n
         f(i,j) = x(i)**2*y(j)
      end do
   end do

   do i = 1,n
      test(i)=integrate(x,f(i,:))
   end do

   test2=integrate(y,test)
   print *, test2
   
  end subroutine test_integrate
    
  pure function integrate(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    real(dp), intent(in)  :: x(:)         !! Variable x
    real(dp), intent(in)  :: y(size(x))   !! Function y(x)
    real(dp)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
  end associate
end function

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

  subroutine correlations(nruns, rhozero,c, cdagger, dt ,filename,t_start, rhovec_zero, H)
    integer, intent(in) :: nruns
    real(dp), intent(in) :: dt,t_start
    complex(dp), intent(in) :: rhovec_zero(:),c(:,:),cdagger(:,:),H(:,:),rhozero(:,:)
    character(len=*), intent(in) :: filename

    complex(dp) :: rho(2,2), rhovec(4),P_zero(4,4), P(4,4), tmp(4,4),tmpvec(4),id(2,2),tmp2(2,2)
    real(dp) :: G1, G2, g2_norm

   
    id = identity_matrix(2)
    P_zero =kronecker(id,id)
    
    open(unit=2, file=filename, action="write")
    do j=1,nruns
       
       if (j*dt .gt. t_start) then
          !rho = reshape ( (/ 0.5_dp,0._dp,0._dp,0.5_dp/),(/2,2/) )
          P = matmul( expm(j*dt-t_start,liouvillian(c,H)), P_zero)
          g1 = real(trace(matmul(cdagger,c)*rho))

          tmpvec = matmul(kronecker(conjg(c),c),rhovec)
          tmpvec = matmul(P,tmpvec)
          tmp2(1,1) = tmpvec(1)
          tmp2(2,1) = tmpvec(2)
          tmp2(1,2) = tmpvec(3)
          tmp2(2,2) = tmpvec(4)
          G2 = real(trace(matmul(cdagger,c)*tmp2))       
          g2_norm = G2/g1**2
          write(2,'(E22.7,A1,E22.7)') j*dt-t_start, char(9), g2_norm
       endif

       rhovec = matmul( expm(j*dt,liouvillian(c,H)), rhovec_zero)
       rho(1,1) = rhovec(1)
       rho(2,1) = rhovec(2)
       rho(1,2) = rhovec(3)
       rho(2,2) = rhovec(4)

    end do
    close(2)
call print_matrix_real(real(rho))
     end subroutine correlations

  subroutine wigner(nruns,gamma, c, dt ,filename, rhovec_zero, H,tau,time,x,p,W)
    integer, intent(in) :: nruns
    real(dp), intent(in) :: dt,gamma,tau,time,x,p
    complex(dp), intent(in) :: rhovec_zero(:),c(:,:),H(:,:)
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: W

    complex(dp) :: rho(2,2), rhovec_t(4)
    real(dp) :: G1, G2, g2_norm,tr,lim
    integer ::n
    real(dp), allocatable :: lambda_1(:), lambda_2(:),result(:,:)

    n=1000
    allocate(lambda_1(n),lambda_2(n),result(n,n))
    lim = 6.0
    
    call linspace(-lim,lim,n,lambda_1)
    call linspace(-lim,lim,n,lambda_2)
   
    open(unit=2, file=filename, action="write")

    
    result=integrand(c,rhovec_zero,gamma,H,tau,time,lambda_1(1), lambda_2(1),x,p) 
   
         
      !    write(2,'(E22.7,A1,E22.7)') j*dt-t_start, char(9), g2_norm
    


    close(2)
  end subroutine wigner

  function integrand(c,rhovec_zero,gamma, H,tau,time,lambda_1, lambda_2,x,p) result(y)
    real(dp), intent(in) :: lambda_1, lambda_2,x,p,tau,time,gamma
    complex(dp), intent(in) :: rhovec_zero(:),c(:,:),H(:,:)
    real(dp) :: tr,y,e
    
    complex(dp) :: rhovec_t(4),L(4,4),tmp(4,4),tst(4),tr_matrix(2,2)
    L = identity_matrix(4) 
    

     
    rhovec_t = matmul( expm(time,liouvillian(c,H)), rhovec_zero)

     tmp = expm(tau,liouvillian(c,H) - L*(lambda_1**2+lambda_2**2)/2._dp &
          & + superM(gamma,c,lambda_1,lambda_2))

     tst = matmul(tmp,rhovec_t)
     
      tr_matrix(1,1) = tst(1)
       tr_matrix(2,1) = tst(2)
       tr_matrix(1,2) = tst(3)
       tr_matrix(2,2) = tst(4)

       tr = trace(tr_matrix)
     
        e=exp(2._dp*i*(p*lambda_1 - x*lambda_2))   
        y = tr*e/pi**2
    
  end function integrand

     
  pure function trace(matrix) result(y)
    complex(dp), intent(in) :: matrix(:,:)
    real(dp) :: y
    integer :: n
    n = size(matrix,1)
    y = 0

    do i=1,n
       y = y + real(matrix(i,i))
    enddo
  end function trace

  function superH(r,rho) result(y)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::y(:,:)

    integer :: n
    n=size(r,1)

    allocate(y(n,n))

    y = matmul(r,rho) + matmul(rho,dagger(r)) - trace( matmul(r,rho) + matmul(rho,dagger(r)) )*rho
  end function superH

  pure function superG(r,rho) result(y)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::y(:,:)
    real(dp) :: tracevalue

    integer :: n
    n=size(r,1)

    allocate(y(n,n))

    tracevalue = real( trace( matmul(matmul(r,rho),dagger(r)) ) )

    if (tracevalue /= 0) then
       y = matmul(matmul(r,rho),dagger(r))/tracevalue - rho
    else
       y = 0!- rho
    end if

  end function superG

  pure function superD(c,rho) result(y)

    complex(dp), intent(in) :: c(:,:), rho(:,:)
    complex(dp), allocatable ::y(:,:)

    complex(dp), allocatable :: cdagger(:,:)
    integer :: n
    n=size(c,1)

    allocate(y(n,n),cdagger(n,n))

    cdagger = conjg(transpose(c))

    y = matmul(matmul(c,rho),cdagger) - 0.5_dp*matmul(matmul(cdagger,c),rho) - 0.5_dp*matmul(matmul(rho,cdagger),c)

  end function superD

  function superM(gamma,c,lambda_1,lambda_2) result(y)
    real(dp), intent(in) :: gamma,lambda_1,lambda_2
    complex(dp), intent(in) :: c(:,:)
    complex(dp), allocatable :: y(:,:)
    integer :: n
    complex(dp), dimension(2,2) :: identity
    complex(dp) :: i
    i = complex(0,1)
    n = size(c,1)
    
    allocate(y(n**2,n**2))
    identity = reshape ( (/1,0,0,1/),(/2,2/) )

    y= sqrt(gamma)*(lambda_1+i*lambda_2)*kronecker(conjg(c),identity) - &
         & sqrt(gamma)*(lambda_1-i*lambda_2)*kronecker(identity,c)
    
  end function superM


  pure function dagger(matrix) result(y)
    complex(dp), intent(in) :: matrix(:,:)
    complex(dp), allocatable :: y(:,:)
    integer :: n
    n = size(matrix,1)
    allocate(y(n,n))

    y = conjg(transpose(matrix))
  end function dagger

  pure function identity_matrix(n) result(y)
    integer, intent(in) :: n
    real(dp), dimension(n,n) :: y

    y = 0
    do i=1,n
       y(i,i) = 1
    end do
  end function identity_matrix

pure subroutine linspace(d1,d2,n,grid)
implicit none

integer, intent(in) :: n ! number of points
double precision, intent(in) :: d1, d2 ! endpoints
double precision, dimension(n), intent(out) :: grid

integer :: indxi


grid(1) = d1
do indxi= 0,n-2
   grid(indxi+1) = d1+(dble(indxi)*(d2-d1))/dble(n-1)
end do
grid(n) = d2

!matlab
!grid = [d1+(0:n-2)*(d2-d1)/(floor(n)-1) d2];
end subroutine



  function liouvillian(c, H_eff) result(y)
    complex(dp), intent(in) :: H_eff(:,:), c(:,:)
    complex(dp) :: i = complex(0,1)
    complex(dp), dimension(2,2) :: identity
    complex(dp), allocatable :: y(:,:)
    integer :: n
    n = size(c,1)
    allocate(y(n**2,n**2))

    identity = reshape ( (/1,0,0,1/),(/2,2/) )

    y = -i*kronecker(identity,H_eff) + i*kronecker(conjg(H_eff),identity) +&
         &kronecker(conjg(c),c)

  end function liouvillian

  function kronecker(A,B) result(y)
    complex(dp), dimension (:,:), intent(in)  :: A, B
    complex(dp), dimension (:,:), allocatable :: y
   integer :: i = 0, j = 0, m = 0, n = 0, p = 0, q = 0

   allocate(y(size(A,1)*size(B,1),size(A,2)*size(B,2)))
   y=0
   do i = 1,size(A,1)
      do j = 1,size(A,2)
      y(1+size(B,1)*(i-1):size(B,1)+size(B,1)*(i-1),1+size(B,2)*(j-1):size(B,2)+size(B,2)*(j-1))=A(i,j)*B
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

end module module1