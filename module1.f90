module module1

  integer, parameter :: dp = selected_real_kind(15, 307)

contains

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


  pure function trace(matrix)
    complex(dp), intent(in) :: matrix(:,:)
    complex(dp) :: trace
    integer :: n 
    n = size(matrix,1)
    trace = 0

    do i=1,n
       trace = trace + matrix(i,i)
    enddo
  end function trace

  pure function superH(r,rho)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::superH(:,:) 

    integer :: n
    n=size(r,1)

    allocate(superH(n,n))

    superH = matmul(r,rho) + matmul(rho,dagger(r)) - trace( matmul(r,rho) + matmul(rho,dagger(r)) )*identity_matrix(n)

  end function superH

  pure function superG(r,rho)

    complex(dp), intent(in) :: r(:,:), rho(:,:)
    complex(dp), allocatable ::superG(:,:) 

    integer :: n
    n=size(r,1)

    allocate(superG(n,n))

    superG = matmul(matmul(r,rho),dagger(r))/trace( matmul(matmul(r,rho),dagger(r)) ) - rho

  end function superG

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

  pure function delta_rho(rho,H,c,cdagger,delta_t,dN)
    complex(dp), intent(in) :: rho(:,:), H(:,:), c(:,:), cdagger(:,:)
    real(dp), intent(in) :: dN, delta_t
    complex(dp), allocatable ::delta_rho(:,:) 

    complex(dp) :: i 
    integer :: n
    i = complex(0,1)
    n = size(rho,1)
    allocate(delta_rho(n,n))

    delta_rho = dN*superG(c,rho)-delta_t*superH(i*H + matmul(cdagger,c)/2,rho)

  end function delta_rho

  function liouvillian(c,cdagger, H_eff) 
    complex(dp), intent(in) :: H_eff(:,:), c(:,:), cdagger(:,:)
    complex(dp) :: i = complex(0,1)
    complex(dp), dimension(2,2) :: identity
    complex(dp), allocatable :: liouvillian(:,:)
    real(dp) :: Omega = 20
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
    integer :: i = 0, j = 0, k = 0, l = 0
    integer :: m = 0, n = 0, p = 0, q = 0

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
    double precision t
    integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
    complex*16       H(ldh,m), wsp(lwsp)

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
    integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
    double precision hnorm
    complex*16 cp, cq, scale, scale2, ZERO, ONE

    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
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


