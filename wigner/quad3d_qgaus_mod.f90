module quad3d_qgaus_mod
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: lim = 1000
  
  private
  !hide all names from the outside,
  public quad3d_qgaus
  ! except quad3d itself.
  real(dp) :: xsav,ysav
  interface
     !user-supplied functions.
     function func(x,y,z) result(r)
       !the three-dimensional function to be integrated.
       use module1
       real(dp), intent(in) :: x,y, z
       real(dp) :: r

       r = x**2 + y**2 + z**2
     end function func
     
     function y1(x) result(r)
       use module1
       real(dp), intent(in) :: x
       real(dp) :: r
       
     end function y1
     
     function y2(x) result(r)
       use module1
       real(dp), intent(in) :: x
       real(dp) :: r
        r = 2._dp
     end function y2
     
     function z1(x,y) result(r)
       use module1
       real(dp), intent(in) :: x,y
       real(dp) :: r
       r = 2._dp
     end function z1
     
     function z2(x,y) result(r)
       use module1
       real(dp), intent(in) :: x,y
       real(dp) :: r
        r = 2._dp
     end function z2
     
  end interface


  !the routine quad3d qgaus returns as ss the integral of a user-supplied function func over a three-dimensional region specified by the limits x1 , x2 , and by the user-supplied functions y1 , y2 , z1 ,and z2 , as defined in (4.6.2). integration is performed by calling qgaus recursively.

contains
  function H(x)
    ! this is H of eq. (4.6.5).
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: H
    integer :: i
    do i=1,size(x)
       xsav=x(i)
       H(i)=qgaus(g,y1(xsav),y2(xsav))
    end do
  end function h
  
  function G(y)
    ! this is G of eq. (4.6.4).
    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(size(y)) :: G
    integer :: j
    do j=1,size(y)
       ysav=y(j)
       G(j)=qgaus(f,z1(xsav,ysav),z2(xsav,ysav))
    end do
  end function g
  
  function f(z)! the integrand f ( x,y,z ) evaluated at fixed x and y .
    real(dp), dimension(:), intent(in) :: z
    real(dp), dimension(size(z)) :: f
    f=func(xsav,ysav,z)
  end function f
  
  recursive function qgaus(func,a,b)
    real(dp), intent(in) :: a,b
    real(dp) :: qgaus
    interface
       function func(x)
         use module1
         real(dp), dimension(:), intent(in) :: x
         real(dp), dimension(size(x)) :: func
       end function func
    end interface
    real(dp) :: xm,xr
    real(dp), dimension(5) :: dx, w = (/ 0.2955242247_sp,0.2692667193_sp,&
         0.2190863625_sp,0.1494513491_sp,0.0666713443_sp /),&
         x = (/ 0.1488743389_sp,0.4333953941_sp,0.6794095682_sp,&
         0.8650633666_sp,0.9739065285_sp /)
    xm=0.5_sp*(b+a)
    xr=0.5_sp*(b-a)
    dx(:)=xr*x(:)
    qgaus=xr*sum(w(:)*(func(xm+dx)+func(xm-dx)))
  end function qgaus

  subroutine quad3d_qgaus(x1,x2,ss)
    real(dp), intent(in) :: x1,x2
    real(dp), intent(out) :: ss
    ss=qgaus(h,x1,x2)
  end subroutine quad3d_qgaus

end module quad3d_qgaus_mod
