
! Author: Mohammad Asif Zaman
! Original version: April, 2012


program nonlinear_solver

implicit none

integer, parameter :: max_iter = 270, data_type = 10
real(kind = data_type) ::  err_lim = 1e-14, err_store = 0.1, eps = 1e-19
real(kind = data_type), dimension(max_iter) :: x = 13, y = 13, fx = 13, fy = 13, dfx = 13, dfy = 13
integer :: n, iter_store, function_select, method_select


print *, eps

print *, "Enter test function number (1 - 11): "
read *, function_select
print *

print *, "Enter method number (1 - 4): "
print *, "1. NM"
print *, "2. DNM"
print *, "3. JM"
print *, "4. KJ"


read *, method_select
print * , "================================================"
print *

! Initial values for the test functions
select case (function_select)
	case(1) ;	x(1) = 1
	case(2) ;	x(1) = 0
	case(3) ;	x(1) = -1.5
	case(4) ;	x(1) = 2
	case(5) ;	x(1) = 3.5
	case(6) ;	x(1) = 3
	case(7) ;	x(1) = 1
	case(8) ;	x(1) = 1
	case(9); x(1) = 1
	case(10); x(1) = 3.5
	case(11); x(1) = 2
end select






Do n = 1, max_iter

	
	if (n > 1) then
	if (abs(x(n) - x(n-1)) < err_lim) then
		err_store = abs(x(n) - x(n-1))
		iter_store = n - 1
		exit
	end if
	end if 
	
	
	
	select case (method_select)
	
		case(1)
			! Newton's Method Start
			call test_function(function_select, x(n), fx(n), dfx(n)) 
			y(n) = x(n) - fx(n)/(dfx(n) + eps) 
			x(n+1) = y(n)
			!Neton's Method end
			
		case(2)
			! Double Newton's Method Start
			call test_function(function_select, x(n), fx(n), dfx(n)) 
			y(n) = x(n) - fx(n)/(dfx(n) + eps) 
			call test_function(function_select, y(n), fy(n), dfy(n)) 
			x(n + 1) = y(n) - fy(n)/(dfy(n) + eps) 
			! Double Neton's Method end
		
		case(3)
			! Jarratt Method Start
			call test_function(function_select, x(n), fx(n), dfx(n)) 
			y(n) = x(n) - (2/3) * fx(n) / dfx(n)
			call test_function(function_select, y(n), fy(n), dfy(n)) 
			x(n + 1) = x(n) - (1 - (3/2)*(dfy(n) - dfx(n))/(3*dfy(n) - dfx(n)))*(fx(n)/dfx(n))
			! Jarratt Method End
			
		case(4)
			! Kou Jesheng's Method Start
			call test_function(function_select, x(n), fx(n), dfx(n)) 
			y(n) = x(n) - fx(n)/(dfx(n) + eps) 
			call test_function(function_select, y(n), fy(n), dfy(n)) 
			x(n+1) = x(n) - (fx(n)**2 + fy(n)**2)/( dfx(n)*(fx(n) - fy(n)) + eps)
			! Kog Jesheng's Method End
		
					
			
	end select


end do




Do n = 1, max_iter
	!if (fx(n) == 0) exit			! The exit condition,
	!print *, x(n),  "  " ,fx(n)
end do


print *
Print *, "================================================"
print *, "Summary Result"
Print *, "================================================"

Print *, "Test function number       = ", function_select
print *, "Method number              = ", method_select
Print *, "================================================"
print *
print *, "Number of iterations       = ", iter_store
print *, "Defined Error limit        = ", err_lim
print *, "Final value of error       = ", err_store 
print *, "Calculated root            = ", x(iter_store) 
print *, "f(x)                       = ", fx(iter_store) 
print *
Print *, "================================================"







 contains 
  
	
  ! Test functions
  subroutine test_function(function_select, x, fx, dfx) 
    integer, intent(in) :: function_select
    !integer, parameter, intent(in) ::data_type
    real(kind = data_type) , intent(in) :: x
    real(kind = data_type) , intent(out) :: fx, dfx
    
    select case (function_select)
	case(1);    fx = x**3 + 4*x**2- 10;    dfx = 3*x**2 + 8*x
	case(2);    fx = x**2 - exp(x) - 3*x + 2;    dfx = 2*x - exp(x) - 3
	case(3)  
			  fx = x*exp(x**2) - (sin(x))**2 + 3*cos(x) + 5
 			  dfx = x*exp(x**2)*2*x + exp(x**2) - 2*sin(x)*cos(x) - 3*sin(x)
	case(4)
			  fx = exp(x)*sin(x) + log(x**2 + 1)
			  dfx = exp(x)*cos(x) + sin(x)*exp(x) + (1/(x**2 + 1))* 2*x
				
	case(5);    fx = (x - 1)**2 - 1; dfx = 2*(x - 1);
	case(6);    fx = (x - 1)**3 - 2;    dfx = 3*(x - 1)**2
	case(7);    fx = 10*x*exp(-x**2) - 1; dfx = 10*x*exp(-x**2)*(-2*x) + 10*exp(-x**2)
	case(8);    fx = cos(x) - x;    dfx = -sin(x) - 1
	case(9);    fx = (sin(x))**2 - x**2 + 1;    dfx = 2*sin(x)*cos(x) - 2*x
	case(10);  fx = exp(x**2 + 7*x - 30) - 1;    dfx = exp(x**2 + 7*x - 30)*(2*x + 7)
	case(11);  fx = (x + 2)*exp(x) - 1; dfx = (x + 2)*exp(x) + exp(x)	
	
    end select
    
  end subroutine test_function
  
end program nonlinear_solver