!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************

!**********************************************************************************************************
!       dynamicalfrictionforce.F90:  Computes the dynamical friction force according to 
!                                    Chandrasekhar (1943) by assuming a constant log lambda term
!**********************************************************************************************************

#include "choice.def"
subroutine dynamicalfrictionforce 

!**********************************************************************************************************
        use parameters, only : nr, M_BH, G, pi


	use variables, only : r
	use variables, only : rho, mass
        use variables, only : F_grid
	use variables, only : loglambda, velocityintegral




!**********************************************************************************************************
	implicit none
	integer :: i

!**********************************************************************************************************
	print *, 'Calculating dynamical friction matrix'
!**************************** Computes linear r mesh *************************************

	do i = 1,nr+1
		F_grid(i,:,:)=(-4D0*pi*(G**2D0)*(M_BH**2D0)*loglambda*rho(i,:,:)*velocityintegral(i,:,:))/(G*mass(i,:,:)/r(i))
	enddo

!**************************** Computes log r mesh ****************************************

	open(10, file='F_data.txt', status='replace')
	do i = 1, nr+1
		write(10,*) F_grid(i,1,1)
	enddo
	close(10)


return
!**********************************************************************************************************
end subroutine dynamicalfrictionforce
!**********************************************************************************************************
