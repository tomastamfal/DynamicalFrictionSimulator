!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************


!**********************************************************************************************************
!       integration.F90:  Integrates the mass  
!**********************************************************************************************************


#include "choice.def"
subroutine integration
	use parameters, only : pi, nr
	use variables, only  : r, mass, rho

	implicit none
	real(8)		:: y(1:nr+1)	
	integer 	:: ii,jj 
	!SPhercial integration

	do jj = 1,nr+1
		y(jj) = rho(jj,1,1) *r(jj)**2D0	
	enddo
	mass(1,:,:)	= 0.5D0 * y(1)*r(1)	


	!Proper integration:
	do ii = 1,nr
		mass(ii+1,:,:)	= mass(ii,:,:) + 0.5D0 * (y(ii) + y(ii+1))*(r(ii+1)-r(ii))		
  	enddo

	mass(:,:,:) = 4D0*pi*mass(:,:,:)

end subroutine integration

