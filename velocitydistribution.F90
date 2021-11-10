!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************



!**********************************************************************************************************
!            Create_mesh.F90:  Computes the mesh of the computing grid and outputs it
!**********************************************************************************************************

#include "choice.def"
subroutine velocitydistribution

!**********************************************************************************************************
        use parameters, only : pi
        use parameters, only : G, sigma

	use variables, only : r
	use variables, only : nr, nphi, ntheta
	use variables, only : mass
	use variables, only : velocityintegral

#ifdef TEST   
        use parameters, only : M_sun
        use variables, only : rho
#endif!     TEST  

!**********************************************************************************************************
	implicit none
	integer :: i
        real(8) :: X(1:nr+1, 1:nphi+1, 1:ntheta+1) 
!**********************************************************************************************************
	print *, 'Calculating velocity distribution '



	do i =1,nr+1
                X(i,:,:) = DSQRT(G*mass(i,:,:)/(2D0*r(i))) / sigma                             
                velocityintegral(i,:,:) = DERF(X(i,:,:)) - ( ( (2D0*X(i,:,:))/DSQRT(pi) ) * DEXP(-X(i,:,:)**2D0) )
        enddo
#ifdef TEST
            velocityintegral(:,:,:) = 1
#endif!     TEST     



return
!**********************************************************************************************************
end subroutine velocitydistribution
!**********************************************************************************************************
