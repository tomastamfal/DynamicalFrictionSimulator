!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************

!**********************************************************************************************************
!            Create_mesh.F90:  Computes the mesh of the computing grid and outputs it
!**********************************************************************************************************

#include "choice.def"
subroutine mesh

!**********************************************************************************************************
	use parameters, only : pi, rin, rout
        use parameters, only : nr, ntheta, nphi 

	use variables, only : r, phi, ntheta
	use variables, only : dr, dphi, dtheta
	use variables, only : r_center, phi_center, theta_center
        use variables, only : theta



!**********************************************************************************************************
	implicit none
	integer :: i, j, k

!**********************************************************************************************************
!print *, 'Creating mesh'

!**************************** Computes linear r mesh *************************************
#ifdef LIN_MESH
	r(1)	= rin
	dr(1)	= (rout - rin) / dble(nr)
	r_center(1) = r(1) + 0.5D0*dr(1)
	do i = 2, nr
		dr(i) 		= dr(1)
		r(i) 		= r(i-1) +  dr(i)  
		r_center(i) 	= r(i)   + 0.5D0*dr(i)
	enddo
	r(nr+1) = r(nr) +  dr(nr) 

#endif!     LIN_MESH

!**************************** Computes log r mesh ****************************************
#ifdef LOG_MESH
	do i = ini, nr
		r(i) 		= rin * ratio**(nx*dble(my_i_mpi-1) + dble(i) - dble(ini+nghost) + 0.5D0)
		dr(i) 		= r(i) * (dsqrt(ratio) - 1.D0/dsqrt(ratio))
		r_center(i) 	= r(i)*dsqrt(ratio)
	enddo  
#endif LOG_MESH
!**************************** Computes phi mesh ********************************************
	phi(1)	= 0D0
	dphi(1)	= 2D0*pi / dble(nphi)
	phi_center(1) = phi(1) + 0.5D0*dphi(1)
	do i = 2, nphi
		dphi(i) 	= dphi(1)
		phi(i) 		= phi(i-1) + dphi(i)
		phi_center(i) 	= phi(i)   + 0.5D0*dphi(i)		
	enddo
	phi(nphi+1) = phi(nphi) +  dphi(nphi)

!**************************** Computes theta mesh ********************************************
	theta(1)	= 0D0
	dtheta(1)	= 1D0*pi / dble(ntheta)
	theta_center(1) = theta(1) + 0.5D0*dtheta(1)
	do i = 2, ntheta
		dtheta(i) 		= dtheta(1)
		theta(i) 		= theta(i-1) + dtheta(i)
		theta_center(i) 	= theta(i)   + 0.5D0*dtheta(i)		
	enddo
	theta(ntheta+1) = theta(ntheta) + dtheta(ntheta)

#ifdef LOG_MESH

!**************************** Computes dimensions of cells *******************************
	do i = ini, x_end         
		cell_vol(i) = 0.5D0 * (x_r(i) + x_r(i-1)) * dx(i) * dy   ! Volume of the cell
	enddo

#endif LOG_MESH
!***************** Print the mesh in r.data and theta.data *******************************
	open(10, file='r_grid.txt', status='replace')
	do i = 1, nr+1
		write(10,*) r(i)
	enddo
	close(10)

	open(20, file='phi_data.txt', status='replace')
	do j = 1, nphi+1
		write(20,*) phi(j)
	enddo
	close(20)

	open(30, file='theta_data.txt', status='replace')
	do k = 1, ntheta+1
		write(30,*) theta(k)
	enddo
	close(30)


return
!**********************************************************************************************************
end subroutine mesh
!**********************************************************************************************************
