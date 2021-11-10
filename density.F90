!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************

!**********************************************************************************************************
!       density.F90:  Computes the density on various meshes and stores the mass as a txt file
!**********************************************************************************************************

#include "choice.def"
subroutine density

!**********************************************************************************************************    
        use parameters, only : nr 


	use variables, only : r
        use variables, only : rho, mass


#ifdef ALPHABETAGAMMA
        use parameters, only : rho_s 
        use parameters, only : alpha, beta, gamma
        use parameters, only : rs
	!use parameters, only : rvir, rs, rdecay, c, epsilon
#endif!     ALPHABETAGAMMA

#ifdef HERNQUIST
        use parameters, only : rho_s 
        use parameters, only : rs
#endif!     HERNQUIST

#ifdef TEST
        use parameters, only : rho_s 
        use parameters, only : alpha, beta, gamma
        use parameters, only : rs
#endif!     TEST


#ifdef NFW
        use parameters, only : rho_s 
        use parameters, only : alpha, beta, gamma
        use parameters, only : rs
#endif!     NFW

#ifdef SIS
        use parameters, only : rho_s 
#endif!     SIS



!**********************************************************************************************************
	implicit none
	integer :: i 

!**********************************************************************************************************
	print *, 'Distributing the density profile'

!**************************** Computes linear r mesh *************************************
#ifdef HERNQUIST
	print *, 'Choosen model: HERNQUIST '
        print *,  rs

	do i = 1, nr+1
                rho(i,:,:) = rho_s * (rs/r(i)) * (1D0/( r(i) +rs)**3)
                !rho(i,:,:) = rho_s / ( ((r(i)/rs)**gamma) * ((1+(r(i)/rs)**alpha)**((beta-gamma)/alpha)) )
	enddo

#endif!     HERNQUIST

#ifdef TEST
	print *, 'Choosen model: TEST'
	do i = 1, nr+1
                rho(i,:,:) = rho_s / ( ((r(i)/rs)**gamma) * ((1+(r(i)/rs)**alpha)**((beta-gamma)/alpha)) )
	enddo

#endif!     TEST

#ifdef NFW
	print *, 'Choosen model: NFW '
	do i = 1, nr+1
                rho(i,:,:) = rho_s / ( ((r(i)/rs)**gamma) * ((1+(r(i)/rs)**alpha)**((beta-gamma)/alpha)) )
	enddo

#endif!     NFW


#ifdef ALPHABETAGAMMA
	print *, 'Choosen model: ALPHA-BETA-GAMMA'
        print *, 'ALPHA'
        print *, alpha        
        print *, 'BETA'
        print *, beta
        print *, 'GAMMA'
        print *, gamma


	do i = 1, nr+1
		!do j = 1, nphi
			!do k = 1, ntheta
				rho(i,:,:) = rho_s / ( ((r(i)/rs)**gamma) * ((1+(r(i)/rs)**alpha)**((beta-gamma)/alpha)) )
			!enddo
		!enddo
	enddo
				!if (r(i) .le. rvir) then
				!	rho(i,:,:) = rho_s / (((r(i)/rs)**gamma) * (1+(r(i)/rs)**alpha)**((beta-gamma)/alpha))
				!else
				!	rho(i,:,:) = (rho_s / ((c**gamma) *(1+c**alpha)**((beta-gamma)/alpha)) ) &
				!		  *((r(i)/rvir)**epsilon) *DEXP((rvir-r(i))/rdecay)
				!endif

#endif!     ALPHABETAGAMMA

#ifdef SIS
	
	!print *, 'Choosen model: ALPHA-BETA-GAMMA    

	do i = 1, nr+1
                rho(i,:,:) = rho_s / (r(i)**2D0) ! Msol / kpc^3
	enddo
#endif!     SIS

	open(11, file='density_data.txt', status='replace')
	do i = 1, nr+1
		write(11,*) rho(i,1,1)
	enddo
	close(11)

	print *, 'Calculating the mass'

	call integration


!***************** Print the mesh in r.data and theta.data *******************************
	open(10, file='mass_data.txt', status='replace')
	do i = 1, nr+1
		write(10,*) mass(i,1,1)
	enddo
	close(10)

return
!**********************************************************************************************************
end subroutine density
!**********************************************************************************************************
