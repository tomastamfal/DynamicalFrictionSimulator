!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************


!**********************************************************************************************************
!            Parameters.F90:  Sets different physical parameters and normalization
!**********************************************************************************************************

#include "choice.def"
module parameters

!**********************************************************************************************************
	implicit none    
!**************************** Grid parameters ********************************************
	integer, parameter :: nr 	= 1D4 ! resolution in x direction
	integer, parameter :: nphi 	= 1D2   ! resolution in y direction
	integer, parameter :: ntheta 	= 10

        real(8), parameter :: kpc_to_m	= 1D0 !3.086D19
	real(8), parameter :: Myr_to_s	= 1D0 !3.154D13
        real(8), parameter :: kms_kpcGyr= 1.02268944D0!
	real(8), parameter :: G		= 4.49975332435D-6 ! kpc^3 / (Gyr^2 Msol)

        real(8), parameter :: M_sun 	= 1D0 !1.989D30          ! Msol
        real(8), parameter :: rin	= 0.1D0 * kpc_to_m

        real(8), parameter :: pi  	= 3.14159265358979323846D0
	real(8), parameter :: r_stop 	= 1D0 * kpc_to_m


#ifdef ALPHABETAGAMMA
        real(8), parameter :: sigma     = 20D0 * kms_kpcGyr
        real(8), parameter :: c		= 20D0
        real(8), parameter :: M_BH	= M_sun * 1D5 

        real(8), parameter :: rvir 	= 71D0 * kpc_to_m
        real(8), parameter :: rdecay    = 71D0 * kpc_to_m
        real(8), parameter :: rout 	= 100D0 * kpc_to_m
        real(8), parameter :: r_start	= 800D0 * kpc_to_m


        real(8), parameter :: alpha 	= 1D0 
        real(8), parameter :: beta 	= 3D0 
	real(8), parameter :: gamma 	= 1D0

        !MODEL
        real(8), parameter :: rs 	= (rvir/c)  
	real(8), parameter :: rho_crit	= 135.950245067
        real(8), parameter :: DELTA_VIR = 103.5 
        real(8), parameter :: delta_c   = DELTA_VIR * c**3.0 /( 3.0*( dlog(1+c) - (c/(1+c)) ) )
        real(8), parameter :: rho_s	= delta_c*rho_crit 

        real(8), parameter :: epsilon	= (-gamma-beta*c**alpha)/(1+c**alpha) + (rvir/rdecay)
#endif! ALPHABETAGAMMA


#ifdef HERNQUIST
        real(8), parameter :: sigma     = 200D0  * kms_kpcGyr
        real(8), parameter :: M_BH	= 1D11   * M_sun

        real(8), parameter :: rvir 	= 310D0  * kpc_to_m
        real(8), parameter :: rout 	= 130D0  * kpc_to_m
        real(8), parameter :: r_start	= 122D0  * kpc_to_m

        real(8), parameter :: c		= 12D0
        real(8), parameter :: r_h 	= (0.6082 - 0.1843*DLOG10(c) - 0.1011*DLOG10(c)**2 + 0.03918*DLOG10(c)**3)*rvir
        real(8), parameter :: rs 	= (r_h/(1D0+DSQRT(2D0)))

        real(8), parameter :: M_galaxy	= 1e12 * M_sun
        real(8), parameter :: rho_s	= M_galaxy/(2D0*pi) 

#endif! HERNQUIST

#ifdef TEST
        real(8), parameter :: sigma     = 0D0 !DUMMY INDEX ...
        real(8), parameter :: M_BH	= 5D5   * M_sun
        real(8), parameter :: M_galaxy	= 1D10  * M_sun

        real(8), parameter :: rvir 	= 0.0001D0  * kpc_to_m
        real(8), parameter :: rout 	= 0.286D0 * kpc_to_m
        real(8), parameter :: r_start	= 0.285D0 * kpc_to_m

        !MODEL
        real(8), parameter :: alpha 	= 1D0 
        real(8), parameter :: beta 	= 4D0 
	real(8), parameter :: gamma 	= 1D0

        real(8), parameter :: rs 	= 0.7D0 *kpc_to_m
        real(8), parameter :: rho_s	= M_galaxy/( 2D0*pi) 
#endif! TEST



#ifdef NFW
        
        real(8), parameter :: sigma     = 230D0  * kms_kpcGyr
        real(8), parameter :: c		= 12D0
        real(8), parameter :: M_BH	= M_sun * 1D11
        real(8), parameter :: M_galaxy	= M_sun * 1D12

        real(8), parameter :: rvir 	= 310D0 *kpc_to_m
        real(8), parameter :: rdecay    = 310D0  * (20D0**(1D0/3D0))*kpc_to_m
        real(8), parameter :: rout 	= 130D0  * kpc_to_m
        real(8), parameter :: r_start	= 122D0 * kpc_to_m

        !MODEL
        real(8), parameter :: alpha 	= 1D0 
        real(8), parameter :: beta 	= 3D0 
	real(8), parameter :: gamma 	= 1D0

        real(8), parameter :: rs 	= (rvir/c)  
	real(8), parameter :: rho_crit	= 135.950245067
        real(8), parameter :: DELTA_VIR = 103.5 
        real(8), parameter :: delta_c   = DELTA_VIR * c**3.0 /( 3.0*( dlog(1+c) - (c/(1+c)) ) )
        real(8), parameter :: rho_s	= delta_c*rho_crit 

        real(8), parameter :: epsilon	= (-gamma-beta*c**alpha)/(1+c**alpha) + (rvir/rdecay)
#endif! NFW



#ifdef SIS
        real(8), parameter :: sigma     = 200D0 * kms_kpcGyr    ! kpc/Gyr
        real(8), parameter :: M_BH	= 1D8   * M_sun         ! Solar mass
        real(8), parameter :: rout 	= 5.1D0 * kpc_to_m      ! kpc
        real(8), parameter :: r_start	= 5.0D0 * kpc_to_m      ! kpc

        !MODEL
        real(8), parameter :: vc        = DSQRT(2D0)*sigma      ! kpc/Gyr
        real(8), parameter :: rho_s	= vc**2D0 / (4*pi*G)    ! kpc^2 / kpc^3 (Gyr^2 Msol^1)/Gyr^2 = Msol / kpc (REMEMBER: rho = rho_s/r^2)
#endif! SIS







!**********************************************************************************************************
end module parameters
!**********************************************************************************************************
