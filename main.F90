!*********************************************************************************************************
!**                                     Developed by Tomas Tamfal 					**
!**                                       Version 0.1 April 2017	                               	**
!*********************************************************************************************************


#include "choice.def"

program stopping_time
!**********************************************************************************************************
	use parameters, only : r_stop, r_start
        use parameters, only : nr

	use variables, only : r
	use variables, only : r_position, pos
        use variables, only : F_grid, L_grid

	implicit none
	real(8) :: pos_min
	real(8) :: test(1:nr+1)
	real(8) :: time

	write(*,*) " "
	write(*,*) " " 
	write(*,*) '------------------------------------'
	write(*,*) " "                                  
        write(*,*) "    ,---,        ,---,.  .--.--.    "
        write(*,*) "  .'  .' `\    ,'  .' | /  /    '.  "
        write(*,*) " ,---.'     \ ,---.'   ||  :  /`. / "
        write(*,*) " |   |  .`\  ||   |   .';  |  |--`  " 
        write(*,*) " :   : |  '  |:   :  :  |  :  ;_    " 
        write(*,*) " |   ' '  ;  ::   |  |-, \  \    `. " 
        write(*,*) " '   | ;  .  ||   :  ;/|  `----.   \" 
        write(*,*) " |   | :  |  '|   |   .'  __ \  \  |" 
        write(*,*) " '   : | /  ; '   :  '   /  /`--'  /" 
        write(*,*) " |   | '` ,/  |   |  |  '--'.     / " 
        write(*,*) " ;   :  .'    |   :  \    `--'---'  "
        write(*,*) " |   ,.'      |   | ,'              "
        write(*,*) " '---'        `----'                " 
        write(*,*) " "                        
	write(*,*) "------------------------------------"
	write(*,*) " Dynamical Friction Simulator"
	write(*,*) "------------------------------------"
	write(*,*) " "
	write(*,*) " "


       	!**************************** Creating mesh **********************************************
	call mesh
	!**************************** Creating density and mass **********************************
	call density
	!**************************** Creating density and mass **********************************
	call angularmomentum
	!**************************** Creating density and mass **********************************
	call lnlambda
	!**************************** Creating density and mass **********************************
	call velocitydistribution
	!**************************** Creating density and mass **********************************
	call dynamicalfrictionforce
	

	!**************************** Find starting grid point ***********************************
	test(:) 	= ABS(r(:) - r_start)
	pos 		= MINLOC(test, DIM=1)
	r_position 	= r(pos)

	time = 0D0
	test(:) 	= ABS(r(:) - r_stop)
	pos_min 	= MINLOC(test, DIM=1)

	!**************************** Find starting grid point ***********************************

	!**************************** Main loop **************************************************
	
	open(20, file='time_data.txt', status='replace')
	write(20,*) time
        close(20)

	open(30, file='r_data.txt', status='replace')
	write(30,*) r(pos)
        close(30)

	IF (pos_min == 1D0) THEN
		write(*,*) 'Inner boundary to small'
		stop
	END IF
	
	write(*,*) "Starting calculation"
	
	do while (pos .GE. pos_min)
		time = time +  ABS(L_grid(pos,1,1)-L_grid(pos-1,1,1))/ABS((r(pos)*F_grid(pos,1,1)))
		pos  = pos - 1
		
		open(10, file='time_data.txt',status='old', position='append')
		write(10,*) time                
                close(10)

		open(40, file='r_data.txt',status='old', position='append')
		write(40,*) r(pos)
                close(40)


	enddo
	write(*,*) '======================='
        write(*,*) 'Stopping time in Gyr:'        
	write(*,*) time 
	write(*,*) '======================='
	write(*,*) 'Done!'
	write(*,*) ' ' 
	write(*,*) ' '


end program stopping_time
