	PROGRAM sort
	IMPLICIT NONE
	INTEGER ::i,j,k,n,nres
	REAL :: mode1,mode2


	OPEN(unit=24,file='holder.dat',status='old')
	OPEN(unit=25,file='old_mode_new_mode.dat',status='unknown')
	OPEN(unit=26,file='protname.txt',status='old')
	READ(26,*)
	READ(26,*)nres
	CLOSE(26)
	n=nres-1
	WRITE(*,*)'Number of modes : ',n
	DO i=1,n
		READ(24,*)mode1,mode2
		WRITE(25,*)INT(mode1),INT(mode2)
	END DO
	CLOSE(25)
	CLOSE(24)
	STOP
	END PROGRAM sort
		
