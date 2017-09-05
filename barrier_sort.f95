	PROGRAM sort
	IMPLICIT NONE
	INTEGER ::i,j,k,n,nres,oldmode,newmode
	REAL :: mode1,mode2
	REAL,ALLOCATABLE :: bar(:),newbar(:)


	OPEN(unit=24,file='fmad_mp_60.dat',status='old')
	OPEN(unit=25,file='fmad_sorted.dat',status='unknown')
	OPEN(unit=26,file='protname.txt',status='old')
	OPEN(unit=27,file='old_mode_new_mode.dat',status='old')
	READ(26,*)
	READ(26,*)nres
	CLOSE(26)
	n=nres-1
	WRITE(*,*)'Number of modes : ',n
	ALLOCATE(bar(n))
	ALLOCATE(newbar(n))
	DO i=1,n
		READ(24,*)bar(i)
	END DO
	DO i=1,n
		READ(27,*)oldmode,newmode
		newbar(i)=bar(oldmode)
		WRITE(25,*)newbar(i)
	END DO
	CLOSE(25)
	CLOSE(24)
	CLOSE(27)
	DEALLOCATE(bar)
	STOP
	END PROGRAM sort
		
