	program sort
	INTEGER ::n,i,bottom,mp,nres
	REAL :: bar
	REAL,ALLOCATABLE :: fmad(:,:)

	OPEN(unit=24,file='protname.txt',status='old')
	READ(24,*)
	READ(24,*)nres
	WRITE(*,*)nres
	CLOSE(24)

	n=nres-1
	WRITE(*,*)'Number of modes :',n	
	ALLOCATE(fmad(n,2))

	OPEN(unit=24,file='fmad_mp_60.dat',status='old')
	DO i=1,n
		READ(24,*)bar
		fmad(1,i)=i
		fmad(2,i)=bar
	END DO
	
	STOP
	END PROGRAM sort
