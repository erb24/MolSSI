	program inputreader
	IMPLICIT NONE
	integer :: nfrs,n,nmol,nstates
	CHARACTER(64) :: protname
	open(unit=5,file="protname.txt",status='old')
	read(5,*)protname
	read(5,*)n
	read(5,*)nfrs
	close(5)
	protname=adjustl(protname)
	nstates=2
	!open(unit=5,file="nmol.dat",status='old')
	!read(5,*)nmol
	!close(5)
	!nfrs=10000 !Testing
	WRITE(*,*)n,nfrs
	call umatrix(n,nfrs,nstates,protname)
	!call writer(n,nfrs,nstates,protname)
	STOP
	End program inputreader

	subroutine umatrix(n,nfrs,nstates,protname)
	IMPLICIT NONE
	INTEGER :: i,j,k,a,n,nfrs,nstates,count1,count2,nmol,state,xcount,ycount,zcount,ncount(nstates)
	DOUBLE PRECISION :: thresh,umat1(3*n,3*n),umat2(3*n,3*n),rmsd(nfrs),prob(nstates,nstates),C(3*n,3*n),rhs1(3*n,3*n)
	DOUBLE PRECISION :: rx1(n,nfrs),ry1(n,nfrs),rz1(n,nfrs),rx2(n,nfrs),ry2(n,nfrs),rz2(n,nfrs),biglist(3*n,nfrs)
	DOUBLE PRECISION :: lx1(n,nfrs),ly1(n,nfrs),lz1(n,nfrs),lx2(n,nfrs),ly2(n,nfrs),lz2(n,nfrs),bigrlist(nstates,3*n,nfrs)
	DOUBLE PRECISION :: lavm1(n),lavm2(n),lavmsq1(n),lavmsq2(n),lmag1(n),lmag2(n),rx(n,nfrs),ry(n,nfrs),rz(n,nfrs)
	DOUBLE PRECISION :: rxref(n,nfrs),ryref(n,nfrs),rzref(n,nfrs),fk(nstates),avgx(n),avgy(n),avgz(n),avglist(3*n)
	DOUBLE PRECISION :: rcomx(nfrs),rcomy(nfrs),rcomz(nfrs),sumx(n),sumy(n),sumz(n),summa(n),p1,p2,cjam(nstates,3*n,3*n)
	DOUBLE PRECISION :: bigavglist(nstates,3*n)
	CHARACTER(64) ::protname
	
	rx=0.0
	ry=0.0
	rz=0.0
	rx1=0.0
	ry1=0.0
	rz1=0.0
	rx2=0.0
	ry2=0.0
	rz2=0.0
	lavm1=0.0
	lavm2=0.0
	lmag1=0.0
	lmag2=0.0
	lx1=0.0
	lx2=0.0
	ly1=0.0
	ly2=0.0
	lz1=0.0
	lz2=0.0
	umat1=0.0
	umat2=0.0
	rxref=0.0
	ryref=0.0
	rzref=0.0
	avgx=0.0
	avgy=0.0
	avgz=0.0
	rcomx=0.0
	rcomy=0.0
	rcomz=0.0
	sumz=0.0
	sumy=0.0
	sumx=0.0
	summa=0.0
	C=0.0
	biglist=0.0
	cjam=0.0
	ncount=0
	bigrlist=0.0
	avglist=0.0
	rhs1=0.0
	bigavglist=0.0
	protname=adjustl(protname)
	WRITE(*,*)protname
	rmsd=0.0
        prob=0.0

	!read from trajectory
	OPEN(unit=10,file='rmsd.dat',status='unknown')
	open(unit=11,file=trim(protname)//'.g96',status='old')
	OPEN(unit=12,file='traj1.dat',status='unknown')
	OPEN(unit=13,file='traj2.dat',status='unknown')
	OPEN(unit=14,file='frames_1',status='unknown')
	OPEN(unit=15,file='frames_2',status='unknown')
	OPEN(unit=16,file='nframes1',status='unknown')
	OPEN(unit=17,file='nframes2',status='unknown')
	WRITE(14,'(A)')'[ State 1 ]'
	WRITE(15,'(A)')'[ State 2 ]'
	!skip first 7,now read and calculate stuff
	do i=1,7
	  read(11,*)
	end do

	do k=1,nfrs
          IF (MOD(k,1000) .EQ. 0)WRITE(*,*)'READING FRAME ',k
	  do j=1,n
	    read(11,*)rx(j,k),ry(j,k),rz(j,k)
	  end do
	  IF(k .EQ. 1)THEN !Reference structure
	    DO i=1,n
	      rxref(i,1)=rx(i,1)
	      ryref(i,1)=ry(i,1)
	      rzref(i,1)=rz(i,1)
	    END DO
	  END IF
	  !Calculate the rmsd
	  DO i=1,n
	    rmsd(k)=rmsd(k)+(SQRT(rx(i,k)**2+ry(i,k)**2+rz(i,k)**2)-SQRT(rxref(i,1)**2+ryref(i,1)**2+rzref(i,1)**2))**2
          end do
	  rmsd(k)=rmsd(k)/REAL(n)
	  rmsd(k)=SQRT(rmsd(k))
          !WRITE(*,*)rmsd(k)
	  WRITE(10,*)1.0*k,rmsd(k)
	  !skip 8 lines
          IF(k .EQ. nfrs)THEN
            do j=1,4
              READ(11,*)
            END DO
          ELSE
	    do j=1,8
	      read(11,*)
	    end do
          END IF
	  !come out of time loop
	end do
	CLOSE(10)
	CLOSE(11)
	rxref=0.0
	ryref=0.0
	rzref=0.0

	ncount(1)=0
	ncount(2)=0
	thresh=0.30
        state=1 !Start in the folded state
	DO k=1,nfrs
	  IF(MOD(k,1000) .EQ. 0)WRITE(*,*)'FRAME ',k
	  IF(rmsd(k) .LE. thresh)THEN
	    ncount(1)=ncount(1)+1
	    WRITE(14,*)k
            IF(state .EQ. 1)prob(1,1)=prob(1,1)+1 !Using Lay's notation
            IF(state .EQ. 2)prob(1,2)=prob(1,2)+1 
	    DO i=1,n
	      rx1(i,ncount(1))=rx(i,k)
	      ry1(i,ncount(1))=ry(i,k)
	      rz1(i,ncount(1))=rz(i,k)
	    END DO
            state=1 !Previous state
	  ELSE IF(rmsd(k) .GT. thresh)THEN
	    WRITE(*,*)'THRESHOLD'
	    ncount(2)=ncount(2)+1
            IF(state .EQ. 1)prob(2,1)=prob(2,1)+1 !Using Lay's notation
            IF(state .EQ. 2)prob(2,2)=prob(2,2)+1
            state=2 !Previous state
	    WRITE(15,*)k
	    DO i=1,n
	      rx2(i,ncount(2))=rx(i,k)
	      ry2(i,ncount(2))=ry(i,k)
	      rz2(i,ncount(2))=rz(i,k)
	    END DO
	    !Calculate bond lengths
	  END IF
	END DO
	WRITE(*,*)ncount(1),ncount(2)
        WRITE(16,*)ncount(1)
        WRITE(17,*)ncount(2)
	CLOSE(12)
	CLOSE(13)
	CLOSE(14)
	CLOSE(15)
        CLOSE(16)
        CLOSE(17)

	!Find weight of each state
	!Do i=1,nstates
	  fk(1)=REAL(ncount(1))/REAL(ncount(1)+ncount(2))
	  fk(2)=REAL(ncount(2))/REAL(ncount(1)+ncount(2))
	  p1=fk(1)
	  p2=fk(2)
	  WRITE(*,*)p1,p2

	!Calculate the C matrix for each basin
	!This calculation will set up the second term on the RHS in the JAM model
	DO k=1,ncount(1)
	    do j=1,n
	      sumx(j)=sumx(j)+rx1(j,k)
	      sumy(j)=sumy(j)+ry1(j,k)
	      sumz(j)=sumz(j)+rz1(j,k)
	      rcomx(k)=rcomx(k)+rx1(j,k)
	      rcomy(k)=rcomy(k)+ry1(j,k)
	      rcomz(k)=rcomz(k)+rz1(j,k)
	    end do
	    rcomx(k)=rcomx(k)/REAL(n)
	    rcomy(k)=rcomy(k)/REAL(n)
	    rcomz(k)=rcomz(k)/REAL(n)
	    DO j=1,n !Calculate radial vectors wrt the com of the molecule, which should remove rotations
	      rx1(j,k)=rx1(j,k)-rcomx(k)
	      ry1(j,k)=ry1(j,k)-rcomy(k)
	      rz1(j,k)=rz1(j,k)-rcomz(k)
	    END DO
	    DO i=1,3*n,3
	      bigrlist(1,i,k)=rx1(i,k)
	      bigrlist(1,i+1,k)=ry1(i,k)
	      bigrlist(1,i+2,k)=rz1(i,k)
	    END DO
	END DO
	!Calculate averages and generate the covariance matrix
	avgx=0.0
	avgy=0.0
	avgz=0.0

	OPEN(unit=24,file='com_avg_coor_1',status='unknown')
	DO i=1,n ! Average coordinates of each residue
	  avgx(i)=sumx(i)/REAL(ncount(1))
	  avgy(i)=sumy(i)/REAL(ncount(1))
	  avgz(i)=sumz(i)/REAL(ncount(1))
	  WRITE(*,*)avgx(i),avgy(i),avgz(i)
	  WRITE(24,*)avgx(i),avgy(i),avgz(i)
	END DO
	CLOSE(24)
	DO i=1,3*n,3
	  bigavglist(1,i)=avgx(i)
	  bigavglist(1,i+1)=avgy(i)
	  bigavglist(1,i+2)=avgz(i)
	END DO
	rxref=0.0
	ryref=0.0
	rzref=0.0
	DO k=1,ncount(1) !Subtract reference values (here the average) from each coordinate
	  DO i=1,n
	    rxref(i,k)=rx1(i,k)-avgx(i)
	    ryref(i,k)=ry1(i,k)-avgy(i)
	    rzref(i,k)=rz1(i,k)-avgz(i)
	  END DO
	END DO
	!Make big list of all augmented coordinates
	xcount=0
	ycount=0
	zcount=0
	WRITE(*,*)'Making the big list of coordinates...'
	OPEN(unit=24,file='com_biglist_1.dat',status='unknown')
	DO k=1,ncount(1)
	  DO i=1,3*n
	    IF(MOD(i,3) .EQ. 1)THEN
	      xcount=xcount+1
	      biglist(i,k)=rxref(xcount,k)
	    ELSE IF(MOD(i,3) .EQ. 2)THEN
	      ycount=ycount+1
	      biglist(i,k)=ryref(ycount,k)	
	    ELSE IF(MOD(i,3) .EQ. 0)THEN
	      zcount=zcount+1
	      biglist(i,k)=rzref(zcount,k)	 
	    END IF
	    WRITE(24,*)biglist(i,k)
	  END DO
	  !Reset counters after each time step
	  xcount=0
	  ycount=0
	  zcount=0
	END DO  
	CLOSE(24) 
	CLOSE(11)
	DO k=1,ncount(1)
	  DO i=1,3*n
	    DO j=1,3*n
	      umat1(i,j)=umat1(i,j)+biglist(i,k)*biglist(j,k)
	    END DO
	  END DO
	END DO
	DO i=1,3*n
	  DO j=1,3*n
	    umat1(i,j)=umat1(i,j)/REAL(ncount(1))
	    !WRITE(*,*)umat1(i,j)
	  END DO
	END DO

	!Do the second basin
	sumx=0.0
	sumy=0.0
	sumz=0.0
	rxref=0.0
	ryref=0.0
	rzref=0.0
	rcomx=0.0
	rcomy=0.0
	rcomz=0.0
	!This calculation will set up the second term on the RHS in the JAM model
	DO k=1,ncount(2)
	    do j=1,n
	      sumx(j)=sumx(j)+rx2(j,k)
	      sumy(j)=sumy(j)+ry2(j,k)
	      sumz(j)=sumz(j)+rz2(j,k)
	      rcomx(k)=rcomx(k)+rx2(j,k)
	      rcomy(k)=rcomy(k)+ry2(j,k)
	      rcomz(k)=rcomz(k)+rz2(j,k)
	    end do
	    rcomx(k)=rcomx(k)/REAL(n)
	    rcomy(k)=rcomy(k)/REAL(n)
	    rcomz(k)=rcomz(k)/REAL(n)
	    DO j=1,n !Calculate radial vectors wrt the com of the molecule, which should remove rotations
	      rx2(j,k)=rx2(j,k)-rcomx(k)
	      ry2(j,k)=ry2(j,k)-rcomy(k)
	      rz2(j,k)=rz2(j,k)-rcomz(k)
	    END DO
	    DO i=1,3*n,3
	      bigrlist(2,i,k)=rx2(i,k)
	      bigrlist(2,i+1,k)=ry2(i,k)
	      bigrlist(2,i+2,k)=rz2(i,k)
	    END DO
	END DO
	!Calculate averages and generate the covariance matrix
	avgx=0.0
	avgy=0.0
	avgz=0.0

	OPEN(unit=24,file='com_avg_coor_2',status='unknown')
	DO i=1,n ! Average coordinates of each residue
	  avgx(i)=sumx(i)/REAL(ncount(2))
	  avgy(i)=sumy(i)/REAL(ncount(2))
	  avgz(i)=sumz(i)/REAL(ncount(2))
	  WRITE(*,*)avgx(i),avgy(i),avgz(i)
	  WRITE(24,*)avgx(i),avgy(i),avgz(i)
	END DO
	CLOSE(24)
	DO i=1,3*n,3
	  bigavglist(2,i)=avgx(i)
	  bigavglist(2,i+1)=avgy(i)
	  bigavglist(2,i+2)=avgz(i)
	END DO
	rxref=0.0
	ryref=0.0
	rzref=0.0
	DO k=1,ncount(2) !Subtract reference values (here the average) from each coordinate
	  DO i=1,n
	    rxref(i,k)=rx2(i,k)-avgx(i)
	    ryref(i,k)=ry2(i,k)-avgy(i)
	    rzref(i,k)=rz2(i,k)-avgz(i)
	  END DO
	END DO
	!Make big list of all augmented coordinates
	xcount=0
	ycount=0
	zcount=0
	WRITE(*,*)'Making the big list of coordinates...'
	OPEN(unit=24,file='com_biglist_2.dat',status='unknown')
	DO k=1,ncount(2)
	  DO i=1,3*n
	    IF(MOD(i,3) .EQ. 1)THEN
	      xcount=xcount+1
	      biglist(i,k)=rxref(xcount,k)
	    ELSE IF(MOD(i,3) .EQ. 2)THEN
	      ycount=ycount+1
	      biglist(i,k)=ryref(ycount,k)	
	    ELSE IF(MOD(i,3) .EQ. 0)THEN
	      zcount=zcount+1
	      biglist(i,k)=rzref(zcount,k)	 
	    END IF
	    WRITE(24,*)biglist(i,k)
	  END DO
	  !Reset counters after each time step
	  xcount=0
	  ycount=0
	  zcount=0
	END DO  
	CLOSE(24) 
	CLOSE(11)

	DO k=1,ncount(2)
	  DO i=1,3*n
	    DO j=1,3*n
	      umat2(i,j)=umat2(i,j)+biglist(i,k)*biglist(j,k)
	    END DO
	  END DO
	END DO
	DO i=1,3*n
	  DO j=1,3*n
	    umat2(i,j)=umat2(i,j)/REAL(ncount(2))
	    !WRITE(*,*)umat2(i,j)
	  END DO
	END DO
	!Weight the terms
	DO i=1,3*n
	  Do j=1,3*n
	    umat1(i,j)=p1*umat1(i,j)
	    umat2(i,j)=p2*umat2(i,j)
	  END DO
	END DO

	!Set up the trajectory for covariance analysis
	  !read from trajectory (.g96)
	  !open(unit=11,file=trim(protname)//'.g96',status='old')
	  !OPEN(unit=13,file='com_traj',status='unknown')
	  !skip first 7,now read and calculate stuff
	  !do i=1,7
	  !  read(11,*)
	  !end do
	  !sumx=0.0
	  !sumy=0.0
	  !sumz=0.0
	  !summa=0.0
	  !do k=1,nfrs
	  !  IF(MOD(k,500) .EQ. 0)WRITE(*,*)'Reading frame ',k
	  !  do j=1,n
	  !    read(11,*)rx(j,k),ry(j,k),rz(j,k)
	  !    sumx(j)=sumx(j)+rx(j,k)
	  !    sumy(j)=sumy(j)+ry(j,k)
	  !    sumz(j)=sumz(j)+rz(j,k)
	  !    rcomx(k)=rcomx(k)+rx(j,k)
	  !    rcomy(k)=rcomy(k)+ry(j,k)
	  !    rcomz(k)=rcomz(k)+rz(j,k)
	  !  end do
	  !  rcomx(k)=rcomx(k)/REAL(n)
	  !  rcomy(k)=rcomy(k)/REAL(n)
	  !  rcomz(k)=rcomz(k)/REAL(n)
	  !  WRITE(13,*)rcomx(k),rcomy(k),rcomz(k)
	  !  DO j=1,n !Calculate radial vectors wrt the com of the molecule, which should remove rotations
	  !    rx(j,k)=rx(j,k)-rcomx(k)
	  !    ry(j,k)=ry(j,k)-rcomy(k)
	  !    rz(j,k)=rz(j,k)-rcomz(k)
	  !  END DO
	    !do j=1,n-1
	    !  lx(j,k)=rx(j+1,k)-rx(j,k)
	    !  ly(j,k)=ry(j+1,k)-ry(j,k)
	    !  lz(j,k)=rz(j+1,k)-rz(j,k)
	    !end do
 
  	  !skip 8 lines
	  !  do j=1,8
	  !    read(11,*)
	  !  end do
	  !come out of time loop
	  !end do
	  !CLOSE(11)
	  !CLOSE(13)

	!Calculate averages and generate the covariance matrix
	!avgx=0.0
	!avgy=0.0
	!avgz=0.0

	!OPEN(unit=24,file='com_avg_coor',status='unknown')
	!DO i=1,n ! Average coordinates of each residue
	!  avgx(i)=sumx(i)/REAL(nfrs)
	!  avgy(i)=sumy(i)/REAL(nfrs)
	!  avgz(i)=sumz(i)/REAL(nfrs)
	!  WRITE(*,*)avgx(i),avgy(i),avgz(i)
	!  WRITE(24,*)avgx(i),avgy(i),avgz(i)
	!END DO
	!CLOSE(24)
	!rxref=0.0
	!ryref=0.0
	!rzref=0.0
	!DO k=1,nfrs !Subtract reference values (here the average) from each coordinate
	!  DO i=1,n
	!    rxref(i,k)=rx(i,k)-avgx(i)
	!    ryref(i,k)=ry(i,k)-avgy(i)
	!    rzref(i,k)=rz(i,k)-avgz(i)
	!  END DO
	!END DO
	!Make big list of all augmented coordinates
	!xcount=0
	!ycount=0
	!zcount=0
	!WRITE(*,*)'Making the big list of coordinates...'
	!OPEN(unit=24,file='com_biglist.dat',status='unknown')
	!DO k=1,nfrs
	!  DO i=1,3*n
	!    IF(MOD(i,3) .EQ. 1)THEN
	!      xcount=xcount+1
	!      biglist(i,k)=rxref(xcount,k)
	!    ELSE IF(MOD(i,3) .EQ. 2)THEN
	!      ycount=ycount+1
	!      biglist(i,k)=ryref(ycount,k)	
	!    ELSE IF(MOD(i,3) .EQ. 0)THEN
	!      zcount=zcount+1
	!      biglist(i,k)=rzref(zcount,k)	 
	!    END IF
	!    WRITE(24,*)biglist(i,k)
	!  END DO
	!  !Reset counters after each time step
	!  xcount=0
	!  ycount=0
	!  zcount=0
	!END DO  
	!CLOSE(24) 
	!CLOSE(11)

	OPEN(unit=13,file='C1.dat',status='unknown')
	OPEN(unit=14,file='C2.dat',status='unknown')
	DO i=1,3*n
	  DO j=1,3*n
	    WRITE(13,*)umat1(i,j)
	    WRITE(14,*)umat2(i,j)
	    WRITE(*,*)umat1(i,j),umat2(i,j)
	  END DO
	END DO
	CLOSE(13)
	CLOSE(14)

	!Now, let's move on to the first term on the RHS:

	sumx=0.0
	sumy=0.0
	sumz=0.0
	rxref=0.0
	ryref=0.0
	rzref=0.0
	rcomx=0.0
	rcomy=0.0
	rcomz=0.0
	DO k=1,nfrs
	    do j=1,n
	      sumx(j)=sumx(j)+rx(j,k)
	      sumy(j)=sumy(j)+ry(j,k)
	      sumz(j)=sumz(j)+rz(j,k)
	      rcomx(k)=rcomx(k)+rx(j,k)
	      rcomy(k)=rcomy(k)+ry(j,k)
	      rcomz(k)=rcomz(k)+rz(j,k)
	    end do
	    rcomx(k)=rcomx(k)/REAL(n)
	    rcomy(k)=rcomy(k)/REAL(n)
	    rcomz(k)=rcomz(k)/REAL(n)
	    DO j=1,n !Move to the centre of the box
	      rx2(j,k)=rx(j,k)-rcomx(k)
	      ry2(j,k)=ry(j,k)-rcomy(k)
	      rz2(j,k)=rz(j,k)-rcomz(k)
	    END DO
	END DO

	!Calculate averages and generate the covariance matrix
	avgx=0.0
	avgy=0.0
	avgz=0.0
	OPEN(unit=24,file='com_avg_coor',status='unknown')
	DO i=1,n !Average coordinates of each residue
	  avgx(i)=sumx(i)/REAL(nfrs)
	  avgy(i)=sumy(i)/REAL(nfrs)
	  avgz(i)=sumz(i)/REAL(nfrs)
	  WRITE(*,*)avgx(i),avgy(i),avgz(i)
	  WRITE(24,*)avgx(i),avgy(i),avgz(i)
	END DO
	DO i=1,3*n,3
	  avglist(i)=avgx(i)
	  avglist(i+1)=avgy(i)
	  avglist(i+2)=avgz(i)
	END DO
	CLOSE(24)
	rxref=0.0
	ryref=0.0
	rzref=0.0
	DO k=1,nfrs !Subtract reference values (here the average) from each coordinate
	  DO i=1,n
	    rxref(i,k)=rx(i,k)-avgx(i)
	    ryref(i,k)=ry(i,k)-avgy(i)
	    rzref(i,k)=rz(i,k)-avgz(i)
	  END DO
	END DO
	!Make big list of all augmented coordinates
	xcount=0
	ycount=0
	zcount=0
	!Don't need to make this biglist right now
	!WRITE(*,*)'Making the big list of coordinates...'
	!OPEN(unit=24,file='com_biglist.dat',status='unknown')
	!DO k=1,nfrs
	!  DO i=1,3*n
	!    IF(MOD(i,3) .EQ. 1)THEN
	!      xcount=xcount+1
	!      biglist(i,k)=rxref(xcount,k)
	!    ELSE IF(MOD(i,3) .EQ. 2)THEN
	!      ycount=ycount+1
	!      biglist(i,k)=ryref(ycount,k)	
	!    ELSE IF(MOD(i,3) .EQ. 0)THEN
	!      zcount=zcount+1
	!      biglist(i,k)=rzref(zcount,k)	 
	!    END IF
	!    WRITE(24,*)biglist(i,k)
	!  END DO
	!  !Reset counters after each time step
	!  xcount=0
	!  ycount=0
	!  zcount=0
	!END DO  
	!CLOSE(24) 
	CLOSE(11)
	!Now the contribution from jumping among minima
	DO a=1,nstates
	  !DO k=1,ncount(a)
	    DO i=1,3*n
	      DO j=1,3*n
	        cjam(a,i,j)=cjam(a,i,j)+(bigavglist(a,i)-avglist(i))*(bigavglist(a,j)-avglist(j))
	      END DO
	    END DO
	  !END DO
	  DO i=1,3*n
	    DO j=1,3*n
	      cjam(a,i,j)=fk(a)*cjam(a,i,j)
	      rhs1(i,j)=rhs1(i,j)+cjam(a,i,j) !Build the weighted sum
	    END DO
	  END DO
	END DO
	OPEN(unit=24,file='CJAM.dat',status='unknown')
	DO i=1,3*n
	  DO j=1,3*n
	    WRITE(24,*)rhs1(i,j)
	  END DO
	END DO
	   	
	end subroutine

	subroutine writer(n,nfrs,nstates,protname)
	integer i,j,k,imin,jmin,a,nfe,ir,nbins,rot,imad,imadmid,nmol,nstates
	real, dimension(n) :: rx,ry,rz,lavm,sigfe,fricorr,pvol
	real, dimension(n) :: lx,ly,lz,lmag,avfe,avfesq,fenorm
	real, dimension(n,n) :: sigij,rij,qinvm,qm
	real, dimension(n,nfrs) :: xix,xiy,xiz,dipcorr,xim,theta,phi
	real dotpij,um,rrij,bl,hrtheta,hrphi,Rb,T,r,dr,delphi
	integer itheta,iphi
	character(32)protname
	character(16)aa,ii,cbins,bb
	real hisp,hismax,delha,rdeg,degr,hnorm(n),x,y,z,emad
	real femax,pi,delr
        nbins=60
	nmol=1
	Rb=.00198 !(boltzmanns constant in kcal/mol*K)
	open(unit=10,file='temp')
	read(10,*)T
	close(10)
	felim=0.0
	femin=0.0
	feminp=0.0
	sigfe=0.0
	fricorr=0.0
	avfe=0.0
	avfesq=0.0
	rij=0.0
	hisp=100.0
	hismax=0.0
	xix=0.0
	xiy=0.0
	xiz=0.0
	xim=0.0
	qinvm=0.0
	qm=0.0
	dipcorr=0.0	
	xim=0.0
	theta=0.0
	phi=0.0
	hisang=0.0
	hrtheta=0.0
	hrphi=0.0
	itheta=0
	iphi=0
	pi=3.1415927
	delha=(2.0*360.0)/real(2*nbins)
	degr=((2.0*pi)/360.0) !deg to rad
	rdeg=1.0/degr !rad to deg
	pvol=0.0
	fenorm=0.0
	r=0.0
	dr=5.0/(real(nfrs))
	ir=0
	dc=1.0/real(nfrs)
	fmax=0.0
	fmaxp=0.0
	fmin=0.0
	fminp=200.0
	fbg=0.0
	fmad=0.0
	fmadord=0.0

	delr=delha*degr
	write(*,*)"delr",delr
	DO ir=1,2 !Two states
        WRITE(aa,*)ir
        aa=adjustl(aa)
	hnorm=0.0
        qinvm=0.0
	open(unit=21,file='QINVmatrix'//TRIM(aa),status='old')
	do i=1,n-nmol
	  do j=1,n-nmol
	    read(21,*)qinvm(i,j)
	  end do
	end do
	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)
	rx=0.0
	ry=0.0
	rz=0.0
	lx=0.0
	ly=0.0
	lz=0.0
	lmag=0.0
	lavm=0.0
	dotpij=0.0
	sigij=0.0
	um=0.0
	rrij=0.0
	rij=0.0
	bl=0.0
	imin=0
	jmin=0

	!read from trajectory
	do a=1,10!n-nmol !mode loop; only over the first ten modes for now
	  WRITE(*,*)'Mode ',a
	  WRITE(bb,*)a
	  bb=adjustl(bb)
	  open(unit=11,file='traj'//trim(aa)//'.dat',status='old')
	  open(unit=12,file='nframes'//trim(aa),status='old')
          READ(12,*)nfrs
          WRITE(*,*)nfrs
          CLOSE(12)
	  OPEN(unit=100+a,file='theta_phi_'//TRIM(aa)//'_'//TRIM(bb)//'.dat',status='unknown')
	  do k=1,nfrs
	    !IF(MOD(k,500) .EQ. 0)WRITE(*,*)'Reading frame number ',k
	    do j=1,n-1
	      read(11,*)lx(j),ly(j),lz(j)
	    end do
	    do j=1,n-nmol !residue loop
	      xix(a,k)=qinvm(a,j)*lx(j)+xix(a,k)
	      xiy(a,k)=qinvm(a,j)*ly(j)+xiy(a,k)
	      xiz(a,k)=qinvm(a,j)*lz(j)+xiz(a,k)
	    end do
	    xim(a,k)=(xix(a,k)**2+xiy(a,k)**2+xiz(a,k)**2)**.5
	    !calculate theta, phi
	    theta(a,k)=acos(xiz(a,k)/xim(a,k))
	    phi(a,k)=atan(xiy(a,k)/xix(a,k))
            if(xix(a,k).lt.0.0)phi(a,k)=phi(a,k)+pi
            theta(a,k)=theta(a,k)*rdeg
            if(phi(a,k).lt.0.0)phi(a,k)=phi(a,k)+2.0*pi
            phi(a,k)=phi(a,k)*rdeg
	    WRITE(100+a,*)theta(a,k),phi(a,k),xim(a,k) !Write the mode magnitude, just in case
	  end do !End loop over time
	CLOSE(100+a)
	CLOSE(200+a)
	CLOSE(11)
	end do !End loop over modes
        end do !End loop over states

        end subroutine
