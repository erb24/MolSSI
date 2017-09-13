      program inputreader
	  IMPLICIT NONE
	  integer nfrs,n
      character(32) :: protname
	  open(unit=5,file="protname.txt",status='old')
	  read(5,'(A)')protname
	  read(5,*)n
	  read(5,*)nfrs
	  close(5)
	  !nfrs=1000 !Testing
	  write(*,*)n,nfrs
	  call contact(n,nfrs,protname)
	  STOP
	  End program inputreader

      subroutine contact(n,nfrs,protname)
      implicit none
      integer :: n,nfrs,i,j,k,nres,counter
      DOUBLE PRECISION :: cmap(n,n),rx,ry,rz,r(n,nfrs),cut
      Character(32) :: protname
      protname=adjustl(protname)
      cmap=0.0
      r=0.0
      open(unit=24,file=TRIM(protname)//'.g96',status='old')
      DO i=1,7
        READ(24,*)
      END DO
      DO k=1,nfrs
        IF(MOD(k,500) .EQ. 0)WRITE(*,*)'Reading frame ',k
        DO i=1,n
          READ(24,*)rx,ry,rz
          r(i,k)=SQRT(rx**2+ry**2+rz**2)
        END DO
        DO i=1,8
          READ(24,*)
        END DO
      END DO
      CLOSE(24)
      DO k=1,nfrs
        IF(MOD(k,500) .EQ. 0)WRITE(*,*)'Writing frame ',k
        DO i=1,n
          DO j=1,n
            cut=ABS(r(i,k)-r(j,k))
            IF(cut .LE. 1.0)cmap(i,j)=cmap(i,j)+1.0 !Cutoff distance of 10 Angstroms
            !cmap(i,j)=cmap(i,j)+(r(i,k)-r(j,k))
          END DO
        END DO
      END DO
      !Normalize
      open(unit=24,file='cmap.dat',status='unknown')  
      DO i=1,n
        DO j=1,n
          cmap(i,j)=cmap(i,j)/REAL(nfrs)
          WRITE(24,*)cmap(i,j)
        END DO
      END DO
      CLOSE(24)
      END SUBROUTINE contact
