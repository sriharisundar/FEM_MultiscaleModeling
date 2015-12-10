!PROGRAM readmesh
MODULE FUNCTIONS
	
	TYPE :: node
		integer :: number
		double precision :: x,y,dispx,dispy,fx,fy
	END Type node

	TYPE :: element
		integer :: number
	  double precision,pointer :: node1(:),node2(:),node3(:),node4(:)
	  double precision,dimension(8,8):: Kele
    double precision :: stress(4,3),strain(4,3)
	END TYPE element

CONTAINS

 SUBROUTINE readinput(filename,nnode,nele,nodes,elements,nfixed,nforce,ndisp,nodes_fixed,nodes_force,nodes_disp)
  IMPLICIT NONE
	character(100):: filename,c
 	integer, INTENT(OUT):: nnode,nele,nfixed,nforce,ndisp
	integer,allocatable :: nodes_fixed(:),nodes_force(:),nodes_disp(:)
 	integer :: num,eof,I,ndof,p,q
 	integer,dimension(4) :: nodenum
 	double precision :: x,y,fx,fy,dx,dy
 	double precision,allocatable,target,INTENT(OUT) :: nodes(:,:)	!num,x,y,fx,fy,dx,dy,xfix,yfix
	type(element),allocatable,INTENT(OUT) :: elements(:)

	open(unit=20,file=filename,action='read')
	
	read(20,*) nnode , nele

	ndof=nnode*2

	ALLOCATE(nodes(nnode,9), elements(nele))

	read(20,*) 

	DO I= 1,nnode
		read(20,*,IOSTAT=eof) num,x,y
		nodes(I,1)=num
		nodes(I,2)=x
		nodes(I,3)=y
!	   write(*,*) nodes(I,:)
  END DO

	DO I=1,nele
		read(20,*) num,nodenum
    elements(I)%number=num
		elements(I)%node1=>nodes(nodenum(1),:)
		elements(I)%node2=>nodes(nodenum(2),:)
		elements(I)%node3=>nodes(nodenum(3),:)
		elements(I)%node4=>nodes(nodenum(4),:)
	END DO

	read(20,*) nfixed
  ALLOCATE(nodes_fixed(nfixed))

	DO I=1,nfixed
		read(20,*) nodes_fixed(I),p,q
    nodes(nodes_fixed(I),8:9)= (/p,q/)
	END DO
  !write(*,*) 'nodes_disp'

	read(20,*) nforce
  ALLOCATE(nodes_force(nforce))

	DO I=1,nforce
		read(20,*) nodes_force(I),fx,fy
		nodes(nodes_force(I),4)=fx
		nodes(nodes_force(I),5)=fy
	END DO

	read(20,*) ndisp
  ALLOCATE(nodes_disp(ndisp))

	DO I=1,ndisp
    read(20,*) nodes_disp(I),dx,dy
    !write(*,*) dx,dy
		nodes(nodes_disp(I),6)=dx
		nodes(nodes_disp(I),7)=dy
    !write(*,*) nodes(nodes_force(I),6)
	END DO

	close(20)
  RETURN
 END SUBROUTINE readinput

 SUBROUTINE elestiff(ele)
 	    IMPLICIT NONE
 	    TYPE(element) :: ele
      DOUBLE PRECISION::GMAT(4,8),AMAT(3,4),S(4),T(4),J(4)
      DOUBLE PRECISION:: B(3,8), BT(8,3),D(3,3),K(8,8)
      DOUBLE PRECISION:: KI(8,8), GP(4,4)
      DOUBLE PRECISION:: GP_W(4), C(3,8), FORCE(8,1), DSTRAN(3,1), DSTRESS(3,1)
      DOUBLE PRECISION:: E,NU,DETJ
      DOUBLE PRECISION:: COORDS(4,2),DISP(8)
      INTEGER:: I,J1,K1,NGPT,NDOFEL,NNODE
      INTEGER,PARAMETER:: FOUR=4.D0

      E=210
      NU=0.3
      NGPT=FOUR
      NNODE=FOUR
      NDOFEL=NNODE*2
!	    OPEN(UNIT=22,FILE= 'debug')

      CALL INITIALISE_COORDS (COORDS,ele)
      CALL INITIALISE_DISP (DISP,ele)
      CALL INITIALISE_GP (GP,GP_W,4,2)
      CALL INITIALISE_D (D,E,NU,3,3)
!      CALL INITIALISE_ZERO(RHS,NDOFEL,1)
!
	  !write(*,*) ele%number
	  !write(*,*) COORDS
      DO I=1,NGPT
! CALCULATE S AND T - PASS XHI,ETA FOR INFO
        CALL CALCULATE_ST(S,T,GP(I,1),GP(I,2),NNODE)
!        WRITE(121,*) I
!        WRITE(121,*) S(:NNODE)
!        WRITE(121,*) T(:NNODE)
! CALCULATE GMAT
        CALL CALCULATE_G(GMAT,S,T,NDOFEL,NNODE)
!        WRITE(120,*) I
!        DO J1=1,4
!          WRITE(120,*) GMAT(J1,:16)
!        END DO
! CALCULATE AMAT
        CALL CALCULATE_A(COORDS,AMAT,S,T,NNODE,DETJ)
!        WRITE(120,*) I d
!        DO J1=1,3
!          WRITE(*,*) AMAT(J1,:4)
!        END DO
! MULTIPLY G AND A TO GET B
        CALL MATRIXMULTIPLY(B,AMAT,GMAT,3,NDOFEL,4)
! MULTIPLY D AND B
        CALL MATRIXMULTIPLY(C,D,B,3,NDOFEL,3)
! CALCULATE B TRANSPOSE
        CALL MATRIXTRANSPOSE(B,BT,3,NDOFEL)
! MULTIPLY BT AND C
        CALL MATRIXMULTIPLY(KI,BT,C,NDOFEL,NDOFEL,3)
!
!        DO J1 = 1 , 3
!          STRAN(J1) = SVARS(J1+(I-1)*3)
!          STRESS(J1,1)= SVARS(12+J1+(I-1)*3)          
!        END DO
        CALL MATRIXMULTIPLY(DSTRAN,B,DISP,3,1,NDOFEL)
        ele%strain(I,:)=DSTRAN(:,1)
!        DO J1=1,3
!          STRAN(J1) = STRAN(J1) + DSTRAN(J1,1)
!        END DO
        CALL MATRIXMULTIPLY(DSTRESS,D,DSTRAN,3,1,3)
        ele%stress(I,:)=DSTRESS(:,1)
! 
!        DO J1=1,3
!          STRESS(J1,1) = STRESS(J1,1) + DSTRESS(J1,1)
!        END DO
!        CALL MATRIXMULTIPLY(FORCE,BT,DSTRESS,NDOFEL,1,3)
!        do J1=1,8
!      	  	write(22,*) ( KI(J1,K1), K1=1,8 )
!   	   enddo

        DO J1=1,NDOFEL
          DO K1=1,NDOFEL
            ele%Kele(J1,K1) = ele%Kele(J1,K1) + KI(J1,K1)*GP_W(I)*DETJ
          END DO 
        END DO

!        DO J1=1,3
!          SVARS(J1+(I-1)*3) = STRAN(J1)
!          SVARS(12+J1+(I-1)*3) = STRESS(J1,1)
!        END DO
      END DO
!        WRITE(120,*) JELEM,KINC
!        WRITE(120,*) 'STRAIN'
!        WRITE(120,*) SVARS(1:12)
!        WRITE(120,*) 'STRESS'
!        WRITE(120,*) SVARS(13:24)
!      END IF
      END

 SUBROUTINE INITIALISE_COORDS(COORDS,ele)
  TYPE(element) :: ele
  DOUBLE PRECISION:: COORDS(4,2)
  COORDS(1,1)=ele%node1(2)
  COORDS(1,2)=ele%node1(3)
  COORDS(2,1)=ele%node2(2)
  COORDS(2,2)=ele%node2(3)
  COORDS(3,1)=ele%node3(2)
  COORDS(3,2)=ele%node3(3)
  COORDS(4,1)=ele%node4(2)
  COORDS(4,2)=ele%node4(3)
  RETURN
 END

 SUBROUTINE INITIALISE_DISP(DISP,ele)
  TYPE(element) :: ele
  DOUBLE PRECISION:: DISP(8)
  DISP(1)=ele%node1(6)
  DISP(2)=ele%node1(7)
  DISP(3)=ele%node2(6)
  DISP(4)=ele%node2(7)
  DISP(5)=ele%node3(6)
  DISP(6)=ele%node3(7)
  DISP(7)=ele%node4(6)
  DISP(8)=ele%node4(7)
  RETURN
 END

 SUBROUTINE INITIALISE_ZERO(A,n,m)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(n,m)
 INTEGER :: n,m,i,j
 DO i = 1, n
  DO j = 1, m
     A(i,j) = 0.0
  END DO
 END DO
 RETURN
 END

! REFER ENWISTLE PG 108-109

 SUBROUTINE INITIALISE_GP(A,B,n,m)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(n,m),B(n)
 INTEGER :: n,m      
 CALL INITIALISE_ZERO(A,n,m)
 A(1,1)=-0.577350269189626
 A(1,2)=-0.577350269189626
 A(2,1)=0.577350269189626
 A(2,2)=-0.577350269189626
 A(3,1)=0.577350269189626
 A(3,2)=0.577350269189626
 A(4,1)=-0.577350269189626
 A(4,2)=0.577350269189626
 B=(/1,1,1,1/)
 RETURN
 END

 SUBROUTINE INITIALISE_D(A,E,NU,n,m)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(n,m),E,NU,E1,E2
 INTEGER :: n,m
 CALL INITIALISE_ZERO(A,n,m)
 E1 = E/(1-NU**2)
 A(1,1)=E1
 A(2,2)=E1
 A(1,2)=E1*NU
 A(2,1)=A(1,2)
 A(3,3)=E1*(1-NU)/2
 RETURN
 END
!
! N1=0.25*(1-X)(1-E)
! N2=0.25*(1+X)(1-E)
! N3=0.25*(1+X)(1+E)
! N4=0.25*(1-X)(1+E)
!
! S ARE XHI DERIVATIVE, T ARE XHI DERIVATIVE
 
 SUBROUTINE CALCULATE_ST(A,B,X,E,n)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(n),B(n),X,E
 INTEGER :: n
 A(1)=-1*0.25*(1-E)
 A(2)=0.25*(1-E)
 A(3)=0.25*(1+E)
 A(4)=-1*0.25*(1+E)
 B(1)=-1*0.25*(1-X)
 B(2)=-1*0.25*(1+X)
 B(3)=0.25*(1+X)
 B(4)=0.25*(1-X)
 RETURN
 END
!
!
!
 
 SUBROUTINE CALCULATE_G(A,B,C,n,m)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(4,n),B(m),C(m)
 INTEGER :: n,m,i
 CALL INITIALISE_ZERO(A,4,n)
 DO i=1,n,2
   A(1,i)=B((i+1)/2)
   A(2,i)=C((i+1)/2)
 END DO
 DO i=2,n,2
   A(3,i)=B(i/2)
   A(4,i)=C(i/2)
 END DO
 RETURN
 END
!
!
! J11=A(1) J12=A(2) J21=A(3) J(22)=A(4)
 
 SUBROUTINE CALCULATE_A(D,A,B,C,n,DETJ)
 IMPLICIT NONE
 DOUBLE PRECISION :: D(n,2),A(3,4),B(n),C(n),J(4),DETJ
 INTEGER :: n,I,J1
 CALL INITIALISE_ZERO(A,3,4)
! DO J1=1,2
!     WRITE(121,*) D(J1,:8)
! END DO
!   WRITE(121,*) B(:8)
!   WRITE(121,*) C(:8)
 DO I=1,n
   J(1)=J(1)+B(I)*D(I,1)
   J(2)=J(2)+B(I)*D(I,2)
   J(3)=J(3)+C(I)*D(I,1)
   J(4)=J(4)+C(I)*D(I,2)
 END DO
!
 DETJ=J(1)*J(4)-J(2)*J(3)

!DO I=1,2
!WRITE(*,*) D(:,I)   
!WRITE(*,*) '' 
!END DO 
!WRITE(*,*) B    
!WRITE(*,*) ''  
!WRITE(*,*) C    
!WRITE(*,*) ''  
!WRITE(*,*) J    
!WRITE(*,*) ''  
!WRITE(*,*) DETJ    
!WRITE(*,*) ''  

!
 A(1,1)=J(4)/DETJ
 A(1,2)=-1*J(2)/DETJ
 A(2,3)=-1*J(3)/DETJ
 A(2,4)=J(1)/DETJ
 A(3,1)=-1*J(3)/DETJ
 A(3,2)=J(1)/DETJ
 A(3,3)=J(4)/DETJ
 A(3,4)=-1*J(2)/DETJ
 RETURN
 END
!
!
 
 SUBROUTINE MATRIXMULTIPLY(A,B,C,n,m,l)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(n,m),B(n,l),C(l,m)
 INTEGER :: n,m,l,i,j,k
 CALL INITIALISE_ZERO(A,n,m)
 DO i = 1, n
   DO j = 1, m
     DO k = 1, l
       A(i,j) = A(i,j) + B(i,k) * C(k,j)
     END DO
   END DO
 END DO
 RETURN
 END
!
!
 
 SUBROUTINE MATRIXTRANSPOSE(A,B,n,m)
 IMPLICIT NONE
 DOUBLE PRECISION :: A(n,m), B(m,n)
 INTEGER :: i,j,n,m
 CALL INITIALISE_ZERO(B,m,n)
 do i = 1, n
    do j = 1, m
       B(j,i) = A(i,j)
    end do
 end do
 RETURN
 END
 

SUBROUTINE assemble(globalstiff,force,nodes,elements,nnode,nele,nfixed,nforce,ndisp,nodes_fixed,nodes_force,nodes_disp)
  
  IMPLICIT NONE
  character(100):: filename
  integer, INTENT(IN):: nnode,nele,nfixed,nforce,ndisp
  integer :: num,eof,I,J,K,node1,node2,ndof,n
  integer,dimension(4) :: nodenum
  double precision :: x,y
  double precision,allocatable :: force(:)
  double precision,INTENT(IN) :: nodes(:,:)
  integer,allocatable :: nodes_fixed(:),nodes_force(:),nodes_disp(:)
  type(element),allocatable,INTENT(IN) :: elements(:)
  double precision,allocatable,INTENT(INOUT) :: globalstiff(:,:)

  ndof=nnode*2
  globalstiff=0
  force=0

  DO I=1,nele

    nodenum(1)=elements(I)%node1(1)
    nodenum(2)=elements(I)%node2(1)
    nodenum(3)=elements(I)%node3(1)
    nodenum(4)=elements(I)%node4(1)

!    write(*,*) nodenum

    DO J=1,8,2
      DO K=1,8,2
        node1=nodenum((J+1)/2)
        node2=nodenum((K+1)/2)
        !write(*,*) node1,node2
        globalstiff(node1*2-1:node1*2,node2*2-1:node2*2)=globalstiff(node1*2-1:node1*2,node2*2-1:node2*2) &
        	+elements(I)%Kele(J:J+1,K:K+1)
      END DO
    END DO
  END DO

  DO I=1,nforce
    force((nodes_force(I)*2)-1)=nodes(nodes_force(I),4)
    force((nodes_force(I)*2))=nodes(nodes_force(I),5)
  END DO

  DO I=1,ndisp    
    n=nodes_disp(I)
    force((n*2)-1)=nodes(n,6)*(1e20+globalstiff((n*2)-1,(n*2)-1))
    force((n*2))=nodes(n,7)*(1e20+globalstiff((n*2),(n*2)))
    globalstiff((n*2)-1,(n*2)-1)=1e20+globalstiff((n*2)-1,(n*2)-1)
    globalstiff((n*2),(n*2))=1e20+globalstiff((n*2),(n*2))
  END DO

  RETURN
 END SUBROUTINE

 SUBROUTINE findnewcoords(nodes,nnode,disp)
  double precision,INTENT(INOUT) :: nodes(:,:),disp(:)
  integer :: nnode,I

  DO I=1,nnode
    nodes(I,2)=nodes(I,2)+disp((I*2)-1)
    nodes(I,3)=nodes(I,3)+disp((I*2))
    nodes(I,6)=disp(I*2-1)
    nodes(I,7)=disp(I*2)
  END DO

  END SUBROUTINE
END MODULE FUNCTIONS
