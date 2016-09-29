C 8 NODE 4 GAUSS POINT ELEMENT BY SSH
      SUBROUTINE UEL (RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     & PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     & DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     & PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     & NJPRO, PERIOD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*),
     & SVARS(*), ENERGY(8), COORDS(MCRD, NNODE), U(NDOFEL),
     & DU(MLVARX,*), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(*),
     & JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*), DDLMAG(MDLOAD,*),
     & PREDEF(2, NPREDF, NNODE), LFLAGS(*), JPROPS(*)
      PARAMETER(FOUR=4.D0)
C
      DOUBLE PRECISION GMAT(4,NDOFEL),AMAT(3,4),S(NNODE),T(NNODE),J(4),
     & B(3,NDOFEL), BT(NDOFEL,3),D(3,3),K(NDOFEL,NDOFEL),
     & KI(NDOFEL,NDOFEL), STRAN(3), STRESS(3,1), GP(4,MCRD),
     & GP_W(4), C(3,NDOFEL), FORCE(NDOFEL,1), DSTRAN(3,1), DSTRESS(3,1)
C
C
C
      CHARACTER(512) FILENAME1,FILENAME2
      CHARACTER(512) JOBDIR      
      DOUBLE PRECISION E,NU
      E=PROPS(1)
      NU=PROPS(2)
      NGPT=FOUR
      CALL K_INITIALISE_ZERO(AMATRX,16,16)
      CALL K_INITIALISE_ZERO(RHS,16,1)
      CALL K_INITIALISE_GP (GP,GP_W,4,2)
      CALL K_INITIALISE_D (D,E,NU,3,3)
C
C
      CALL GETOUTDIR(JOBDIR,LENJOBDIR) 
      FILENAME1=JOBDIR(:LENJOBDIR)//'/uel_2ele.120'     
      FILENAME2=JOBDIR(:LENJOBDIR)//'/uel_2ele.121'
      IF (JELEM .eq. 1 .and. KINC==1) THEN
        OPEN(120,file=FILENAME1)
C        OPEN(121,file=FILENAME2)
      ELSE
        OPEN(120,file=FILENAME1,Access = 'append')
C        OPEN(121,file=FILENAME2,Access = 'append')
      END IF
C
C      WRITE(121,*) JELEM,KINC
C      WRITE(121,*) U(:NDOFEL)
      DO I=1,NGPT
C CALCULATE S AND T - PASS XHI,ETA FOR INFO
        CALL K_CALCULATE_ST(S,T,GP(I,1),GP(I,2),NNODE)
c        WRITE(121,*) I
c        WRITE(121,*) S(:NNODE)
c        WRITE(121,*) T(:NNODE)
C CALCULATE GMAT
        CALL K_CALCULATE_G(GMAT,S,T,NDOFEL,NNODE)
c        WRITE(120,*) I
C        DO J1=1,4
C          WRITE(120,*) GMAT(J1,:16)
C        END DO
C CALCULATE AMAT
C        WRITE(120,*) MCRD    
C        DO J1=1,2
C            WRITE(120,*) COORDS(J1,:8)
C        END DO
        CALL K_CALCULATE_A(COORDS,AMAT,S,T,NNODE,MCRD)
C        WRITE(120,*) I d
C        DO J1=1,3
C          WRITE(120,*) AMAT(J1,:4)
C        END DO
C MULTIPLY G AND A TO GET B
        CALL K_MATRIXMULTIPLY(B,AMAT,GMAT,3,NDOFEL,4)
C MULTIPLY D AND B
        CALL K_MATRIXMULTIPLY(C,D,B,3,NDOFEL,3)
C CALCULATE B TRANSPOSE
        CALL K_MATRIXTRANSPOSE(B,BT,3,NDOFEL)
C MULTIPLY BT AND C
        CALL K_MATRIXMULTIPLY(KI,BT,C,NDOFEL,NDOFEL,3)
C
        DO J1 = 1 , 3
          STRAN(J1) = SVARS(J1+(I-1)*3)
          STRESS(J1,1)= SVARS(12+J1+(I-1)*3)          
        END DO
        CALL K_MATRIXMULTIPLY(DSTRAN,B,DU,3,1,NDOFEL)
        DO J1=1,3
          STRAN(J1) = STRAN(J1) + DSTRAN(J1,1)
        END DO
C        IF (KINC .eq. 1) THEN
C          WRITE(*,*) D(:3,:3), DSTRAN(:3,1)
C        END IF
        CALL K_MATRIXMULTIPLY(DSTRESS,D,DSTRAN,3,1,3)
C        IF (KINC .eq. 1) WRITE(*,*) DSTRESS(:3,1)
        DO J1=1,3
          STRESS(J1,1) = STRESS(J1,1) + DSTRESS(J1,1)
        END DO
        CALL K_MATRIXMULTIPLY(FORCE,BT,DSTRESS,NDOFEL,1,3)
        DO J1=1,NDOFEL
          RHS(J1,1) = RHS(J1,1) - FORCE(J1,1)*GP_W(I)
          DO K1=1,NDOFEL
            AMATRX(J1,K1) = AMATRX(J1,K1) + KI(J1,K1)*GP_W(I)
          END DO 
        END DO
        DO J1=1,3
          SVARS(J1+(I-1)*3) = STRAN(J1)
          SVARS(12+J1+(I-1)*3) = STRESS(J1,1)
        END DO
      END DO
C      IF (.not. mod(KINC,100)) THEN
        WRITE(120,*) JELEM,KINC
        WRITE(120,*) 'STRAIN'
        WRITE(120,*) SVARS(1:12)
        WRITE(120,*) 'STRESS'
        WRITE(120,*) SVARS(13:24)
C      END IF
      close(120)
      close(121) 
      END
C 
      SUBROUTINE K_INITIALISE_ZERO(A,n,m)
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
C
C REFER ENWISTLE PG 136-137
C
      SUBROUTINE K_INITIALISE_GP(A,B,n,m)
      IMPLICIT NONE
      DOUBLE PRECISION :: A(n,m),B(n)
      INTEGER :: n,m      
      CALL K_INITIALISE_ZERO(A,n,m)
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
C
      SUBROUTINE K_INITIALISE_D(A,E,NU,n,m)
      IMPLICIT NONE
      DOUBLE PRECISION :: A(n,m),E,NU,E1,E2
      INTEGER :: n,m
      CALL K_INITIALISE_ZERO(A,n,m)
      E1 = E/(1-NU**2)
      A(1,1)=E1
      A(2,2)=E1
      A(1,2)=E1*NU
      A(2,1)=A(1,2)
      A(3,3)=E1*(1-NU)/2
      RETURN
      END
C
C N1=-(1-XHI)(1-ETA)(1+XHI+ETA)/4
C N2=-(1+XHI)(1-ETA)(1-XHI+ETA)/4
C N3=-(1+XHI)(1+ETA)(1-XHI-ETA)/4
C N4=-(1-XHI)(1+ETA)(1+XHI-ETA)/4
C N5=(1-XHI^2)(1-ETA)/2
C N6=(1+XHI)(1-ETA^2)/2
C N7=(1-XHI^2)(1+ETA)/2
C N8=(1-XHI)(1-ETA^2)/2
C
C S ARE XHI DERIVATIVE, T ARE XHI DERIVATIVE
      SUBROUTINE K_CALCULATE_ST(A,B,X,E,n)
      IMPLICIT NONE
      DOUBLE PRECISION :: A(n),B(n),X,E
      INTEGER :: n
      A(1)=0.25*(1-E)*(2*X+E)
      A(2)=0.25*(1-E)*(2*X-E)
      A(3)=0.25*(1+E)*(2*X+E)
      A(4)=0.25*(1+E)*(2*X-E)
      A(5)=(-1)*X*(1-E)
      A(6)=(0.5)*(1-E**2)
      A(7)=(-1)*X*(1+E)
      A(8)=(-0.5)*(1-E**2)
      B(1)=0.25*(1-X)*(2*E+X)
      B(2)=0.25*(1+X)*(2*E-X)
      B(3)=0.25*(1+X)*(2*E+X)
      B(4)=0.25*(1-X)*(2*E-X)
      B(5)=(-0.5)*(1-X**2)
      B(6)=(-1)*E*(1+X)
      B(7)=(0.5)*(1-X**2)
      B(8)=(-1)*E*(1-X)
      RETURN
      END
C
C
C
      SUBROUTINE K_CALCULATE_G(A,B,C,n,m)
      IMPLICIT NONE
      DOUBLE PRECISION :: A(4,n),B(m),C(m)
      INTEGER :: n,m,i
      CALL K_INITIALISE_ZERO(A,4,n)
      DO i=1,n,2
        A(1,i)=B((i+1)/2)
        A(2,i)=C((i+1)/2)
      END DO
C
      DO i=2,n,2
        A(3,i)=B(i/2)
        A(4,i)=C(i/2)
      END DO
      RETURN
      END
C
C
C J11=A(1) J12=A(2) J21=A(3) J(22)=A(4)
      SUBROUTINE K_CALCULATE_A(D,A,B,C,n,m)
      IMPLICIT NONE
      DOUBLE PRECISION :: D(2,8),A(3,4),B(n),C(n),J(4),DETJ
      INTEGER :: n,m,I,J1
      CALL K_INITIALISE_ZERO(A,3,4)
C      WRITE(121,*) 'BLAAAH'    
C      DO J1=1,2
C          WRITE(121,*) D(J1,:8)
C      END DO
C        WRITE(121,*) B(:8)
C        WRITE(121,*) C(:8)
      DO I=1,n
        J(1)=J(1)+B(I)*D(1,I)
        J(2)=J(2)+B(I)*D(2,I)
        J(3)=J(3)+C(I)*D(1,I)
        J(4)=J(4)+C(I)*D(2,I)
      END DO
C
      DETJ=J(1)*J(4)-J(2)*J(3)
C      WRITE(121,*) J(:4)
C      WRITE(121,*) DETJ
C
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
C
C
      SUBROUTINE K_MATRIXMULTIPLY(A,B,C,n,m,l)
      IMPLICIT NONE
      DOUBLE PRECISION :: A(n,m),B(n,l),C(l,m)
      INTEGER :: n,m,l,i,j,k
      CALL K_INITIALISE_ZERO(A,n,m)
      DO i = 1, n
        DO j = 1, m
          DO k = 1, l
            A(i,j) = A(i,j) + B(i,k) * C(k,j)
          END DO
        END DO
      END DO
      RETURN
      END
C
C
      SUBROUTINE K_MATRIXTRANSPOSE(A,B,n,m)
      IMPLICIT NONE
      DOUBLE PRECISION :: A(n,m), B(m,n)
      INTEGER :: i,j,n,m
      CALL K_INITIALISE_ZERO(B,m,n)
      do i = 1, n
         do j = 1, m
            B(j,i) = A(i,j)
         end do
      end do
      RETURN
      END
