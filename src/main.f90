PROGRAM platehole2d
	USE functions
	
	IMPLICIT NONE

	integer :: nnode,nele,ndof,i,j,k,ok,nfixed,nforce,ndisp
	integer,allocatable :: pivot(:),nodes_fixed(:),nodes_force(:),nodes_disp(:)
 	double precision,allocatable,target:: nodes(:,:)
	type(element),allocatable :: elements(:)
	double precision,allocatable :: globalstiff(:,:)
	double precision,allocatable :: force(:)
	character(100):: filename

	
	!call getarg(1,filename)
	filename='platehole_3.in'
	CALL readinput(filename,nnode,nele,nodes,elements,nfixed,nforce,ndisp,nodes_fixed,nodes_force,nodes_disp)
	

	ndof=nnode*2
	ALLOCATE(globalstiff(ndof,ndof))
	ALLOCATE(force(ndof))
	ALLOCATE(pivot(ndof))
	globalstiff=0.0
	
	OPEN(UNIT=20,FILE= filename(:LEN_TRIM(filename)-3)//'_coords_disp.out')
	OPEN(UNIT=21,FILE= filename(:LEN_TRIM(filename)-3)//'_stress_strain.out')

	!write(*,*) nnode,nele
	!DO i=1,nnode
	!write(*,*) nodes(i,:)
	!write(*,*) ' '
	!END DO
	!write(*,*) elements(1)%node2(1)

	do k=1,nele
		CALL elestiff(elements(k))
		!do i=1,8
    		!write(20,*) ( elements(k)%Kele(i,j), j=1,8 )
		!enddo
		!write(20,*) ''
	enddo


	CALL assemble(globalstiff,force,nodes,elements,nnode,nele,nfixed,nforce,ndisp,nodes_fixed,nodes_force,nodes_disp)
	
	!do i=1,ndof
	    !write(20,*) ( globalstiff(i,j), j=1,ndof )
	!enddo

!	write(*,*) f
	
	write(*,*) nodes_fixed
	write(*,*) nodes_force
	write(*,*) nodes_disp
	
	call DGESV(ndof, 1, globalstiff, ndof, pivot, force, ndof, ok)

	write(*,*) ok

	CALL findnewcoords(nodes,nnode,force)

	DO I=1,nnode	
		write(20,'(E20.7E2, ",", E20.7E2, ",", E20.7E2, ",", E20.7E2, ",", E20.7E2)') nodes(I,1:3),nodes(I,6:7)
	END DO

	write(20,*) 

	DO k=1,nele
		CALL elestiff(elements(k))
		write(21,*) 'element ',k
		DO i=1,4
			write(21,*) i,elements(k)%stress(i,:),elements(k)%strain(i,:)
		END DO
	END DO


	DEALLOCATE(globalstiff,nodes,elements,pivot,force)
	CLOSE(20)
END PROGRAM platehole2d