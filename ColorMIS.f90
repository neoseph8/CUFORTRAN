!! Starting with a symetric, sparse matrix with values W(i,j), i, j having length n
!! This program will create a maximal independent set from the data. It will also solve the independent color scheme for the data


   
	
	
	
	
	
Program MISGraph
use ifport
USE OMP_LIB

	implicit none
 
Type :: Graph
	Integer :: ActiveFlag                       !The status of the node
	Integer :: Charisma                         !The number of connections the node has
	Integer, allocatable, dimension (:) :: Edge !The vertices for neighboring nodes.
	
end type








!!Variable Declarations
Integer :: status, io , iunit     !Status value for array allocation
Integer :: i, j, k, L, iNcount      !Counting integers
Type(Graph), allocatable, dimension(:) :: Node
Integer :: n, nnz                                       !Size of the Array being inputted and number of non-zero values
real, allocatable, dimension(:,:) :: W              !Weight values for the input matrix
Logical Color, VState                               !Color asks if all Nodes have a color assigned. VState Asks if all nodes are accounted for in generating th MIS
real :: Wavelength                                  !Randomly generated wavelength of visible light, in nanometers
Real :: my_rando                                     !Place holder for randomly generated number
Integer :: Start, Finish                            !Start and End time place holders
real ::  CountRate, Sparcity				!Countrate for the clock, desired Sparcity of the fabricated matrices, and  
logical :: debugFlag
character ( 20) :: FileName, CharN
real, dimension (:), allocatable :: ColorV						!Color Vector. 1000 should be more than enough

!End Variable Declaration
debugFlag = 0

Sparcity = 0.95     !Value between 0 and 1. The larger the value, the fewer non-zero values the matrix will have




do iNCount = 1,5
	n = 10 ** iNCount
Allocate (W(n,n), STAT = status)
write(*,*)"W Allocation Status: ",Status
Allocate (Node(n), STAT = status)  
write(*,*)"Node Allocation Status: ",Status
Allocate(ColorV(n), STAT = status)
write(*,*)"ColorV Allocation Status: ",Status





L = 1 !Number of Colors





!!Input W matrix

iunit = 99
write(CharN, *) n
FileName = trim("MatrixIN"//trim(CharN)//".in")

open(unit = iunit, file = FileName, status = 'replace', access = 'sequential', form = 'formatted', action = 'write', IOSTAT = io)

write(iunit,*) "Input File Created for calculating MIS color scheme for ", n,"x",n," Symmetric, Sparse Matrix"
write(iunit,*) "Note that upper tirangular values are eschewed"
write(iunit,*) "Row, Column, Value"
 ! READ(iunit,*) nrow , ncol, nelem

!  n = nrow
  


!initialize W values


!Create Symmetric, positively weighted matrix


!!$OMP PARALLEL DO
do i = 1,n
	do j = i,n
		W(i,j)=100.0*(rand()-(sparcity))
		if (W(i,j) < 0) then
			W(i,j) = 0.0
		end if
		W(j,i) = W(i,j)
		if (W(i,j)>0) write(iunit, *) i, j, W(i,j)
		
    end do
    W(i,i) = 100.0
end do

 !!$OMP END PARALLEL DO

!!Read values
!do k = 1,nelem
!	 Read(iunit,*) i, j, W(i,j)
!	 W(j,i) = W(i,j)
!
!end do


close (unit=iunit)


  



call system_clock(Start, CountRate)



!!$OMP PARALLEL DO
do i=1,n        !Par for
	Node(i)%ActiveFlag = 0      !Set initial neutral value for active flag
	Node(i)%Charisma = count (W(i,:)>0)-1       !Count Non Zeroes for Charisma score sans center value
	Allocate(Node(i)%Edge(Node(i)%Charisma))
	k=1
	do j = 1,n
		if (W(i,j) > 0 .AND. (i/=j)) then
			
			Node(i)%Edge(k) = j  !Add another neighbor
			k = k+1
		end if
		
	end do
  
end do

 ! !$OMP END PARALLEL DO



	
	
	

!!Start Color Scheme While Loop
call srand(86753)
Color = 1
do while (Color)
!call RANDOM_NUMBER (Wavelength)
	
	
ColorV(L) = 400 + 400 * RAND(0)

if (debugFlag) write(*,*) "Color for this MIS (nm): ", ColorV(L)





!! Start VState Loop. VState will come up false once every vertice is accounted for

VState = 1
do while (VState)





!!ActiveFlag Values:
!!-2 Removed from set of vertices
!!-1 Removed from proposal for MIS
!!0 Neutral
!!1 Proposed for MIS
!!>=2 Assigned color

!!Set Active flag = 1 for randomly selected Nodes

	!$OMP PARALLEL DO
	
do i=1,n
   if (ABS(Node(i)%ActiveFlag)<2) then
	  !call RANDOM_NUMBER (my_rando)
	   my_rando = rand(0)
	   
	   if (debugFlag) then
			write(*,*) "Random Number and Current CHA derivative: ", my_rando, "  ", 1/(2*real(Node(i)%Charisma))
	   endif
	   
	  if (my_rando < 1/(2*real(Node(i)%Charisma))    .or. Node(i)%Charisma==0) then!This Node has been selected for consideration!
		  
		  Node(i)%ActiveFlag = 1
	  endif
	   
   endif
   
end do

!$OMP END PARALLEL DO

if (debugFlag) then
	
write(*,*) "1 - Neutral Nodes (0) selected for consideration (1)"
do i=1,n
		
	Write(*,*) Node(i)%ActiveFlag
	
end do

endif

	
	!$OMP PARALLEL DO
	

!!Set Active Flag = -1 for contested nodes
do i=1,n
	if ((Node(i)%ActiveFlag) == 1) then  !Only handles Nodes that are preselected for MIS
		do k=1,Node(i)%Charisma
			if (Node( Node(i)%Edge(k) )%ActiveFlag == 1)  then !Contended Edge
				if (Node( Node(i)%Edge(k) )%Charisma < Node(i)%Charisma  ) then
					Node( Node(i)%Edge(k) )%ActiveFlag = -1
				else
					Node(i)%ActiveFlag = -1
				end if
			end if
		end do
	end if
end do
!$OMP END PARALLEL DO

if (debugFlag) then

write(*,*) "2 - Contested Nodes (1) removed from consideration (-1)"
do i=1,n
		
	Write(*,*) Node(i)%ActiveFlag
	
end do    
		
endif







!!Set active flag = -2 for border nodes removed from the base set, and 2 or greater for nodes declared part of the MIS
!$OMP PARALLEL DO
   
do i = 1,n
	if (Node(i)%ActiveFlag == 1) then   !Shift Active Flag to Current Color
		Node(i)%ActiveFlag = ColorV(L)
		
		do k = 1,Node(i)%Charisma
			if (Node(Node(i)%Edge(k))%ActiveFLag<2) then  !Only remove Nodes that border and are not already removed
				Node( Node(i)%Edge(k) )%ActiveFLag = -2  !Remove Node from consideration 
			endif
		end do
		
	end if
	
end do
 !$OMP END PARALLEL DO

if (debugFlag) then

Write(*,*) "3 - Colors Assigned (2 -> ", ColorV(L)
	do i=1,n
		
	Write(*,*)  Node(i)%ActiveFlag
	
end do    
	  
endif

	  
!!Check for Vertice state doneness (ie, all active flags should be -2 or 2+)      
!!Restore unremoved nodes to neutral (set active flag to 0)      
VState = 0      
	  

!$OMP PARALLEL DO
	
do i = 1,n
   if (Node(i)%ActiveFlag == 0 ) then
	   VState = 1
   else if (Node(i)%ActiveFlag == -1 ) then
	   VState = 1
	   Node(i)%ActiveFlag = 0
   end if
	
end do
!$OMP END PARALLEL DO

if (debugFlag) then

Write(*,*) "4 - Reset Discarded Nodes (-1) to neutral (0)"
	do i=1,n
		
	Write(*,*)  Node(i)%ActiveFlag
	
end do    
	
endif

!Write(*,*) "Node 3 Status"
! 
!        
!    Write(*,*)"Active Flag: ",  Node(3)%ActiveFlag, " and Charisma Score " ,  Node(3)%Charisma
!    





!!End VSTate While loop
	end do 
	
	

Color = 0
!!Reset Node Active Flags set to -2 to 0 for next color.

!$OMP PARALLEL DO
	
do i = 1,n
	if (Node(i)%ActiveFlag == -2) then
		Node(i)%ActiveFlag = 0      !reset flags
		Color = 1
	else if (Node(i)%ActiveFlag <400) then
		Color = 1
		
	end if
end do
!$OMP END PARALLEL DO

!
L = L+1   !Measure the number of colors for a MIS set
if (L>1000) then
    write(*,*) "1000 is not enough"
endif
!



!!Check for color completeness



!!End Color While Loop


	end do
	
	

!!Output useful data

 call system_clock(Finish, CountRate)
	
write(*,*) "Time to complete Parallel MIS (seconds): ",(Finish - Start)/CountRate


iunit = 25
FileName = trim("MatrixOUT"//trim(CharN)//".out")

open(unit = iunit, file = FileName, status = 'replace', access = 'sequential', form = 'formatted', action = 'write', IOSTAT = io)

write (iunit, *) "Maximally Independent Colored Sets for ", n, "x", n, "Matrix. ", L-1, "Colors used."

do i=1,L-1
    write(iunit, *) "MIS: ", i,"Wavelength (nm): ", ColorV(i)
end do

write(iunit, *) "Time to Completion(seconds): ",(Finish - Start)/CountRate


write (iunit, *) "Vertice and Visible Light Wavelength (nm)"

do i = 1,n
	write (iunit, *) i, Node(i)%ActiveFlag
end do


close(unit = iunit)


	Deallocate (Node)
	Deallocate (W)
    Deallocate (ColorV)
	
	
	
	
	
end do

	
	
End program