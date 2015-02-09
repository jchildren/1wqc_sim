!*********************************************
!
!Module containing the variable cz_operation
!creation and application subroutine
!
!In modular form to make use of multi-level allocatable arrays
!that are standard in Fortran 95
!
!*********************************************

module cz_op_module

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! DEPENDENCIES !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Use the kronecker products module that contains linear alegbra
  !routines necessary for module to function
  use kronecker_module


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!! GLOBAL VARIABLES !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Precision of numerical values
  !integer        :: dp=selected_real_kind(15, 300)			!IEEE 754 Double Precision

contains

  !*********************************************
  !
  !Forms the cz_operator between control and target qubits
  !in state vector of n qubits
  !
  !Applies cz_operation to state vector and returns output
  !
  !*********************************************

  subroutine cz_operation(state_vector, n, ctrl, trgt)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    !Loop integers
    integer :: i, j

    !Qubit number, input from function
    integer :: n

    !Control and target qubit integers, input from function
    !Determines matrix multiplication proceedure for final cz_matrix
    integer :: ctrl, trgt	

    !Constant size matrices
    complex(kind=dp), dimension(:, :), Allocatable :: identity							!2x2 identity matrix
    complex(kind=dp), dimension(:, :), Allocatable :: z_matrix							!2x2 z pauli matrix
  
    !Matrices dependant upon integer n for sizing
    !Two dimensional array of size 2**n by 2**n
    complex(kind=dp), dimension(2**n, 2**n) :: cz_matrix	            !Final cz_matrix

    !State vector array for which size is dependant upon the array input in function
    complex(kind=dp), dimension(:), Allocatable :: state_vector 		!wavefunction arrays

    !Work matrices used in calculation of cz_matrix
    !Allocatable as constantly resized in loops for kronecker products
    complex(kind=dp), dimension(:, :), Allocatable :: work_matrix	    !Primary variable work matrix
    complex(kind=dp), dimension(:, :), Allocatable :: out_matrix	    !Loop output variable work matrix


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Allocate(identity(2, 2), z_matrix(2,2))

    !Typical 2x2 matrix initialisation
    identity(:,:)  = 0.0_dp		!Identity matrix values
    identity(1, 1) = 1.0_dp		!Identity matrix values
    identity(2, 2) = 1.0_dp		!Identity matrix values
  
    z_matrix(:,:)  = 0.0_dp		!z matrix values
    z_matrix(1, 1) = 1.0_dp		!z matrix values
    z_matrix(2, 2) = -1.0_dp    !z matrix values

  
    !Initialises all values of CZ matrix to double precision zero
    !As this matrix is made of several additions this step is
    !necessary to minmise error in calculation
    cz_matrix(:,:) = 0.0_dp


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!! MATRIX GENERATION !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    !Loops over values 1 to 4 with J
    !Reflective of the unitary operator which has 4 terms
    do j = 1, 4

      !Allocates the secondary work matrix into the initial size
      !for Kronecker product multiplication
      Allocate(out_matrix(2, 2))
  
      !Performs a check to see if the first matrix in the order of
      !multiplication corresponds to a control or target 
      !if they do, assigns the appropriate matrix dependant
      !upon the value of j.
      !Otherwise assigns the identity matrix
      if((trgt.eq.1).and.((j.eq.2).or.(j.eq.4))) then
        out_matrix = z_matrix
      elseif((ctrl.eq.1).and.((j.eq.3).or.(j.eq.4))) then
        out_matrix = z_matrix
      else
        out_matrix = identity
      end if

      !Loops over each individual qubit
      !This will produce a matrix of appropriate size for each
      !term of the unitary operator
      do i = 2, n
      
        !Assigns a size value to the work matrix dependent
        !on the step in the loop, thus allowing it
        !to contain the appropriate size of matrix at this step
        Allocate(work_matrix(2**(i-1), 2**(i-1)))
      
        !Assigns the work matrix the value of the output matrix
        !Takes the value from the output of last step
        !Allowing output matrix to be deallocated
        work_matrix = out_matrix
      
        !Deallocate output matrix
        Deallocate(out_matrix)
       
        !Allocates new size to the output matrix
        !New size is appropriate for storage of
        !Kronecker product between work matrix and a 2x2 matrix
        Allocate(out_matrix(2**i, 2**i))
      
        !Determines whether the next operator of the multiplication
        !Will be a control or target qubit, depending on the term of U
        !Assigns the appropriate value if so for kronecker products
        !Otherwise uses the identity matrix for the kronecker product
        if((ctrl.eq.i).and.((j.eq.3).or.(j.eq.4))) then
	      call kronecker_product_complex(work_matrix, z_matrix, out_matrix)
        elseif((trgt.eq.i).and.((j.eq.2).or.(j.eq.4))) then
	      call kronecker_product_complex(work_matrix, z_matrix, out_matrix)
        else
	      call kronecker_product_complex(work_matrix, identity, out_matrix)
        end if
      
        !Deallocates the work matrix ready
        !For allocation in next loop
        Deallocate(work_matrix)
      
      !Ends do loop for the jth term of the Operator
      end do
    
      !Determines which term of operator the loop is on
      !If j is four, removes output from cz_matrix
      !otherwise adds output to cz_matrix
      !Reflective of signs of unitary operator
      if(j.eq.4) then
        cz_matrix = cz_matrix - out_matrix
      else
        cz_matrix = cz_matrix + out_matrix
      end if
    
      !Deallocates out matrix for next loop of j
      Deallocate(out_matrix)
  
    end do
  
  
    !Multiplies every term of the cz_matrix by half
    !This is the first constant of the operator
     cz_matrix = 0.5_dp * cz_matrix

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!! OPERATOR APPLICATION !!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !call for double precision complex matrix vector multiply from BLAS
    !Applies cz_matrix operator to state vector and outputs it
    call ZGEMV( 'N', 2**n, 2**n, 1.0_dp, cz_matrix, 2**n, &
    state_vector, 1, 0.0_dp, state_vector, 1)


   !Returns calculated state vectors to main program
   return
  
  end subroutine cz_operation

end module
