!*********************************************
!
!Module containing routines for kronecker products
!and vector decomposition
!
!In modular form to make use of multi-level allocatable arrays
!that are standard in Fortran 95
!
!*********************************************
module kronecker_module
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!! GLOBAL VARIABLES !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Precision of numerical values
  integer, parameter		 :: dp=selected_real_kind(15, 300)		!IEEE 754 Double Precision
  
  !Error reporting
  integer, 			save :: ierr 										!Error integer

contains

  
  !*********************************************
  !
  !Performs Kronecker product of two matrices A and B producing matrix P
  !n and m correspond to the dimensions of A and x and y correspond to the dimensions of B
  !dimensions of P are assumed to be n*x and m*y
  !
  !*********************************************
  
  subroutine kronecker_product_complex(A, B, P)
    implicit none
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !loop integers
    integer :: i, j, k, l

    !array size integers
    integer :: n, m, x, y

    !arrays for matrix storage
    complex(kind=dp), dimension(:, :), Allocatable :: A, B, P
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !calls the size of each dimension of matrix A 
    !and assigns output to integer values
    n = SIZE(A, 1)
    m = SIZE(A, 2)
    
    !calls the size of each dimension of matrix B
    !and assigns output to integer values
    x = SIZE(B, 1)
    y = SIZE(B, 2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! ARRAY SIZE CHECK !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Checks to see if array sizes match for kronecker
    !product, aborts program if mismatch and prints error message
	if(SIZE(P, 1).eq.(n*x)) then
	  continue
	else if(SIZE(P, 2).eq.(m*y)) then
	  continue
	else
	  stop 'Array size mismatch in kronecker product routine'
	end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, n
      do j = 1, m
    
	    do k = 1, x
	      do l = 1, y
	        
            !Computes kronecker product values and
            !assigns them to matrix P
	        P((i-1)*x + k, (j-1)*y + l) = A(i, j) * B(k, l)

	      end do
	     end do

      end do
    end do
    
    return
  
  end subroutine kronecker_product_complex
 
  !*********************************************
  !
  !Performs Kronecker product of two vectors A and B producing vector P
  !n is the array size of of A and x is the array size of B
  !P is allocated a size of n*x
  !
  !*********************************************

  subroutine kronecker_product_complex_vector(A, B, P)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !loop integers
    integer :: i, j
    
    !array size integers
    integer :: n, x

    !arrays for matrix storage
    complex(kind=dp), dimension(:), Allocatable :: A, B, P
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !calls the size of vector A
    !and assigns output to integer value n
    n = SIZE(A)
    
    !calls the size of vector B
    !and assigns output to integer value x
    x = SIZE(B)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! ARRAY SIZE CHECK !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Checks to see if array sizes match for kronecker
    !product, aborts program if mismatch and prints error message
	if(SIZE(P, 1).eq.(n*x)) then
	  continue
	else
	  stop 'Array size mismatch in kronecker product routine'
	end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, n
      do j = 1, x
         
         !Computes kronecker product values and
         !assigns them to vector P
	     P((i-1)*x + j) = A(i) * B(j)

      end do
    end do
    
    
    return
  
  end subroutine kronecker_product_complex_vector
  
  
    !*********************************************
  !
  !Performs Kronecker product of two vectors A and B producing vector P
  !n is the array size of of A and x is the array size of B
  !P is allocated a size of n*x
  !
  !*********************************************

  subroutine outer_product_complex_vector(A, B, P)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !loop integers
    integer :: i, j
    
    !array size integers
    integer :: n, x

    !arrays for matrix storage
    complex(kind=dp), dimension(:), Allocatable :: A, B
    complex(kind=dp), dimension(:, :), Allocatable :: P
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !calls the size of vector A
    !and assigns output to integer value n
    n = SIZE(A)
    
    !calls the size of vector B
    !and assigns output to integer value x
    x = SIZE(B)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! ARRAY SIZE CHECK !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Checks to see if array sizes match for outer
    !product, aborts program if mismatch and prints error message
	if((SIZE(P, 1).eq.(n)).and.(SIZE(P, 2).eq.(x))) then
	  continue
	else
	  stop 'Array size mismatch in between vectors and product matrix'
	end if
	
	if(x == n) then
	  continue
	else
	  stop 'Array size mismatch between vectors in outer product routine'
	end if
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, n
      do j = 1, x
         
         !Computes outer product values and
         !assigns them to matrix P
         P(i, j) = A(i) * CONJG(B(j))

      end do
    end do
    
    
    return
  
  end subroutine outer_product_complex_vector
  
  !*********************************************
  !
  !Performs decomposition of vector A by performing inverse vectorisation
  !then taking a rank 1 decomposition
  !
  !Outputs vectors C and D
  !*********************************************
  
  subroutine rank_decomposition_complex(A, C, D)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    !loop integers
    integer :: i, j

    !array size integers
    integer :: m, n, o

    !magnitude variable for normalisation
    real(kind=dp) :: magnitude

    !input and output vectors
    complex(kind=dp), dimension(:), Allocatable :: A, C, D

    !matrix for decomposition calculator
    complex(kind=dp), dimension(:,:), Allocatable :: B
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !calls the size of vector A
    !and assigns output to integer value m
    m = SIZE(A)

    !calls the size of vector C
    !and assigns output to integer value n
    n = SIZE(C)

    !calls the size of vector D
    !and assigns output to integer value o
    o = SIZE(D)
    
    !Allocates B with dimensions from size of input vectors
    Allocate(B(n, o))
    

    !Perform an inverse of the vec() operator by reshaping
    !vector A into a matrix of dimensions n x o 
    B = TRANSPOSE(RESHAPE(A, (/ o, n /)))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Finds the first non zero row from matrix B
    !and assigns it to vector C
    do i = 1, o
      do j = 1, n
	if(B(j, i) /= 0.0_dp) then
	  C(j) = B(j, i)
	  exit
	else
	  continue
	end if
      end do
    end do

    
    !Calculate vector magnitude of C
    do i = 1, n
      magnitude = magnitude + C(i)*C(i)
    end do    
    
    magnitude = sqrt(magnitude)

    
    !Divide components of C by magnitude to normalise
    C = C / magnitude
    
    !Find the value of D required for matrix B to be formed
    !Reversing the kronecker product
    do i = 1, o
      do j = 1, n
	if(C(j) /= 0.0_dp) then
	  D(i) = (1.0_dp/(C(j))) * B(j, i)
	  exit
	else
	  continue
	end if
      end do
    end do
    
    Deallocate(B)
    
    return
    
    
  end subroutine
  
    !*********************************************
  !
  !Performs decomposition of vector A by performing inverse vectorisation
  !then taking a rank 1 decomposition
  !
  !Outputs vectors C and D
  !*********************************************
  
  subroutine known_rank_decomposition_complex(A, C, D)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    !loop integers
    integer :: i, j

    !array size integers
    integer :: m, n, o

    !input and output vectors
    complex(kind=dp), dimension(:), Allocatable :: A, C, D

    !matrix for decomposition calculator
    complex(kind=dp), dimension(:,:), Allocatable :: B
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !calls the size of vector A
    !and assigns output to integer value m
    m = SIZE(A)

    !calls the size of vector C
    !and assigns output to integer value n
    n = SIZE(C)

    !calls the size of vector D
    !and assigns output to integer value o
    o = SIZE(D)
    
    !Allocates B with dimensions from size of input vectors
    Allocate(B(n, o))
    

    !Perform an inverse of the vec() operator by reshaping
    !vector A into a matrix of dimensions n x o 
    B = TRANSPOSE(RESHAPE(A, (/ o, n /)))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !Find the value of D required for matrix B to be formed
    !Reversing the kronecker product
    do i = 1, o
      do j = 1, n
	if(C(j) /= 0.0_dp) then
	  D(i) = (1.0_dp/(C(j))) * B(j, i)
	  exit
	else
	  continue
	end if
      end do
    end do
    
    Deallocate(B)
    
    return
    
    
  end subroutine

  !*********************************************
  !
  !Performs decomposition of vector A by performing inverse vectorisation
  !then taking a single value decomposition and summing values of matrices
  !
  !Outputs vectors C and D
  !
  !*********************************************  

  subroutine sv_decomposition_complex(A, C, D, U, S, VT)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !loop variables
    integer :: i, j

    !array size integers
    integer :: m, n, o

    !Input and output vectors
    complex(kind=dp), dimension(:), Allocatable :: A, C, D

    !Matrix for computation
    complex(kind=dp), dimension(:,:), Allocatable :: B
    
    !Output matrices and vectors of svd
    real(kind=dp), dimension(:), Allocatable :: S
    complex(kind=dp), dimension(:,:), Allocatable :: U, VT
    
    !variables required zgesvd
    integer :: LWORK, LDA, LDU, LDVT
    complex(kind=dp), dimension(:), Allocatable :: WORK
    real(kind=dp), dimension (:), Allocatable :: RWORK
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !calls the size of vector A
    !and assigns output to integer value m
    m = SIZE(A)

    !calls the size of vector C
    !and assigns output to integer value n
    n = SIZE(C)

    !calls the size of vector D
    !and assigns output to integer value o
    o = SIZE(D)
    
    !Assigns values to work variables
    LWORK = MAX(1, 4*(2*MIN(n, o) + MAX(n, o)))
    LDA = MAX(1, n)
    LDU = MAX(1, n)
    LDVT = MAX(1, o)
    
    !Allocates matrix for decomposition size n x o
    Allocate(B(n, o))

    !Allocates outputs of SVD
    Allocate(S(MIN(n, o)))
    Allocate(U(LDU, n))
    Allocate(VT(LDVT, o))
    
    !Allocates work arrays
    Allocate(WORK(MAX(1, LWORK)))
    Allocate(RWORK(5*MIN(n, o)))
    
    !Perform an inverse of the vec() operator by reshaping
    !vector A into a matrix of dimensions n x o 
    B = TRANSPOSE(RESHAPE(A, (/ o, n /)))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Calls the zgesvd LAPACK function to find SVD
    !stops program and prints error in case of failure
    call zgesvd('A', 'A', n, o, B, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, ierr)
      if(ierr.eq.0) then
	continue
      elseif(ierr.gt.0) then
	stop 'ZBDSQR did not converge'
      elseif(ierr.lt.0) then
	stop 'argument had an illegal value'
      end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!! OUTPUT !!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Recombines values of svd to form vectors
    do i = 1, n
      C(i) = SUM(U(:,i))
    end do
    
    do i = 1, o
      D(i) = S(i)*SUM(VT(i,:))
    end do

    !Deallocates values no longer needed
    Deallocate(B, WORK, RWORK)
    
    return
    
    
  end subroutine
  
  !*********************************************
  !
  !Finds a 2 dimensional state vector for a specific qubit
  !in larger state vector
  !
  !Calls for decompositions to find appropriate qubit then recombines state vector
  !*********************************************
  
  subroutine get_vector_complex(A, E, G, m, trgt, U, S, VT)
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer :: trgt
    integer :: m, n, o, p
    complex(kind=dp), dimension(:), Allocatable :: A, E, G
    complex(kind=dp), dimension(:), Allocatable :: C, D, F
    
    !Output matrices and vectors of svd
    real(kind=dp), dimension(:), Allocatable :: S
    complex(kind=dp), dimension(:,:), Allocatable :: U, VT
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Assigns vector sizes depending on subroutine input
    if(trgt.eq.1) then
      n = 2
      o = 2**(m - 1)
    elseif(trgt.eq.m) then
      o = 2
      n = 2**(m - 1)
    else
      o = 2**(m - trgt)
      n = 2**(trgt)
      p = 2**(trgt-1)
    end if
    
    Allocate(D(o), C(n))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    !Calls for svd of vector A
    call sv_decomposition_complex(A, C, D, U, S, VT)

    !Assigns values of decomposition based on target qubit
    !If target is not at either end of state vector, performs
    !another decomposition to find it
    if(trgt.eq.1) then
      E = C
      G = D
    elseif(trgt.eq.m) then
      G = C
      E = D
    else
      !Allocate(F(p))    
    

    
      !call sv_decomposition_complex(D, E, F, U, S, VT)

      
      !call kronecker_product_complex_vector(C, F, G)
      
      !Deallocate(F)
    end if


    
    Deallocate(C, D)

    
    return
    
  end subroutine
  
end module
