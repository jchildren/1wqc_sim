!*********************************************
!
!Module containing routines for simulating measurement
!proceedures. Includes general and pauli basis
!Measurements
!
!Also contains feed forward subroutine
!
!In modular form to make use of multi-level allocatable arrays
!that are standard in Fortran 95
!
!*********************************************

module measurement_module
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! DEPENDENCIES !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Use the kronecker products module that contains linear alegbra
  !routines necessary for module to function
  use kronecker_module
  
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!! GLOBAL VARIABLES !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !Precision of numerical values
  !integer,          parameter :: dp=selected_real_kind(15, 300)		!IEEE 754 Double Precision
  
  !Data format value
  !Character, two integers then three double precision numbers
  !Spacing of 3 between entries
  character(len=70), save :: data_format = '(A4, 2I4, 3E25.16E3)'
  
  !Value of pi to be used by program
  !Used to convert phase angles
  real(kind=dp), parameter :: pi = 3.1415926535897932384626433832795
  
contains
  
  !*********************************************
  !
  !Performs generalised basis measurement on target qubit
  !in state vector
  !
  !Measurement phase read through arguments of subroutine
  !Writes measurement data to file
  !*********************************************

  subroutine phi_measurement(state_out, state_in, n, m_phase, trgt)
    implicit none
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !input arguments
    integer :: trgt 								!target qubit
    integer :: n    								!number of qubits
    real(kind=dp) :: m_phase						!measurement phase value

    !states for calculation    
    complex(kind=dp), dimension(:), Allocatable :: state_in, state_out

    !measurement variables
    integer :: m_result								!measurement outcome (0 or 1)
    complex(kind=dp) :: m_value						!measurement inner product value (braket)
    complex(kind=dp), dimension(2) :: m_vector		!measurement vector
    
    !variables for function calls
    real(kind=dp) :: rand_num						!random number storage
    complex(kind=dp) :: ZDOTC						!complex dot product subroutine
    
    !work vector for subroutine calculation
    complex(kind=dp), dimension(:), Allocatable :: work_vector
    
    !common matrices
    complex(kind=dp), dimension(2, 2) :: x_matrix	!pauli x
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    x_matrix(:, :) = 0.0_dp		!Pauli x
    x_matrix(1, 2) = 1.0_dp		!Pauli x
    x_matrix(2, 1) = 1.0_dp		!Pauli x

    m_vector(1) = EXP(CMPLX(0.0_dp, (- m_phase) / 2.0_dp, kind=dp))
    m_vector(2) = EXP(CMPLX(0.0_dp, (m_phase) / 2.0_dp, kind=dp))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! MEASUREMENT OUTCOME !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Calls for random number from intrinsic subroutine
    !rounds output to nearest integer (0 or 1)
    call RANDOM_NUMBER(rand_num)
    m_result = NINT(rand_num)
    
    !adjusts measurement vector based on measurement
    !result. multiplies by pauli x if result = 1
    if(m_result.eq.0) then
      continue
    elseif(m_result.eq.1) then
      m_vector = MATMUL(x_matrix, m_vector)
    else
      print *, 'Error with random measurement results'
      stop
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Gets the state vector of the target qubit and outputs it to work vector
    call get_vector_complex(state_in, work_vector, state_out, n, trgt)
    
    !Calculates inner product between measurement result and target qubit
    !Uses either BLAS ZDOTU or intrinsic function depending on preprocessor
#ifdef lblas

    m_value = ZDOTC(2**n, m_vector, 1, work_vector, 1)

#else

    m_value = DOT_PRODUCT(m_vector, work_vector)

#endif

    !Multiplies result of dot product with rest of state vector
    state_out = m_value * state_out
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FILE I/O !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Prints data about measurements to file for use in
    !feed forward subroutine later
    open(100, file='measurements.dat', access='append', iostat=ierr)
      if (ierr/=0) stop 'Error in opening file measurements.dat'
    
    write(100, data_format) 'R', trgt, m_result, m_phase, REALPART(m_value), IMAGPART(m_value)
    
    close(100, iostat=ierr)
      if (ierr/=0) stop 'Error in closing file measurements.dat'
    
    return
    
  end subroutine phi_measurement

  !*********************************************
  !
  !Perfoms measurements on target qubits in state vector
  !in the pauli basis of choice
  !
  !Writes measurement data to file
  !*********************************************

  subroutine pauli_measurement(state_out, state_in, n, m_type, trgt)
    implicit none
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !input arguments
    integer :: trgt 								!target qubit
    integer :: n    								!number of qubits
    real(kind=dp) :: m_phase						!measurement phase value

    !states for calculation    
    complex(kind=dp), dimension(:), Allocatable :: state_in, state_out

    !measurement variables
    integer :: m_result								!measurement outcome (0 or 1)
    complex(kind=dp) :: m_value						!measurement inner product value (braket)
    character(len=1) :: m_type						!measurement time string
    complex(kind=dp), dimension(2) :: m_vector		!measurement vector

    !variables for function calls
    real(kind=dp) :: rand_num						!random number storage
    complex(kind=dp) :: ZDOTC						!complex dot product subroutine
    
    !work vector for subroutine calculation
    complex(kind=dp), dimension(:), Allocatable :: work_vector

    !common vectors
    complex(kind=dp), dimension(2) :: x_evector					!x pauli eigenvector
    complex(kind=dp), dimension(2) :: z_evector					!z pauli eigenvector
    complex(kind=dp), dimension(2) :: y_evector					!y pauli eigenvector

    !common matrices
    complex(kind=dp), dimension(2, 2) :: x_matrix
    complex(kind=dp), dimension(2, 2) :: z_matrix
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Basis vectors
    x_evector(1) = (1.0_dp / SQRT(2.0_dp)) * 1.0_dp		!Pauli x eigenvector
    x_evector(2) = (1.0_dp / SQRT(2.0_dp)) * 1.0_dp		!Pauli x eigenvector
    
    z_evector(1) = (1.0_dp / SQRT(2.0_dp)) * 1.0_dp		!Pauli z eigenvector
    z_evector(2) = (1.0_dp / SQRT(2.0_dp)) * 0.0_dp		!Pauli z eigenvector
    
    y_evector(1) = (1.0_dp / SQRT(2.0_dp)) * (1.0_dp, 0.0_dp)  !Pauli y eigenvector
    y_evector(2) = (1.0_dp / SQRT(2.0_dp)) * (0.0_dp, 1.0_dp)  !Pauli y eigenvector
    
    !pauli matrices
    x_matrix(1, 2) = 1.0_dp		!Pauli x
    x_matrix(1, 1) = 0.0_dp		!Pauli x
    x_matrix(2, 2) = 0.0_dp		!Pauli x
    x_matrix(2, 1) = 1.0_dp		!Pauli x
    
    z_matrix(1, 1) = 1.0_dp		!Pauli z
    z_matrix(1, 2) = 0.0_dp		!Pauli z
    z_matrix(2, 1) = 0.0_dp		!Pauli z
    z_matrix(2, 2) = -1.0_dp	!Pauli z

    !Selects appropriate measurement vector and measurement
    !phase based on input arguments
    if(m_type == 'X') then
      m_vector = x_evector
      m_phase = 0.0_dp
    elseif(m_type == 'Z') then
      m_vector = z_evector
      m_phase = 0.0_dp
    elseif(m_type == 'Y') then
      m_vector = y_evector
      m_phase = pi / 2.0_dp
    else
      print *, 'Error selecting measurement matrix check subroutine input'
      stop
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! MEASUREMENT OUTCOME !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Calls for random number from intrinsic subroutine
    !rounds output to nearest integer (0 or 1)
    call RANDOM_NUMBER(rand_num)
    m_result = 0!NINT(rand_num)
    
    !adjusts measurement vector based on measurement
    !result. multiplies by pauli x if result = 1 for z
    !multiplies by pauli z if result = 1 for x and y.
    if(m_result.eq.0) then
      continue
    elseif((m_result.eq.1).and.(m_type.eq.'Z')) then
      m_vector = MATMUL(x_matrix, m_vector)
    elseif((m_result.eq.1).and.(m_type.eq.'X')) then
      m_vector = MATMUL(z_matrix, m_vector)
    elseif((m_result.eq.1).and.(m_type.eq.'Y')) then
      m_vector = MATMUL(z_matrix, m_vector)
    else
      print *, 'Error with random measurement results'
      stop 
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Gets the state vector of the target qubit and outputs it to work vector
    call get_vector_complex(state_in, work_vector, state_out, n, trgt)
    

    !Calculates inner product between measurement result and target qubit
    !Uses either BLAS ZDOTU or intrinsic function depending on preprocessor
#ifdef lblas

    m_value = ZDOTC(2**n, m_vector, 1, work_vector, 1)

#else

    m_value = DOT_PRODUCT(m_vector, work_vector)

#endif
    

    !Multiplies result of dot product with rest of state vector
    state_out = m_value * state_out
    


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FILE I/O !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Prints data about measurements to file for use in
    !feed forward subroutine later
    open(100, file='measurements.dat', access='append', iostat=ierr)
      if (ierr/=0) stop 'Error in opening file measurements.dat'
    
    write(100, data_format) m_type, trgt, m_result, m_phase, REALPART(m_value), IMAGPART(m_value)
    
    close(100, iostat=ierr)
      if (ierr/=0) stop 'Error in closing file measurements.dat'
    
    return
  
  end subroutine pauli_measurement
  
  !*********************************************
  !
  !Reads measurement data from file and performs calculations
  !to obtain output states for qubits
  !
  !
  !*********************************************
  
  subroutine feed_forward(state_vector, n, m_number, trgt)
    implicit none
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!loop integers
	integer :: i

    !input arguments
    integer :: trgt, n, m_number
    
    !measurement variables
    integer :: m_target
    integer :: m_result
    real(kind=dp) :: m_phase, m_value_re, m_value_im
    character(len=1) :: m_type
    complex(kind=dp) :: m_value
    
    !state vectors for processing
    complex(kind=dp), dimension(:), Allocatable :: state_vector
    complex(kind=dp), dimension(:), Allocatable :: trgt_vector
    
    !common matrices
    complex(kind=dp), dimension(2, 2) :: x_matrix
    complex(kind=dp), dimension(2, 2) :: h_matrix

    !calculated matrices
    complex(kind=dp), dimension(2, 2) :: rz_matrix
    complex(kind=dp), dimension(2, 2) :: h_rz_matrix
    
    !work vectors
    complex(kind=dp), dimension(:), Allocatable :: work_vector
    complex(kind=dp), dimension(:), Allocatable :: top_vector, bot_vector

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!! INITIALISATION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    h_matrix(1, 1) = (1.0_dp / SQRT(2.0_dp))		!Hadamard matrix
    h_matrix(1, 2) = (1.0_dp / SQRT(2.0_dp))		!Hadamard matrix
    h_matrix(2, 1) = (1.0_dp / SQRT(2.0_dp))		!Hadamard matrix
    h_matrix(2, 2) = -(1.0_dp / SQRT(2.0_dp))		!Hadamard matrix
    
    x_matrix(1, 2) = 1.0_dp							!Pauli x
    x_matrix(1, 1) = 0.0_dp							!Pauli x
    x_matrix(2, 2) = 0.0_dp							!Pauli x
    x_matrix(2, 1) = 1.0_dp							!Pauli x

    Allocate(trgt_vector(2))
    
    !Open measurements file for reading by main function
    open(100, file='measurements.dat', access='sequential', iostat=ierr)
      if (ierr/=0) stop 'Error in opening file measurements.dat'
    
    !determines number of qubits, if qubit no > 1, finds the target qubit state vector
    if(n.eq.1) then
      trgt_vector = state_vector
    else
      Allocate(work_vector(2**(n-1)))
      call get_vector_complex(state_vector, trgt_vector, work_vector, n, trgt)
    end if

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !loops over the number of measurements performed
    do i = 1, m_number

      read(100, data_format) m_type, m_target, m_result, m_phase, m_value_re, m_value_im
      
      !forms rotation matrix from phase information
      rz_matrix(:,:) = 0.0_dp
      rz_matrix(1,1) = EXP(CMPLX(0.0_dp, (- m_phase) / 2.0_dp, kind=dp))
      rz_matrix(2,2) = EXP(CMPLX(0.0_dp, (m_phase) / 2.0_dp, kind=dp))
      
      !forms matrix that is product of rotation and hadamard matrices
      h_rz_matrix = MATMUL(h_matrix, rz_matrix)
      
      !computes output state based on measurement information and stores in target vector
      if(m_result.eq.0) then
	    trgt_vector = EXP(CMPLX(0.0_dp, (- m_phase) / 2.0_dp, kind=dp)) * MATMUL(h_rz_matrix, trgt_vector)
      elseif(m_result.eq.1) then
	    trgt_vector = EXP(CMPLX(0.0_dp, (- m_phase) / 2.0_dp, kind=dp)) * MATMUL(MATMUL(x_matrix, h_rz_matrix), trgt_vector)
      endif
      
    end do
    
    !close measurement information file as it is no longer needed
    close(100, iostat=ierr)
      if (ierr/=0) stop 'Error in closing file measurements.dat'
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! STATE VECTOR RECONSTRUCTION !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !checks to see if state vector needs to be reconstructed
    !(i.e. if n is not 1 then must be recombine)
    if(n.eq.1) then
      !if n = 1 no recombination necessary
      state_vector = trgt_vector
    else
      !allocates work vectors based on target and number of qubits
      Allocate(top_vector(2**trgt), bot_vector(2**(n-trgt)))
    
      !breaks state vector apart for recombination
      call sv_decomposition_complex(work_vector, top_vector, bot_vector)
  
      !resizes work vector for next calculation
      Deallocate(work_vector)
      Allocate(work_vector(2**trgt))

      !recombines state vectors
      call kronecker_product_complex_vector(top_vector, trgt_vector, work_vector)
      call kronecker_product_complex_vector(work_vector, bot_vector, state_vector)
      
      Deallocate(top_vector, bot_vector, work_vector)

    end if
    
    Deallocate(trgt_vector)
  
    return
  
  end subroutine

  !*********************************************
  !Finds "fidelity" of two states
  !
  !
  !Can be compiled with or without BLAS
  !
  !
  !*********************************************
  function fidelity(state_one, state_two)
  	implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Function variables
    integer :: n                    									!integer size of state vectors
    real(kind=dp) :: fidelity       									!Fidelity variable for function output
    complex(kind=dp) :: ZDOTC											!ZDOT variable for blas library function

    !Input state vectors
    complex(kind=dp), dimension(:), Allocatable :: state_one		!First input state vector
    complex(kind=dp), dimension(:), Allocatable :: state_two		!Second input state vector

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!! ARRAY SIZE CHECK !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print *, SIZE(state_one), SIZE(state_two)

    !Checks to see if array sizes of two state vectors match
    !If they do assigns value to n for inner product calculation
    !Otherwise stops program and prints error message
    if(SIZE(state_one).eq.SIZE(state_two)) then
      n = SIZE(state_one)
    else
      stop 'Array size mismatch in fidelity function'
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!! MAIN FUNCTION !!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!If BLAS is enabled when compiling (-lblas), this code segment
!will be compiled into final program using BLAS functions
#ifdef lblas
    
    fidelity = (abs(ZDOTC(n, state_one,1, state_two, 1)))**2
    
!If BLAS is not enabled when compiling, this section of intrinsic functions
!will be used instead.
#else
    
  	fidelity = (abs(DOT_PRODUCT(state_one, state_two)))**2
    
!ends the preprocessor if statement
#endif
    
  	!returns value of fidelity from subroutine to main program
  	return
 
  end function fidelity

end module measurement_module
