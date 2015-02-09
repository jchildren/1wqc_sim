!*********************************************
!Program for simulation of measurements on a chain of
!cluster states
!
!Relies on several external modules and BLAS and LAPACK to function
!
!
!*********************************************
program single_chain

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! DEPENDENCIES !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Use the kronecker products module that contains linear alegbra
  !routines necessary for module to function
  use kronecker_module

  !Uses the module for the generation of control z operators
  !Necessary for creation of cluster states
  use cz_op_module

  !Use the measurement module, containing subroutines handling
  !The measurement of single qubits in state vectors
  use measurement_module
  
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!! FUNCTIONALITY VARIABLES !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Variables required for program functionality

  !Precision of numerical values
  !integer,         :: dp=selected_real_kind(15, 300)		!IEEE 754 Double Precision
  
  !Loop integers
  integer :: i
  integer :: m_number

  !measurement time string
  character(len=1) :: m_type
  
  !variables for file i/o
  real(kind=dp) :: fid_out
  real(kind=dp) :: m_phase

  !input and output state vectors
  complex(kind=dp), dimension(:), Allocatable :: state_in, state_out
  
  !Output matrices and vectors of svd
    real(kind=dp), dimension(:), Allocatable :: S
    complex(kind=dp), dimension(:,:), Allocatable :: U, VT
    
        !work vector for subroutine calculation
    complex(kind=dp), dimension(:), Allocatable :: work_vector
    
    
        !variables for function calls
    complex(kind=dp) :: ZDOTC							!complex dot product subroutine
    
    
        integer :: m_result, trgt								!measurement outcome (0 or 1)
    complex(kind=dp), dimension(:), Allocatable :: m_vector					!measurement vector
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! INPUT  VARIABLES !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  !Fixed variables governing program behaviour

  integer :: chain_length = 13

  complex(kind=dp), dimension(:), Allocatable :: init_state
  complex(kind=dp), dimension(:), Allocatable :: plus_state

  !Must be allocated so it can be used with kronecker product subroutine
  Allocate(plus_state(2), init_state(2))

  !init_state = (/ (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)/)
  init_state = (1.0_dp / SQRT(2.0_dp)) * (/ 1.0_dp, 1.0_dp/)
  plus_state = (1.0_dp / SQRT(2.0_dp)) * (/ 1.0_dp, 1.0_dp/)  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! I/O SETUP !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Open the measurements output file so that it can be reset
  !'replace' status removes file
  open(100, file='measurements.dat', status='replace', iostat=ierr)
    if (ierr/=0) stop 'Error in opening file measurements.dat'
  
  !Close the file again as it is not needed until
  !the measurement subroutine is called
  close(100)
  
  !Opens a file of measurement instructions to be read from
  !contains measurement type and phase information
  open(200, file='m_instructions.dat', iostat=ierr)
   if (ierr/=0) stop 'Error in opening file m_instructions.dat'
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! PROGRAM SETUP !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Initialised the seed for the fortran intrinsic random number
  !generator, currently unused
  !call init_random_seed
  
  !Allocate the initial state vector to the size of a single qubit
  Allocate(state_in(2))
  
  Allocate(m_vector(2))
  
  !Initialises the first state as input initial state
  state_in = init_state
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! MAIN SIMULATION !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  !Repeats the process depending on the length of qubits
  do i = 1, (chain_length - 1)
    
    
    !Allocates memory for next state, will always be four in this case
    Allocate(state_out(2**(i+1)))
    
    !Calls the subroutine to perform a kronecker product between the input state
    !and a plus state qubit, forming a cluster state
    call kronecker_product_complex_vector(state_in, plus_state, state_out)
    
    !Deallocates and resizes the input state for shifting of variables
    Deallocate(state_in)
    Allocate(state_in(2**(i+1)))
    
    state_in = state_out
    
    Deallocate(state_out)
    
  end do
  
  do i=1, (chain_length - 1)
  
    !forms a cz operator of size 4x4 and applies it to input state
    call cz_operation(state_in, chain_length, i, i+1)

  end do
  
  call multi_measurement(state_in, chain_length, chain_length - 1)
  
  
    !Prints data about measurements to file for use in
    !feed forward subroutine later
  open(100, file='measurements.dat', access='direct', recl=40, iostat=ierr, form='formatted')
    if (ierr/=0) stop 'Error in opening file measurements.dat'
    
    
  do i = 1, chain_length - 1

    read(100, fmt=data_format, rec=(chain_length-1)) m_type, trgt, m_result, m_phase
      
    call measurement_type(m_vector, m_type, m_phase, m_result)
    
    Allocate(state_out(2**(chain_length - i)))
  
    call known_rank_decomposition_complex(state_in, m_vector, state_out)
   
    Deallocate(state_in)
    
    if(i < (chain_length - 1)) then
      Allocate(state_in(2**(chain_length - i)))
      state_in = state_out
      Deallocate(state_out)
    end if
  
  end do
  
  close(100, iostat=ierr)
    if (ierr/=0) stop 'Error in closing file measurements.dat'
  
  close(200)
  
  print *, CMPLX(state_out, kind=4)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! FEED FORWARD !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Calls the feed forward subroutine that determines the desired
  !information state. 
  call feed_forward(state_out, 1, chain_length - 1, 1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!! DATA OUTPUT !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Prints the input state to standard output
  print *, CMPLX(state_out, kind=4)
  
  !Opens the output fidelity file, reports error and aborts program on failure
  open(100, file="fideity.dat",iostat=ierr)
    if (ierr/=0) stop 'Error in opening file fidelity.dat'
 
  fid_out = fidelity(state_out, init_state)
  
  !Writes output of fidelity function to file
  write(100, *) i, fid_out
  
  close(100)
  
  Deallocate(init_state, plus_state)
  
end program single_chain
