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

  !measurement time string
  character(len=1) :: m_type
  integer :: m_number
  
  !variables for file i/o
  real(kind=dp) :: fid_out
  real(kind=dp) :: m_phase

  !input and output state vectors
  complex(kind=dp), dimension(:), Allocatable :: state_in, state_out

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! INPUT  VARIABLES !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  !Fixed variables governing program behaviour

  integer :: chain_length = 1000000

  complex(kind=dp), dimension(:), Allocatable :: init_state
  complex(kind=dp), dimension(:), Allocatable :: plus_state

  !Must be allocated so it can be used with kronecker product subroutine
  Allocate(plus_state(2), init_state(2))

  init_state = (/ (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)/)
  !init_state = (1.0_dp / SQRT(2.0_dp)) * (/ 1.0_dp, 1.0_dp/)
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
  
  !Initialises the first state as input initial state
  state_in = init_state
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!! MAIN SIMULATION !!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Repeats the process depending on the length of qubits
  do i = 1, (chain_length - 1)
    
    !Allocates memory for next state, will always be four in this case
    Allocate(state_out(4))
    
    !Calls the subroutine to perform a kronecker product between the input state
    !and a plus state qubit, forming a cluster state
    call kronecker_product_complex_vector(state_in, plus_state, state_out)
    
    !Deallocates and resizes the input state for shifting of variables
    Deallocate(state_in)
    Allocate(state_in(4))

    !reassign input state to value of output state
    state_in = state_out
    
    !forms a cz operator of size 4x4 and applies it to input state
    call cz_operation(state_in, 2, 1, 2)

    !Adds error to input state
    !call add_error()

    !Deallocates and resizes output state for shifting of variables
    Deallocate(state_out)
    Allocate(state_out(2))
    
    !Reads measurement instructions from file
    !READ(200, *) m_number, m_type, m_phase
    
    m_number = i
    m_phase = 0.0_dp
    
    call general_measurement(state_out, state_in, 2, 'X', m_phase, 1, i)
    
    !Resizes state input to fit new size of output state
    Deallocate(state_in)
    Allocate(state_in(2))
    
    state_in = state_out
    
    Deallocate(state_out)
    
  end do
  
  close(200)
  
  !print *, CMPLX(state_in, kind=4)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!! FEED FORWARD !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Calls the feed forward subroutine that determines the desired
  !information state. 
  call feed_forward(state_in, 1, i-1, 1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!! DATA OUTPUT !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Prints the input state to standard output
  !print *, CMPLX(state_in, kind=4)
  
  !Opens the output fidelity file, reports error and aborts program on failure
  open(100, file="fideity.dat",iostat=ierr)
    if (ierr/=0) stop 'Error in opening file fidelity.dat'
 
  fid_out = fidelity(state_in, init_state)
  
  !Writes output of fidelity function to file
  write(100, *) i, fid_out
  
  close(100)
  
  Deallocate(state_in, init_state, plus_state)
  
end program single_chain
