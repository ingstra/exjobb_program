program test_random_seed
  implicit none
  integer, allocatable :: seed(:)
  integer :: n
CALL RANDOM_SEED()
  call random_seed(size = n)
  allocate(seed(n))

  call random_seed(get=seed)
  write (*, *) seed
end program test_random_seed
