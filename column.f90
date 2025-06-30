
program column  
    implicit none
    include 'mpif.h'
    integer, parameter :: dp = selected_real_kind(15)
    integer, parameter :: bigM = 20, bigN = 20, T=80
    integer :: N 
    integer :: i 
    integer :: ierr, rank, numproc, new_comm
    real(dp) :: start_time
    integer :: l_rank, r_rank
    integer, allocatable :: domain(:, :)

    call INTIALIZE(ierr, rank, numproc, new_comm, domain, &
    l_rank, r_rank, bigM, bigN, N, T)

    call print_full_grid(domain, N, bigM, bigN, ierr, rank, numproc, new_comm, t=0)
    if (rank==0)    start_time = MPI_WTIME()

    do i = 1, T
        call seance(domain, bigM, N, ierr, rank, l_rank, r_rank) ! cause you talk to ghost cells
        ! Game of life rules
        domain(1:bigM, 1:N) = domain(0:bigM-1, 0:N-1) + domain(0:bigM-1, 1:N) + domain(0:bigM-1, 2:N+1) + & ! 3x3 convolution
                        domain(1:bigM, 0:N-1) + 9*domain(1:bigM, 1:N) +   domain(1:bigM, 2:N+1) + &
                        domain(2:bigM+1, 0:N-1) + domain(2:bigM+1, 1:N) + domain(2:bigM+1, 2:N+1)
        where((domain==3) .or. (domain==11) .or. (domain==12)) domain = 18 ! evil logic       
        domain = domain/18 ! integer division for more evil logic

        if ((i==20) .or. (i==40)) then
            call print_full_grid(domain, N, bigM, bigN, ierr, rank, numproc, new_comm, t=i)
        end if
    end do
    if (rank==0)    print *, 'time elapsed = ', MPI_WTIME()-start_time
    call print_full_grid(domain, N, bigM, bigN, ierr, rank, numproc, new_comm, t=T)

    call FINALIZE(ierr, domain)

contains

    subroutine intialize(ierr, rank, numproc, NEW_COMM, domain, &
        l_rank, r_rank, bigM, bigN, N, T)

        integer, intent(in) :: T, bigM, bigN
        integer, intent(out) :: ierr, rank, numproc, NEW_COMM
        integer, intent(out) :: N, l_rank, r_rank
        integer, allocatable, intent(out), dimension(:, :) :: domain
        integer :: dummynumproc, full_domain(bigM, bigN), i, j

        ! MPI initialize
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, dummynumproc, ierr )

        numproc = dummynumproc
        do while((mod(bigN, numproc) .ne. 0))
            numproc = numproc-1
        end do
            ! make a comm split. This helps with gathering for printing the whole grid
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, merge(1, 0, (rank < numproc)), rank, NEW_COMM, ierr)
        ! stop all unused processors
        if (numproc <= rank) then
            print *, 'stopping ', rank
            call MPI_Finalize(ierr)
            stop 
        end if
        N = bigN/numproc
        if (rank==0) print *, 'numproc = ', numproc, 'N = ', N
        print *, 'running with ', rank

    
        ! initial conditions
        allocate(domain(0:bigM+1, 0:N+1))
        domain = 0  
        full_domain = 0
        full_domain(1:4, 1:4) = transpose(reshape([0, 0, 1, 0, & 
                                                   1, 0, 1, 0, &
                                                   0, 1, 1, 0, &
                                                   0, 0, 0, 0 ], (/4, 4/)))

        domain(1:bigM, 1:N) = full_domain( 1:bigM, (rank*N+1):((rank+1)*N)   )

        ! if (rank==0) then
        ! do i = 1, bigM
        !     do j = 1, N
        !         write(*,'(I1)', advance='no') domain(i, j)
        !     end do
        !     write(*, *) ''
        ! end do
        ! end if  

        ! Set up boundary communcation
        l_rank = mod(numproc + rank-1, numproc) 
        r_rank = mod(rank+1, numproc) 
    end subroutine intialize

    subroutine finalize(ierr, domain)
        integer, intent(out) :: ierr
        integer, allocatable, intent(inout) :: domain(:, :)
        deallocate(domain)
        call MPI_FINALIZE( ierr)
    end subroutine finalize

    subroutine seance(domain, bigM, N, ierr, rank, l_rank, r_rank)
        integer, intent(in) :: bigM, N, rank, l_rank, r_rank
        integer, intent(inout) :: ierr, domain(0:bigM+1, 0:N+1)
        integer :: stat(MPI_STATUS_SIZE), requests(4), i, j
        ! send left,    recv right
        call MPI_Isend(domain(1:bigM, 1),   bigM, MPI_INTEGER, l_rank,  rank*10+3,    MPI_COMM_WORLD, requests(1), ierr)
        call MPI_Irecv(domain(1:bigM, N+1), bigM, MPI_INTEGER, r_rank,  r_rank*10+3,  MPI_COMM_WORLD, requests(2), ierr)

        ! send right,   recv left
        call MPI_Isend(domain(1:bigM, N),   bigM, MPI_INTEGER, r_rank,  rank*10+4,    MPI_COMM_WORLD, requests(3), ierr)
        call MPI_Irecv(domain(1:bigM, 0),   bigM, MPI_INTEGER, l_rank,  l_rank*10+4,  MPI_COMM_WORLD, requests(4), ierr)

        call MPI_Waitall(4, requests, MPI_STATUSES_IGNORE, ierr)
        
        !up
        domain(0, 0:N+1) = domain(bigM, 0:N+1)
        !down
        domain(bigM+1, 0:N+1) = domain(1, 0:N+1)



    end subroutine seance

    subroutine print_full_grid(domain, N, bigM, bigN, ierr, rank, numproc, new_comm, t)
        integer, intent(in) :: domain(0:bigM+1, 0:N+1), N, bigM, bigN, ierr, rank, numproc, new_comm, t
        integer :: i, j
        integer :: full_domain(bigM, bigN)
        integer, allocatable :: recv_buffer(:, :, :)
        allocate(recv_buffer(bigM, N, numproc))

    
        call MPI_Gather(domain(1:bigM, 1:N), bigM*N, MPI_INTEGER, recv_buffer, bigM*N, MPI_INTEGER, 0, new_comm, ierr)

        if (rank == 0) then
            do i = 1, numproc
                full_domain( 1:bigM, ((i-1)*N+1):((i)*N)   ) = recv_buffer(:, :, i)
            end do

            print *, "Grid at t = ", t
            do i = 1, bigM
                do j = 1, bigN
                    write(*,'(I1)', advance='no') full_domain(i, j)
                end do
                write(*, *) ''
            end do
        end if
        deallocate(recv_buffer)
    end subroutine print_full_grid

end program
