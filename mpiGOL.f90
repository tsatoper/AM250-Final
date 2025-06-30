
program mpiGOL ! mpif90 mpiGOL.f90 -o mpiGOL.exe | mpirun -np X mpiGOL.exe  , X can be any positive value
    implicit none
    include 'mpif.h'
    integer, parameter :: dp = selected_real_kind(15)
    integer, parameter :: bigN = 20, T= 80! These can be any positive value, grid is (bigN x bigN)
    integer :: N 
    integer :: i 
    integer :: ierr, rank, numproc, sqrtNP, row_type, new_comm
    real(dp) :: start_time
    integer :: l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank
    integer, allocatable :: domain(:, :)

    call INTIALIZE(ierr, rank, numproc, sqrtNP, row_type, new_comm, domain, &
    l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank, bigN, N, T)

    call print_full_grid(domain, N, bigN, ierr, rank, numproc, new_comm, sqrtNP, t=0)

    if (rank==0)    start_time = MPI_WTIME()

    do i = 1, T
        call seance(domain, N, ierr, rank, row_type, l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank) ! cause you talk to ghost cells
        ! Game of life rules
        domain(1:N, 1:N) = domain(0:N-1, 0:N-1) + domain(0:N-1, 1:N) + domain(0:N-1, 2:N+1) + & ! 3x3 convolution
                        domain(1:N, 0:N-1) + 9*domain(1:N, 1:N) +   domain(1:N, 2:N+1) + &
                        domain(2:N+1, 0:N-1) + domain(2:N+1, 1:N) + domain(2:N+1, 2:N+1)
        where((domain==3) .or. (domain==11) .or. (domain==12)) domain = 18 ! evil logic       
        domain = domain/18 ! integer division for more evil logic

        if ((i==20) .or. (i==40)) then
            call print_full_grid(domain, N, bigN, ierr, rank, numproc, new_comm, sqrtNP, t=i)
        end if
    end do
    if (rank==0)    print *, 'time elapsed = ', MPI_WTIME()-start_time
    call print_full_grid(domain, N, bigN, ierr, rank, numproc, new_comm, sqrtNP, t=T)

    call FINALIZE(ierr, domain)

contains

    subroutine intialize(ierr, rank, numproc, sqrtNP, row_type, NEW_COMM, domain, &
        l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank, bigN, N, T)
        integer, intent(in) :: T, bigN
        integer, intent(out) :: ierr, rank, numproc, sqrtNP, row_type, NEW_COMM
        integer, intent(out) :: N, l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank
        integer, allocatable, intent(out), dimension(:, :) :: domain
        integer :: dummynumproc, full_domain(bigN, bigN), i, j

        ! MPI initialize
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, dummynumproc, ierr )

        ! Set up sizes and squares, reduce numproc to perfect square s.t. sqrt(numproc) divides bigN
        numproc = dummynumproc
        sqrtNP = int(sqrt(numproc*1.0_dp))
        do while((int(sqrtNP)**2 .ne. numproc) .or. (mod(bigN, sqrtNP) .ne. 0))
            numproc = numproc-1
            sqrtNP = int(sqrt(numproc*1.0_dp))
        end do
            ! make a comm split. This helps with gathering for printing the whole grid
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, merge(1, 0, (rank < numproc)), rank, NEW_COMM, ierr)
        ! stop all unused processors
        if (numproc <= rank) then
            print *, 'stopping ', rank
            call MPI_Finalize(ierr)
            stop 
        end if
        N = bigN/sqrtNP
        if (rank==0)    print *, 'numproc = ', numproc, 'sqrtNP = ', sqrtNP, 'N = ', N, 'T = ', T
        ! print *, 'running with ', rank

        ! make noncontinugous derived type
        call MPI_TYPE_VECTOR(N, 1, N+2, MPI_INTEGER, row_type, ierr)
        call MPI_TYPE_COMMIT(row_type, ierr)
    
        ! initial conditions    
        allocate(domain(0:N+1, 0:N+1))
        domain = 0
        full_domain = 0
        full_domain(1:4, 1:4) = transpose(reshape([0, 0, 1, 0, & 
                                                   1, 0, 1, 0, &
                                                   0, 1, 1, 0, &
                                                   0, 0, 0, 0 ], (/4, 4/)))
        do i = 0, numproc-1
            if (rank == i) domain(1:N, 1:N) = full_domain(  N*(i/sqrtNP)+[(j, j=1,N)], N*mod(i, sqrtNP)+[(j, j=1,N)]   )
        end do

        ! Set up boundary communcation
        l_rank = sqrtNP*(rank/sqrtNP) + mod(sqrtNP + rank-1, sqrtNP) 
        r_rank = sqrtNP*(rank/sqrtNP) + mod(rank+1, sqrtNP) 
        u_rank = mod(numproc+rank-sqrtNP, numproc) 
        d_rank = mod(rank+sqrtNP, numproc)

        lu_rank = mod(numproc+l_rank-sqrtNP, numproc) 
        ru_rank = mod(numproc+r_rank-sqrtNP, numproc) 
        ld_rank = mod(l_rank+sqrtNP, numproc)
        rd_rank = mod(r_rank+sqrtNP, numproc)   
    end subroutine intialize

    subroutine finalize(ierr, domain)
        integer, intent(out) :: ierr
        integer, allocatable, intent(inout) :: domain(:, :)
        deallocate(domain)
        call MPI_FINALIZE( ierr)
    end subroutine finalize

    subroutine seance(domain, N, ierr, rank, row_type, l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank)
        integer, intent(in) :: N, rank, row_type, l_rank, r_rank, u_rank, d_rank, lu_rank, ru_rank, ld_rank, rd_rank
        integer, intent(inout) :: ierr, domain(0:N+1, 0:N+1)
        integer :: stat(MPI_STATUS_SIZE), requests(16)

        ! send up, recv down
        call MPI_Isend(domain(1, 1),     1, row_type,    u_rank,  rank*10+1,    MPI_COMM_WORLD, requests(1), ierr) 
        call MPI_Irecv(domain(N+1, 1),   1, row_type,    d_rank,  d_rank*10+1,  MPI_COMM_WORLD, requests(9), ierr)
        ! send down,    recv up
        call MPI_Isend(domain(N, 1),     1, row_type,    d_rank,  rank*10+2,    MPI_COMM_WORLD, requests(2), ierr) 
        call MPI_Irecv(domain(0, 1),     1, row_type,    u_rank,  u_rank*10+2,  MPI_COMM_WORLD, requests(10),  ierr) 

        ! send left,    recv right
        call MPI_Isend(domain(1:N, 1),   N, MPI_INTEGER, l_rank,  rank*10+3,    MPI_COMM_WORLD, requests(3), ierr)
        call MPI_Irecv(domain(1:N, N+1), N, MPI_INTEGER, r_rank,  r_rank*10+3,  MPI_COMM_WORLD, requests(11), ierr)

        ! send right,   recv left
        call MPI_Isend(domain(1:N, N),   N, MPI_INTEGER, r_rank,  rank*10+4,    MPI_COMM_WORLD, requests(4), ierr)
        call MPI_Irecv(domain(1:N, 0),   N, MPI_INTEGER, l_rank,  l_rank*10+4,  MPI_COMM_WORLD, requests(12), ierr)

        ! send left up, recv right down
        call MPI_Isend(domain(1, 1),     1, MPI_INTEGER, lu_rank, rank*10+5,    MPI_COMM_WORLD, requests(5), ierr)
        call MPI_Irecv(domain(N+1, N+1), 1, MPI_INTEGER, rd_rank, rd_rank*10+5, MPI_COMM_WORLD, requests(13), ierr)

        ! send right up, recv left down
        call MPI_Isend(domain(1, N),     1, MPI_INTEGER, ru_rank, rank*10+6,    MPI_COMM_WORLD, requests(6), ierr)
        call MPI_Irecv(domain(N+1, 0),   1, MPI_INTEGER, ld_rank, ld_rank*10+6, MPI_COMM_WORLD, requests(14), ierr)

        ! send left down, recv left up
        call MPI_Isend(domain(N, 1),     1, MPI_INTEGER, ld_rank, rank*10+7,    MPI_COMM_WORLD, requests(7), ierr)
        call MPI_Irecv(domain(0, N+1),   1, MPI_INTEGER, ru_rank, ru_rank*10+7, MPI_COMM_WORLD, requests(15), ierr)

        ! send right down, recv left up
        call MPI_Isend(domain(N, N),     1, MPI_INTEGER, rd_rank, rank*10+8,    MPI_COMM_WORLD, requests(8), ierr)
        call MPI_Irecv(domain(0, 0),     1, MPI_INTEGER, lu_rank, lu_rank*10+8, MPI_COMM_WORLD, requests(16), ierr)

        call MPI_Waitall(16, requests, MPI_STATUSES_IGNORE, ierr)
    end subroutine seance

    subroutine print_full_grid(domain, N, bigN, ierr, rank, numproc, new_comm, sqrtNP, t)
        integer, intent(in) :: domain(0:N+1, 0:N+1), N, bigN, ierr, rank, numproc, new_comm, sqrtNP, t
        integer :: i, j
        integer :: full_domain(bigN, bigN)
        integer, allocatable :: recv_buffer(:, :, :)
        allocate(recv_buffer(N, N, numproc))

        call MPI_Gather(domain(1:N, 1:N), N**2, MPI_INTEGER, recv_buffer, N**2, MPI_INTEGER, 0, new_comm, ierr)

        if (rank == 0) then
            do i = 0, numproc-1
                full_domain(    N*(i/sqrtNP)+[(j, j=1,N)], N*mod(i, sqrtNP)+[(j, j=1,N)]   ) = recv_buffer(:, :, i+1)
            end do

            print *, "Grid at t = ", t
            do i = 1, bigN
                do j = 1, bigN
                    write(*,'(I1)', advance='no') full_domain(i, j)
                end do
                write(*, *) ''
            end do
        end if
        deallocate(recv_buffer)
    end subroutine print_full_grid

end program
