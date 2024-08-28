module mpijobs_mod
  use mpi_f08
  implicit none
  private
  public :: mpijobs_init, mpijobs_finalize, wait_for_job
  public :: report_complete, dispatch_and_collect, UNASSIGNED, ROOT

  type(MPI_Comm), parameter :: comm = MPI_COMM_WORLD
  integer, parameter :: ROOT = 0
  integer, parameter :: TAG_COMPLETE = 11, TAG_ASSIGN = 12
  integer, parameter :: UNASSIGNED = -1

  integer, parameter :: DEBUG = 0
    !! DEBUG = 0, 1

contains

  subroutine mpijobs_init(rank, np)
    integer, intent(out) :: rank, np
    integer :: ierr
    call mpi_init(ierr)
    if (ierr /= 0) error stop 'mpi_init error'
    call mpi_comm_rank(comm, rank, ierr)
    if (ierr /= 0) error stop 'mpi_comm_rank error'
    call mpi_comm_size(comm, np, ierr)
    if (ierr /= 0) error stop 'mpi_comm_size error'
  end subroutine


  subroutine mpijobs_finalize()
    integer :: ierr
    call mpi_finalize(ierr)
    if (ierr /= 0) error stop 'mpi_finalize error'
  end subroutine


  subroutine wait_for_job(rank, jobindex, got_job)
!
! If worker then wait until receiving a new job from root.
! No more jobs is signalized by "jobindex=0" on output.
!
! If root then just check if mew job message has been send.
! If no message then return with "got_job==.false."
!
    integer, intent(in)  :: rank
    integer, intent(out) :: jobindex
    logical, intent(out) :: got_job
    type(MPI_STATUS) :: status
    logical :: flag
    integer :: ierr

    ! Avoid deadlocking the root if no message is in queue
    if (rank == ROOT) then
      call mpi_iprobe(ROOT, TAG_ASSIGN, comm, flag, status, ierr)
      if (.not. flag) then
        got_job = .false.
        return
      end if
    end if

    call mpi_recv(jobindex, 1, MPI_INTEGER, MPI_ANY_SOURCE, TAG_ASSIGN, &
      comm, status, ierr)
    if (ierr /= 0) error stop 'mpi_recv error'
    if (DEBUG > 0) print &
      '("Process ",i0," received job index ",i0,".")', rank, jobindex
    got_job = .true.
  end subroutine


  subroutine report_complete(rank, jobindex, able_to_report)
!
! Send completed job message to root
!
    integer, intent(in) :: jobindex
    logical, intent(out) :: able_to_report

    type(MPI_REQUEST), save :: send_req = MPI_REQUEST_NULL
    type(MPI_STATUS) :: status
    integer :: rank, ierr
    logical :: flag

    if (rank==ROOT) then
      ! non-blocking send for ROOT, send is not made if the send from
      ! the last call has not been yet received
      if (send_req /= MPI_REQUEST_NULL) then
        call mpi_test(send_req, flag, status, ierr)
        if (.not. flag) then
          able_to_report = .false.
          return
        end if
      end if
      able_to_report = .true.

      if (jobindex /= 0) then
        call mpi_issend(jobindex, 1, MPI_INTEGER, ROOT, TAG_COMPLETE, &
          comm, send_req, ierr)
        if (ierr /= 0) error stop 'mpi_issend error'
        if (DEBUG > 0) print &
          '("Root sended message complete ",i0," to root")', jobindex
      end if

    else
      ! blocking send for workers
      call mpi_ssend(jobindex, 1, MPI_INTEGER, ROOT, TAG_COMPLETE, &
        comm, ierr)
      if (ierr /= 0) error stop 'mpi_ssend error'
      if (DEBUG > 0) print &
        '("Process ",i0," sended message complete ",i0," to root")', rank, jobindex
    end if
  end subroutine


  subroutine dispatch_and_collect(rank, np, jobs, is_end)
    integer, intent(in) :: rank, np
    integer, intent(inout) :: jobs(:)
    logical, intent(out) :: is_end
!
! Receive "job_complete" messages and send "assign_job" messages
! Array "jobs" keeps track to which processes have been jobs assigned,
! Set it to "UNASSIGNED" initially.
! Root will send an index of jobs array or "0" if no more jobs
! are left
!
    integer :: ierr
    logical :: flag
    type(MPI_STATUS) :: status
    logical, save, allocatable :: is_working(:), is_finished(:)
    type(MPI_REQUEST), save, allocatable :: send_req(:)
    integer :: message, i, current_job_index

    is_end = .false.
    if (rank /= ROOT) return

    ! root node only
    if (.not. allocated(send_req)) &
      allocate(send_req(0:np-1), source=MPI_REQUEST_NULL)
    if (.not. allocated(is_working)) &
      allocate(is_working(0:np-1), source=.false.)
    if (.not. allocated(is_finished)) &
      allocate(is_finished(0:np-1), source=.false.)

    ! process all pending "completed" messages
    RECEIVE_LOOP: do
      call mpi_iprobe(MPI_ANY_SOURCE, TAG_COMPLETE, comm, flag, status, ierr)
      if (ierr /= 0) error stop 'mpi_iprobe error'
      if (.not. flag) exit RECEIVE_LOOP

      ! mark node as not_working after receiving "completed" message
      call mpi_recv(message, 1, MPI_INTEGER, &
        MPI_ANY_SOURCE, TAG_COMPLETE, comm, status, ierr)
      if (ierr /=0) error stop 'mpi_recv error'
      associate(src => status%MPI_SOURCE)
        if (DEBUG > 0) print &
          '("Root received complete ",i0," from ",i0)', message, src
        ! verify message is correct
        if (.not. is_working(src)) &
          error stop 'bogus complete message, process is not working'
        if (findloc(jobs, src, back=.true., dim=1) /= message) then
          print *, 'message / source ', message, src
          print *, 'jobs ', jobs
          error stop 'bogus complete message, wrong jobindex'
        end if
        is_working(src) = .false.
      end associate
    end do RECEIVE_LOOP

    ! clear requests from a previous run
    do i=0, np-1
      if (send_req(i)==MPI_REQUEST_NULL) cycle
      if (i==ROOT) then
        call mpi_test(send_req(i), flag, status, ierr)
        if (.not. flag .and. DEBUG>0) &
          print *, 'send to root hanging during this call'
      else
        call mpi_wait(send_req(i), status, ierr)
      end if
      if (ierr /= 0) error stop 'mpi_wait/test error'
    end do

    ! send to free processes jobs
    JOB_LOOP: do
      current_job_index = findloc(jobs, UNASSIGNED, dim=1)

      do i=0, np-1
        ! skip busy or finished processes
        if (is_working(i) .or. is_finished(i)) cycle
        if (send_req(i) /= MPI_REQUEST_NULL) cycle

        ! send assigning message and exit
        call mpi_issend(current_job_index, 1, MPI_INTEGER, i, TAG_ASSIGN, &
            comm, send_req(i), ierr)
        if (ierr /= 0) error stop 'mpi_issend error'
        is_working(i) = .true.
        if (current_job_index /= 0) then
          jobs(current_job_index) = i
          if (DEBUG > 0) print &
            '("Task no ",i0," assigned to process ",i0,".")', current_job_index, i
        else
          is_finished(i) = .true.
print '("No more tasks send to process ",i0,".")', i
        end if
        exit
      end do
      ! all processes are now busy exit
      if (i==np) exit JOB_LOOP
    end do JOB_LOOP
    is_end = all(is_finished)

  end subroutine
end module
