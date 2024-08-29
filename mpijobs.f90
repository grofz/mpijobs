module mpijobs_mod
  use mpi_f08
  implicit none
  private
  public :: UNASSIGNED, ROOT

  integer, parameter :: ROOT = 0
  integer, parameter :: TAG_COMPLETED = 11, TAG_ASSIGNED = 12
  integer, parameter :: UNASSIGNED = -1

  ! Set this to switch on/off printing send/received messages
  integer, parameter :: DEBUG = 0
    !! DEBUG = 0, 1

  integer, parameter, public :: &
    S_WAITING   = 0, & ! waiting for a new job
    S_NEWJOB    = 1, & ! received a new job (or signal no more jobs)
    S_WORKING   = 2, & ! is working on the job
    S_COMPLETED = 3    ! job completed, has to send completed message

  type, public :: mpijobs_t
    private
    type(MPI_Comm) :: comm
    integer :: rank=-1, np=-1
    integer :: jobindex=1 ! index to first unassigned job in jobs array
    type(MPI_Request) :: completed_send_req
    type(MPI_Request), allocatable :: assigned_send_reqs(:)
    logical, allocatable :: is_working(:), is_finished(:)
    integer :: state = S_WAITING
  contains
    procedure :: myrank, getstate, set_working, set_completed
    procedure :: init => mpijobs_init
    procedure :: finalize => mpijobs_finalize
    procedure :: receive_job
    procedure :: send_completed
    procedure :: dispatcher
  end type

contains

  pure function myrank(this)
    class(mpijobs_t), intent(in) :: this
    integer myrank
    myrank = this%rank
  end function


  pure function getstate(this)
    class(mpijobs_t), intent(in) :: this
    integer :: getstate
    getstate = this%state
  end function


  pure subroutine set_working(this)
    class(mpijobs_t) ,intent(inout) :: this
    if (this%state /= S_NEWJOB) error stop 'should not call set_working'
    this%state = S_WORKING
  end subroutine


  pure subroutine set_completed(this)
    class(mpijobs_t) ,intent(inout) :: this
    if (this%state /= S_WORKING) error stop 'should not call set_completed'
    this%state = S_COMPLETED
  end subroutine


  subroutine mpijobs_init(this)
    class(mpijobs_t), intent(inout) :: this
    integer :: ierr

    this%comm = MPI_COMM_WORLD
    call mpi_init(ierr)
    if (ierr /= 0) error stop 'mpi_init error'

    call mpi_comm_rank(this%comm, this%rank, ierr)
    if (ierr /= 0) error stop 'mpi_comm_rank error'

    call mpi_comm_size(this%comm, this%np, ierr)
    if (ierr /= 0) error stop 'mpi_comm_size error'

    this%jobindex = 1 ! first unassigned job in the index
    this%completed_send_req = MPI_REQUEST_NULL
    if (allocated(this%assigned_send_reqs)) deallocate(this%assigned_send_reqs)
    allocate(this%assigned_send_reqs(0:this%np-1), source=MPI_REQUEST_NULL)
    if (allocated(this%is_working)) deallocate(this%is_working)
    allocate(this%is_working(0:this%np-1), source=.false.)
    if (allocated(this%is_finished)) deallocate(this%is_finished)
    allocate(this%is_finished(0:this%np-1), source=.false.)
  end subroutine


  subroutine mpijobs_finalize(this)
    class(mpijobs_t), intent(inout) :: this
    integer :: ierr, i
    type(MPI_STATUS) :: status

    do i=lbound(this%assigned_send_reqs,1),ubound(this%assigned_send_reqs,1)
      if (this%assigned_send_reqs(i) == MPI_REQUEST_NULL) cycle
 print *, 'finalize - waiting to clear send request', i
      call mpi_wait(this%assigned_send_reqs(i), status, ierr)
      if (ierr /= 0) error stop 'mpi_wait error'
    end do
    if (this%completed_send_req /= MPI_REQUEST_NULL) then
 print *, 'finalize - waiting to clear send request'
      call mpi_wait(this%completed_send_req, status, ierr)
      if (ierr /= 0) error stop 'mpi_wait error'
    end if

    deallocate(this%assigned_send_reqs, this%is_working, this%is_finished)

    call mpi_finalize(ierr)
    if (ierr /= 0) error stop 'mpi_finalize error'
  end subroutine


  subroutine receive_job(this, jobindex)
!
! Receive TAG_ASSIGNED messages from root.
! "jobindex==0" on output signalizes no more jobs can be assigned.
!
! For root process: "got_job==.false." on output means that no message
! was received.
! Other processes: This is a blocking subroutine and "got_job" is always
! true on return.
!
    class(mpijobs_t), intent(inout) :: this
    integer, intent(out) :: jobindex

    type(MPI_STATUS) :: status
    logical :: flag
    integer :: ierr

    ! Is in the correct state?
    if (this%state /= S_WAITING) return

    ! Avoid deadlocking the root if no message is waiting
    if (this%rank == ROOT) then
      call mpi_iprobe(ROOT, TAG_ASSIGNED, this%comm, flag, status, ierr)
      if (ierr /= 0) error stop 'mpi_iprobe error'
      if (.not. flag) then
        return
      end if
    end if

    call mpi_recv(jobindex, 1, MPI_INTEGER, &
      ROOT, TAG_ASSIGNED, this%comm, status, ierr)
    if (ierr /= 0) error stop 'mpi_recv error'
    if (DEBUG > 0) print &
      '("Process ",i0," received job index ",i0,".")', this%rank, jobindex
    this%state = S_NEWJOB
  end subroutine receive_job


  subroutine send_completed(this, jobindex)
!
! Send "completed" message to root.
!
! For root process: "able_to_send=.false." signalizes that the send from the
! previous call has not yet been received, and nothing has been sent during
! this call.
! Other processes: This is a blocking subroutine and "able_to_send" is
! always true on return.
!
    class(mpijobs_t), intent(inout) :: this
    integer, intent(in) :: jobindex

    type(MPI_STATUS) :: status
    integer :: ierr
    logical :: flag

    if (this%state /= S_COMPLETED) return

    if (this%rank==ROOT) then
      ! non-blocking send for ROOT, send is not made if the send from
      ! the last call has not been yet received
      if (this%completed_send_req /= MPI_REQUEST_NULL) then
        call mpi_test(this%completed_send_req, flag, status, ierr)
        if (ierr /= 0) error stop 'mpi_test error'
        if (.not. flag) then
          return
        end if
      end if

      if (jobindex /= 0) then
        call mpi_issend(jobindex, 1, MPI_INTEGER, ROOT, TAG_COMPLETED, &
          this%comm, this%completed_send_req, ierr)
        if (ierr /= 0) error stop 'mpi_issend error'
        if (DEBUG > 0) print &
          '("Process ",i0," sended message complete ",i0," to root")', this%rank, jobindex
      end if

    else
      ! blocking send for workers
      call mpi_ssend(jobindex, 1, MPI_INTEGER, ROOT, TAG_COMPLETED, &
        this%comm, ierr)
      if (ierr /= 0) error stop 'mpi_ssend error'
      if (DEBUG > 0) print &
        '("Process ",i0," sended message complete ",i0," to root")', this%rank, jobindex
    end if
    this%state = S_WAITING
  end subroutine send_completed


  subroutine dispatcher(this, jobs, is_end)
    class(mpijobs_t), intent(inout) :: this
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
    integer :: message, i

    ! root node only
    is_end = .false.
    if (this%rank /= ROOT) return

    ! process all pending "completed" messages
    RECEIVE_LOOP: do
      call mpi_iprobe(MPI_ANY_SOURCE, TAG_COMPLETED, this%comm, flag, &
        status, ierr)
      if (ierr /= 0) error stop 'mpi_iprobe error'

      if (.not. flag) exit RECEIVE_LOOP

      ! mark node as not_working after receiving "completed" message
      call mpi_recv(message, 1, MPI_INTEGER, &
        MPI_ANY_SOURCE, TAG_COMPLETED, this%comm, status, ierr)
      if (ierr /=0) error stop 'mpi_recv error'

      if (DEBUG > 0) print &
        '("Root received completed ",i0," from ",i0)', message, status%MPI_SOURCE
      ! verify message is correct
      if (.not. this%is_working(status%MPI_SOURCE)) &
        error stop 'bogus complete message, process is not working'
      if (findloc(jobs, status%MPI_SOURCE, back=.true., dim=1) /= message) then
        print *, 'message / source ', message, status%MPI_SOURCE
        print *, 'jobs ', jobs
        error stop 'bogus complete message, wrong jobindex'
      end if
      this%is_working(status%MPI_SOURCE) = .false.
    end do RECEIVE_LOOP

    ! clear requests from a previous run
    do i=0, this%np-1
      if (this%assigned_send_reqs(i)==MPI_REQUEST_NULL) cycle
      if (i==ROOT) then
        call mpi_test(this%assigned_send_reqs(i), flag, status, ierr)
        if (.not. flag .and. DEBUG>0) &
          print *, 'send to root hanging during this call'
      else
        call mpi_wait(this%assigned_send_reqs(i), status, ierr)
      end if
      if (ierr /= 0) error stop 'mpi_wait/test error'
    end do

    ! send to free processes jobs
    JOB_LOOP: do
      if (this%jobindex > size(jobs)) this%jobindex = 0
      do i=0, this%np-1
        ! skip busy or finished processes, or processes with incomplete sends
        if (this%is_working(i) .or. this%is_finished(i) .or. &
            this%assigned_send_reqs(i) /= MPI_REQUEST_NULL) cycle

        ! send assigning message and exit
        call mpi_issend(this%jobindex, 1, MPI_INTEGER, i, TAG_ASSIGNED, &
            this%comm, this%assigned_send_reqs(i), ierr)
        if (ierr /= 0) error stop 'mpi_issend error'

        this%is_working(i) = .true.

        if (this%jobindex /= 0) then
          jobs(this%jobindex) = i
          if (DEBUG > 0) print &
            '("Task no ",i0," assigned to process ",i0,".")', this%jobindex, i
          this%jobindex = this%jobindex+1
        else
          this%is_finished(i) = .true.
        end if
        exit
      end do
      ! all processes are now busy, exit
      if (i == this%np) exit JOB_LOOP
    end do JOB_LOOP

    is_end = all(this%is_finished) .and. &
             all(request_null(this%assigned_send_reqs))

  end subroutine dispatcher


  impure elemental function request_null(a)
    type(MPI_REQUEST), intent(in) :: a
    logical request_null
    request_null = a==MPI_REQUEST_NULL
  end function request_null

end module
