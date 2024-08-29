!
! This is a simple scheduler of jobs among processes
!
! See "example" code in main.f90 how to use this module
!
! (c)2024 Z. Grof (grofz@vscht.cz)
!

module mpijobs_mod
  use mpi_f08
  implicit none
  private

  ! Set this to switch on/off printing send/received messages
  integer, parameter :: DEBUG = 1
    !! DEBUG = 0, 1-warn, 2-info

  ! process states used to control flow in the event loop "state"
  !
  !    ------------------------------------
  !    |                                  |
  !    v                                  |
  ! WAITING --> NEWJOB --> WORKING --> COMPLETED
  !    |
  !    -------> FINISHED
  !
  integer, parameter, public :: &
    S_WAITING   = 0, & ! waiting for a new job (initial state)
    S_NEWJOB    = 1, & ! received a new job
    S_WORKING   = 2, & ! is working on the job
    S_COMPLETED = 3, & ! job completed, has to send a "completed" message
    S_FINISHED  = 4    ! received signal that there will be no more jobs

  type, public :: mpijobs_t
    private
    type(MPI_Comm) :: comm
    integer :: rank=-1, np=-1
    integer :: jobindex=1 ! index to first unassigned job in jobs array
    type(MPI_Request) :: completed_send_req
    type(MPI_Request), allocatable :: assigned_send_reqs(:)
    integer, allocatable :: pstates(:)
    integer :: state = S_WAITING
  contains
    procedure :: myrank, getstate, set_working, set_completed, amiroot
    procedure :: init => mpijobs_init
    procedure :: finalize => mpijobs_finalize
    procedure :: receive_job
    procedure :: send_completed
    procedure :: dispatcher
  end type

  integer, parameter :: ROOT = 0
  integer, parameter :: TAG_COMPLETED = 11, TAG_ASSIGNED = 12
  integer, parameter :: UNASSIGNED = -1     ! taken values are 0..np-1
  integer, parameter :: NULL_JOBINDEX = 0   ! taken values are 1..size(jobs)

  ! process states used in dispatcher's "pstates" array
  !
  !  ------------------a--------------
  !  | ----b-----                    |
  !  | |        |                    |
  !  v v        |                    |
  ! FREE --> ASSI_SENT --> WORKING --|
  !    |
  !    |---> FINI_SENT --> FINISHED
  !
  integer, parameter :: &
    P_FREE          = 0, & ! available processes (initial state)
    P_ASSIGNED_SENT = 1, & ! "assign job" message not yet completed
    P_WORKING       = 2, & ! "assign job" message completed
    P_FINISHED_SENT = 3, & ! "finished" message not yet completed
    P_FINISHED      = 4    ! "finished" message completed (final state)

contains

  pure function amiroot(this)
    class(mpijobs_t), intent(in) :: this
    logical amiroot
    amiroot = this%rank == ROOT
  end function


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


  subroutine mpijobs_init(this, jobs)
    class(mpijobs_t), intent(inout) :: this
    integer, intent(inout) :: jobs(:)
    integer :: ierr

    this%comm = MPI_COMM_WORLD
    call mpi_init(ierr)
    if (ierr /= 0) error stop 'mpi_init error'

    call mpi_comm_rank(this%comm, this%rank, ierr)
    if (ierr /= 0) error stop 'mpi_comm_rank error'

    call mpi_comm_size(this%comm, this%np, ierr)
    if (ierr /= 0) error stop 'mpi_comm_size error'

    this%completed_send_req = MPI_REQUEST_NULL

    if (this%rank==ROOT) then
      jobs = UNASSIGNED
      this%jobindex = 1 ! first unassigned job in the index
      if (allocated(this%assigned_send_reqs)) deallocate(this%assigned_send_reqs)
      allocate(this%assigned_send_reqs(0:this%np-1), source=MPI_REQUEST_NULL)
      if (allocated(this%pstates)) deallocate(this%pstates)
      allocate(this%pstates(0:this%np-1), source=P_FREE)
    end if
  end subroutine mpijobs_init


  subroutine mpijobs_finalize(this)
    class(mpijobs_t), intent(inout) :: this
    integer :: ierr, i
    type(MPI_STATUS) :: status

    if (this%completed_send_req /= MPI_REQUEST_NULL) then
      if (DEBUG > 0) print &
        '("finalized rank ",i0," waiting to clear request")',this%rank
      call mpi_wait(this%completed_send_req, status, ierr)
      if (ierr /= 0) error stop 'mpi_wait error'
    end if

    if (this%rank==ROOT) then
      do i=lbound(this%assigned_send_reqs,1),ubound(this%assigned_send_reqs,1)
        if (this%assigned_send_reqs(i) == MPI_REQUEST_NULL) cycle
        if (DEBUG > 0) print &
          '("finalized root waiting to clear request to rank ",i0,".")', i
        call mpi_wait(this%assigned_send_reqs(i), status, ierr)
        if (ierr /= 0) error stop 'mpi_wait error'
      end do

      deallocate(this%assigned_send_reqs, this%pstates)
    end if

    call mpi_finalize(ierr)
    if (ierr /= 0) error stop 'mpi_finalize error'
  end subroutine mpijobs_finalize


  subroutine receive_job(this, jobindex)
!
! Receive TAG_ASSIGNED messages from root.
! State is changed from S_WAITING to S_NEWJOB or S_FINISHED
! if a message is waiting and was received.
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
    if (DEBUG > 1) print &
      '("Process ",i0," received job index ",i0,".")', this%rank, jobindex
    if (jobindex==NULL_JOBINDEX) then
      this%state = S_FINISHED
    else
      this%state = S_NEWJOB
    end if
  end subroutine receive_job


  subroutine send_completed(this, jobindex)
!
! Send "completed" message to root.
! State is changed from S_COMPLETED to S_WAITING if "completed" message
! has been send. If last message was not yet received, nothing is done
!
    class(mpijobs_t), intent(inout) :: this
    integer, intent(in) :: jobindex

    type(MPI_STATUS) :: status
    integer :: ierr
    logical :: flag

    ! Is in the correct state?
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

      call mpi_issend(jobindex, 1, MPI_INTEGER, ROOT, TAG_COMPLETED, &
        this%comm, this%completed_send_req, ierr)
      if (ierr /= 0) error stop 'mpi_issend error'
      if (DEBUG > 1) print &
        '("Process ",i0," sended message complete ",i0," to root")', this%rank, jobindex

    else
      ! blocking send for workers
      call mpi_ssend(jobindex, 1, MPI_INTEGER, ROOT, TAG_COMPLETED, &
        this%comm, ierr)
      if (ierr /= 0) error stop 'mpi_ssend error'
      if (DEBUG > 1) print &
        '("Process ",i0," sended message complete ",i0," to root")', this%rank, jobindex
    end if
    this%state = S_WAITING
  end subroutine send_completed


  subroutine dispatcher(this, jobs, is_end)
    class(mpijobs_t), intent(inout) :: this
    integer, intent(inout) :: jobs(:)
    logical, intent(out) :: is_end
!
! Receive "completed" messages and send "assigned" messages.
! Array "jobs" keeps track to which processes have been jobs assigned,
! Must be set to "UNASSIGNED" in "init" procedure.
! Root will send an index of jobs array or NULL_JOBINDEX if no more jobs
! are left.
!
    integer :: ierr
    logical :: flag
    type(MPI_STATUS) :: status
    integer :: message, i

    ! root node only
    is_end = .false.
    if (this%rank /= ROOT) return

    ! receive all pending "completed" messages
    RECEIVE_LOOP: do
      call mpi_iprobe(MPI_ANY_SOURCE, TAG_COMPLETED, this%comm, flag, &
        status, ierr)
      if (ierr /= 0) error stop 'mpi_iprobe error'

      if (.not. flag) exit RECEIVE_LOOP

      ! mark node as not_working after receiving "completed" message
      call mpi_recv(message, 1, MPI_INTEGER, &
        MPI_ANY_SOURCE, TAG_COMPLETED, this%comm, status, ierr)
      if (ierr /=0) error stop 'mpi_recv error'

      if (DEBUG > 1) print &
        '("Root received completed ",i0," from ",i0)', message, status%MPI_SOURCE

      ! verify message is correct
      if (message == NULL_JOBINDEX) &
        error stop 'bogus complete message, NULL_JOBINDEX received'
      if (findloc(jobs, status%MPI_SOURCE, back=.true., dim=1) /= message) then
        print '("message ",i0,"  source ",i0)', message, status%MPI_SOURCE
        print '("jobs array: ",/,10(i0,1x))', jobs
        error stop 'bogus message, received an unexpected value of jobindex'
      end if
      if (this%pstates(status%MPI_SOURCE) == P_WORKING) then
        ! this is expected (branch a)
        continue
      else if (this%pstates(status%MPI_SOURCE) == P_ASSIGNED_SENT) then
        ! this is rare, but possible (branch b)
        ! job was completed before root cleared send request
        ! must clear it now
        call mpi_test(this%assigned_send_reqs(status%MPI_SOURCE), &
          flag, status, ierr)
        if (ierr /= 0) error stop 'mpi_test error'
        if (.not. flag) &
          error stop 'bogus complete message, process did not receive assign'
      else
        error stop 'bogus complete message, process is not working'
      end if

      this%pstates(status%MPI_SOURCE) = P_FREE
    end do RECEIVE_LOOP

    ! clear requests from a previous run
    do i=0, this%np-1
      if (this%assigned_send_reqs(i)==MPI_REQUEST_NULL) cycle

      ! assert we are in "P_ASSIGNED/FINISHED_SEND" state
      if (this%pstates(i)/=P_ASSIGNED_SENT .and. &
          this%pstates(i)/=P_FINISHED_SENT) &
          error stop 'unexpect state in this context'

      if (i==ROOT) then
        call mpi_test(this%assigned_send_reqs(i), flag, status, ierr)
      else
        call mpi_wait(this%assigned_send_reqs(i), status, ierr)
        flag = .true.
      end if
      if (.not. flag) then
        if (DEBUG>0) print &
          '("Send assigned to rank ",i0," is still hanging")',i
      else
        ! send completed, update state
        select case(this%pstates(i))
        case(P_ASSIGNED_SENT)
          this%pstates(i) = P_WORKING
        case(P_FINISHED_SENT)
          this%pstates(i) = P_FINISHED
        case default
          error stop 'unexpect state in this context'
        end select
      end if
      if (ierr /= 0) error stop 'mpi_wait/test error'
    end do

    ! send unassigned jobs to free processes
    JOBS_LOOP: do
      if (this%jobindex > size(jobs)) this%jobindex = NULL_JOBINDEX
      PROCS: do i=0, this%np-1
        ! skip busy or finished processes, or processes with incomplete sends
        if (this%pstates(i) /= P_FREE) cycle PROCS

        ! defeensive check
        if (this%assigned_send_reqs(i) /= MPI_REQUEST_NULL) &
          error stop 'assigned_send_reqs not null where ecpected to be null'

        ! send assigning message and exit
        call mpi_issend(this%jobindex, 1, MPI_INTEGER, i, TAG_ASSIGNED, &
            this%comm, this%assigned_send_reqs(i), ierr)
        if (ierr /= 0) error stop 'mpi_issend error'

        if (this%jobindex /= NULL_JOBINDEX) then
          jobs(this%jobindex) = i
          if (DEBUG > 1) print &
            '("Task no ",i0," assigned to process ",i0,".")', this%jobindex, i
          this%jobindex = this%jobindex+1
          this%pstates(i) = P_ASSIGNED_SENT
        else
          this%pstates(i) = P_FINISHED_SENT
        end if
        exit PROCS
      end do PROCS
      ! we looped over all processes, exit
      if (i == this%np) exit JOBS_LOOP
    end do JOBS_LOOP

    is_end = all(this%pstates==P_FINISHED)

  end subroutine dispatcher

  ! WE ARE NOT USING THIS AT THE MOMENT
  ! because == operator is not elemental for MPI_REQUEST objects,
  ! this is an in-house fix...
  impure elemental function request_null(a)
    type(MPI_REQUEST), intent(in) :: a
    logical request_null
    request_null = a==MPI_REQUEST_NULL
  end function request_null

end module mpijobs_mod
