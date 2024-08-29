program main
  use mpijobs_mod
  implicit none

  type(mpijobs_t) :: this
  integer :: jobindex
  logical :: is_end
  integer :: jobs(30)

  jobs = UNASSIGNED

  ! initialize is compulsory
  call this%init()
  print *, 'hello from rank ', this%myrank()

  do
    call this%dispatcher(jobs, is_end)
    if (is_end) exit

    call this%receive_job(jobindex)

    if (this%getstate()==S_NEWJOB) then
      if (jobindex==0) then
        if (this%myrank() /= ROOT) exit
      else
        call this%set_working
      end if
    end if

    if (this%getstate()==S_WORKING) then
      print *, 'hello from ',this%myrank(),' doiing job ', jobindex
      call this%set_completed()
    end if

    call this%send_completed(jobindex)
  end do

  ! report who did what
  if (this%myrank()==ROOT) print *, jobs

  ! clear non-blocking send_req
  call this%finalize()

end program
