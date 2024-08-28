program main
  use mpijobs_mod
  implicit none

  integer :: rank, np, jobindex
  logical :: able_to_report, is_end, got_new_job
  integer :: jobs(30)

  jobs = UNASSIGNED
  call mpijobs_init(rank, np)
  print *, 'hello from rank ', rank
  able_to_report = .true.

  do
    call dispatch_and_collect(rank, np, jobs, is_end)
    if (is_end .and. rank==ROOT) exit

    if (able_to_report) call wait_for_job(rank, jobindex, got_new_job)

    if (jobindex==0) then
      if (rank /= ROOT) exit
    else
      print *, 'hello from ',rank,' doiing job ', jobindex
      if (got_new_job) &
        call report_complete(rank, jobindex, able_to_report)
    end if
  end do

  ! just clear non-blocking send_req
  if (rank==ROOT) &
    call report_complete(rank, jobindex, able_to_report)
  if (.not. able_to_report) error stop 'not able to finish last call'
  if (rank==ROOT) print *, jobs

  call mpijobs_finalize()

end program
