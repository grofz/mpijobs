! =============================================
! This is an example of using MPIJOBS scheduler
! =============================================

program main
  use mpijobs_mod, only : mpijobs_t, S_NEWJOB, S_WORKING, S_FINISHED
  implicit none

  type(mpijobs_t) :: mp
  integer :: jobindex, jobs(30)
  logical :: is_end

  ! initialize is compulsory
  call mp%init(jobs) ! jobs are set to unassigned on root

  print *, 'hello from ', mp%myrank()

  MAINLOOP: do
    ! root is sending "assigned" and receiving "completed" messages here
    call mp%dispatcher(jobs, is_end)
    if (is_end) exit MAINLOOP ! root process exits here

    ! receive new job from root
    call mp%receive_job(jobindex)

    ! non-root processes exit here
    if (mp%getstate()==S_FINISHED .and. .not. mp%amiroot()) exit MAINLOOP

    ! job initialization
    if (mp%getstate()==S_NEWJOB) then
      ! ... HERE ADD JOB INITIALIZATION CODE ...
      call mp%set_working
    end if

    ! job processing
    if (mp%getstate()==S_WORKING) then
      ! ... HERE ADD JOB RUNNING CODE ...
      print *, 'hello from ',mp%myrank(),' doiing job ', jobindex
      ! when job is done, call "set_completed"
      call mp%set_completed()
    end if

    ! send job completed message to root
    call mp%send_completed(jobindex)
  end do MAINLOOP

  ! who did what? (only root has this information)
  if (mp%amiroot()) print *, jobs

  ! clear non-blocking send_req and mpi
  call mp%finalize()

end program
