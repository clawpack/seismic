!==================
subroutine setprob
!==================

    use fault_module, only: load_fault

    implicit none

    character*12 fname
    integer iunit

!
    iunit = 7
    fname = 'fault.data'

    call load_fault(fname)

    return
end
