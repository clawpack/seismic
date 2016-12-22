!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault

    implicit none

    character*12 :: fname

    fname = 'fault.data'

    call load_fault(fname)

    return
    end
