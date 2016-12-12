!   ==================
    subroutine setprob
!   ==================

    use fault_module, only: load_fault

    implicit none

    call load_fault('fault.data')

    return
    end
