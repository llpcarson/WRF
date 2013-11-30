RECURSIVE SUBROUTINE THREADED_INTEGRATION_DRIVER(theDomainId)
    USE MODULE_SYNC_DOMAINS     !Required to retrieve the domain grid pointer given the domain ID.
    USE MODULE_INTEGRATE_NESTS  !Simply a copy of MODULE_INTEGRATE, which conatins the recursive INTEGRATE subroutine.

    IMPLICIT NONE

    !Variables
    INTEGER,INTENT(IN)      ::  theDomainId
    TYPE(DOMAIN),POINTER    ::  grid_ptr
    INCLUDE "mpif.h"

    !Get the grid pointer for the corresponding domain ID
    CALL RETRIEVE_DOMAIN_POINTER(theDomainId, grid_ptr) !References to all domains are stored into MODULE_SYNC_DOMAINS when they are spawned in the INTEGRATE routine.
    IF(theDomainId .NE. grid_ptr%id) THEN
        CALL WRF_ERROR_FATAL('In frame/THREADED_INTEGRATION_DRIVER.F, could not retrieve domain grid pointer.')
    ENDIF

    !Invoke the recursive INTEGRATE subroutine
    !IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid_ptr%id,' START P-THREAD INTEGRATION'
    !WRITE(*,*) '@@@@@ DOMAIN ',grid_ptr%id,' START P-THREAD INTEGRATION ON ',mpi_wtime()
    CALL INTEGRATE(grid_ptr)
    !IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid_ptr%id,' END P-THREAD INTEGRATION'
    !WRITE(*,*) '@@@@@ DOMAIN ',grid_ptr%id,' END P-THREAD INTEGRATION ON ',mpi_wtime()

END SUBROUTINE
