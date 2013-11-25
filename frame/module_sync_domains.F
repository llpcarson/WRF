! This module contains routines for allowing domains to sequentially synchronize on model tasks such as integration, forcing, feedback, and halo updates.
! It automatically figures out if the model is configured using a single domain, double-nested (i.e. multiple single nesting level domains), or triple nested
! (i.e. multiple 2-level nested domains). This module assumes that the outermost domain has ID "1", but makes no assumption about the IDs of other domains in
! any levels, as long as the domain IDs include all values from 1 to "grid%max_dom". While trasversing through all domains within a nesting level, this module
! trasverses through the domains in monotonic order, from the lowest domain ID to the highest domain ID within a nesting level. It labels the domains having
! the lowest/highest IDs within a nesting level as the first/last domains in that level for trasversal purpose. Thus, domain IDs can be distributed in any
! fashion among the nesting levels. This module will ensure that all domains within a nesting level will be visited at least once during the traversal process.
!
! Author: Thiago.Quirino@noaa.gov, 2013-03-12, DOC/NOAA/AOML/HRD, (305)361-4337 / (305)409-9587.
!
MODULE MODULE_SYNC_DOMAINS
    USE MODULE_DRIVER_CONSTANTS                                                                     !Needed to access the value of "max_domains"
    USE MODULE_DOMAIN                                                                               !Needed to access TYPE(DOMAIN)
    IMPLICIT NONE

    !An array of this TYPE will hold pointers to all model domain grids.
    TYPE DOMAIN_POINTER
        TYPE(DOMAIN),POINTER :: grid_ptr
    END TYPE DOMAIN_POINTER

    TYPE(DOMAIN_POINTER),SAVE,DIMENSION(1:max_domains),PRIVATE  ::  modelDomainsPointerList         !Holds a pointer to every model domain grid object.
    INTEGER,SAVE,PRIVATE                                        ::  currDomainWithPrivilege         !ID of current domain/grid with synchronized task privilege. This variable must be updated atomically.
    INTEGER,SAVE,PRIVATE                                        ::  currNestingLevelWithPrivilege   !Nesting level (0, 1, or 2) of the current domain/grid with synchronized task privilege. This variable must be updated atomically.
    INTEGER,SAVE,PRIVATE                                        ::  totalDomainsInSimulation        !Holds the total number of domains in the simulation (i.e., the value of grid%max_dom).
    INTEGER,SAVE,PRIVATE                                        ::  totalStormsInSimulation         !Holds the total number of storms in the simulation that are being tracked by high-resolution moving nests.
    INTEGER,SAVE,PRIVATE                                        ::  maxSimulationNestingLevel       !The max nesting level (0, 1, or 2) of the current model configuration.
    INTEGER,SAVE,DIMENSION(1:max_domains),PRIVATE               ::  domainNestingLevel              !Holds the nesting level (0, 1, or 2) of each domain (IDs ranging from 1 to max_domains).
    INTEGER,SAVE,DIMENSION(0:2,1:max_domains),PRIVATE           ::  domainsInEachNestingLevel       !Holds a list of the IDs of all domains in each nesting level. The first domain in each nesting level have sequence index "1".
    INTEGER,SAVE,DIMENSION(1:max_domains),PRIVATE               ::  domainIndexInNestingLevel       !Holds the sequence index (1 to max_domains) of each domain within their nesting level.
    INTEGER,SAVE,DIMENSION(0:2),PRIVATE                         ::  numDomainsInEachNestingLevel    !Holds the number of domains in each nesting level, that is, the number of slots in the "domainsInEachNestingLevel" that actually hold valid domain IDs.
    INTEGER,SAVE,DIMENSION(0:2),PRIVATE                         ::  privilegedDomainInNestingLevel  !For each nesting level, holds the index into the "domainsInEachNestingLevel" array of the current domain with synchronized task privilege.
    INTEGER,SAVE,DIMENSION(0:2),PRIVATE                         ::  idOfFirstDomainInLevel          !Holds the ID of the first domain in each nesting level. Needed for circular synchronized task privilege assignment.
    INTEGER,SAVE,DIMENSION(0:2),PRIVATE                         ::  idOfLastDomainInLevel           !Holds the ID of the last domain in each nesting level. Needed for circular synchronized task privilege assignment.
    INTEGER,SAVE,DIMENSION(0:2,1:max_domains),PRIVATE           ::  nestingLevelDomainBarrier       !This array is initialize to 0. It is used by the NEST_LEVEL_BARRIER routine to create an execution barrier until all domains from a nesting level call the routine.
    INTEGER,SAVE,DIMENSION(1:max_domains)                       ::  domainCommunicatorList          !Holds a different MPI communicators for each model domain
    INTEGER,SAVE,DIMENSION(1:max_domains)                       ::  domainXCommunicatorList         !Holds a different MPI communicators along x-direction of mesh for each model domain
    INTEGER,SAVE,DIMENSION(1:max_domains)                       ::  domainYCommunicatorList         !Holds a different MPI communicators along y-direction of mesh for each model domain
    LOGICAL,PARAMETER                                           ::  isIntegrationVerbose = .FALSE.  !Whether to display debug information.

    CONTAINS

    !This subroutine is used to initialize this module.
    SUBROUTINE INIT_MODULE_SYNC_DOMAINS
        USE MODULE_DM,ONLY: local_communicator,local_communicator_x,local_communicator_y    !These MPI communicators will be duplicated for each domain

        IMPLICIT NONE
        INCLUDE 'mpif.h'                                                                    !Required to invoke MPI_COMM_DUP to create MPI compute node communicators for each model domain
        INTEGER                         ::  domainId, parentId, level, nestid,ierr
        INTEGER,DIMENSION(0:2)          ::  lastDomainIndexAssignedInLevel

        WRITE(*,*) '@@@@@ INITIALIZING MODULE_SYNC_DOMAINS.'

        !Get the number of domains in the simulation and store it in the global static variable
        CALL nl_get_max_dom(1,totalDomainsInSimulation)

        !Initialize all static global SAVE variables
        currDomainWithPrivilege = 1         !Domain D01 has the initial task privilege
        currNestingLevelWithPrivilege = 0   !Domain D01 is always assigned to nesting level 0        
        domainNestingLevel = 0              !This simulation must have at least 1 outer-most domain.
        maxSimulationNestingLevel = 0       !This simulation must have at least 1 outer-most domain.
        domainsInEachNestingLevel = 0       !This count will be updated as we iterate through every domain.
        domainIndexInNestingLevel = 1       !This variable will be updated as we iterate through every domain.
        numDomainsInEachNestingLevel = 0    !This variable will be updated as we iterate through every domain.
        privilegedDomainInNestingLevel = 1  !A suitable initial value since the outer-most domain is initially assigned the task privilege
        idOfFirstDomainInLevel = 1          !This variable will be updated as we iterate through every domain.
        idOfLastDomainInLevel = 1           !This variable will be updated as we iterate through every domain.
        nestingLevelDomainBarrier = 1       !Used by the NEST_LEVEL_BARRIER routine. Reset to zero by the first domain entering the routine.

        !Initialize the array used to hold the index of the last domain assigned to different nesting levels
        lastDomainIndexAssignedInLevel = 0  !No domain has been assigned to a spot in a level yet. Valid indexing starts at 1. Variable is incremented before being used.

        !Iterate through every domain and extract their properties
        DO domainId=1,totalDomainsInSimulation
            !Get the parent domain ID and determine the nesting level of this domain
            CALL nl_get_parent_id(domainId,parentId)    !For nests, their parent ID must be >=1
            IF(parentId < 1) THEN
                domainNestingLevel(domainId) = 0    !This is the D01 domain
            ELSE IF(parentId == 1) THEN  
                domainNestingLevel(domainId) = 1    !This is D01's direct child nest
            ELSE
                domainNestingLevel(domainId) = 2    !This nest is a child of one of D01's child nests
            ENDIF
                
            !Update the count of domains in the respective nesting level
            numDomainsInEachNestingLevel(domainNestingLevel(domainId)) = numDomainsInEachNestingLevel(domainNestingLevel(domainId)) + 1

            !Store the ID of this domain in the next available spot in the list of domains belonging to the corresponding nesting level
            lastDomainIndexAssignedInLevel(domainNestingLevel(domainId)) = lastDomainIndexAssignedInLevel(domainNestingLevel(domainId)) + 1 !Increment the domain index variable for this nesting level to point to the next available spot in the array
            domainsInEachNestingLevel(domainNestingLevel(domainId),lastDomainIndexAssignedInLevel(domainNestingLevel(domainId))) = domainId !Store the domain ID in the next available slot for this nesting level
            domainIndexInNestingLevel(domainId) = lastDomainIndexAssignedInLevel(domainNestingLevel(domainId))                              !Store this domain's index within its nesting level

            !Update the ID of the first domain in the corresponding nesting level if it has not been updated before
            IF(lastDomainIndexAssignedInLevel(domainNestingLevel(domainId))==1) idOfFirstDomainInLevel(domainNestingLevel(domainId)) = domainId

            !Update the ID of the last domain in the correspoding nesting level
            idOfLastDomainInLevel(domainNestingLevel(domainId)) = domainId
        ENDDO

        !Determine the max nesting level in this simulation
        maxSimulationNestingLevel = MAXVAL(domainNestingLevel)

        !Determine the total number of storms in the simulation that are being tracked by high-resolution nests
        totalStormsInSimulation = numDomainsInEachNestingLevel(maxSimulationNestingLevel)

        !Print information about every domain found
        WRITE(*,*) '@@@@@ Total domains: ',totalDomainsInSimulation
        WRITE(*,*) '@@@@@ Total storms:  ',totalStormsInSimulation
        WRITE(*,*) '@@@@@ Max level:     ',maxSimulationNestingLevel
        DO domainId=1,totalDomainsInSimulation
            WRITE(*,*) '@@@@@ Domain ',domainId,': ',domainNestingLevel(domainId),'/',domainIndexInNestingLevel(domainId)
        ENDDO

        !Print information about every valid nesting level
        DO level=0,maxSimulationNestingLevel
            WRITE(*,*) '@@@@@ Level ',level,': ',numDomainsInEachNestingLevel(level),'/',idOfFirstDomainInLevel(level),'/',idOfLastDomainInLevel(level)
            WRITE(*,*) '@@@@@ ',domainsInEachNestingLevel(level,1:numDomainsInEachNestingLevel(level))
        ENDDO

        !Let's instantiate MPI communicators for each model domain
        !Initialize communicator arrays
        domainCommunicatorList=-999 ; domainXCommunicatorList=-999 ; domainYCommunicatorList=-999

        !Duplicate the compute node communicators
        WRITE(*,*) '@@@@@ MPI ORIG: ',local_communicator,'/',local_communicator_x,'/',local_communicator_y
        DO nestid=1,totalDomainsInSimulation
            !Duplicate the compute nodes communicator
            CALL mpi_comm_dup(local_communicator,domainCommunicatorList(nestid),ierr)
            IF(ierr .NE. 0) CALL WRF_ERROR_FATAL("Failed to duplicate compute node MPI communicator for each domain in RECURSIVE SUBROUTINE INTEGRATE.")

            !Duplicate the compute nodes communicator in the x-direction of the mesh
            CALL mpi_comm_dup(local_communicator_x,domainXCommunicatorList(nestid),ierr)
            IF(ierr .NE. 0) CALL WRF_ERROR_FATAL("Failed to duplicate compute node MPI communicator for each domain in RECURSIVE SUBROUTINE INTEGRATE.")

            !Duplicate the compute nodes communicator in the y-direction of the mesh
            CALL mpi_comm_dup(local_communicator_y,domainYCommunicatorList(nestid),ierr)
            IF(ierr .NE. 0) CALL WRF_ERROR_FATAL("Failed to duplicate compute node MPI communicator for each domain in RECURSIVE SUBROUTINE INTEGRATE.")

            !Added by Thiago to test multiple threads using serialized access to a single communicator
            !domainCommunicatorList(nestid)=local_communicator
            !domainXCommunicatorList(nestid)=local_communicator_x
            !domainYCommunicatorList(nestid)=local_communicator_y

            WRITE(*,*) '@@@@@ MPIDUPP ',nestid,': ',domainCommunicatorList(nestid),'/',domainXCommunicatorList(nestid),'/',domainYCommunicatorList(nestid)
            CALL MPI_BARRIER(domainCommunicatorList(nestid),ierr)
            IF(ierr .NE. 0) CALL WRF_ERROR_FATAL("Failed to duplicate compute node MPI communicator for each domain in RECURSIVE SUBROUTINE INTEGRATE.")
            CALL MPI_BARRIER(domainXCommunicatorList(nestid),ierr)
            IF(ierr .NE. 0) CALL WRF_ERROR_FATAL("Failed to duplicate compute node MPI communicator for each domain in RECURSIVE SUBROUTINE INTEGRATE.")
            CALL MPI_BARRIER(domainYCommunicatorList(nestid),ierr)
            IF(ierr .NE. 0) CALL WRF_ERROR_FATAL("Failed to duplicate compute node MPI communicator for each domain in RECURSIVE SUBROUTINE INTEGRATE.")
        ENDDO

        !Initialize the "C" parallel integration p-threads driver. The current implementation does not use a p-thread pool. Instead, it spawn new p-threads every time the nests are integrated.
        CALL init_pthread_sync_interface(totalDomainsInSimulation,totalStormsInSimulation,ierr)
        IF(ierr .NE. 0) THEN
            CALL WRF_ERROR_FATAL('Could not initialize the P-thread interface driver.')
        ENDIF

        WRITE(*,*) '@@@@@ END MODULE_SYNC_DOMAINS INITIALIZATION.'
    END SUBROUTINE


    !This subroutine returns a pointer to the TYPE(DOMAIN) grid of the specified domain
    SUBROUTINE RETRIEVE_DOMAIN_POINTER(theDomainId,grid_ptr)
        USE MODULE_DOMAIN
        IMPLICIT NONE
        INTEGER,INTENT(IN)               :: theDomainId
        TYPE(DOMAIN),POINTER,INTENT(INOUT) :: grid_ptr

        IF(.NOT. ASSOCIATED(modelDomainsPointerList(theDomainId)%grid_ptr)) THEN
            CALL WRF_ERROR_FATAL('Attempt to retrieve an invalid model domain pointer.')
        ENDIF
        grid_ptr => modelDomainsPointerList(theDomainId)%grid_ptr
    END SUBROUTINE


    !This subroutine sets the pointer to the TYPE(DOMAIN) grid of the specified domain
    !TODO: Call WRF_FATAL if the GRID_PTR argument is not associated
    SUBROUTINE STORE_DOMAIN_POINTER(grid_ptr)
        USE MODULE_DOMAIN
        IMPLICIT NONE
        TYPE(DOMAIN),POINTER,INTENT(IN) :: grid_ptr

        IF(.NOT. ASSOCIATED(grid_ptr)) THEN
            CALL WRF_ERROR_FATAL('Attempt to store a non-associated model domain pointer.')
        ENDIF

        IF(.NOT. ASSOCIATED(modelDomainsPointerList(grid_ptr%id)%grid_ptr)) THEN
            modelDomainsPointerList(grid_ptr%id)%grid_ptr => grid_ptr
            IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',grid_ptr%id,' STORED ITS GRID POINTER'
        ELSE IF(isIntegrationVerbose) THEN
            WRITE(*,*) '@@@@@ DOMAIN ',grid_ptr%id,' GRID POINTER IGNORED'
        ENDIF
    END SUBROUTINE


    !This function returns the unique MPI communicator for the specified domain
    INTEGER FUNCTION GET_DOMAIN_MPI_COMMUNICATOR(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId     !The ID of the domain
        GET_DOMAIN_MPI_COMMUNICATOR = domainCommunicatorList(theDomainId)
    END FUNCTION


    !This function returns the unique MPI communicator along x-direction of mesh for the specified domain
    INTEGER FUNCTION GET_DOMAIN_MPI_COMMUNICATOR_X(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId     !The ID of the domain
        GET_DOMAIN_MPI_COMMUNICATOR_X = domainXCommunicatorList(theDomainId)
    END FUNCTION


    !This function returns the unique MPI communicator along x-direction of mesh for the specified domain
    INTEGER FUNCTION GET_DOMAIN_MPI_COMMUNICATOR_Y(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId     !The ID of the domain
        GET_DOMAIN_MPI_COMMUNICATOR_Y = domainYCommunicatorList(theDomainId)
    END FUNCTION


    !This function returns the ID (from "1" to "totalDomainsInSimulation") of the domain currently having the privilege to perform synchronized tasks
    INTEGER FUNCTION GET_DOMAIN_WITH_PRIVILEGE
        IMPLICIT NONE
        !!$OMP CRITICAL (module_sync_domains_give_privilege)
        GET_DOMAIN_WITH_PRIVILEGE = currDomainWithPrivilege
        !!$OMP END CRITICAL (module_sync_domains_give_privilege)
    END FUNCTION


    !This function returns the nesting level (0, 1, or 2) whose domains currently have the privilege to perform synchronized tasks
    INTEGER FUNCTION GET_NEST_LEVEL_WITH_PRIVILEGE
        IMPLICIT NONE
        !!$OMP CRITICAL (module_sync_domains_give_privilege)
        GET_NEST_LEVEL_WITH_PRIVILEGE = currNestingLevelWithPrivilege
        !!$OMP END CRITICAL (module_sync_domains_give_privilege)
    END FUNCTION


    !This function returns the simulation nesting level (0 for single-domain, 1 for double-nested, or 2 for triple-nested) in the simulation
    INTEGER FUNCTION GET_MAX_SIMULATION_NEST_LEVEL
        IMPLICIT NONE
        GET_MAX_SIMULATION_NEST_LEVEL = maxSimulationNestingLevel
    END FUNCTION


    !This function returns the number of storms in the simulation that are being tracked by high-resolution nests
    INTEGER FUNCTION GET_TOTAL_STORMS_IN_SIMULATION
        IMPLICIT NONE
        GET_TOTAL_STORMS_IN_SIMULATION = totalStormsInSimulation
    END FUNCTION


    !This function returns whether debug information should be displayed for the parallel integration
    LOGICAL FUNCTION IS_INTEGRATION_VERBOSE
        IMPLICIT NONE
        IS_INTEGRATION_VERBOSE = isIntegrationVerbose
    END FUNCTION


    !This function returns the nesting level ("0" for D01, "1" for intermediary domains, or "2" for inner-most domains) of a specified domain
    !TODO: Check if the specified domain ID is withing valid range and call WRF_FATAL if not
    INTEGER FUNCTION GET_DOMAIN_NESTING_LEVEL(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId     !The ID of the domain
        GET_DOMAIN_NESTING_LEVEL = domainNestingLevel(theDomainId)
    END FUNCTION


    !This function returns the ID of the first domain in a specified nesting level
    !TODO: Check if the specified nesting level is within valid range and call WRF_FATAL if not
    INTEGER FUNCTION GET_FIRST_DOMAIN_IN_A_LEVEL(theNestingLevel)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theNestingLevel  !The nesting level of interest

        !Determine the ID of the first domain in the specified level
        GET_FIRST_DOMAIN_IN_A_LEVEL = idOfFirstDomainInLevel(theNestingLevel)
    END FUNCTION


    !This function returns the ID of the first domain in the same nesting level of a specified domain
    !TODO: Check if the specified domain ID is within valid range and call WRF_FATAL if not
    INTEGER FUNCTION GET_FIRST_DOMAIN_IN_MY_LEVEL(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId  !The ID of the domain
        GET_FIRST_DOMAIN_IN_MY_LEVEL = idOfFirstDomainInLevel(domainNestingLevel(theDomainId))
    END FUNCTION


    !This function returns .TRUE. if the specified domain is the first domain in the specified domain's nesting level
    !TODO: Check if the specified domain ID is within valid range and call WRF_FATAL if not
    LOGICAL FUNCTION AM_I_FIRST_DOMAIN_IN_MY_LEVEL(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId  !The ID of the domain
        AM_I_FIRST_DOMAIN_IN_MY_LEVEL = (theDomainId == idOfFirstDomainInLevel(domainNestingLevel(theDomainId)))
    END FUNCTION


    !This function returns the ID of the last domain in a specified nesting level
    !TODO: Check if the specified nesting level is within valid range and call WRF_FATAL if not
    INTEGER FUNCTION GET_LAST_DOMAIN_IN_A_LEVEL(theNestingLevel)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theNestingLevel  !The nesting level of interest
        
        !Determine the ID of the last domain in the specified nesting level
        GET_LAST_DOMAIN_IN_A_LEVEL = idOfLastDomainInLevel(theNestingLevel)
    END FUNCTION


    !This function returns the ID of the last domain in the same nesting level  of a specified domain
    !TODO: Check if the specified domain ID is withing valid range and call WRF_FATAL if not
    INTEGER FUNCTION GET_LAST_DOMAIN_IN_MY_LEVEL(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId  !The ID of the domain
        GET_LAST_DOMAIN_IN_MY_LEVEL = idOfLastDomainInLevel(domainNestingLevel(theDomainId))
    END FUNCTION


    !This function returns .TRUE. if the specified domain is the last domain in the specified domain's nesting level
    !TODO: Check if the specified domain ID is within valid range and call WRF_FATAL if not
    LOGICAL FUNCTION AM_I_LAST_DOMAIN_IN_MY_LEVEL(theDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId  !The ID of the domain
        AM_I_LAST_DOMAIN_IN_MY_LEVEL = (theDomainId == idOfLastDomainInLevel(domainNestingLevel(theDomainId)))
    END FUNCTION


    !This subroutine blocks and waits until the specified domain receives the privilege to perform synchronized tasks
    !TODO: Check if the specified domain ID is withing valid range and call WRF_FATAL if not
    SUBROUTINE WAIT_FOR_DOMAIN_PRIVILEGE(theDomainId, optionalSleepTimeMillisecs)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: theDomainId                  !The ID of the domain
        INTEGER,INTENT(IN),OPTIONAL :: optionalSleepTimeMillisecs   !The optional sleep time, in milliseconds
        INTEGER                     :: theActualSleepTime           !Will hold the actual sleep time for each iteration of the wait loop
        INTEGER                     :: ierr                         !Used to check errors from using a p-thread MUTEX

        !Determine the sleep time for the wait loop
        theActualSleepTime = 500 !Default is check iterations every 10ms
        IF(PRESENT(optionalSleepTimeMillisecs)) THEN
            IF(optionalSleepTimeMillisecs>0) theActualSleepTime = optionalSleepTimeMillisecs
        ENDIF

        !Wait until the specified domain has task privilege
        CALL block_till_domain_privilege(theDomainId,currDomainWithPrivilege,ierr)
        IF(ierr .NE. 0) THEN
            CALL WRF_ERROR_FATAL('An error occurred while block-waiting for domain privilege.')
        ENDIF
        !DO WHILE(currDomainWithPrivilege .NE. theDomainId)
        !    CALL SLEEPQQ(theActualSleepTime)
        !ENDDO
        !!DO WHILE(.TRUE.)
        !!    IF(currDomainWithPrivilege == theDomainId) EXIT
        !!ENDDO
    END SUBROUTINE


    !This subroutine blocks and waits until the specified nesting level receives the privilege for its domains to perform synchronized tasks
    !TODO: Check if the specified nesting level is within valid range and call WRF_FATAL if not
    SUBROUTINE WAIT_FOR_NEST_LEVEL_PRIVILEGE(theNestingLevel, optionalSleepTimeMillisecs)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: theNestingLevel              !The nesting level of interest
        INTEGER,INTENT(IN),OPTIONAL :: optionalSleepTimeMillisecs   !The optional sleep time, in milliseconds
        INTEGER                     :: theActualSleepTime           !Will hold the actual sleep time for each iteration of the wait loop

        !Determine the sleep time for the wait loop
        theActualSleepTime = 1 !Default is check iterations every 10ms
        IF(PRESENT(optionalSleepTimeMillisecs)) THEN
            IF(optionalSleepTimeMillisecs>0) theActualSleepTime = optionalSleepTimeMillisecs
        ENDIF

        !Wait until the specified nesting level has task privilege
        DO WHILE(currNestingLevelWithPrivilege .NE. theNestingLevel)
            CALL SLEEPQQ(theActualSleepTime)
        ENDDO
        !DO WHILE(.TRUE.)
        !    IF(currNestingLevelWithPrivilege == theNestingLevel) EXIT
        !ENDDO
    END SUBROUTINE


    !This subroutine blocks and waits until all of the domains within a nesting level invoke this routine. It can only be invoked
    !by domains belonging to the currently privileged nesting level. Otherwise, WRF_ERROR_FATAL is called.
    SUBROUTINE NEST_LEVEL_BARRIER(myDomainId, optionalSleepTimeMillisecs)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: myDomainId                   !The ID of a domain in the specified nesting level
        INTEGER,INTENT(IN),OPTIONAL :: optionalSleepTimeMillisecs   !The optional sleep time, in milliseconds
        INTEGER                     :: theActualSleepTime           !Will hold the actual sleep time for each iteration of the wait loop
		INTEGER						:: theSum						!Holds the sum of the number of domains inside the barrier

		!Check if the specified domain is in the currently privileged nesting level
        IF(domainNestingLevel(myDomainId) .NE. currNestingLevelWithPrivilege ) THEN
            CALL WRF_ERROR_FATAL("Invalid nesting level barrier request by an unprivileged domain.")
        ENDIF

        !Determine the sleep time for the wait loop
        theActualSleepTime = 5  !Default is check iterations every 10ms
        IF(PRESENT(optionalSleepTimeMillisecs)) THEN
            IF(optionalSleepTimeMillisecs>0) theActualSleepTime = optionalSleepTimeMillisecs
        ENDIF

        !Set the correspoding element in the "nestingLevelDomainBarrier" to 1. The first domain to enter here resets all flags to 0.
		!!$OMP CRITICAL (module_sync_domains_level_barrier)
		theSum=SUM(nestingLevelDomainBarrier(domainNestingLevel(myDomainId),1:numDomainsInEachNestingLevel(domainNestingLevel(myDomainId))))
		IF(theSum == numDomainsInEachNestingLevel(domainNestingLevel(myDomainId))) THEN
			!Reset the flags
			nestingLevelDomainBarrier(domainNestingLevel(myDomainId),:) = 0
            IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' RESET ALL NEST LEVEL BARRIER FLAGS'
        ELSE IF(nestingLevelDomainBarrier(domainNestingLevel(myDomainId),domainIndexInNestingLevel(myDomainId)) == 1) THEN	!There is a duplicate call by the same domain
			CALL WRF_ERROR_FATAL("Invalid nesting level barrier request by duplicate domain.")
        ENDIF
		nestingLevelDomainBarrier(domainNestingLevel(myDomainId),domainIndexInNestingLevel(myDomainId)) = 1
        IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' SET ITS NEST LEVEL BARRIER FLAG'
		!!$OMP END CRITICAL (module_sync_domains_level_barrier)

        !Wait until all domains are locked in this subroutine
        DO WHILE(SUM(nestingLevelDomainBarrier(domainNestingLevel(myDomainId),1:numDomainsInEachNestingLevel(domainNestingLevel(myDomainId)))) < numDomainsInEachNestingLevel(domainNestingLevel(myDomainId)))
            nestingLevelDomainBarrier(domainNestingLevel(myDomainId),domainIndexInNestingLevel(myDomainId)) = 1 !Keep setting the flag to 1 in case another thread resets it to 0 by mistake. Remove this once made thread-safe.
            CALL SLEEPQQ(theActualSleepTime)
        ENDDO
        !DO WHILE(.TRUE.)
        !    IF(SUM(nestingLevelDomainBarrier(domainNestingLevel(myDomainId),1:numDomainsInEachNestingLevel(domainNestingLevel(myDomainId)))) == numDomainsInEachNestingLevel(domainNestingLevel(myDomainId))) EXIT
        !ENDDO
        IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' LEAVING NEST LEVEL BARRIER'
    END SUBROUTINE


    !This subroutine blocks and waits until all of the domains within a nesting level invoke this routine. It can only be invoked
    !by domains belonging to the currently privileged nesting level. Otherwise, WRF_ERROR_FATAL is called.
    SUBROUTINE NEST_LEVEL_BARRIER_OMP(myDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: myDomainId                   !The ID of a domain in the specified nesting level

        !Check if the specified domain is in the currently privileged nesting level
        IF(domainNestingLevel(myDomainId) .NE. currNestingLevelWithPrivilege ) THEN
            CALL WRF_ERROR_FATAL("Invalid nesting level barrier request by an unprivileged domain.")
        ENDIF

        !Invoke an OpenMP BARRIER command so that all threads will synchronize at this point
        IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' ENTERED OMP NEST LEVEL BARRIER'
        !!$OMP BARRIER
        IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' EXITED OMP NEST LEVEL BARRIER'
    END SUBROUTINE


    !This subroutine blocks and waits until all of the domains within a nesting level invoke this routine. It can only be invoked
    !by domains belonging to the currently privileged nesting level. Otherwise, WRF_ERROR_FATAL is called.
    SUBROUTINE NEST_LEVEL_BARRIER_PTHREAD(myDomainId)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: myDomainId                   !The ID of a domain in the specified nesting level
        INTEGER                     :: ierr                         !Used to determine p-thread barrier errors

        !Check if the specified domain is in the currently privileged nesting level
        IF(domainNestingLevel(myDomainId) .NE. currNestingLevelWithPrivilege ) THEN
            CALL WRF_ERROR_FATAL("Invalid nesting level barrier request by an unprivileged domain.")
        ENDIF

        !Domain D01 does not need to wait at the barrier
        IF(myDomainId==1) RETURN

        !Invoke an P-thread BARRIER command so that all threads will synchronize at this point
        IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' ENTERED PTHREAD NEST LEVEL BARRIER'
        CALL integration_pthread_barrier(ierr)
        IF(ierr .NE. 0) THEN
            CALL WRF_ERROR_FATAL('An error occurred in p-thread barrier in NEST_LEVEL_BARRIER_PTHREAD.')
        ENDIF
        IF(isIntegrationVerbose) WRITE(*,*) '@@@@@ DOMAIN ',myDomainId,' EXITED PTHREAD NEST LEVEL BARRIER'
    END SUBROUTINE


    !This subroutine allows a specific domain to receive task privilege, which also switches which nesting level has privilege.
    !TODO: Check if the specified domain ID is withing valid range and call WRF_FATAL if not
    SUBROUTINE ASSIGN_PRIVILEGE_TO_DOMAIN(theDomainId)
		USE MODULE_WRF_ERROR
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theDomainId          !The ID of the domain who should receive task privilege
        INTEGER             :: ierr                 !Used to check errors from using p-thread mutex

        !Check if the domain ID actually exists and call WRF_FATAL if not.
        IF(theDomainId > totalDomainsInSimulation) THEN
             CALL WRF_ERROR_FATAL("Invalid request to assign synchronized task privilege to an invalid nesting level.")
        ENDIF

        !!$OMP CRITICAL (module_sync_domains_give_privilege)
        currNestingLevelWithPrivilege = domainNestingLevel(theDomainId)                                         !Assign privilege to the specified domain's respective nesting level
        privilegedDomainInNestingLevel(currNestingLevelWithPrivilege) = domainIndexInNestingLevel(theDomainId)  !Update the privileged domain index pointer for the spceified nesting level
        currDomainWithPrivilege = theDomainId                                                                   !Assign privilege to the specified domain
        !!$OMP END CRITICAL (module_sync_domains_give_privilege)

        !Forces any waiting thread to check if their domain was assigned privilege
        CALL force_check_domain_privilege(theDomainId,ierr)
        IF(ierr .NE. 0) THEN
            CALL WRF_ERROR_FATAL('An error occurred while forcing a domain privilege check in ASSIGN_PRIVILEGE_TO_DOMAIN.')
        ENDIF
    END SUBROUTINE


    !This subroutine allows a specific nesting level to receive task privilege, starting with the first domain in that nesting level.
    !This routine calls WRF_FATAL if the nesting level receiving task privilege is greater than "maxSimulationNestingLevel".
    SUBROUTINE ASSIGN_PRIVILEGE_TO_NEST_LEVEL(theNestingLevel)
		USE MODULE_WRF_ERROR
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: theNestingLevel  !The nesting level whose first domain should receive task privilege
        INTEGER             :: ierr             !Used to check errors from using p-thread mutex

        !Check if the nesting level actually exists and call WRF_FATAL if not.
        IF(theNestingLevel > maxSimulationNestingLevel) THEN
             CALL WRF_ERROR_FATAL("Invalid request to assign synchronized task privilege to an invalid nesting level.")
        ENDIF

        !!$OMP CRITICAL (module_sync_domains_give_privilege)
        currNestingLevelWithPrivilege = theNestingLevel                     !Assign privilege to the nesting level
        privilegedDomainInNestingLevel(theNestingLevel) = 1                 !Update the privileged domain index pointer for the specified nesting level
        currDomainWithPrivilege = idOfFirstDomainInLevel(theNestingLevel)   !Assign privilege to the first domain in the specified nesting level
        !!$OMP END CRITICAL (module_sync_domains_give_privilege)

        !Forces any waiting thread to check if their domain was assigned privilege
        CALL force_check_domain_privilege(idOfFirstDomainInLevel(theNestingLevel),ierr)
        IF(ierr .NE. 0) THEN
            CALL WRF_ERROR_FATAL('An error occurred while forcing a domain privilege check in ASSIGN_PRIVILEGE_TO_NEST_LEVEL.')
        ENDIF
    END SUBROUTINE


    !This subroutine gives task privilege to the next domain in the same nesting level in a circular fashion.
    !If the optional argument "isMoveToNextNestingLevel" is true and "myDomainId" corresponds to the ID of the
    !last domain in the currently privileged nesting level, then privilege will be trasnferred to the first
    !domain in the next nesting level. This routine calls WRF_FATAL if the nesting level receiving task
    !privilege is greater than "maxSimulationNestingLevel". This routine also call WRF_FATAL if the ID of the
    !domain invoking this routine is the ID of the domain currently having the task privilege
    SUBROUTINE CYCLE_DOMAIN_PRIVILEGE(myDomainId,isMoveToNextNestingLevel)
        USE MODULE_WRF_ERROR
		IMPLICIT NONE
        INTEGER,INTENT(IN)          :: myDomainId               !The ID of the requesting domain
        LOGICAL,INTENT(IN),OPTIONAL :: isMoveToNextNestingLevel !Whether to transfer privilege from last domain in a level to first domain in the same (false) or next (true) level
        INTEGER                     :: nestingLevelIncrement    !Holds the nesting level increment that determines which nesting level will receive privilege
        INTEGER                     :: nextLevel,nextDomain     !Holds the next nesting level or next domain to be assigned privilege
        INTEGER                     :: ierr                     !Used to check errors from using p-thread mutex

        !!$OMP CRITICAL (module_sync_domains_give_privilege)
        !Check if the specified domain currently has the task privilege
        IF(myDomainId .NE. currDomainWithPrivilege .OR. domainNestingLevel(myDomainId) .NE. currNestingLevelWithPrivilege ) THEN
            CALL WRF_ERROR_FATAL("Invalid request to assign synchronized task privilege by an unprivileged domain.")
        ENDIF

        !The last domain in the current nesting level will either transfer privilege (1) the first domain within its nesting level
        !(if "isMoveToNextNestingLevel" is missing or set to .FALSE.) or (2) transfer privilege to first domain in the next nesting
        !level (if "isMoveToNextNestingLevel" is available and set to .TRUE.). Otherwise, the subsequent domain within the same
        !nesting level will receive the synchronized task privilege.
        IF(myDomainId == idOfLastDomainInLevel(domainNestingLevel(myDomainId))) THEN
            !Determine whether or not to move to the next nesting level
            nextLevel = currNestingLevelWithPrivilege   !Assume we will remain in the current nesting level
            IF(PRESENT(isMoveToNextNestingLevel)) THEN
                IF(isMoveToNextNestingLevel) THEN
                    !Check if the nesting level actually exists and call WRF_FATAL if not.
                    nextLevel = currNestingLevelWithPrivilege + 1   !Increment the nesting level
                    IF(nextLevel > maxSimulationNestingLevel) THEN
                        CALL WRF_ERROR_FATAL("Invalid request to assign synchronized task privilege to an invalid nesting level.")
                    ENDIF
                ENDIF
            ENDIF

            !Assign privilege to first domain in the chosen (i.e., either the current or subsequent) nesting level
            currNestingLevelWithPrivilege = nextLevel                     !Assign privilege to the nesting level
            privilegedDomainInNestingLevel(nextLevel) = 1                 !Update the privileged domain index pointer for the specified nesting level
            currDomainWithPrivilege = idOfFirstDomainInLevel(nextLevel)   !Assign privilege to the first domain in the specified nesting level
        ELSE
            !Determine the next domain within the currently privileged nesting level that will receive privilege
            privilegedDomainInNestingLevel(currNestingLevelWithPrivilege) = privilegedDomainInNestingLevel(currNestingLevelWithPrivilege) + 1                   !Increment the index pointing to the privilege domain within the current nesting level
            currDomainWithPrivilege = domainsInEachNestingLevel(currNestingLevelWithPrivilege,privilegedDomainInNestingLevel(currNestingLevelWithPrivilege))    !Get the ID of the domain that has been assigned privilege
        ENDIF
        !!$OMP END CRITICAL (module_sync_domains_give_privilege)

        !Forces any waiting thread to check if their domain was assigned privilege
        CALL force_check_domain_privilege(currDomainWithPrivilege,ierr)
        IF(ierr .NE. 0) THEN
            CALL WRF_ERROR_FATAL('An error occurred while forcing a domain privilege check in CYCLE_DOMAIN_PRIVILEGE.')
        ENDIF
    END SUBROUTINE

END MODULE MODULE_SYNC_DOMAINS
