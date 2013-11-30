!WRF:DRIVER_LAYER:INTEGRATION
!

MODULE module_integrate

CONTAINS

RECURSIVE SUBROUTINE integrate ( grid )



   USE module_domain
   USE module_driver_constants
   USE module_nesting
   USE module_configure
   USE module_timing
   USE module_utility
   USE module_sync_domains	!Added by Thiago to allow domain-synchronized execution of various parts of the model.

   IMPLICIT NONE

   !  Input data.

   TYPE(domain) , POINTER :: grid

! module_integrate:integrate
! <DESCRIPTION> 
! This is a driver-level routine that controls the integration of a
! domain and subdomains rooted at the domain. 
! 
! The integrate routine takes a domain pointed to by the argument
! <em>grid</em> and advances the domain and its associated nests from the
! grid's current time, stored within grid%domain_clock, to a given time
! forward in the simulation, stored as grid%stop_subtime. The
! stop_subtime value is arbitrary and does not have to be the same as
! time that the domain finished integrating.  The simulation stop time
! for the grid is known to the grid's clock (grid%domain_clock) and that
! is checked with a call to domain_clockisstoptime prior to beginning the
! loop over time period that is specified by the
! current time/stop_subtime interval.
! 
! The clock, the simulation stop time for the domain, and other timing
! aspects for the grid are set up in the routine
! (<a href="setup_timekeeping.html">setup_timekeeping</a>) at the time
! that the domain is initialized.
! The lower-level time library and the type declarations for the times
! and time intervals used are defined either in 
! external/esmf_time_f90/module_utility.F90 or in 
! external/io_esmf/module_utility.F90 depending on a build-time decision to 
! incorporate either the embedded ESMF subset implementation contained in 
! external/esmf_time_f90 or to use a site-specific installation of the ESMF 
! library.  This decision is made during the configuration step of the WRF 
! build process.  Note that arithmetic and comparison is performed on these 
! data types using F90 operator overloading, also defined in that library.
! 
! This routine is the lowest level of the WRF Driver Layer and for the most
! part the WRF routines that are called from here are in the topmost level
! of the Mediation Layer.  Mediation layer routines typically are not 
! defined in modules. Therefore, the routines that this routine calls
! have explicit interfaces specified in an interface block in this routine.
!
! As part of the Driver Layer, this routine is intended to be non model-specific
! and so a minimum of WRF-specific logic is coded at this level. Rather, there
! are a number of calls to mediation layer routines that contain this logic, some
! of which are merely stubs in WRF Mediation Layer that sits below this routine
! in the call tree.  The routines that integrate calls in WRF are defined in
! share/mediation_integrate.F.
! 
! Flow of control
! 
! 1. Check to see that the domain is not finished 
! by testing the value returned by domain_clockisstoptime for the
! domain.
! 
! 2. <a href=model_to_grid_config_rec.html>Model_to_grid_config_rec</a> is called to load the local config_flags
! structure with the configuration information for the grid stored
! in model_config_rec and indexed by the grid's unique integer id. These
! structures are defined in frame/module_configure.F.
! 
! 3. The current time of the domain is retrieved from the domain's clock
! using domain_get_current_time.  
! 
! 4. Iterate forward while the current time is less than the stop subtime.
! 
! 4.a. Start timing for this iteration (only on node zero in distributed-memory runs)
! 
! 4.b. Call <a href=med_setup_step.html>med_setup_step</a> to allow the mediation layer to 
! do anything that's needed to call the solver for this domain.  In WRF this means setting
! the indices into the 4D tracer arrays for the domain.
! 
! 4.c. Check for any nests that need to be started up at this time.  This is done 
! calling the logical function <a href=nests_to_open.html>nests_to_open</a> (defined in 
! frame/module_nesting.F) which returns true and the index into the current domain's list
! of children to use for the nest when one needs to be started.
! 
! 4.c.1  Call <a href=alloc_and_configure_domain.html>alloc_and_configure_domain</a> to allocate
! the new nest and link it as a child of this grid.
! 
! 4.c.2  Call <a href=setup_Timekeeping.html>setup_Timekeeping</a> for the nest.
! 
! 4.c.3  Initialize the nest's arrays by calling <a href=med_nest_initial.html>med_nest_initial</a>. This will
! either interpolate data from this grid onto the nest, read it in from a file, or both. In a restart run, this
! is also where the nest reads in its restart file.
! 
! 4.d  If a nest was opened above, check for and resolve overlaps (this is not implemented in WRF 2.0, which
! supports multiple nests on the same level but does not support overlapping).
! 
! 4.e  Give the mediation layer an opportunity to do something before the solver is called by
! calling <a href=med_before_solve_io.html>med_before_solve_io</a>. In WRF this is the point at which history and
! restart data is output.
! 
! 4.f  Call <a href=solve_interface.html>solve_interface</a> which calls the solver that advance the domain 
! one time step, then advance the domain's clock by calling domain_clockadvance.  
! The enclosing WHILE loop around this section is for handling other domains 
! with which this domain may overlap.  It is not active in WRF 2.0 and only 
! executes one trip.  
! 
! 4.g Call med_calc_model_time and med_after_solve_io, which are stubs in WRF.
! 
! 4.h Iterate over the children of this domain (<tt>DO kid = 1, max_nests</tt>) and check each child pointer to see
! if it is associated (and therefore, active).
! 
! 4.h.1  Force the nested domain boundaries from this domain by calling <a href=med_nest_force.html>med_nest_force</a>.
! 
! 4.h.2  Setup the time period over which the nest is to run. Sine the current grid has been advanced one time step
! and the nest has not, the start for the nest is this grid's current time minus one time step.  The nest's stop_subtime
! is the current time, bringing the nest up the same time level as this grid, its parent.
! 
! 4.h.3  Recursively call this routine, integrate, to advance the nest's time.  Since it is recursive, this will
! also advance all the domains who are nests of this nest and so on.  In other words, when this call returns, all
! the domains rooted at the nest will be at the current time.
! 
! 4.h.4  Feedback data from the nested domain back onto this domain by calling <a href=med_nest_feedback.html>med_nest_feedback</a>.
! 
! 4.i  Write the time to compute this grid and its subtree. This marks the end of the loop begun at step 4, above.
! 
! 5. Give the mediation layer an opportunity to do I/O at the end of the sequence of steps that brought the
! grid up to stop_subtime with a call to <a href=med_last_solve_io.html>med_last_solve_io</a>.  In WRF, this 
! is used to generate the final history and/or restart output when the domain reaches the end of it's integration.
! There is logic here to make sure this occurs correctly on a nest, since the nest may finish before its parent.
! </DESCRIPTION>

   !  Local data.

   CHARACTER*32                           :: outname, rstname
   TYPE(domain) , POINTER                 :: grid_ptr , new_nest
   TYPE(domain)                           :: intermediate_grid
   INTEGER                                :: step
   INTEGER                                :: nestid , kid
   LOGICAL                                :: a_nest_was_opened
   INTEGER                                :: fid , rid
   LOGICAL                                :: lbc_opened
   REAL                                   :: time, btime, bfrq
   CHARACTER*256                          :: message, message2,message3
   TYPE (grid_config_rec_type)            :: config_flags
   LOGICAL , EXTERNAL                     :: wrf_dm_on_monitor
   INTEGER                                :: idum1 , idum2 , ierr , open_status
   LOGICAL                                :: should_do_last_io

   ! interface
   INTERFACE
       ! mediation-supplied solver
     SUBROUTINE solve_interface ( grid )
       USE module_domain
       TYPE (domain) grid
     END SUBROUTINE solve_interface
       ! mediation-supplied routine to allow driver to pass time information
       ! down to mediation/model layer
     SUBROUTINE med_calc_model_time ( grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain) grid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_calc_model_time
       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! perform I/O before the call to the solve routine
     SUBROUTINE med_before_solve_io ( grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain) grid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_before_solve_io
       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! perform I/O after the call to the solve routine
     SUBROUTINE med_after_solve_io ( grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain) grid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_after_solve_io
       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! perform I/O to initialize a new nest
     SUBROUTINE med_pre_nest_initial ( parent , newid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain), POINTER ::  parent
       INTEGER, INTENT(IN)    ::  newid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_pre_nest_initial
     SUBROUTINE med_nest_initial ( parent , grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain), POINTER ::  grid , parent
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_nest_initial
       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! provide parent->nest forcing
     SUBROUTINE med_nest_force ( parent , grid )
       USE module_domain
       USE module_configure
       TYPE (domain), POINTER ::  grid, parent
     END SUBROUTINE med_nest_force

#ifdef MOVE_NESTS
     SUBROUTINE med_nest_move ( parent , grid )
       USE module_domain
       USE module_configure
       TYPE (domain), POINTER ::  grid , parent
     END SUBROUTINE med_nest_move
#endif

       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! provide parent->nest feedback
     SUBROUTINE med_nest_feedback ( parent , grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain), POINTER ::  grid , parent
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_nest_feedback

       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! perform I/O prior to the close of this call to integrate
     SUBROUTINE med_last_solve_io ( grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain) grid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_last_solve_io
       ! mediation-supplied routine that gives mediation layer opportunity to 
       ! perform setup before iteration over steps in this call to integrate
     SUBROUTINE med_setup_step ( grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain) grid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_setup_step
       ! mediation-supplied routine that gives mediation layer opportunity to
       ! perform setup before iteration over steps in this call to integrate
     SUBROUTINE med_endup_step ( grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain) grid
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_endup_step
       ! mediation-supplied routine that intializes the nest from the grid
       ! by interpolation

     SUBROUTINE Setup_Timekeeping( grid )
       USE module_domain
       TYPE(domain), POINTER :: grid
     END SUBROUTINE

     SUBROUTINE dfi_accumulate( grid )
       USE module_domain
       TYPE(domain), POINTER :: grid
     END SUBROUTINE

   END INTERFACE

!Added by Thiago: Forces only the first call to this INTEGRATE routine to initialize the module responsible for synchronizing
!domain access to different non-thread safe model components when running in threaded integration mode.
#if defined(HRD_MULTIPLE_STORMS) && defined(HRD_THREADED_INTEGRATION)
   IF(grid%id == 1) THEN
       !Initialize the module used to synchronize model tasks among domains in the same nesting level
       CALL init_module_sync_domains()

       !Reset all of domain D01's MPI communicators
       grid%communicator   = GET_DOMAIN_MPI_COMMUNICATOR(1)     !Initially, D01's communicator was the same as "MODULE_DM:local_communicator", now it carries its own unique communicator.
       grid%communicator_x = GET_DOMAIN_MPI_COMMUNICATOR_X(1)   !Set the MPI communicator for the x-direction of the mesh
       grid%communicator_y = GET_DOMAIN_MPI_COMMUNICATOR_Y(1)   !Set the MPI communicator for the y-direction of the mesh

       !Store the D01 domain reference for global access
       CALL STORE_DOMAIN_POINTER(grid)

       !Initialize the parallel integration P-threads driver, a C-based bridge. The current implementation spawns new P-thread in each time-step.
       CALL init_threaded_integration(grid%max_dom,ierr)
       IF(ierr .NE. 0) THEN
           CALL WRF_ERROR_FATAL('Could not initialize the parallel integration p-threads driver.')
       ENDIF
   ENDIF
#endif

   ! This allows us to reference the current grid from anywhere beneath 
   ! this point for debugging purposes.  
   CALL set_current_grid_ptr( grid )

   IF ( .NOT. domain_clockisstoptime( grid ) ) THEN
      CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )
      IF ( config_flags%grid_allowed ) THEN
         CALL domain_clockprint ( 150, grid, 'DEBUG:  top of integrate(),' )
         DO WHILE ( .NOT. domain_clockisstopsubtime(grid) )
            !Added by Thiago for debug purposes
            IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ TOP OF MAIN INTEGRATION LOOP FOR DOMAIN ',grid%id

            IF ( wrf_dm_on_monitor() ) THEN
#if defined(HRD_THREADED_INTEGRATION)
               IF(grid%id == 1) CALL start_timing   !Condition added by Thiago: Only measure timing for D01 till MODULE_TIMING is made thread-safe
#else
               CALL start_timing
#endif
            END IF
            CALL med_setup_step ( grid , config_flags )

            CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id) !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            a_nest_was_opened = .false.
            ! for each nest whose time has come...
            DO WHILE ( nests_to_open( grid , nestid , kid ) )
               !Added by Thiago to catch OpenMP nest bug
               IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' opened nest=',nestid,'/kid=',kid

               ! nestid is index into model_config_rec (module_configure) of the grid
               ! to be opened; kid is index into an open slot in grid's list of children
               a_nest_was_opened = .true.
               CALL med_pre_nest_initial ( grid , nestid , config_flags )
               CALL alloc_and_configure_domain ( domain_id  = nestid ,   &
                                                 grid       = new_nest , &
                                                 parent     = grid ,     &
                                                 kid        = kid        )
               CALL Setup_Timekeeping (new_nest)
               CALL med_nest_initial ( grid , new_nest , config_flags )
               IF ( grid%dfi_stage == DFI_STARTFWD ) THEN
                  CALL wrf_dfi_startfwd_init(new_nest)
               ENDIF

#if defined(HRD_MULTIPLE_STORMS) && defined(HRD_THREADED_INTEGRATION)
               !Added by Thiago to allow each grid to use its own MPI communicator
               new_nest%communicator   = GET_DOMAIN_MPI_COMMUNICATOR(nestid)    !Set the unique communicator for this domain
               new_nest%communicator_x = GET_DOMAIN_MPI_COMMUNICATOR_X(nestid)  !Set the unique x-direction communicator for this domain. AFAIK, this is needed only by routine VT4_MOVE (dyn_nmm/module_STATS_FOR_MOVE.F)
               new_nest%communicator_y = GET_DOMAIN_MPI_COMMUNICATOR_Y(nestid)  !Set the unique y-direction communicator for this domain. AFAIK, this is needed only by routine VT4_MOVE (dyn_nmm/module_STATS_FOR_MOVE.F)

               !Store reference to this nest's domain pointer for global access through MODULE_SYNC_DOMAINS.
               CALL STORE_DOMAIN_POINTER(new_nest)  !Added by Thiago
#endif
            END DO
            IF ( a_nest_was_opened ) THEN
               CALL set_overlaps ( grid )   ! find overlapping and set pointers
            END IF
            CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)    !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel

            ! Accumulation calculation for DFI
            CALL dfi_accumulate ( grid )

            CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id) !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' calling med_before_solve_io'
            CALL med_before_solve_io ( grid , config_flags )
            CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)    !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel

            grid_ptr => grid
            DO WHILE ( ASSOCIATED( grid_ptr ) )
               CALL set_current_grid_ptr( grid_ptr )
               CALL wrf_debug( 100 , 'module_integrate: calling solve interface ' )

               !Wait until the last domain in current nesting level gives I/O priviledge back to the first domain in the current nesting level.
               !At this point, all domains in the same nesting level have completed the previous task and are synced/ready to integrate in parallel.
               CALL NEST_LEVEL_BARRIER_PTHREAD(grid_ptr%id) ; IF(AM_I_FIRST_DOMAIN_IN_MY_LEVEL(grid_ptr%id) .AND. grid_ptr%id > 1 .AND. IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ STARTING PARALLEL INTEGRATION'
               IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' calling solve_interface'
               CALL solve_interface ( grid_ptr ) 
               CALL NEST_LEVEL_BARRIER_PTHREAD(grid_ptr%id) ; IF(AM_I_LAST_DOMAIN_IN_MY_LEVEL(grid_ptr%id) .AND. grid_ptr%id > 1 .AND. IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ COMPLETED PARALLEL INTEGRATION'

               CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id)  !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
               CALL domain_clockadvance ( grid_ptr )
               CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)     !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel

               CALL wrf_debug( 100 , 'module_integrate: back from solve interface ' )
               ! print lots of time-related information for testing
               ! switch this on with namelist variable self_test_domain
               CALL domain_time_test( grid_ptr, 'domain_clockadvance' )
               grid_ptr => grid_ptr%sibling
            END DO
            CALL set_current_grid_ptr( grid )
            CALL med_calc_model_time ( grid , config_flags )

            CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id) !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' calling med_after_solve_io'
            CALL med_after_solve_io ( grid , config_flags )
            CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)    !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel

#if !defined(HRD_MULTIPLE_STORMS) || !defined(HRD_THREADED_INTEGRATION)
            !Thiago comment: This is the original recursive and sequential child nest integration loop that performs forcing, integration, and feedback sequentially for each nest.
            grid_ptr => grid
            DO WHILE ( ASSOCIATED( grid_ptr ) )
               DO kid = 1, max_nests
                 IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                   CALL set_current_grid_ptr( grid_ptr%nests(kid)%ptr )
                   ! Recursive -- advance nests from previous time level to this time level.
                   CALL wrf_debug( 100 , 'module_integrate: calling med_nest_force ' )
                   CALL med_nest_force ( grid_ptr , grid_ptr%nests(kid)%ptr )
                   CALL wrf_debug( 100 , 'module_integrate: back from med_nest_force ' )
                   grid_ptr%nests(kid)%ptr%start_subtime = &
                     domain_get_current_time(grid) - domain_get_time_step(grid)
                   grid_ptr%nests(kid)%ptr%stop_subtime = &
                     domain_get_current_time(grid)
                   CALL integrate ( grid_ptr%nests(kid)%ptr )
                   CALL wrf_debug( 100 , 'module_integrate: back from recursive call to integrate ' )
                   IF ( .NOT. ( domain_clockisstoptime(head_grid              ) .OR. &
                                domain_clockisstoptime(grid                   ) .OR. &
                                domain_clockisstoptime(grid_ptr%nests(kid)%ptr) ) )  THEN
                     CALL wrf_debug( 100 , 'module_integrate: calling med_nest_feedback ' )
                     CALL med_nest_feedback ( grid_ptr , grid_ptr%nests(kid)%ptr , config_flags )
                     CALL wrf_debug( 100 , 'module_integrate: back from med_nest_feedback ' )
                   END IF
#ifdef MOVE_NESTS
                   IF ( .NOT. domain_clockisstoptime( head_grid ) ) THEN
                     CALL med_nest_move ( grid_ptr , grid_ptr%nests(kid)%ptr )
                   ENDIF
#endif
                 END IF
               END DO
               grid_ptr => grid_ptr%sibling
            END DO
            CALL set_current_grid_ptr( grid )
#else
            !Added by Thiago: The following 3 loops perform the same operation as the original sequential loop above, but allows nested domains to be integrated in parallel, while forcing and feedback are still executed sequentially.
            CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id) !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            grid_ptr => grid
            DO WHILE ( ASSOCIATED( grid_ptr ) )
               DO kid = 1, max_nests
                 IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                   IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' forcing nest=',grid_ptr%nests(kid)%ptr%id,'/kid=',kid
                   CALL set_current_grid_ptr( grid_ptr%nests(kid)%ptr )
                   CALL wrf_debug( 100 , 'module_integrate: calling med_nest_force ' )
                   CALL med_nest_force ( grid_ptr , grid_ptr%nests(kid)%ptr )
                   CALL wrf_debug( 100 , 'module_integrate: back from med_nest_force ' )
                 END IF
               END DO
               grid_ptr => grid_ptr%sibling
            END DO
            CALL set_current_grid_ptr( grid )

            !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            IF(GET_DOMAIN_NESTING_LEVEL(grid%id) < GET_MAX_SIMULATION_NEST_LEVEL()) THEN
                CALL CYCLE_DOMAIN_PRIVILEGE(grid%id,.TRUE.) !The last domain in current nesting level transfer priviledge the first domain in the next nesting level
            ELSE
                CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)        !The current domain transfers priviledge to the next domain in its nesting level
            ENDIF

            !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            !Integrate D01's children nests in parallel
            grid_ptr => grid
            DO WHILE ( ASSOCIATED( grid_ptr ) )
               DO kid = 1, max_nests
                 IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                   IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' integrating nest=',grid_ptr%nests(kid)%ptr%id,'/kid=',kid
                   CALL set_current_grid_ptr( grid_ptr%nests(kid)%ptr )
                   ! Recursive -- advance nests from previous time level to this time level.
                   grid_ptr%nests(kid)%ptr%start_subtime = &
                     domain_get_current_time(grid) - domain_get_time_step(grid)
                   grid_ptr%nests(kid)%ptr%stop_subtime = &
                     domain_get_current_time(grid)

                   !Added by Thiago: Domain D01 spawns threads to invoke INTEGRATE for multiple nests. Otherwise, all domains invoke INTEGRATE directly.
                   IF(grid%num_nests == 1) THEN
                       CALL integrate ( grid_ptr%nests(kid)%ptr )
                       CALL wrf_debug( 100 , 'module_integrate: back from recursive call to integrate ' )
                   ELSE
                       !Only D01 should reach this code since it is the only domain possibly having multiple children nests.
                       CALL integrate_domain_by_thread(grid_ptr%nests(kid)%ptr%id,ierr) !Pass the nest ID to the "C" function because the grid object itself cannot be passed.
                       IF(ierr .NE. 0) THEN
                           CALL WRF_ERROR_FATAL('D01 failed to spawn p-threads to integrate child nest.')
                       ENDIF
                   ENDIF
                 END IF
               END DO
               grid_ptr => grid_ptr%sibling
            END DO
            !Added by Thiago to wait until all of D01 children' nests have been integrated
            IF(grid%num_nests > 1) THEN
                CALL wait_for_integration_done() !Wait for parallel integration to complete (i.e., all p-threads have joined).
                CALL wrf_debug( 100 , 'module_integrate: back from recursive call to integrate ' )
            ENDIF
            CALL set_current_grid_ptr( grid )

            !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            !The last domain in current nesting level transfer priviledge the first domain within its nesting level for the next iteration
            IF(GET_DOMAIN_NESTING_LEVEL(grid%id) < GET_MAX_SIMULATION_NEST_LEVEL() .AND. AM_I_LAST_DOMAIN_IN_MY_LEVEL(grid%id)) CALL ASSIGN_PRIVILEGE_TO_NEST_LEVEL(GET_DOMAIN_NESTING_LEVEL(grid%id))

            !Synchronize the feedback and moving steps of every domain from their children nests
            CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id) !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            grid_ptr => grid
            DO WHILE ( ASSOCIATED( grid_ptr ) )
               DO kid = 1, max_nests
                 IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                   IF(IS_INTEGRATION_VERBOSE()) WRITE(*,*) '@@@@@ DOMAIN ',grid%id,' feeding back nest=',grid_ptr%nests(kid)%ptr%id,'/kid=',kid
                   CALL set_current_grid_ptr( grid_ptr%nests(kid)%ptr )
                   IF(.NOT. (domain_clockisstoptime(head_grid) .OR. domain_clockisstoptime(grid) .OR. domain_clockisstoptime(grid_ptr%nests(kid)%ptr))) THEN
                       CALL wrf_debug( 100 , 'module_integrate: calling med_nest_feedback ' )
                       CALL med_nest_feedback ( grid_ptr , grid_ptr%nests(kid)%ptr , config_flags )
                       CALL wrf_debug( 100 , 'module_integrate: back from med_nest_feedback ' )
                   ENDIF
#ifdef MOVE_NESTS
                   IF ( .NOT. domain_clockisstoptime( head_grid ) ) THEN
                     CALL med_nest_move ( grid_ptr , grid_ptr%nests(kid)%ptr )
                   ENDIF
#endif
                 END IF
               END DO
               grid_ptr => grid_ptr%sibling
            END DO
            CALL set_current_grid_ptr( grid )
            CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)    !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
#endif

            !  Report on the timing for a single time step.
            IF ( wrf_dm_on_monitor() ) THEN
#if defined(HRD_THREADED_INTEGRATION)
              IF(grid%id == 1) THEN   !Condition added by Thiago: Only measure timing for D01 till MODULE_TIMING is made thread-safe
#endif
               CALL domain_clock_get ( grid, current_timestr=message2 )
#if (EM_CORE == 1)
               if (config_flags%use_adaptive_time_step) then
                  WRITE ( message , FMT = '("main (dt=",F6.2,"): time ",A," on domain ",I3)' ) grid%dt, TRIM(message2), grid%id
               else
                  WRITE ( message , FMT = '("main: time ",A," on domain ",I3)' ) TRIM(message2), grid%id
               endif
#else
                  WRITE ( message , FMT = '("main: time ",A," on domain ",I3)' ) TRIM(message2), grid%id
#endif
               CALL end_timing ( TRIM(message) )
#if defined(HRD_THREADED_INTEGRATION)
              ENDIF
#endif
            END IF
            CALL med_endup_step ( grid , config_flags )
         END DO

         ! Accumulation calculation for DFI
         CALL dfi_accumulate ( grid )

         ! Avoid double writes on nests if this is not really the last time;
         ! Do check for write if the parent domain is ending.
         IF ( grid%id .EQ. 1 ) THEN               ! head_grid
            CALL med_last_solve_io ( grid , config_flags )
         ELSE
! zip up the tree and see if any ancestor is at its stop time
            should_do_last_io = domain_clockisstoptime( head_grid )
            grid_ptr => grid 
            DO WHILE ( grid_ptr%id .NE. 1 )
               IF ( domain_clockisstoptime( grid_ptr ) ) THEN
                  should_do_last_io = .TRUE. 
               END IF
               grid_ptr => grid_ptr%parents(1)%ptr
            ENDDO

            CALL WAIT_FOR_DOMAIN_PRIVILEGE(grid%id) !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
            IF ( should_do_last_io ) THEN
               grid_ptr => grid
               CALL med_nest_feedback ( grid_ptr%parents(1)%ptr, grid , config_flags )
               CALL med_last_solve_io ( grid , config_flags )
            ENDIF
            CALL CYCLE_DOMAIN_PRIVILEGE(grid%id)    !Added by Thiago to synchronize I/O access by the high-res nest-pairs integrated in parallel
         ENDIF
      ENDIF
   END IF

END SUBROUTINE integrate

!Added by Thiago: copied from solve_nmm, this function provides an alternative to the non-thread safe module_timing.
REAL*8 FUNCTION HIRESTIMER()
    implicit none
#if defined(OLD_TIMERS)
    INTEGER :: IC,IR
    CALL SYSTEM_CLOCK(COUNT=IC,COUNT_RATE=IR)
    HIRESTIMER=REAL(IC)/REAL(IR)*1000.0
#else
    call hires_timer(HIRESTIMER)
    HIRESTIMER=HIRESTIMER*1000
#endif
END FUNCTION HIRESTIMER

END MODULE module_integrate

