/**
    This code is part of "Fortran-C-Fortran" bridge designed to allow WRF's recursive INTEGRATE subroutine to
    spawn P-threads that integrate multiple nest-pairs in parallel. Each thread handles the integrattion of
    nests tracking a specific storm.
**/

//Required header files
#include <pthread.h>
#include <stdio.h>


//Prototype of the Fortran routine responsible for invoking WRFS's recursive INTEGRATE subroutine.
extern void threaded_integration_driver_(int *domainId);


//Some global variables required for parallel p-thread based integration.
#define MAX_DOMAINS_IN_PARALLEL 40
static pthread_t integrationPThreads[MAX_DOMAINS_IN_PARALLEL];                      //Pointer to the threads used to integrate model domains in parallel
static int integrationPThreadsActive[MAX_DOMAINS_IN_PARALLEL];                      //Indicates whether a certain thread is active, 0 or 1.
static pthread_attr_t integrationPThreadsAttrs;                                     //The attributes of the threads, i.e., 128MB stack size for WRF-NMM.


//This function initializes the thread variables required for parallel integration. It should be called once by WRF.
void init_hrd_pthread_interface_(int* max_dom, int* ierr)
{
    //Local variables
    int i;

    //Ensure the number of domains is within valid range
    if(*max_dom > MAX_DOMAINS_IN_PARALLEL)
    {
        *ierr=99999999;
        return;
    }

    //Set all thread activity flags to 0
    for(i=0;i<MAX_DOMAINS_IN_PARALLEL;i++)
        integrationPThreadsActive[i] = 0;

    //Set the desired p-thread attributes
    *ierr=pthread_attr_init(&integrationPThreadsAttrs);                      //Initialize the p-thread attribute object.
    if(*ierr!=0) return;

    *ierr=pthread_attr_setstacksize(&integrationPThreadsAttrs, 134217728);   //Set the p-thread stack size to 128MB, because the SOLVE_NMM call tree uses lots of temporary memory.
    if(*ierr!=0) return;
}


//This function is invoked by a worker p-thread. It invokes a Fortran routine that then calls WRFS's recursive INTEGRATE subroutine for the specified domain.
void* pthread_integrate_domain_worker(void* domainIdVoidPtr)
{
    //Call the Fortran routine that invokes WRFS's recursive INTEGRATE subroutine
    threaded_integration_driver_((int *)domainIdVoidPtr);

    //Exit the thread
    pthread_exit(NULL);
}


//This function is invoked by the Fortran INTEGRATE subroutine to spawn threads that will integrate multiple domains in parallel.
void pthread_integrate_domain_(int *domainId, int *ierr)
{
    //threaded_integration_driver_(domainId);
    pthread_t* thread = &integrationPThreads[(*domainId)-1];                                                //Get the corresponding thread. Note:  C indexing starts at 0, but lowest domain ID is 1.
    *ierr = pthread_create(thread, &integrationPThreadsAttrs, pthread_integrate_domain_worker, domainId);   //Start the thread
    if(*ierr == 0)
    {
        integrationPThreadsActive[(*domainId)-1] = 1;   //Store thread activity information
        //printf("Activating thread %d,%lu,%lu,%lu.\n",((*domainId)-1),*thread,thread,&integrationPThreads[(*domainId)-1]);
    }
}


//This function waits till all active threads have completed their tasks.
void wait_for_integration_done_()
{
    int i;
    pthread_t thread;

    for(i=0;i<MAX_DOMAINS_IN_PARALLEL;i++)
    {
        if(integrationPThreadsActive[i]>0)
        {
            thread = integrationPThreads[i];
            //printf("Deactivating thread %d,%lu,%lu.\n",i,thread,integrationPThreads[i]);
            pthread_join(thread, NULL);
            integrationPThreadsActive[i]=0; //Reset the thread activity flag
        }
    }
}
