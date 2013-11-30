#if defined(HRD_MULTIPLE_STORMS) && defined(HRD_THREADED_INTEGRATION)
/**
    This code is part of "Fortran-C-Fortran" bridge designed to allow MODULE_SYNC_DOMAINS to use MUTEXES and CONDITION VARIABLES
    to schedule domain task privilege access instead of using the Fortran SLEEP subroutine to create a polling wait loop.
**/

//Required header files
#include <pthread.h>
#include <stdio.h>


//Some global variables required for parallel p-thread based synchronization.
#define MAX_DOMAINS_IN_PARALLEL 40
static pthread_mutex_t domainTaskPrivilegeSyncMutex[MAX_DOMAINS_IN_PARALLEL];       //Allows MODULE_SYNC_DOMAINS to use MUTEXES and CONDITION VARIABLES to schedule..
static pthread_cond_t domainTaskPrivilegeSyncCondition[MAX_DOMAINS_IN_PARALLEL];    //...domain task privilege instead of using the Fortran SLEEP subroutine.

static pthread_barrier_t domainTaskPrivilegeSyncBarrier;                            //Used to create a barrier for all integration threads


//This function initializes the p-thread mutex variables required for synchronization. It should be called once by WRF.
void init_pthread_sync_interface_(int* max_dom, int* totalStorms, int* ierr)
{
    //Local variables
    int i, count;

    //Ensure the number of domains is within valid range
    *ierr=0;    //Initialize the error flag to 0
    if(*max_dom < 1 && *max_dom > MAX_DOMAINS_IN_PARALLEL)
    {
        *ierr=99999999;
        return;
    }

    //Initialize the MUTEX and CONDITION used by MODULE_SYNC_DOMAIN, only enough to cover the number of domains in the current simulation (conserve memory).
    for(i=0;i<*max_dom;i++) //Note:  C indexing starts at 0, but lowest domain ID is 1.
    {
        *ierr=pthread_mutex_init(&domainTaskPrivilegeSyncMutex[i], NULL);
        if(*ierr!=0) return;

        *ierr=pthread_cond_init (&domainTaskPrivilegeSyncCondition[i], NULL);
        if(*ierr!=0) return;
    }

    //Initialize the synchronization barrier
    //printf("@@@@@ Initializing p-thread barrier.\n");
    count = (*totalStorms) < 1 ? 1 : *totalStorms;  //A value of 1 causes the P-thread barrier to return immediately
    *ierr=pthread_barrier_init(&domainTaskPrivilegeSyncBarrier, NULL, count);
    if(*ierr!=0) return;
}


//This function creates a barrier for the integration threads
void integration_pthread_barrier_(int* ierr)
{
    //Wait for all domains to reach the barrier
    int rc = pthread_barrier_wait(&domainTaskPrivilegeSyncBarrier);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
        *ierr=-1;
    else
        *ierr=0;
}


//This function blocks until the domain with ID "theDomainId" has the privilege to perform syncrhonized tasks.
void block_till_domain_privilege_(int* theDomainId, int* currDomainWithPrivilege, int* ierr)
{
    //Lock the mutex
    //*ierr=pthread_mutex_lock(&domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);  //Lock the corresponding mutex for this domain. Note: C indexing starts at 0, but lowest domain ID is 1.
    *ierr=pthread_mutex_lock(&domainTaskPrivilegeSyncMutex[0]);  //Lock the corresponding mutex for this domain. Note: C indexing starts at 0, but lowest domain ID is 1.
    if(*ierr!=0) return;

    //Block until the domain has the privilege
    while (*theDomainId != *currDomainWithPrivilege)
    {
        //printf("@@@@@ WAITING FOR DOMAIN %d PRIVILEGE\n",*theDomainId);
        //*ierr=pthread_cond_wait(&domainTaskPrivilegeSyncCondition[(*theDomainId)-1], &domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);
        *ierr=pthread_cond_wait(&domainTaskPrivilegeSyncCondition[0], &domainTaskPrivilegeSyncMutex[0]);
        if(*ierr!=0)
        {
            //pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);
            pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[0]);
            return;
        }
    }
    //printf("@@@@@ DOMAIN %d HAS PRIVILEGE\n",*theDomainId);

    //Unlock the mutex
    //*ierr=pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);
    *ierr=pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[0]);
    if(*ierr!=0) return;
}


//This function awakes all threads waiting for the specified domain task privilege
void force_check_domain_privilege_(int* theDomainId, int* ierr)
{
    //Lock the mutex
    //*ierr=pthread_mutex_lock(&domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);
    *ierr=pthread_mutex_lock(&domainTaskPrivilegeSyncMutex[0]);
    if(*ierr!=0) return;

    //Awaken any/all threads waiting on this condition
    //*ierr=pthread_cond_broadcast(&domainTaskPrivilegeSyncCondition[(*theDomainId)-1]);
    *ierr=pthread_cond_broadcast(&domainTaskPrivilegeSyncCondition[0]);
    if(*ierr!=0)
    {
        //pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);
        pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[0]);
        return;
    }

    //Unlock the mutex
    //printf("@@@@@ AWAKING THREADS WAITING FOR DOMAIN %d PRIVILEGE\n",*theDomainId);
    //*ierr=pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[(*theDomainId)-1]);
    *ierr=pthread_mutex_unlock(&domainTaskPrivilegeSyncMutex[0]);
    if(*ierr!=0) return;
}

#else
void dummyuselessfunction_(){}
#endif
