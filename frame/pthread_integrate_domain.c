#if defined(HRD_MULTIPLE_STORMS) && defined(HRD_THREADED_INTEGRATION)
/**
    This code is part of "Fortran-C-Fortran" bridge designed to allow WRF's recursive INTEGRATE subroutine to
    spawn P-threads that integrate multiple nest-pairs in parallel. Each thread handles the integrattion of
    nests tracking a specific storm.
**/

//Required header files and defitions
#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>


//Prototype of the Fortran routine responsible for invoking WRFS's recursive INTEGRATE subroutine.
extern void threaded_integration_driver_(int *domainId);


//Some global variables required for parallel p-thread based integration.
#define MAX_DOMAINS_IN_PARALLEL 40
static pthread_t integrationPThreads[MAX_DOMAINS_IN_PARALLEL];                      //Pointer to the threads used to integrate model domains in parallel
static int integrationPThreadsActive[MAX_DOMAINS_IN_PARALLEL];                      //Indicates whether a certain thread is active, 0 or 1.
static int integrationPThreadsAffinity[MAX_DOMAINS_IN_PARALLEL];                    //Points to the cores where each domains should be integrated. A negative value means no set affinity.
static int mainProcessThreadAffinity;                                               //Points to the core where the main program process should be executed. A negative value means no set affinity.
static pthread_attr_t integrationPThreadsAttrs;                                     //The attributes of the threads, i.e., 128MB stack size for WRF-NMM.


//This function removes spaces from a string.
void removeTheWhitespaces(char *s)
{
    char *p = s;
    if(!*s) return;
    while(*s)
    {
        if(isspace(*s))
        {
            s++;
            continue;
        }
        *p++ = *s++;
    }
    *p = '\0';
}


//This subroutine indicates whether all characters in a string are numeric.
int isTheStringNumeric(const char* str)
{
    if(!*str) return 0; //String is NULL
    while(*str)
    {
        if(!isdigit(*str)) return 0;
        str++;
    }
    return 1;
}


//This function initializes the thread variables required for parallel integration. It should be called once by WRF.
//Argument max_dom: total number of domains in simulation
//Argument num_storms: total number of storms integrated in high resolution using moving nests
//Argument first_domains: array carrying a unique sequence of domains IDs representing the first moving nesting level.
//                        This array should have *num_storms elements.
//Argument mpi_node_rank: the local node rank of this process, starting at 0.
//
//This function will attempt to parse the environmental variable HRD_THREADED_INTEGRATION_AFFINITY. This variable should
//contain a colon ":" separated list of sets of comma-separated core IDs. The number of colon-separated entries should equal
//to the number of MPI ranks (e.g. WRF tasks) spawned in each node. Each of these entries describe the IDs of the cores
//into which the main WRF task and the integration p-threads should be bound. Here is an example assuming 12 core nodes and
//3 WRF tasks spawned per node, each  task spawning 4 pthreads that each integrate a different storm in parallel:
//
//  export HRD_THREADED_INTEGRATION_AFFINITY="0,0,1,2,3:4,4,5,6,7:8,8,9,10,11"
//
//The statement above indicates that the first WRF task in the node (local rank 0) can use cores "0,0,1,2,3" to bind its processes.
//Thus, the WRF task will bind its main WRF process and the pthread for storm 1 to the core having ID 0 (the first two comma-separated
//entries are 0) and will bind the pthreads integrating storms 2, 3, and 4 to cores 0, 1, 2, and 3, respectively. The main WRF task
//and one of the pthreads can share the same core because they don't perform work simultaneously. The same idea applies to the second
//and third WRF tasks in the node (local ranks 1 and 2), which will use the cores specified respectively by the entries "4,4,5,6,7" and
//"8,8,9,10,11" to bind their main WRF processes and corresponding integration pthreads to cores.
//
//Any of the comma-separted entries can be left blank or given a negative value to indicate that no specific core affinity should be
//automatically set for any process.
void init_threaded_integration_(int* max_dom, int *num_storms, int *first_domains, int *mpi_node_rank, int* ierr)
{
    //Local variables
    int i, myGlobalRank, myLocalRank, entryFound;
    char *envThreadAffinity, *affinitySet, *coreId;

    //Ensure the number of domains is within valid range
    *ierr=0;    //Default value indicates no error
    if(*max_dom > MAX_DOMAINS_IN_PARALLEL)
    {
        *ierr=99999999;
        return;
    }

    //Set all thread activity flags to 0, and core affinity to a negative value.
    mainProcessThreadAffinity = -1;
    for(i=0;i<MAX_DOMAINS_IN_PARALLEL;i++)
    {
        integrationPThreadsActive[i] = 0;       //Indicates there are no active pthreads
        integrationPThreadsAffinity[i] = -1;    //Indicates no specific core affinity
    }

    //Set the desired p-thread attributes
    *ierr=pthread_attr_init(&integrationPThreadsAttrs);                      //Initialize the p-thread attribute object.
    if(*ierr!=0) return;

    *ierr=pthread_attr_setstacksize(&integrationPThreadsAttrs, 134217728);   //Set the p-thread stack size to 128MB, because the SOLVE_NMM call tree uses lots of temporary memory.
    if(*ierr!=0) return;

    //Search for the environmental variable HRD_THREADED_INTEGRATION_AFFINITY.
    envThreadAffinity=getenv("HRD_THREADED_INTEGRATION_AFFINITY");
    if(envThreadAffinity != NULL && *num_storms >=0 && *mpi_node_rank >=0)
    {
        //Remove whitespaces
        removeTheWhitespaces(envThreadAffinity);
        if(strlen(envThreadAffinity)<=0) return;    //Ignore an empty string

        //Get the core set entry corresponding to the specified MPI local node rank.
        for(i=0;i<=*mpi_node_rank;i++)
        {
            //Get the next core affinity set
            if(i==0) affinitySet=strtok(envThreadAffinity,":"); //Get the first core set
            else affinitySet=strtok(NULL,":");               //Get a core ID entry
        }

        //Exit if no entries were found for the specified local node rank
        if(affinitySet==NULL) return;   //No affinity mask was set
        printf("-Found affinity mask '%s' for local rank %d.\n",affinitySet,*mpi_node_rank);

        //Process the entry that was found, trying to extract as many core IDs as storms
        for(i=0;i<=*num_storms;i++)
        {
            //Get the next core ID and ensure it is numeric
            if(i==0) coreId=strtok(affinitySet,",");    //Get the first core ID entry
            else coreId=strtok(NULL,",");               //Get a core ID entry

            //Ensure the core ID is valid
            if(coreId==NULL)
            {
                if(i==0) printf("   -Skipping empty core ID entry for main WRF task process.");
                else printf("   -Skipping empty core ID entry for storm %d.\n",i);
                continue;
            }
            else if(!isTheStringNumeric(coreId))
            {
                if(i==0) printf("   -Skipping invalid core ID entry '%s' for main WRF task process.",coreId);
                else printf("   -Skipping invalid core ID entry '%s' for storm %d.\n",coreId,i);
                continue;
            }

            //This is a valid entry. Check if it is the first core ID entry.
            if(i==0)
            {
                printf("   -Assigning core ID %s to the main WRF task process (domain 1).\n",coreId);
                mainProcessThreadAffinity = atoi(coreId);
            }
            else
            {
                printf("   -Assigning core ID %s for storm %d (domain %d).\n",coreId,i,first_domains[i-1]);
                integrationPThreadsAffinity[first_domains[i-1]-1] = atoi(coreId);
            }
        }

        //Assign the main process to run in the specified core
        assign_thread_to_core(mainProcessThreadAffinity);   //Negative core ID values are ignored
    }
}


//This function binds the current thread to a specified core.
int assign_thread_to_core(int coreId)
{
    //Negative core ID values are invalid
    if(coreId<0) return -1;

    //Create the CPU affinity set
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(coreId, &cpuset);

    //Attempt to move the thread to the specified core
    //return pthread_setaffinity_np(pthread_self(),sizeof(cpu_set_t),&cpuset);
    return sched_setaffinity(0,sizeof(cpu_set_t),&cpuset);
}


//This function is invoked by a worker p-thread. It invokes a Fortran routine that then calls WRFS's recursive INTEGRATE subroutine for the specified domain.
void* threaded_integration_worker(void* domainIdVoidPtr)
{
    //Apply thread core affinity, if specified.
    int *domainId = (int *)domainIdVoidPtr;
    assign_thread_to_core(integrationPThreadsAffinity[(*domainId)-1]); //Negative core ID values are ignored

    //Call the Fortran routine that invokes WRFS's recursive INTEGRATE subroutine
    threaded_integration_driver_((int *)domainIdVoidPtr);

    //Exit the thread
    pthread_exit(NULL);
}


//This function is invoked by the Fortran INTEGRATE subroutine to spawn threads that will integrate multiple domains in parallel.
void integrate_domain_by_thread_(int *domainId, int *ierr)
{
    //threaded_integration_driver_(domainId);
    pthread_t* thread = &integrationPThreads[(*domainId)-1];                                                //Get the corresponding thread. Note:  C indexing starts at 0, but lowest domain ID is 1.
    *ierr = pthread_create(thread, &integrationPThreadsAttrs, threaded_integration_worker, domainId);   //Start the thread
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

#else
void dummyuselessfunction(){}
#endif
