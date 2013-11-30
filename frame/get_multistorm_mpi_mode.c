#if defined(HRD_MULTIPLE_STORMS) && defined(HRD_THREADED_INTEGRATION)
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

/*void removeWhitespaces(char *s)
{
    char *p = s;
    if(!*s) return;
    while((*s++))
    {
        if(isspace(*s)) continue;
        *p++ = *s;
    }
    *p = '\0';
}*/
void removeWhitespaces(char *s)
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

void toLowerStr(char *str)
{
    if(!*str) return 0; //String is NULL
    char* pstr = str;
    while (*pstr)
    {
        *pstr = (char)tolower(*pstr);
        pstr++;
    }
}

int startsWithPrefix(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? 0 : strncmp(pre, str, lenpre) == 0;
}

int isStringNumeric(const char* str)
{
    if(!*str) return 0; //String is NULL
    while(*str)
    {
        if(!isdigit(*str)) return 0;
        str++;
    }
    return 1;
}

//This function is invoked by MODULE_DM::SPLIT_COMMUNICATOR() to determine the MPI mode used to run the model
//in threaded multistorm integration mode. It first looks at environmental variable WRF_NMM_MPI_MODE, whose
//value can be one of MPI_THREAD_{SINGLE,FUNNELED,SERIALIZED,MULTIPLE}. If not available, it will then count
//the number of storms from the namelist.input file. It assumes that variables "max_dom" and "parent_id"
//(comma separated list) are listed in the namelist.input file. Upon error, MPI_THREAD_MULTIPLE mode is used.
void get_multistorm_mpi_mode_(int* mpiMode, int* totalStorms, int *ierr)
{
    int max_dom=0;
    const char filename[] = "namelist.input";   //Assumes namelist.input is in the current directory
    char* envMpiMode;
    *ierr = 0;
    *totalStorms = 0;
    *mpiMode = MPI_THREAD_MULTIPLE; //Default MPI mode

    //Search for the environmental varible WRF_NMM_MPI_MODE.
    envMpiMode=getenv("WRF_NMM_MPI_MODE");
    if(envMpiMode != NULL)
    {
        //Remove whitespaces and lower case letters
        printf("WRF_NMM_MPI_MODE=%s\n",envMpiMode);
        removeWhitespaces(envMpiMode);
        printf("WRF_NMM_MPI_MODE=%s\n",envMpiMode);
        toLowerStr(envMpiMode);
        printf("WRF_NMM_MPI_MODE=%s\n",envMpiMode);

        //Translate the MPI mode
        if(startsWithPrefix("mpi_thread_single",envMpiMode))
        {
            *mpiMode = MPI_THREAD_SINGLE;
            return;
        }
        else if(startsWithPrefix("mpi_thread_funneled",envMpiMode))
        {
            *mpiMode = MPI_THREAD_FUNNELED;
            return;
        }
        else if(startsWithPrefix("mpi_thread_serialized",envMpiMode))
        {
            *mpiMode = MPI_THREAD_SERIALIZED;
            return;
        }
        else if(startsWithPrefix("mpi_thread_multiple",envMpiMode))
        {
            *mpiMode = MPI_THREAD_MULTIPLE;
            return;
        }
        else if(strlen(envMpiMode))
            printf("@@@@@ WARNING in get_multistorm_mpi_mode.c: Invalid value of environmental variable WRF_NMM_MPI_MODE.\n");
    }

    //At this stage, the user-specified MPI mode could not be determined from the environment.
    //Let's determine the number of storms from the namelist.input file.
    FILE *file = fopen(filename, "r");
    if ( file )
    {
        int countDomains;   //Used to ensure that the parent ID of all domains are accounted for.
        char line[1024], parentIdLine[1024];
        char* pStr;
        parentIdLine[0]='\0';
        while(fgets(line,sizeof line,file))
        {
            //Remove whitespaces and ensure line is not empty
            removeWhitespaces(line);
            if(strlen(line)<=0) continue;
            toLowerStr(line);

            //Look for max_dom in the namelist
            if(startsWithPrefix("max_dom=",line))
            {
                //Extract the number of domains
                if(sscanf(line,"max_dom=%d%*s",&max_dom)!=1)
                {
                    *ierr = 2;  //The value of max_dom is not numeric
                    return;
                }
            }

            //Look for parent_id in the namelist.
            if(startsWithPrefix("parent_id=",line))
            {
                //Store the line
                strcpy(parentIdLine,line);
            }
        }

        //Close the namelist.input file.
        fclose(file);

        //Count the number of domains having parent ID 1 (D01).
        if(max_dom>0 && parentIdLine[0])
        {
            countDomains=0;
            strtok(parentIdLine,"=");   //Remove "parent_id="
            while((pStr=strtok(NULL,","))!=NULL && (++countDomains)<=max_dom)   //Attempt to process only as many entries as domains.
            {
                //Check if the expected parent domain ID is numeric
                if(!isStringNumeric(pStr))
                {
                    *ierr=4;    //Variable "parent_id" has non-numeric entries.
                    return;
                }
                if(atoi(pStr)==1) (*totalStorms)++;
            }

            //Ensure enough entries were processed
            if(countDomains<max_dom)
            {
                *ierr=5;    //Variable parent_id does not contain at least max_dom comma-separated entries.
                return;
            }
            else if(*totalStorms<=1) *mpiMode = MPI_THREAD_SINGLE;
        }
        else
        {
            *ierr=3;    //Either invalid value of variable max_dom or (comma-separated) parent_id.
            return;
        }
    }
    else
    {
        *ierr=1;    //Cannot locate file "./namelist.input"
        return;
    }
}


/*int main(void)
{
    int totalStorms, mpiMode, ierr;
    get_multistorm_mpi_mode_(&mpiMode,&totalStorms,&ierr);
    printf("ierr=%d totalStorms=%d mpiMode=%d\n",ierr,totalStorms,mpiMode);
}*/
#else
void dummyuselessfunction_(){}
#endif
