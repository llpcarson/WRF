#!/bin/ksh -l

#Search for the "module module_integrate" and "end module" keywords, which coud have spaces and tabs in "frame/module_integrate.F"
startLine=`grep -i "^[\t ]\{0,\}module[\t ]\{1,\}module_integrate" frame/module_integrate.F | awk '(NR==1){print}'`
if(test -z "$startLine") then
    echo "Error: Failed to locate (case-insensitive) sentence 'MODULE MODULE_INTEGRATE' in frame/module_integrate.F. Threaded integration will fail!!!"
    exit 1
fi
endLine=`grep -i "^[\t ]\{0,\}end[\t ]\{1,\}module" frame/module_integrate.F | awk 'END{print}'`
if(test -z "$endLine") then
    echo "Error: Failed to locate (case-insensitive) sentence 'END MODULE' in frame/module_integrate.F. Threaded integration will fail!!!"
    exit 2
fi

#Create the file "frame/module_integrate_nests.F" as a copy of "frame/module_integrate.F" having replaced module begin and end declaration lines.
sed "s:$startLine:MODULE MODULE_INTEGRATE_NESTS:g" frame/module_integrate.F | sed "s:$endLine:END MODULE MODULE_INTEGRATE_NESTS:g" > frame/module_integrate_nests.F
if(! test -s frame/module_integrate_nests.F) then
    echo "Error: Failed to create file 'frame/module_integrate_nests.F'."
    exit 3
fi

#Create the file "frame/pthread_integrate_domain_dummy.c" which corresponds to empty implementations of the routines in "frame/pthread_integrate_domain.c"
printf "%b" "//These are simply dummy declarations of routines used by MODULE_INTEGRATE::INTEGRATE() but never used by MODULE_INTEGRATE_NESTS::INTEGRATE().\n" > frame/pthread_integrate_domain_dummy.c
printf "%b" "void init_threaded_integration_(int* max_dom, int* ierr){}\n" >> frame/pthread_integrate_domain_dummy.c
printf "%b" "void integrate_domain_by_thread_(int *domainId, int *ierr){}\n" >> frame/pthread_integrate_domain_dummy.c
printf "%b" "void wait_for_integration_done_(){}\n" >> frame/pthread_integrate_domain_dummy.c
