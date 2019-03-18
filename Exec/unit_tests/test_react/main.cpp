
#include <new>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

#include <time.h>

#ifdef HAS_DUMPMODEL
#include <DumpModel1d.H>
#endif

#ifdef HAS_XGRAPH
#include <XGraph1d.H>
#endif

#include "Castro_io.H"


extern "C"
{
   void do_burn();
}


std::string inputs_name = "";

int
main (int   argc,
      char* argv[])
{

    //
    // Make sure to catch new failures.
    //
    amrex::Initialize(argc,argv);

    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
	inputs_name = argv[1];
      }
    }

    do_burn();

    amrex::Finalize();

    return 0;
}
