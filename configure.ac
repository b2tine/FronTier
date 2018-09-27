AC_PREREQ(2.57)
AC_INIT([FronTier],[1.0])
AC_CONFIG_SRCDIR([src/intfc/array.h])
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_BUILD 
AC_PROG_AWK

AC_PROG_CC(gcc gcc295 cc)
AC_PROG_CXX(g++)
AC_PROG_F77(gfortran)
AC_LANG(C)


#######################################################
########          TESTS FOR FRONTIER            #######
#######################################################

AC_CHECK_LIB([m],[feenableexcept],[fenv=" -DHAS_FENV "],[fenv=""])
AC_CHECK_LIB([ieee],[main],[],[])

AC_CHECK_HEADER([malloc.h],[malloch='yes'],[malloch='no'])
systemincs=" "
if test ${malloch} = 'no'; then
    AC_CHECK_HEADER(malloc/malloc.h,[mallochd='yes'],[mallochd='no'])
    systemincs=${systemincs}"-I/usr/include/malloc/ "
    if test $mallochd = 'no'; then
        AC_MSG_ERROR([malloc.h cannot be found  please do something about this.])
    fi
fi
 
AC_CHECK_FUNC([malloc],[malloc='yes'],[malloc='no'])
                if test ${malloc} = 'no'; then
                AC_MSG_ERROR([malloc cannot be found  please do something about this.])
                fi

AC_CHECK_FUNC(basename,[basename='yes'],[basename='no'])
AC_CHECK_FUNC(dirname,[dirname='yes'],[dirname='no'])

if test $basename = 'no' && test $dirname = 'no';then
echo "using alternate basename and dirname"
AC_SUBST(name,name)
noname=" -DNO_DIRNAME "
echo ${scriptotherlibsinc}
scriptotherlibsinc=${scriptotherlibsinc}" -L../util/name "
echo ${scriptotherlibsinc}
else
AC_SUBST(name,[])
noname=""
echo not using name.h
fi

#######################################################

#WARNING: The white space on shell-style quoted strings is not trivial.

AC_ARG_WITH([no-debugging],
            [AS_HELP_STRING([--with-no-debugging],[remove -g option when compiling libFrontier])],
            [debug=""],[debug="-g "])

AC_ARG_WITH([deptrack],
            [AS_HELP_STRING([--with-deptrack],[dependency tracking in examples for FronTier files])],
            [checklibs=FT_libs],[checklibs=""])

AC_SUBST(devel,${checklibs})



#NOTE: These options have terrible names, that don't line up with what they do.
basecode="full"
AC_ARG_WITH([adv],
            [AS_HELP_STRING([--with-adv],[Advanced FronTier++ code])],
            [basecode="adv"])

AC_ARG_WITH([base],
            [AS_HELP_STRING([--with-base],[Basic FronTier++ code])],
            [basecode="base"])

AC_ARG_WITH([gd],
            [AS_HELP_STRING([--with-gd=path],[use gd library for 1d post processing])],
            [scriptotherlibsinc="-L"${withval}"/lib "
             scriptotherlibs="-lgd "${scriptotherlibs}
             gdincs="-D__GD__ -I"${withval}"/include "],[gdincs=""])

AC_ARG_WITH([gsl],
            [AS_HELP_STRING([--with-gsl=path],[use gsl library for 1d post processing])],
            [scriptotherlibsinc="-L"${withval}"/lib "
             scriptotherlibs="-lgsl "${scriptotherlibs}
             gslincs="-D__GSL__ -I"${withval}"/include "],[gslincs=""])

AC_ARG_WITH([gpu],
            [AS_HELP_STRING([--with-gpu],[use GPU acceletated algorithms])],
            [gpuincs="-D__GPU__ "],[gpuincs=""])

if test "$gpuincs" = ""; then
    echo "GPU not enabled"
else
    echo "GPU enabled"
fi

AC_SUBST(USEGPU,${USEGPU})        


####################################################
# HDF

AC_ARG_WITH([hdf4],
            [AS_HELP_STRING([--with-hdf4=path],[hdf4 (default is no)])],
            [hdfprefix=${withval}])

AS_IF([test "x$with_hdf4" != xno],
      [hdfincs="-I"${hdfprefix}"/include -DUSE_HDF "
       scriptotherlibsinc=${scriptotherlibsinc}"-L/usr/lib -L"${hdfprefix}"/lib "
       scriptotherlibs=${scriptotherlibs}"-lmfhdf -ldf -ljpeg -lz "])

####################################################
# openmpi

#TODO: add with-openmpi-include and with-openmpi-library

#AC_ARG_WITH([openmpi],
 #           AS_HELP_STRING([--with-openmpi=path],[use openmpi])],
  #          [openmpiprefix=${withval} mpi="yes"])

#AS_IF([test "x$with_openmpi" != xno],
 #    [mpi_incs="-I"${openmpiprefix}"/include -D__MPI__"
  #    scriptotherlibsinc=${scriptotherlibsinc}" -L"${openmpiprefix}"/lib "
   #   scriptotherlibs=${scriptotherlibs}"-lmpi -lmpi_cxx -lpthread "])

####################################################
# mpich

#TODO: add with-mpich-include and with-mpich-library

AC_ARG_WITH([mpich],
            [AS_HELP_STRING([--with-mpich=path],[use mpich])],
            [mpichprefix=${withval} mpi="yes"],[])

AS_IF([test "x$with_mpich" != xno],
      [mpi_incs="-I"${mpichprefix}"/include -D__MPI__ "
       scriptotherlibsinc=${scriptotherlibsinc}"-L"${mpichprefix}"/lib "
       scriptotherlibs=${scriptotherlibs}"-lmpich -lmpichcxx -lpthread "])

####################################################
# PETSc

#Only providing single directory installation since FronTier
#is only compatable with PETSc version 3.1 which is no longer provided
#by most package managers.

AC_ARG_WITH([petsc],
            [AS_HELP_STRING([--with-petsc=path],[needed to build iFluid module])],
            [petsc_dir=${withval}],[petsc_dir=""])

AS_IF([test "x$with_petsc" != xno],
     [AS_IF([test ${mpi} = yes],
           [AC_SUBST(PETSC_INCLUDE,"-I"${petsc_dir}"/include ")
           AC_SUBST(PETSC_LIB,"-L"${petsc_dir}"/lib "${PETSC_LDFLAGS})],
           [AC_MSG_FAILURE([can't use PETSc without mpi])])
     ])

#This needs to be handled more robustly, but for now we will assume
#HYPRE was installed by PETSc, which is the case on our machines.
AC_ARG_WITH([hypre],
            AS_HELP_STRING([--with-hypre=path],[use hypre preconditioners with PETSc]),
            [hypreincs="-D__HYPRE__ "],[hypreincs=""])

AS_IF([test "x$with_hypre" != xno],[echo "HYPRE enabled"])


####################################################
# CGAL

AC_ARG_WITH([cgal],
            [AS_HELP_STRING([--with-cgal=path],[use single directory installed CGAL])],
            [cgal_dir=${withval}],[cgal_dir=""])

AC_ARG_WITH([cgal_include],
            AS_HELP_STRING([--with-cgal-include=path],[specify location of CGAL headers])],
            [cgal_inc=${withval}],[cgal_inc=""])

AC_ARG_WITH(cgal_library,
            [AS_HELP_STRING([--with-cgal-library=path],[specify location of CGAL libs])],
            [cgal_lib=${withval}],[cgal_lib=""])

if test "x$with_cgal" != xno; then
    echo "CGAL enabled"
    echo "cgal_dir = ${cgal_dir}"
    AC_SUBST(CGAL_DIR,${cgal_dir})
    AC_SUBST(CGAL_INCLUDE,"${cgal_dir}/include")
    AC_SUBST(CGAL_LIB,"${cgal_dir}/lib")
    else
        if test "x$with_cgal_include" != xno && test "x$with_cgal_library" != xno; then
            echo "CGAL enabled"
            echo "cgal_inc = ${cgal_inc}"
            echo "cgal_lib = ${cgal_lib}"
            AC_SUBST(CGAL_DIR,${cgal_dir})
            AC_SUBST(CGAL_INCLUDE,${cgal_inc})
            AC_SUBST(CGAL_LIB,${cgal_lib})
        else
            echo "CGAL is not enabled"
        fi
fi


####################################################
# Boost

AC_ARG_WITH(boost,[  --with-boost=path              Support for the BOOST directory],[boost_dir=${withval}],[boost_dir=""])
if test ${boost_dir} = ""; then
        echo "BOOST is not enabled"
else
        echo "BOOST enabled, boost_dir = ${boost_dir}"
        AC_SUBST(BOOST_DIR,${boost_dir})
        AC_SUBST(BOOST_INCLUDE,${boost_dir}"/include ")
        AC_SUBST(BOOST_LIB,${boost_dir}"/lib ")
fi


####################################################
# Cvode (SUNDIALS library)

AC_ARG_WITH([cvode], [ --with-cvode=path                 compile with cvode support], [cvode_dir=${withval} folding="cd folding && \$(MAKE)"], [cvode_dir="" folding=""])
if test ${cvode_dir} = ""; then
    echo "CVODE not enabled, cannot implement implicit solver"
    echo "Folding disabled"
else
    echo "CVODE enabled"
    AC_SUBST(CVODE_DIR,${cvode_dir})
    AC_SUBST(CVODE_INCLUDE,${cvode_dir}"/include ")
    AC_SUBST(CVODE_LIB,${cvode_dir}"/lib ")
fi

####################################################
# Operating System dependent configuration

#linux
if test $build_os = "linux" || test $build_os = "linux-gnu";then
    echo "using linux GCC optimizations when no debug option"
    if test "$debug" = ""; then
        optimization=" -O3 "
    else
        optimization=""
    fi
fi


##############################################################
# Setting C++ as the compiler
###############################################################
AC_SUBST(CC,${CC})
AC_SUBST(CXX,${CXX})
AC_SUBST(F77,${F77})
AC_SUBST(F77_LIBS,${F77_LIBS})
AC_SUBST(libext,"")
cppflagss=""


######################################################################
## The complete varibles as the are sent to their parts of the code.
#####################################################################

# included in all code in: driver 
#                          front 
#                          hyp 
#                          intfc 
#                          tri

scriptftlibincs=" -I.. -I../util "${debug}${optimization}${cppflagss}${hdfincs}${mpi_incs}${itapsd}${fenv}${noname}${hdfincs}${systemincs}${gdincs}${gpuincs}${gslincs}${hypreincs}


#included in all code in: example and testfront 
scriptexampleincs=${debug}${optimization}${mpi_incs}${hdfincs}${systemincs}${gpuincs}${gslincs}${hypreincs}

#included in all code in: util/util
scriptutilincs=" -I../.. -I../ "${debug}${optimization}${cppflagss}${mpi_incs}${hdfincs}${fenv}${noname}${hdfincs}${systemincs}${gpuincs}${gslincs}${hypreincs}

#include in all code in: util/dutil
scriptdutilincs=" -I../.. -I../ "${debug}${optimization}${cppflagss}${mpi_incs}${hdfincs}${fenv}${noname}${hdfincs}${systemincs}${gpuincs}${gslincs}${hypreincs}

flagsss=${FFLAGS}${debug}${optimization}${cppflagss}

headers=`find ./src/util/ ./src/intfc ./src/front -name \*\.*h | perl -e "s/\.\/src\/(.*)\n/include\/FronTier\/\1 /g" -p`

include_dirs="include/FronTier/intfc include/FronTier/front include/FronTier/util"

AC_SUBST(headers, $headers)
AC_SUBST(include_dirs, $include_dirs)

# Substitute these changes in the makefile.ins
AC_SUBST(sn_ellip,$snyesorno)
AC_SUBST(debugg,$debug)
AC_SUBST(tsttyesornoo,$tsttyesorno)
AC_SUBST(FFLAGS,$flagsss)
AC_SUBST(ftlibincs,$scriptftlibincs)
AC_SUBST(utilincs,$scriptutilincs)
AC_SUBST(exampleincs,$scriptexampleincs)
AC_SUBST(dutilincs,$scriptdutilincs)
AC_SUBST(CFLAGS,"")
AC_SUBST(otherlibsinc,$scriptotherlibsinc)
AC_SUBST(otherlibs,$scriptotherlibs)
AC_SUBST(randy,$RANLIB)
AC_SUBST_FILE(ftbuild)
ftbuild=./src/util/ftrules



# Write the makefiles
if test $basecode = "full"; then
filelist='Makefile 
src/Makefile 
src/util/Makefile 
src/front/Makefile 
src/intfc/Makefile 
Collision/Makefile
airfoil/Makefile
folding/Makefile
parachute/Makefile
parachute2d/Makefile
sprModel/Makefile
subsurf/Makefile
crystal/Makefile
example2d/Makefile
example3d/Makefile
iTaps/Makefile
finance/Makefile
melting/Makefile
curvature/Makefile
cell/Makefile
pde/Makefile
MonteCarlo/Makefile
timeSeries/Makefile
iFluid/Makefile
vFluid/Makefile
rgbody/Makefile
climate/Makefile
poisson/Makefile
cFluid/Makefile
wings/Makefile
cim/Makefile
DG/Makefile
inverse/Makefile
solver/Makefile
weno/Makefile
proj/Makefile
devel-deps.inc
iMesh-Defs.inc'
fi

if test $basecode = "adv"; then
filelist='Makefile 
src/Makefile
src/util/Makefile 
src/front/Makefile 
src/intfc/Makefile 
example2d/Makefile
example3d/Makefile
cFluid/Makefile
iFluid/Makefile
solver/Makefile
crystal/Makefile
devel-deps.inc
iMesh-Defs.inc'
fi

if test $basecode = "base"; then
filelist='Makefile 
src/Makefile 
src/util/Makefile 
src/front/Makefile 
src/intfc/Makefile 
example2d/Makefile
example3d/Makefile
cFluid/Makefile
devel-deps.inc
iMesh-Defs.inc'
fi


AC_CONFIG_FILES($filelist)

AC_OUTPUT
mkdir -p lib/
mkdir -p lib/${build_cpu}
mkdir -p src/front/${build_cpu}
mkdir -p src/intfc/${build_cpu}
mkdir -p src/util/${build_cpu}
mkdir -p src/util/d${build_cpu}

echo " Compiling with :"
echo   $ccft
echo   $F77
echo "The util flags are:" $scriptutilincs
echo
echo "The dutil flags are:" $scriptdutilincs
echo
echo "The front flags are:" $scriptftlibincs
echo 
echo "The example/testfront flags are:"$scriptexampleincs
echo ""


