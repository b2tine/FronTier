CXX=mpic++ -std=c++11 -g -rdynamic -fopenmp -pedantic -Wno-long-long -Wno-unused-result
F77=mpif77 -g -rdynamic
F77_LIBS=

libext = 
incs =  # -I/usr/local/pkg/mpich2-3.2/include -D__MPI__  -D__HYPRE__  
libincs =  -L/usr/lib  -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2-3.2/lib 
#libincs =  -L/usr/lib  /usr/lib/libblas.so /usr/lib/x86_64-linux-gnu/libgfortran.so.3 -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2-3.2/lib 
libs = -lgd  -ldf -ljpeg -lz -lmpich -lpthread
#libs = -lgd  /usr/lib/libblas.so /usr/lib/x86_64-linux-gnu/libgfortran.so.3  -ldf -ljpeg -lz  -lmpich -lpthread

CGAL_INCLUDE = -I/usr/include
CGAL_LIB = -L/usr/lib/x86_64-linux-gnu -lCGAL_Core -lCGAL_ImageIO -lCGAL

GMOCK_INCLUDE = -I/usr/local/include -DGTEST_LINKED_AS_SHARED_LIBRARY=1
GMOCK_LIB = -L/usr/local/lib -lgmock -lgmock_main


OBJS = ../HyperSurfElement.o HSEtests.o ../BoundingVolume.o AABBtests.o  ../BVH_Node.o BVH_Nodetests.o ../BVH.o BVHtests.o
HEADERS = ../HyperSurfElement.h ../BoundingVolume.h ../BVH_Node.h ../BVH.h


testAll: main.cpp FT_libs_testdir $(OBJS)
	$(CXX) $< $(OBJS) $(GMOCK_LIB) $(CGAL_INCLUDE) $(CGAL_LIB) -I../../include $(incs) -L../../lib/x86_64 $(libincs) -lFronTier -lm $(libs) -o $@


%.o: %.cpp $(HEADERS)
	$(CXX) -c $< $(CGAL_INCLUDE) $(CGAL_LIB) -I../../include $(incs) -L../../lib/x86_64 $(libincs) -lFronTier -lm $(libs)

../%.o: ../%.cpp
	cd .. && $(CXX) -c $(notdir $<) $(CGAL_INCLUDE) $(CGAL_LIB) -I../include $(incs) -L../lib/x86_64 $(libincs) -lFronTier -lm $(libs)



clean:
	$(RM) testAll
	$(RM) *.o
	cd .. && $(MAKE) clean

tagsfile:
	ctags *.cpp ../../src/*/*.[chf]

-include ../../devel-deps.inc
