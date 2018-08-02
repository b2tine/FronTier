# C++ library for robust collision handling
## Introduction
This is a library implenmented to handle fabric-fabric or fabric-rigid body collision.
The fabric surface is modeled with triangle mesh and its interior dynmaics is simulated with spring model (not included in this library).
For convenience, the data structure SURFACE, TRI, POINT are from FronTier library (a computational fluid dynamic library with interface tracking capability from Stony Brook University).

## Numerical Test
The basic idea of the implementation is originated from the paper "Robust Treatment of Collisions, Contact and Fricition for Cloth Animation", but we modify the method to couple it with FronTier easily. If finding some bugs, please first use the included test cases to verify your modification before commit. A few numerical examples are shown below:

<img style="float: left;" src="http://www.ams.sunysb.edu/~zgao/work/collision/img/fall-sphere.gif" width="230">
<img style="float: left;" src="http://www.ams.sunysb.edu/~zgao/work/collision/img/fall-body.gif" width="230">
<img style="float: left;" src="http://www.ams.sunysb.edu/~zgao/work/collision/img/fall-box.gif" width="230">
<img style="float: center;" src="http://www.ams.sunysb.edu/~zgao/work/collision/img/string_front.gif" width="230">
<img style="float: center;" src="http://www.ams.sunysb.edu/~zgao/work/collision/img/string_rear.gif" width="230">
<img style="float: center;" src="http://www.ams.sunysb.edu/~zgao/work/collision/img/fallstring.gif" width="230">

This project will not be maintained anymore. For more information, please contact via hit.zhenggao@gmail.com

