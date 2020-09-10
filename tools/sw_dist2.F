#include "fintrf.h"
C SW_DIST2 Distance between two lat,lon coordinates in km 
C This is a MEX file for matlab
C
C=================================================================== 
C SW_DIST $Id: sw_dist.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $ 
C Copyright (C) CSIRO, Phil Morgan & Steve Rintoul 1992. 
C 
C USAGE: [dist] = sw_dist(lat,lon) 
C 
C DESCRIPTION: 
C Calculate distance between two positions on glode using the "Plane 
C Sailing" method. Also uses simple geometry to calculate the bearing of 
C the path between position pairs. 
C 
C INPUT: 
C lat = decimal degrees (+ve N, -ve S) [- 90.. +90] 
C lon = decimal degrees (+ve E, -ve W) [-180..+180] 
C 
C OUTPUT: 
C dist = distance between positions in units 
C 
C AUTHOR: Phil Morgan and Steve Rintoul 92-02-10 
C Simplified by B. LE VU (2016) for AMEDA purpose 
C This software is provided "as is" without warranty of any kind. 
C See the file sw_copy.m for conditions of use and licence. 
C 
C REFERENCE: 
C The PLANE SAILING method as descriibed in "CELESTIAL NAVIGATION" 1989 by 
C Dr. P. Gormley. The Australian Antartic Division. 
C================================================================== 
C
C Modifications 
C 99-06-25. Lindsay Pender, Function name change from distance to sw_dist. 
C 99-06-25. Lindsay Pender, Fixed transpose of row vectors. 
C
C CALLER: general purpose 
C CALLEE: none 
C=================================================================== 
C-----------------------------------------------------------------------
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      integer nlhs, nrhs
      mwPointer plhs(*), prhs(*)

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

C     Pointers to input/output mxArrays:
      mwPointer x_ptr1, x_ptr2, y_ptr

C     Array information:
      mwPointer mrows1, ncols1, mrows2, ncols2
      mwSize size1, size2

C     Arguments for computational routine:
      real*8  lat(2), lon(2), dist

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:sw_dist2:nInput',
     +                           'Lat and Lon input required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:sw_dist2:nOutput',
     +                           'Too many output arguments.')
      endif

C     Validate inputs
C     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:sw_dist2:NonNumeric',
     +                           'Lat must be a number.')
      endif
      if(mxIsNumeric(prhs(2)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:sw_dist2:NonNumeric',
     +                           'Lon must be a number.')
      endif

C     Get the size of the input array.
      mrows1 = mxGetM(prhs(1))
      ncols1 = mxGetN(prhs(1))
      size1 = mrows1*ncols1
      mrows2 = mxGetM(prhs(2))
      ncols2 = mxGetN(prhs(2))
      size2 = mrows2*ncols2

C     Create Fortran array from the input argument.
      x_ptr1 = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(x_ptr1,lat,size1)
      x_ptr2 = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(x_ptr2,lon,size2)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(1,1,0)
      y_ptr = mxGetPr(plhs(1))

C     Call the computational subroutine.
      call sw_dist2(dist, lat, lon)

C     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(dist,y_ptr,1)     

      return
      end

C-----------------------------------------------------------------------
C     Computational routine
      subroutine sw_dist2(dist, lat, lon)

      implicit none 

      ! Arguments declarations 
      real*8 dist
      real*8 lat(2), lon(2)
      ! Variable declarations 
      integer DEG2NM
      real*8 DEG2RAD, dep
      real*8 dlat, dlon, NM2KM, pi
      real*8 latrad(2)

      !----------------- 
      ! DEFINE CONSTANTS 
      !----------------- 
      pi = 3.1415926535;
      DEG2RAD = 2*pi/360
      DEG2NM = 60
      NM2KM = 1.8520 ! Defined in Pond & Pickard p303.

      !----------------- 
      ! BEGIN 
      !----------------- 

      dlon = lon(2)-lon(1)
      latrad = abs(lat*DEG2RAD)
      dep = cos( (latrad(2)+latrad(1))/2 ) * dlon
      dlat = lat(2)-lat(1)
      dist = DEG2NM*sqrt(dlat*dlat + dep*dep) ! in n.miles

      dist = dist * NM2KM ! in km

      return
      end subroutine  sw_dist2


