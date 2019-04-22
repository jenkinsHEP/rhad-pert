

      program example
c..
c..   This is the "typical program" as given in Sect.9 of
c..
c..   ***********
c..   "rhad: a program for the evaluation of the hadronic
c..      R-ratio in the perturbative regime of QCD"
c..
c..   by Robert V. Harlander and Matthias Steinhauser
c..   ***********
c..
c..   To compile, copy this file to the source directory of rhad, and say
c..
c..   gmake prog=example
c..   The executable will be named xexample.
c..
      implicit real*8(a-h,m-z)
      implicit integer(i,j)
      implicit character*60(k)
      implicit logical(l)

      include 'common.f'

      sqrts = 2.60d0
c..      sqrts = 3.07d0
c..      sqrts = 3.65d0
      scms  = sqrts*sqrts

      call parameters(scms)

c      alphasmz  = 0.130d0
      do i=1,10000
         alphasmz  = alphasmz - 0.000001d0
         call init(scms)
         rallnew = rhad(scms)
         if(dabs(rallnew-2.18).lt.0.00001d0) then
            goto 500
         endif
      enddo
 500  continue

      lverbose = .true.
      call init(scms)

      rall = rhad(scms)
      ru = ruqrk(scms)
      rd = rdqrk(scms)
      rs = rsqrk(scms)
      rc = rcqrk(scms)
      rb = rbqrk(scms)
      rt = rtqrk(scms)
      rsg = rsinglet(scms)
      rem = rqed(scms)

      print*,'R_had  = ',rall
      print*,'R_u    = ',ru
      print*,'R_d    = ',rd
      print*,'R_s    = ',rs
      print*,'R_c    = ',rc
      print*,'R_b    = ',rb
      print*,'R_t    = ',rt
      print*,'R_sing = ',rsg
      print*,'R_QED  = ',rem

      end
