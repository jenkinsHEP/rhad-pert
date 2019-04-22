      program tabr
c..
c..   This program generates Table 7 of the paper
c..
c..   ***********
c..   "rhad: a program for the evaluation of the hadronic
c..      R-ratio in the perturbative regime of QCD"
c..
c..   by Robert V. Harlander and Matthias Steinhauser
c..   ***********
c..
c..   if the default parameters given in that paper are used.
c..
c..   To compile, copy this file to the source directory of rhad, and say
c..
c..   gmake prog=table_scan
c..
c..   The executable will be named xtable_scan.
c..   Output-file:  table_scan.0.tex.
c..
c..
      implicit real*8(a-h,m-z)
      implicit integer(i,j)
      implicit character*60(k)
      implicit logical(l)

      include 'common.f'

      open(unit=12,file='table_scan.0.tex',status='unknown')

      write(12,102)

      do ii=0,4
         sqrts = 1.8d0 + ii*.4d0
         scms = sqrts*sqrts
         call prepout(scms)
      enddo

      sqrts = 3.73d0
      scms = sqrts*sqrts
      call prepout(scms)

      write(12,104)

c..   ----------------------------------------

      sqrts = 4.8d0
      scms = sqrts*sqrts
      call prepout(scms)

      do ii=0,5
         sqrts = 5.d0 + ii*1.d0
         scms = sqrts*sqrts
         call prepout(scms)
      enddo

      sqrts = 10.52d0
      scms = sqrts*sqrts
      call prepout(scms)

      write(12,104)

c..   ----------------------------------------

      sqrts = 11.2d0
      scms = sqrts*sqrts
      call prepout(scms)

      do ii=0,8
         sqrts = 12.d0 + ii*1.d0
         scms = sqrts*sqrts
         call prepout(scms)
      enddo

      write(12,103)

c..   ----------------------------------------

 102  format('\\documentclass{article}'/,
     &     '\\begin{document}'/,
     &     '\\begin{tabular}{|r|r|r|r|r|}'/,
     &     '\\hline'/,
     &     '$\\sqrt{s}$ (GeV) &',
     &     ' $ \\alpha_s(s) $ & ',
     &     ' $ R_c(s) $ & ',
     &     ' $ R_b(s) $ & ',
     &     ' $ R(s) $ \\\\'/,
     &     '\\hline\\hline')

 103  format('\\end{tabular}'/,
     &     '\\end{document}')

 104  format('\\hline')

      end

c..   ----------------------------------------

      subroutine prepout(scms)

      implicit  real*8(a-h,m-z)
      implicit  integer(i,j)
      implicit  character*60(k)
      implicit  logical(l)

      include 'common.f'

      call parameters(scms)
      call init(scms)
      col1 = dsqrt(scms)
      col2 = pi*api
      col3 = rcqrk(scms)
      col4 = rbqrk(scms)
      col5 = rhad(scms)

      write(12,101) col1,col2,col3,col4,col5
      print*, col1,col2,col3,col4,col5
         
 101  format(f8.2, ' & ', f8.4,' & ', f8.4,' & ', f8.4,
     &     ' & ',f8.4,' \\\\'/,' \\hline ')

      return
      end

