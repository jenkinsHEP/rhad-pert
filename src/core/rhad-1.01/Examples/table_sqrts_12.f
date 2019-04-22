
      program table
c..
c..   This program generates Table 6 of the paper
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
c..   gmake prog=table_sqrts_12
c..
c..   The executable will be named xtable_sqrts_12.
c..   Output-file:  table_sqrts_12.0.tex.
c..
      implicit real*8(a-h,m-z)
      implicit integer(i,j)
      implicit character*60(k)
      implicit logical(l)

      include 'common.f'

      sqrts = 12d0
      scms = sqrts*sqrts
      call parameters(scms)


      open(unit=12,file='table_sqrts_12.0.tex',status='unknown')

      write(12,102)

      do ii=0,4
         iord = ii
         call init(scms)
         colpi = pi*api
         colc = rcqrk(scms)
         colb = rbqrk(scms)
         colrhad = rhad(scms)
         write(12,101) iord,colpi,colc,colb,colrhad
      enddo

      write(12,103) 

 101  format(I2, ' & ', f8.4,' & ',f8.4,' & ',f8.4,' & ',f8.4,
     &     ' \\\\'/,' \\hline ')

 102  format('\\documentclass{article}'/,
     &     '\\begin{document}'/,
     &     '\\begin{tabular}{|c||c|c|c|c|}'/,
     &     '\\multicolumn{5}{c}{$\\sqrt{s} = 12$\\,GeV}\\\\'/,
     &     '\\hline'/,
     &     ' order &',
     &     ' $ \\alpha_s $ & ',
     &     ' $ R_c(s) $ &'/,
     &     ' $ R_b(s) $ &'/,
     &     ' $ R(s) $ \\\\'/,
     &     '\\hline\\hline')

 103  format('\\end{tabular}'/,
     &     '\\end{document}')

      end
