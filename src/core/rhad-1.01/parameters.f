
C-{{{ subroutine parameters:

c..   ------------------------------------------------------------

      subroutine parameters(scms)
c..
c..   User-defined parameters.
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common.f'

c..   verbose mode:
      lverbose = .true.

c..   output unit for parameter list (6 = STDOUT)
      iunit = 6

c..   order or calculation:
      iord  = 4

c..   strong coupling constant at scale mz (5 active flavors):
      alphasmz  = 0.118d0

c..   use MS-bar or pole quark mass?  (.true. == MS-bar mass)
      lmsbar    = .false.

c..   masses
      massc      = 1.65d0        ! charm
      massb      = 4.75d0        ! bottom
      masst      = 175.d0        ! top

c..   renormalization scale:
      mu        = dsqrt(scms)

c..   decoupling scales:
      muc       = 2.d0*massc    ! charm
      mub       = massb         ! bottom
      mut       = masst         ! top

c..   threshold for open quark production
      thrc = 4.8d0
      thrb = 11.2d0
      thrt = 2*masst+10.d0

c..   lower bound of quark threshold region
      thrclow = 3.73d0
      thrblow = 10.52d0
      thrtlow = 2*masst-10.d0

c..   minimum allowed cms energy:
      sqmin     = 1.8d0

c..   some switches
      lmassless = .false.        ! use only massless approximation
      lqed      = .true.        ! QED corrections (.true. = ON)
      lpsup     = .true.       ! include power suppressed terms
      la3m2     = .true.        ! \alpha_s^3 * m^2 terms
      la3m4     = .true.        ! \alpha_s^3 * m^4 terms
      la3sing   = .true.       ! singlet contribution at order \alpha_s^3
      la4m2     = .true.        ! \alpha_s^4 * m^2 terms

      return
      end

C-}}}


