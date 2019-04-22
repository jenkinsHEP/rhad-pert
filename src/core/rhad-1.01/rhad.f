c..
c..   rhad.f
c..
c..   This is the main file for `rhad'.
c..
c..
C-{{{ header:

c..
c..   rhad.f  --  Version 1.01
c..   ------------------------------------------------------------
c..   Program: R(s)
c..   Authors: R. Harlander (CERN) and M. Steinhauser (U Hamburg)
c..   Date   : December 2002
c..   Changes: [March 30, 2009] v1.01
c..            - exact results at alphas^4 included
c..   ------------------------------------------------------------
c..   for documentation, see
c..
c..   ***********
c..   "rhad: a program for the evaluation of the hadronic
c..      R-ratio in the perturbative regime of QCD"
c..
c..   by Robert V. Harlander and Matthias Steinhauser,
c..   Comp. Phys. Comm. 153 (2003) 244,
c..   CERN-TH/2002-379, DESY 02-223, hep-ph/0212294.
c..   ***********
c..

C-}}}
C-{{{ function rhad(scms)

c..   ------------------------------------------------------------

      function rhad(scms)
c..
c..   scms:  center-of-mass energy squared
c..
c..   Returns the hadronic R-ratio, including all factors.
c..
c..   Note: Only photon contribution is included at the moment!
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      rhad = 0.d0

c..   non-singlet:
      rhad = rhad + ruqrk(scms) + rdqrk(scms) + rsqrk(scms)
     &     + rcqrk(scms) + rbqrk(scms) + rtqrk(scms)

c..   singlet:
      rhad = rhad + rsinglet(scms)

c..   QED:
      rhad = rhad + rqed(scms)

      return
      end

C-}}}
C-{{{ function rsinglet:

      function rsinglet(scms)
c..
c..   singlet terms
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if (.not.la3sing) then
         rsinglet = 0.d0
         return
      endif

      rtmp = 0.d0

      rtmp = rtmp
     &     + (vup + vdown + vstrange)**2 * rv3sing(scms,0.d0,0.d0)

      if (lcharm) then
         rtmp = rtmp
     &        + 2*(vup + vdown + vstrange)*vcharm * rv3sing(scms,mc,0
     &        .d0)+ vcharm**2 * rv3sing(scms,mc,mc)
      endif

      if (lbottom) then
         rtmp = rtmp
     &        + 2*(vup + vdown + vstrange)*vbottom * rv3sing(scms,mb
     &        ,0.d0)+ 2*vcharm*vbottom * rv3sing(scms,mb,0.d0)+
     &        vbottom**2 * rv3sing(scms,mb,mb)
      endif

      if (ltop) then
         rtmp = rtmp
     &        + 2*(vup + vdown + vstrange)*vtop * rv3sing(scms,mt,0
     &        .d0)+ 2*vcharm*vtop * rv3sing(scms,mt,0.d0)+ 2*vbottom
     &        *vtop * rv3sing(scms,mt,0.d0)+ vtop**2 * rv3sing(scms
     &        ,mt,mt)
      endif

      rsinglet = api**3 * nc * rtmp

      return
      end

C-}}}
C-{{{ function rqed:

      function rqed(scms)
c..
c..   QED corrections
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if (.not.lqed) then
         rqed = 0.d0
         return
      endif

      charges = vup**4 + vdown**4 + vstrange**4

      if (lcharm) then
         charges = charges + vcharm**4
      endif
      if (lbottom) then
         charges = charges + vbottom**4
      endif
      if (ltop) then
         charges = charges + vtop**4
      endif

      rqed = nc * alpha/pi * charges * r01qed(scms)

      if (iord.gt.1) then
         rqed = rqed*( 1 - 1/3.d0*api )
      endif

      return
      end

C-}}}
C-{{{ function ruqrk:

      function ruqrk(scms)
c..
c..   massless non-singlet contributions
c..   without power-suppressed terms
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      mq = 0.d0
      nl = 3.d0

c..   Born:
      rv0l = rv0(scms,mq)
      rv = rv0l

      if (iord.ge.1) then
c..   1-loop:
         rv1l = r01(scms)

         rv = rv + api*rv1l
      endif


      if (iord.ge.2) then
c..   2-loop:
         rv2l =
     &        + r0cf(scms)      ! cf**2
     &        + r0ca(scms)      ! ca*cf
     &        + nl * r0nl(scms) ! nl
     &        + r0db(scms,mc,lcharm)
     &        + r0db(scms,mb,lbottom)
     &        + r0db(scms,mt,ltop)

         rv = rv + api**2*rv2l
      endif

      if (iord.ge.3) then
c..   3-loop:

         if (.not.lcharm) then
            rv3l = rv3m0(scms)
         elseif (.not.lbottom) then
            rv3l = rv3m0(scms) + delr03(scms,mc)
         elseif (.not.ltop) then
            rv3l = rv3m0(scms) + delr03m2(scms,mc) + delr03(scms,mb)
         else
            rv3l = rv3m0(scms) + delr03(scms,mt)
         endif

         rv = rv + api**3*rv3l
      endif

      if (iord.ge.4) then
c..   4-loop:
         rv4l = rv4m0(scms)

         rv = rv + api**4*rv4l
      endif

      ruqrk = nc * vup**2 * rv

      return
      end

C-}}}
C-{{{ function rdqrk:

      function rdqrk(scms)
c..
c..   massless non-singlet contributions
c..   without power-suppressed terms
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      mq = 0.d0
      nl = 3.d0

c..   Born:
      rv0l = rv0(scms,mq)
      rv = rv0l

      if (iord.ge.1) then
c..   1-loop:
         rv1l = r01(scms)

         rv = rv + api*rv1l
      endif


      if (iord.ge.2) then
c..   2-loop:
         rv2l =
     &        + r0cf(scms)      ! cf**2
     &        + r0ca(scms)      ! ca*cf
     &        + nl * r0nl(scms) ! nl
     &        + r0db(scms,mc,lcharm)   ! c inside u,d,s
     &        + r0db(scms,mb,lbottom)   ! b inside u,d,s
     &        + r0db(scms,mt,ltop)   ! t inside u,d,s

         rv = rv + api**2*rv2l
      endif

      if (iord.ge.3) then
c..   3-loop:

         if (.not.lcharm) then
            rv3l = rv3m0(scms)
         elseif (.not.lbottom) then
            rv3l = rv3m0(scms) + delr03(scms,mc)
         elseif (.not.ltop) then
            rv3l = rv3m0(scms) + delr03m2(scms,mc) + delr03(scms,mb)
         else
            rv3l = rv3m0(scms) + delr03(scms,mt)
         endif

         rv = rv + api**3*rv3l
      endif

      if (iord.ge.4) then
c..   4-loop:
         rv4l = rv4m0(scms)

         rv = rv + api**4*rv4l
      endif

      rdqrk = nc * vdown**2 * rv

      return
      end

C-}}}
C-{{{ function rsqrk:

      function rsqrk(scms)
c..
c..   massless non-singlet contributions
c..   without power-suppressed terms
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      mq = 0.d0
      nl = 3.d0

c..   Born:
      rv0l = rv0(scms,mq)
      rv = rv0l

      if (iord.ge.1) then
c..   1-loop:
         rv1l = r01(scms)

         rv = rv + api*rv1l
      endif


      if (iord.ge.2) then
c..   2-loop:
         rv2l =
     &        + r0cf(scms)      ! cf**2
     &        + r0ca(scms)      ! ca*cf
     &        + nl * r0nl(scms) ! nl
     &        + r0db(scms,mc,lcharm) ! c inside u,d,s
     &        + r0db(scms,mb,lbottom) ! b inside u,d,s
     &        + r0db(scms,mt,ltop) ! t inside u,d,s

         rv = rv + api**2*rv2l
      endif

      if (iord.ge.3) then
c..   3-loop:

         if (.not.lcharm) then
            rv3l = rv3m0(scms)
         elseif (.not.lbottom) then
            rv3l = rv3m0(scms) + delr03(scms,mc)
         elseif (.not.ltop) then
            rv3l = rv3m0(scms) + delr03m2(scms,mc) + delr03(scms,mb)
         else
            rv3l = rv3m0(scms) + delr03(scms,mt)
         endif

         rv = rv + api**3*rv3l
      endif

      if (iord.ge.4) then
c..   4-loop:
         rv4l = rv4m0(scms)

         rv = rv + api**4*rv4l
      endif

      rsqrk = nc * vstrange**2 * rv

      return
      end

C-}}}
C-{{{ function rcqrk:

      function rcqrk(scms)
c..
c..   non-singlet contributions for c-quark production
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if (.not.lcharm) then
         rcqrk = 0.d0
         return
      endif

c..   number of massless flavors:
      nl = 3.d0

c..   Born:
      rv0l = rv0(scms,mc)
      rv = rv0l


      if (iord.ge.1) then
c..   1-loop:
         rv1l = rv1(scms,mc)
         if (lmsbar) then
            rv1l = rv1l + rv1ctms(scms,mc)
         endif

         rv = rv + api*rv1l
      endif


      if (iord.ge.2) then
c..   2-loop:
         rv2l =
     &        + rvcf(scms,mc)   ! cf**2
     &        + rvca(scms,mc)   ! ca*cf
     &        + nl*rvnl(scms,mc) ! cf*tr*nl
     &        + rvct(scms,mc)    ! double bubble (c inside c)
     &        + rvdb(scms,mc,mb,lbottom) ! double bubble (b inside c)
     &        + rvdb(scms,mc,mt,ltop) ! double bubble (t inside c)
         if (lmsbar) then
            rv2l = rv2l + rv2ctms(scms,mc)
         endif

         rv = rv + api**2*rv2l
      endif


      if (iord.ge.3) then
c..   3-loop:

         if (.not.lbottom) then
            rv3l = rv3(scms,mc)
         elseif (.not.ltop) then
            rv3l = rv3m0(scms) + rv3m2(scms,mc) + delr03(scms,mb)
         else
            rv3l = rv3m0(scms) + delr03(scms,mt)
         endif

         rv = rv + api**3*rv3l
      endif


      if (iord.ge.4) then
c..   4-loop:
         rv4l = rv4(scms,mc)

         rv = rv + api**4*rv4l
      endif

      rcqrk = nc * vcharm**2 * rv

      return
      end

C-}}}
C-{{{ function rbqrk:

      function rbqrk(scms)
c..
c..   non-singlet contributions for b-quark production
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if (.not.lbottom) then
         rbqrk = 0.d0
         return
      endif

c..   number of massless flavors:
      nl = 3.d0

c..   Born:
      rv0l = rv0(scms,mb)
      rv = rv0l


      if (iord.ge.1) then
c..   1-loop:
         rv1l = rv1(scms,mb)
         if (lmsbar) then
            rv1l = rv1l + rv1ctms(scms,mb)
         endif

         rv = rv + api*rv1l
      endif


      if (iord.ge.2) then
c..   2-loop:
         rv2l =
     &        + rvcf(scms,mb)   ! cf**2
     &        + rvca(scms,mb)   ! ca*cf
     &        + nl*rvnl(scms,mb) ! cf*tr*nl
     &        + rvct(scms,mb)    ! double bubble (b inside b)
     &        + rvdb(scms,mb,mc,lcharm) ! double bubble (c inside b)
     &        + rvdb(scms,mb,mt,ltop) ! double bubble (t inside b)
         if (lmsbar) then
            rv2l = rv2l + rv2ctms(scms,mb)
         endif

         rv = rv + api**2*rv2l
      endif


      if (iord.ge.3) then
c..   3-loop:

         if (.not.ltop) then
            rv3l = rv3(scms,mb) + delr03m2(scms,mc)
         else
            rv3l = rv3m0(scms) + delr03(scms,mt)
         endif

         rv = rv + api**3*rv3l
      endif


      if (iord.ge.4) then
c..   4-loop:
         rv4l = rv4(scms,mb)

         rv = rv + api**4*rv4l
      endif

      rbqrk = nc * vbottom**2 * rv

      return
      end

C-}}}
C-{{{ function rtqrk:

      function rtqrk(scms)
c..
c..   non-singlet contributions for b-quark production
c..
c..   scms:  center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if (.not.ltop) then
         rtqrk = 0.d0
         return
      endif

c..   number of massless flavors:
      nl = 3.d0

c..   Born:
      rv0l = rv0(scms,mt)
      rv = rv0l


      if (iord.ge.1) then
c..   1-loop:
         rv1l = rv1(scms,mt)
         if (lmsbar) then
            rv1l = rv1l + rv1ctms(scms,mt)
         endif

         rv = rv + api*rv1l
      endif


      if (iord.ge.2) then
c..   2-loop:
         rv2l =
     &        + rvcf(scms,mt)   ! cf**2
     &        + rvca(scms,mt)   ! ca*cf
     &        + nl*rvnl(scms,mt) ! cf*tr*nl
     &        + rvct(scms,mt)    ! double bubble (t inside t)
     &        + rvdb(scms,mt,mc,lcharm) ! double bubble (c inside t)
     &        + rvdb(scms,mt,mb,lbottom)   ! double bubble (b inside t)
         if (lmsbar) then
            rv2l = rv2l + rv2ctms(scms,mt)
         endif

         rv = rv + api**2*rv2l
      endif


      if (iord.ge.3) then
c..   3-loop:
         rv3l = rv3(scms,mt)

         rv = rv + api**3*rv3l
      endif


      if (iord.ge.4) then
c..   4-loop:
         rv4l = rv4(scms,mt)

         rv = rv + api**4*rv4l
      endif

      rtqrk = nc * vtop**2 * rv

      return
      end

C-}}}
C-{{{ subroutine fixpar:

c..   ------------------------------------------------------------

      subroutine fixpar(scms)
c..
c..   Here we define parameters that should not be changed by the user.
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

c..   program version:
      kversion = '1.01'

c..   numerical constants:
      dln2      = dlog(2.d0)
      pi	= 4.d0*datan(1.d0) ! Pi
      zeta2     = pi*pi/6.d0
      zeta3     = 1.2020569031595942854d0   ! \zeta(3)
      zeta4     = 2*zeta2*zeta2/5.d0
      zeta5     = 1.0369277551433699263d0   ! \zeta(5)
      a4        = 0.51747906167389938633d0  !  Li_4(1/2)
      b4        = -1.7628000870737708641d0

c..   initial value of alpha_s is defined at 'infini' flavors
c..   and scale mz:
      infini    = 5
      mz        = 91.1876d0

c..   number of colors:
      nc        = 3.d0

c..   electromagnetic interaction
      alpha	= 1.d0/137.0359895d0 ! electromagnetic coupling
      qel       = dsqrt(4.d0*pi*alpha) ! electron charge

      vdown     =-1.d0/3.d0     ! d-quark charge
      vup       = 2.d0/3.d0     ! u-quark charge
      vstrange  =-1.d0/3.d0     ! s-quark charge
      vcharm    = 2.d0/3.d0     ! c-quark charge
      vbottom   =-1.d0/3.d0     ! b-quark charge
      vtop      = 2.d0/3.d0     ! t-quark charge

c..   color factors:
      cf = 4.d0/3.d0
      ca = 3.d0
      tr = 1.d0/2.d0

c..   buffers for if's:
      xgtbuf = 1.d0 + 1.d-12
      xltbuf = 1.d0 - 1.d-12

c..   use expansion for rvct? (faster than numerical integration)
      lrvctexp = .true.

      return
      end

C-}}}
C-{{{ subroutine init(scms)

      subroutine init(scms)
c..
c..   Prepare everything for the actual calculation.
c..
c..   scms = center-of-mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

c..   load non-user parameters:
      call fixpar(scms)

c..   run some checks on the input:
      call checkin(scms,lpass)
      if (.not.lpass) stop

      lcharm = .false.
      lbottom = .false.
      ltop = .false.
      if (scms.ge.thrc**2 * xltbuf) then
         lcharm = .true.
         if (scms.ge.thrb**2 * xltbuf) then
            lbottom = .true.
            if (scms.ge.thrt**2 * xltbuf) then
               ltop = .true.
            endif
         endif
      endif

c..   determine inffin from scms and thr{c,b,t}
      inffin = 3
      if(ltop)then
         inffin = 6
      elseif(lbottom)then
         inffin = 5
      elseif(lcharm)then
         inffin = 4
      endif

c..   compute alpha_s(mu)
      call rundecalpha(alphasmz,mz,mu,als)
      api = als/pi

      if (lmsbar) then
c..   compute mq(mu) from mq(mq):
         if (lcharm) then
            call rundecmass(massc,4,massc,mu,mc)
         else
            call mms2mos(massc,4,iord,mc)
         endif
         if (lbottom) then
            call rundecmass(massb,5,massb,mu,mb)
         else
            call mms2mos(massb,5,iord,mb)
         endif
         if (ltop) then
            call rundecmass(masst,6,masst,mu,mt)
         else
            call mms2mos(masst,6,iord,mt)
         endif
      else
c..   input: pole quark masses
         mc = massc
         mb = massb
         mt = masst
      endif

c..   fix some inconsistencies:
      if (la3m4) la3m2 = .true.  !  no m^4-terms without m^2-terms
      if (lmassless) then
         la3m2 = .false.
         la3m4 = .false.
         la4m2 = .false.
      endif

c..   print parameters:
      if (lverbose) call printpar(scms)

c..   check if energy region makes sense
      if (((scms.gt.thrclow**2 * xgtbuf).and.
     &     (.not.lcharm)).or.
     &     ((scms.gt.thrblow**2 * xgtbuf).and.
     &     (.not.lbottom)).or.
     &     ((scms.gt.thrtlow**2 * xgtbuf).and.
     &     (.not.ltop))) then
      write(6,2003) dsqrt(scms)
 2003 format(
     &     'WARNING: sqrt(s)=',1f8.3,' GeV not in region where ',/
     &     '         perturbation theory can be applied.'
     &     )
      endif

      end

C-}}}
C-{{{ subroutine checkin(scms,lpass)

      subroutine checkin(scms,lpass)
c..
c..   Perform some sanity checks on the input.
c..   Common block must be filled before calling.
c..
c..   scms = center-of-mass energy squared
c..
c..   Returns
c..   lpass = .true.   if checks are passed,
c..   lpass = .false.  otherwise
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      lpass = .true.

      if ((iord.gt.4).or.(iord.lt.0)) then
         write(6,*) 'iord must be between 0 and 4.'
         write(6,*) 'iord = ',iord
         lpass = .false.
      endif

      if (scms.lt.sqmin*sqmin)then
         write(6,*) '<function checkin>: Energy value too low.'
         write(6,*) '                    Is:',dsqrt(scms)
         write(6,*) '                    Must be > ',sqmin
         lpass = .false.
      endif

      if ( (thrc.lt.2*massc) ) then
         write(6,*) '<function checkin>: thrc < 2*massc:'
         write(6,*) '                    thrc = ',thrc
         write(6,*) '                 2*massc = ',2*massc
         lpass = .false.
      endif

      if ( (thrb.lt.2*massb) ) then
         write(6,*) '<function checkin>: thrb < 2*massb:'
         write(6,*) '                    thrb = ',thrb
         write(6,*) '                 2*massb = ',2*massb
         lpass = .false.
      endif

      if ( (thrt.lt.2*masst) ) then
         write(6,*) '<function checkin>: thrt < 2*masst:'
         write(6,*) '                    thrt = ',thrt
         write(6,*) '                 2*masst = ',2*masst
         lpass = .false.
      endif

      end

C-}}}
C-{{{ subroutine printpar(scms):

      subroutine printpar(scms)
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      write(iunit,3001) kversion
 3001 format(/,
     &     '            rhad.f -- Version ',A4,/
     &     ' by Robert Harlander and Matthias Steinhauser',/
     &     '                December 2002',/)


      write(iunit,3005) iord
 3005 format(
     &     'Order of calculation: ',I3,/
     &     )


      write(iunit,3004) dsqrt(scms),mu,thrc,thrb,thrt
 3004 format(
     &     'Scales:'/
     &     '  sqrt(s)  =  ',1f8.3,' GeV'/
     &     '  mu       =  ',1f8.3,' GeV'/
     &     '  thrc     =  ',1f8.3,' GeV'/
     &     '  thrb     =  ',1f8.3,' GeV'/
     &     '  thrt     =  ',1f8.3,' GeV'/
     &     )


      write(iunit,3014) inffin
 3014 format(
     &     'Number of active flavors:'/
     &     '  nf  =  ',I1/
     &     )

      if (infini-inffin.eq.1) then
         write(iunit,3020) 'b',mub,5,4
      elseif (infini-inffin.eq.2) then
         write(iunit,3020) 'b',mub,5,4
         write(iunit,3020) 'c',muc,4,3
      elseif (infini-inffin.eq.-1) then
         write(iunit,3020) 't',mut,5,6
      endif
 3020 format('Decoupling at'/'  mu',a1,' =',1f8.3,' GeV  (',i1,' -> ',i1
     &     ,')'/)



      write(iunit,3003) 1/alpha,alphasmz,mz,api*pi
 3003 format(
     &     'Coupling constants:',/,
     &     '  alpha_QED    =  1 / (',1f8.3,' )',/
     &     '  alpha_s(Mz)  =  ',1f6.4,
     &     '             [ Mz = ',1f8.4,' GeV ]'/
     &     '  alpha_s(mu)  =  ',1f12.10,/
     &     )


      if (lmsbar) then
         write(iunit,3013)
         if (lcharm) then
            write(iunit,3012) 'c','c',massc,'c',mc
         else
            write(iunit,3002) 'c',mc
         endif
         if (lbottom) then
            write(iunit,3012) 'b','b',massb,'b',mb
         else
            write(iunit,3002) 'b',mb
         endif
         if (ltop) then
            write(iunit,3012) 't','t',masst,'t',mt
         else
            write(iunit,3002) 't',mt
         endif
      else
         write(iunit,3013)
         write(iunit,3002) 'c',mc
         write(iunit,3002) 'b',mb
         write(iunit,3002) 't',mt
      endif

 3013 format('Quark masses:')
 3012 format( '  m_',A1,'(m_',A1,')  =  ',1f6.2,' GeV  --> ',
     &        '  m_',A1,'(mu)   =  ',1f6.2,' GeV')

 3002 format( '  M_',A1,'       =  ',1f6.2,' GeV')

      write(iunit,3006) lmassless,
     &     lpsup,lqed,la3sing,la3m2,la3m4,la4m2
 3006 format(/
     &     'General switches (F=False, T=True):',/
     &     '  only massless terms    :',L2,/
     &     '  power suppressed terms :',L2,/
     &     '  QED corrections        :',L2,/
     &     '  singlet contributions  :',L2,/
     &     '  alphas^3 m^2 included  :',L2,/
     &     '  alphas^3 m^4 included  :',L2,/
     &     '  alphas^4 m^2 included  :',L2
     &     )

C      write(iunit,3007) thrc,thrclow,thrb,thrblow,thrt,thrtlow
C 3007 format(
C     &     'Limits:',/,
C     &     '  thrc     =  ',1f8.3,' GeV',/
C     &     '  thrclow  =  ',1f8.3,' GeV',/
C     &     '  thrb     =  ',1f8.3,' GeV',/
C     &     '  thrblow  =  ',1f8.3,' GeV',/
C     &     '  thrt     =  ',1f8.3,' GeV',/
C     &     '  thrtlow  =  ',1f8.3,' GeV',/
C     &     )

      write(iunit,*) '-- end of parameters --'

      end

C-}}}
