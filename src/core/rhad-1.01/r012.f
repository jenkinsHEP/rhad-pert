c..
c..   r012.f
c..   
c..   This file is part of `rhad'.
c..
c..   routines for Born, 1-loop, 2-loop contributions
c..
C-{{{ function r01qed:

c..   ------------------------------------------------------------

      function r01qed(scms)
c..
c..   massless one-loop QED
c..
c..   scms = center of mass energy squared
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      r01qed = 3.d0/4.d0

      return
      end

C-}}}
C-{{{ function r01:

c..   ------------------------------------------------------------

      function r01(scms)
c..
c..   massless one-loop
c..
c..   scms = center of mass energy squared
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      r01 = cf*3.d0/4.d0

      return
      end

C-}}}
C-{{{ function r0cf:

c..   ------------------------------------------------------------

      function r0cf(scms)
c..
c..   massless two-loop
c..   CF**2 terms
c..
c..   scms = center of mass energy squared
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      r0cf = cf**2 * ( -3.d0/32.d0 )

      return
      end

C-}}}
C-{{{ function r0ca:

c..   ------------------------------------------------------------

      function r0ca(scms)
c..
c..   massless two-loop
c..   CA*CF terms
c..
c..   scms = center of mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlsmu = dlog(scms/mu/mu)

      r0ca = cf*ca * ( 123.d0/32.d0 - 11/4.d0*zeta3 - 11/16.d0*dlsmu )

      return
      end

C-}}}
C-{{{ function r0nl:

c..   ------------------------------------------------------------

      function r0nl(scms)
c..
c..   massless two-loop
c..   double bubble with massless inner quarks
c..
c..   scms = center of mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlsmu = dlog(scms/mu/mu)

      r0nl = cf*tr* ( -11.d0/8.d0 + zeta3 + 1/4.d0*dlsmu )

      return
      end

C-}}}
C-{{{ function r0f:

c..   ------------------------------------------------------------

      function r0f(scms)
c..
c..   massless two-loop
c..   double bubble with inner quark = outer quark
c..
c..   scms = center of mass energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlsmu = dlog(scms/mu/mu)

      r0f = cf*tr * ( - 11.d0/8.d0 + zeta3 + dlsmu/4.d0 )

      return
      end

C-}}}
C-{{{ function rv0:

c..   ------------------------------------------------------------

      function rv0(scms,mq)
c..
c..   tree-level
c..   vector case
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common.f'

      if (lmassless) then
         velo = 1.d0
      else
         velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )
      endif

      rv0 = velo/2.d0*(3.d0-velo**2)

      return
      end

C-}}}
C-{{{ function rv1:

c..   ------------------------------------------------------------

      function rv1(scms,mq)
c..
c..   one-loop
c..   vector case
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      complex*16 spence,csp1,csp2
      real*8 arg
      external spence
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         arg = (1.d0 - velo)/(1.d0 + velo)
         csp1 = spence(arg)
         csp2 = spence(arg**2)
         sp1 = dreal(csp1)
         sp2 = dreal(csp2)

         rtmp = (3.d0 - velo**2)/2.d0*( (1.d0 + velo**2)*(2.d0*sp1 + sp2
     &        +dlog((1.d0 - velo)/(1.d0 + velo))* dlog(8.d0*velo**2/(1
     &        .d0+velo)**3) ) - 2.d0*velo*dlog(8.d0*velo**2/(1.d0+velo)
     &        **3) +( 3.d0*velo - (33.d0 + 22.d0*velo**2 - 7.d0*velo**4)
     &        /8.d0/(3.d0 - velo**2) )*dlog((1.d0-velo)/(1.d0+velo)) +
     &        (15.d0*velo- 9.d0*velo**3)/(4.d0*(3.d0-velo**2)))
      else
         rtmp = 3.d0/4.d0
      endif

      rv1 = cf*rtmp

      return
      end

C-}}}
C-{{{ function rv1ctms:

c..   ------------------------------------------------------------

      function rv1ctms(scms,mq)
c..
c..   one-loop counterterm contribution; OS->MS-bar
c..   vector case
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..
c..   Note, 'rtmp' contains already the colour factor cf 
c..   but not the overall factor nc
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         rtmp = (-0.5*(-1. + velo**2)**2*(4. + 
     -        3.*dlog(mu**2/mq**2)))/velo
      else
         rtmp = 0.d0
      endif

      rv1ctms = rtmp

      return
      end

C-}}}
C-{{{ function rvcf:

c..   ------------------------------------------------------------

      function rvcf(scms,mq)
c..
c..   two-loop
c..   vector case
c..   CF**2 terms
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         z  = 1.d0/(1.d0-velo**2)
         rtmp =  3.d0/4.d0/pi * dimag(-30.92186365499425d0 + 13
     &        .15947253478581d0*dlog(2.d0) + ((2659.144746799566d0 + 994
     &        .9062634178374d0* (-1.d0 + 2.d0/z + (0.d0,2.d0)*dsqrt((-1
     &        .d0 + z)/z**2)) - 2742.9905230464883d0* (-1.d0 + 2.d0/z +
     &        (0.d0,2.d0) *dsqrt((-1.d0 + z)/z**2))**2 - 779
     &        .4006777648011d0 * (-1.d0 + 2.d0 /z + (0.d0,2.d0) *dsqrt((
     &        -1.d0 + z)/z**2))**3 + 562.9840842893414d0* (-1.d0 + 2.d0
     &        /z + (0.d0,2.d0)*dsqrt((-1.d0 + z)/z **2))**4 + 132
     &        .8451536724444d0 * (-1.d0 + 2.d0/z + (0.d0,2.d0) *dsqrt((
     &        -1.d0 + z)/z**2))**5)* (2.d0/z + (0.d0,2.d0)*dsqrt((-1.d0
     &        + z)/z **2))**2)/ ((121.97680950645643d0 + 11
     &        .55904929184703d0* (-1.d0 + 2.d0/z + (0.d0,2.d0)*dsqrt((-1
     &        .d0 + z) /z**2)) - 119 .51160360383739d0* (-1.d0 + 2.d0/z
     &        + (0.d0,2.d0) *dsqrt((-1.d0 + z)/z **2))**2 - 1
     &        .0411385083877673d0* (-1.d0 + 2.d0/z + (0.d0,2.d0)
     &        *dsqrt((-1.d0 + z)/z**2))**3 + 18.465346802086135d0* (-1
     &        .d0 + 2.d0/z + (0.d0,2.d0) *dsqrt((-1.d0 + z)/z **2))**4 +
     &        1.d0*(-1.d0 + 2.d0/z + (0.d0,2.d0) *dsqrt((-1.d0 + z)/z**2
     &        )) **5) * (2.d0 - 2.d0/z - (0.d0,2.d0)*dsqrt((-1.d0 + z)/z
     &        **2))) - (0 .020833333333333332d0*(-1374.d0 - 346.d0*z -
     &        108.d0 *((0.d0,3.141592653589793d0) + dlog((1.d0 - 1.d0
     &        *dsqrt(1.d0 - 1.d0/z))/(1.d0 + dsqrt(1.d0 - 1.d0/z))))**2
     &        - 687.d0 *((0.d0,3.141592653589793d0) + dlog((1.d0 - 1.d0
     &        *dsqrt(1.d0 - 1.d0/z))/(1.d0 + dsqrt(1.d0 - 1.d0/z))))
     &        * dsqrt((-1.d0 + z)/z) - 186.d0 *z*((0.d0,3
     &        .141592653589793d0) + dlog((1.d0 - 1.d0*dsqrt(1.d0 - 1.d0
     &        /z))/(1.d0 + dsqrt(1.d0 - 1.d0/z)))) * dsqrt((-1.d0 + z)/z
     &        )))/z)
         rtmp = rtmp - 4.d0*1/cf*rv1(scms,mq)
      else
         rtmp = -3.d0/32.d0
      endif

      rvcf = cf**2 * rtmp

      return
      end

C-}}}
C-{{{ function rvca:

c..   ------------------------------------------------------------

      function rvca(scms,mq)
c..
c..   two-loop
c..   vector case
c..   CA*CF terms
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      complex*16 spence,csp1,csp2
      complex*16 trilog,ctr1,ctr2,ctr3,ctr4,ctr5
      real*8 arg,argt1,argt2,argt3,argt4,argt5
      external spence
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         z  = 1.d0/(1.d0-velo**2)
         rtmp =  3.d0/4.d0/pi * dimag(1.7695276526121257d0 - 6
     &        .579736267392905d0*dlog(2.d0) + ((53.56935816521245d0 + 23
     &        .42800818317322d0* (-1.d0 + 2.d0/z + (0.d0,2.d0)*dsqrt((-1
     &        .d0 + z)/z **2)) - 28.58891755403711d0* (-1.d0 + 2.d0/z +
     &        (0.d0,2.d0)*dsqrt((-1 .d0 + z)/z**2))**2 - 10
     &        .11103466811403d0* (-1.d0 + 2.d0/z + (0.d0,2.d0)*dsqrt((-1
     &        .d0 + z)/z**2))**3 + 0.4266396586616486d0* (-1.d0 + 2 .d0
     &        /z + (0.d0,2.d0)*dsqrt((-1.d0 + z)/z**2))**4)* (2.d0/z +
     &        (0.d0,2.d0) *dsqrt((-1.d0 + z)/z**2))**2)/ (19
     &        .19224589510177d0 - 1 .767625499084376d0* (-1.d0 + 2.d0/z
     &        + (0.d0,2.d0)*dsqrt((-1.d0 + z)/z **2)) - 13
     &        .2990319119615d0*(-1.d0 + 2.d0/z + (0.d0,2.d0)*dsqrt((-1
     &        .d0 + z)/z**2))** 2 + 1.00664474828887d0* (-1.d0 + 2.d0/z
     &        + (0.d0,2.d0) *dsqrt((-1.d0 + z)/z**2))**3 + 1.d0*(-1.d0 +
     &        2.d0/z + (0.d0,2.d0) *dsqrt((-1.d0 + z)/z**2))**4) - (0
     &        .125d0*(33.d0*((0.d0,3 .141592653589793d0) + dlog((1.d0 -
     &        1.d0*dsqrt(1.d0 - 1.d0/z))/(1.d0 + dsqrt(1.d0 - 1.d0/z))))
     &        - 27.d0*z*((0.d0,3.141592653589793d0) + dlog((1.d0 - 1.d0
     &        *dsqrt(1.d0 - 1.d0/z))/(1.d0 + dsqrt(1.d0 - 1.d0/z)))) - 6
     &        .d0*z**2*((0.d0,3.141592653589793d0) + dlog((1.d0 - 1.d0
     &        *dsqrt(1.d0 - 1.d0/z))/(1.d0 + dsqrt(1.d0 - 1.d0/z)))) -
     &        66.d0*z*dsqrt((-1.d0 + z) /z) + 10.d0*z**2*dsqrt((-1.d0 +
     &        z)/z)))/ (dsqrt((-1.d0 + z)/z)*z **2))

         arg = (1.d0 - velo)/(1.d0 + velo)
         csp1 = spence(arg)
         csp2 = spence(arg**2)
         sp1 = dreal(csp1)
         sp2 = dreal(csp2)
         argt1 = arg
         argt2 = arg**2
         argt3 = 1.d0-arg
         argt4 = (1.d0-arg)*(1.d0+arg)
         argt5 = arg/(1.d0+arg)
         ctr1 = trilog(argt1)
         ctr2 = trilog(argt2)
         ctr3 = trilog(argt3)
         ctr4 = trilog(argt4)
         ctr5 = trilog(argt5)
         tr1 = dreal(ctr1)
         tr2 = dreal(ctr2)
         tr3 = dreal(ctr3)
         tr4 = dreal(ctr4)
         tr5 = dreal(ctr5)

         rtmp = rtmp + 18.70313878328875d0 - 12.182816401997158d0*velo +
     &        9.453046732970753d0*velo**2 + 3.401216578443497d0*velo**3
     &        - 3.972595252513273d0*velo**4 + tr1*(1.375d0 + 0
     &        .9166666666666666d0*velo**2 - 0.4583333333333333d0*velo**4
     &        ) + tr3*(-2.75d0 - 1.8333333333333333d0*velo**2 + 0
     &        .9166666666666666d0*velo**4) + tr2*(-4.125d0 - 2.75d0*velo
     &        **2 + 1.375d0*velo**4) + tr5*(-5.5d0 - 3
     &        .6666666666666665d0*velo**2 + 1.8333333333333333d0*velo**4
     &        ) + tr4*(-6.875d0 - 4.583333333333333d0*velo**2 + 2
     &        .2916666666666665d0*velo**4) + (-8.25d0*velo + 2.75d0*velo
     &        **3)*dlog(2.d0)**2 + (-5.729166666666667d0 - 3
     &        .8194444444444446d0*velo**2 + 1.9097222222222223d0*velo**4
     &        )* dlog((1.d0 - 1.d0*velo)/(1.d0 + velo))**3 + dlog(velo)
     &        **2*(11.d0*velo - 3.6666666666666665d0*velo**3 + (-11.d0 -
     &        7.333333333333333d0*velo**2 + 3.6666666666666665d0*velo**4
     &        )* dlog((1.d0 - 1.d0*velo)/(1.d0 + velo))) + (-0
     &        .11458333333333333d0 - 0.0763888888888889d0*velo**2 + 0
     &        .03819444444444445d0*velo**4)* dlog(0.25d0*(1.d0 - 1.d0
     &        *velo**2))**3 + dlog(2.d0)*(-8.25d0*velo + 2.75d0*velo**3
     &        + (4.125d0 + 2.75d0*velo**2 - 1.375d0*velo**4)* dlog((1.d0
     &        - 1.d0*velo)/(1.d0 + velo))**2 + dlog(velo)*(11.d0*velo -
     &        3.6666666666666665d0*velo**3 + (-24.75d0 - 16.5d0*velo**2
     &        + 8.25d0*velo**4)* dlog((1.d0 - 1.d0*velo)/(1.d0 + velo)))
     &        + dlog((1.d0 - 1.d0*velo)/(1.d0 + velo))* (-3.78125d0 - 2
     &        .5208333333333335d0*velo**2 + 0.8020833333333334d0*velo**4
     &        + (-4.125d0 - 2.75d0*velo**2 + 1.375d0*velo**4)* dlog(0
     &        .25d0*(1.d0 - 1.d0*velo**2)))) + (2.0625d0*velo - 0.6875d0
     &        *velo**3)*dlog(1.d0 - 1.d0*velo**2)**2 + dlog(velo)*(5.5d0
     &        *velo - 1.8333333333333333d0*velo**3 + (-15.125d0 - 10
     &        .083333333333334d0*velo**2 + 5.041666666666667d0*velo**4)
     &        * dlog((1.d0 - 1.d0*velo)/(1.d0 + velo))**2 + (-11.d0*velo
     &        + 3.6666666666666665d0*velo**3)* dlog(1.d0 - 1.d0*velo**2)
     &        + dlog((1.d0 - 1.d0*velo)/(1.d0 + velo))* (3.4375d0 + 2
     &        .2916666666666665d0*velo**2 - 1.1458333333333333d0 *velo
     &        **4 + (15.125d0 + 10.083333333333334d0*velo**2 - 5
     &        .041666666666667d0*velo**4)*dlog(1.d0 - 1.d0*velo**2))) +
     &        sp1*(-1.5416666666666667d0 + 9.972222222222221d0*velo**2 -
     &        3.1527777777777777d0*velo**4 + (-5.5d0 - 3
     &        .6666666666666665d0*velo**2 + 1.8333333333333333d0*velo**4
     &        )* dlog(4.d0) + (2.75d0 + 1.8333333333333333d0*velo**2 - 0
     &        .9166666666666666d0*velo**4)*dlog(1.d0 - 1.d0*velo**2) +
     &        (1.375d0 + 0.9166666666666666d0*velo**2 - 0
     &        .4583333333333333d0*velo**4)* dlog((4.d0*(1.d0 - 1.d0*velo
     &        **2))/velo**4)) + (11.6875d0*velo - 4.8125d0*velo**3)
     &        * dlog((0.5d0*(1.d0 - 1.d0*velo**2))/velo**2) + sp2*(-0
     &        .4270833333333333d0 + 1.5486111111111112d0*velo**2 - 0
     &        .3159722222222222d0*velo**4 + (-2.75d0 - 1
     &        .8333333333333333d0*velo**2 + 0.9166666666666666d0*velo**4
     &        )* dlog(4.d0) + (1.375d0 + 0.9166666666666666d0*velo**2 -
     &        0.4583333333333333d0*velo**4)*dlog(1.d0 - 1.d0*velo**2) +
     &        (2.75d0 + 1.8333333333333333d0*velo**2 - 0
     &        .9166666666666666d0*velo**4)* dlog((0.5d0*(1.d0 - 1.d0
     &        *velo**2))/velo**2)) + (-4.523568683832623d0 - 3
     &        .015712455888415d0*velo**2 + 1.5078562279442076d0*velo**4)
     &        * dlog((0.25d0*(1.d0 - 1.d0*velo**2))/velo) + dlog((1.d0 -
     &        1.d0*velo)/(1.d0 + velo))**2* (-2.75d0 + 1.03125d0/velo -
     &        2.75d0*velo + 2.75d0*velo**2 + 1.2604166666666667d0*velo
     &        **3 - 0.4583333333333333d0*velo**4 + (11.34375d0 + 7
     &        .5625d0*velo**2 - 3.78125d0*velo**4)* dlog(0.25d0*(1.d0 -
     &        1.d0*velo**2)) + (0.6875d0 + 0.4583333333333333d0*velo**2
     &        - 0.22916666666666666d0*velo**4)* dlog(1.d0 - (1.d0*(1.d0
     &        - 1.d0*velo)**2)/(1.d0 + velo)**2)) + (-5
     &        .166666666666667d0*velo + 1.7222222222222223d0*velo**3)
     &        * dlog(1.d0 - (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo)) + (
     &        -2.5833333333333335d0*velo + 0.8611111111111112d0*velo**3)
     &        * dlog(1.d0 + (1.d0 - 1.d0*velo)/(1.d0 + velo)) + dlog(1
     &        .d0 - 1.d0*velo**2)*(1.71875d0*velo - 1.03125d0*velo**3 +
     &        (-5.5d0*velo + 1.8333333333333333d0*velo**3)* dlog(1.d0 -
     &        (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo)) + (-2.75d0*velo +
     &        0.9166666666666666d0*velo**3)* dlog(1.d0 + (1.d0 - 1.d0
     &        *velo)/(1.d0 + velo))) + dlog(mu**2/mq**2)*(1.71875d0*velo
     &        - 1.03125d0*velo**3 + sp1*(2.75d0 + 1.8333333333333333d0
     &        *velo**2 - 0.9166666666666666d0*velo**4) + sp2*(1.375d0 +
     &        0.9166666666666666d0*velo**2 - 0.4583333333333333d0*velo
     &        **4) + (-5.5d0*velo + 1.8333333333333333d0*velo**3)
     &        * dlog(1.d0 - (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo)) + (
     &        -2.75d0*velo + 0.9166666666666666d0*velo**3)* dlog(1.d0 +
     &        (1.d0 - 1.d0*velo)/(1.d0 + velo)) + dlog((1.d0 - 1.d0*velo
     &        )/(1.d0 + velo))* (-1.890625d0 + 4.125d0*velo - 1
     &        .2604166666666667d0*velo**2 - 1.375d0*velo**3 + 0
     &        .4010416666666667d0*velo**4 + (2.75d0 + 1
     &        .8333333333333333d0*velo**2 - 0.9166666666666666d0*velo**4
     &        )* dlog(1.d0 - (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo)) +
     &        (1.375d0 + 0.9166666666666666d0*velo**2 - 0
     &        .4583333333333333d0*velo**4)* dlog(1.d0 + (1.d0 - 1.d0
     &        *velo)/(1.d0 + velo)))) + dlog((1.d0 - 1.d0*velo)/(1.d0 +
     &        velo))* (-4.0415354828340435d0 + 3.875d0*velo + 0
     &        .05564301144397055d0*velo**2 - 1.2916666666666667d0*velo
     &        **3 + 0.5589840498335701d0*velo**4 + (-5.5d0 - 3
     &        .6666666666666665d0*velo**2 + 1.8333333333333333d0*velo**4
     &        )* dlog(0.25d0*(1.d0 - 1.d0*velo**2))**2 + (1.71875d0 - 2
     &        .5208333333333335d0*velo**2 + 0.34375d0*velo**4)* dlog((0
     &        .25d0*(1.d0 - 1.d0*velo**2))/velo**2) + (2
     &        .5833333333333335d0 + 1.7222222222222223d0*velo**2 - 0
     &        .8611111111111112d0*velo**4)* dlog(1.d0 - (1.d0*(1.d0 - 1
     &        .d0*velo))/(1.d0 + velo)) + (1.2916666666666667d0 + 0
     &        .8611111111111112d0*velo**2 - 0.4305555555555556d0*velo**4
     &        )* dlog(1.d0 + (1.d0 - 1.d0*velo)/(1.d0 + velo)) 
     &        + dlog(1.d0 - 1.d0*velo**2)* (-1.890625d0 + 4.125d0*velo -
     &        1.2604166666666667d0*velo**2 - 1.375d0*velo**3 + 0
     &        .4010416666666667d0*velo**4 + (2.75d0 + 1
     &        .8333333333333333d0*velo**2 - 0.9166666666666666d0*velo**4
     &        )* dlog(1.d0 - (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo)) +
     &        (1.375d0 + 0.9166666666666666d0*velo**2 - 0
     &        .4583333333333333d0*velo**4)* dlog(1.d0 + (1.d0 - 1.d0
     &        *velo)/(1.d0 + velo)))) + dlog(4.d0)*(-3.4375d0*velo + 2
     &        .0625d0*velo**3 + (11.d0*velo - 3.6666666666666665d0*velo
     &        **3)* dlog(1.d0 - (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo))
     &        + (5.5d0*velo - 1.8333333333333333d0*velo**3)* dlog(1.d0 +
     &        (1.d0 - 1.d0*velo)/(1.d0 + velo)) + dlog((1.d0 - 1.d0*velo
     &        )/(1.d0 + velo))* (3.78125d0 - 8.25d0*velo + 2
     &        .5208333333333335d0*velo**2 + 2.75d0*velo**3 - 0
     &        .8020833333333334d0*velo**4 + (-5.5d0 - 3
     &        .6666666666666665d0*velo**2 + 1.8333333333333333d0*velo**4
     &        )* dlog(1.d0 - (1.d0*(1.d0 - 1.d0*velo))/(1.d0 + velo)) +
     &        (-2.75d0 - 1.8333333333333333d0*velo**2 + 0
     &        .9166666666666666d0*velo**4)* dlog(1.d0 + (1.d0 - 1.d0
     &        *velo)/(1.d0 + velo))))

      else
         rvca = r0ca(scms)
         return
      endif
         
      rvca = cf*ca * rtmp

      return
      end

C-}}}
C-{{{ function rvnl:

c..   ------------------------------------------------------------

      function rvnl(scms,mq)
c..
c..   two-loop
c..   vector case
c..   double bubble with massless inner loop.
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      complex*16 spence,csp1,csp2
      complex*16 trilog,ctr1,ctr2,ctr3,ctr4,ctr5
      real*8 arg,argt1,argt2,argt3,argt4,argt5
      external spence
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      dlsmu = dlog(scms/mu/mu)

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         arg = (1.d0 - velo)/(1.d0 + velo)
         csp1 = spence(arg)
         csp2 = spence(arg**2)
         sp1 = dreal(csp1)
         sp2 = dreal(csp2)
         argt1 = arg
         argt2 = arg**2
         argt3 = 1.d0-arg
         argt4 = (1.d0-arg)*(1.d0+arg)
         argt5 = arg/(1.d0+arg)
         ctr1 = trilog(argt1)
         ctr2 = trilog(argt2)
         ctr3 = trilog(argt3)
         ctr4 = trilog(argt4)
         ctr5 = trilog(argt5)
         tr1 = dreal(ctr1)
         tr2 = dreal(ctr2)
         tr3 = dreal(ctr3)
         tr4 = dreal(ctr4)
         tr5 = dreal(ctr5)

         rtmp = (velo*(-75 + 29*velo**2))/48.d0 + (pi**2*(-1 + velo)
     &        * (51 - 45*velo - 27*velo**2 + 5*velo**3))/144.d0 + ((sp1
     &        + sp2)*(15 - 6*velo**2 - velo**4))/24.d0 + (sp1*(7 - 22
     &        *velo**2 + 7*velo**4))/8.d0 + (velo*(-3 + velo**2)*(-3
     &        *dlog(2.d0) + 2*dlog(velo)))/3.d0 + ((237 + 62*velo**2 -
     &        59*velo**4)* dlog((1 - velo)/(1 + velo)))/96.d0 + ((1 +
     &        velo)*(-9 + 33*velo - 9*velo**2- 15*velo**3 + 4*velo**4)
     &        * dlog((1 - velo)/(1 + velo))**2)/(24.d0*velo) + (velo*(-3
     &        + velo**2)*(2*dlog(2.d0) - 4*dlog(velo) + dlog(1 - velo**2
     &        ))* (-6*dlog(2.d0) - 4*dlog(velo) + 3*dlog(1 - velo**2)))
     &        /12.d0 + dlog((1 - velo)/(1 + velo))* (((33 + 22*velo**2 -
     &        7*velo**4)*dlog(2.d0))/24.d0 + (5*(-3 + velo**2)*(1 + velo
     &        **2)*dlog(velo))/12.d0 + ((-15 + 22*velo**2 - 3*velo**4)
     &        * dlog((1 - velo**2)/(4.d0*velo**2)))/ 24.d0) + (velo*(-17
     &        + 7*velo**2)* dlog((1 - velo**2)/(2.d0*velo**2)))/ 4.d0 -
     &        ((1.6666666666666667d0 - 2*dlog(4.d0) + dlog(mu**2/mq**2)
     &        + dlog(1 - velo**2))* ((3*velo*(5 - 3*velo**2))/8.d0 - ((1
     &        - velo)*(33 - 39*velo - 17*velo**2 + 7*velo**3)* dlog((1 -
     &        velo)/(1 + velo)))/16.d0 + ((3 - velo**2)*(-2*velo* (2
     &        *dlog(1 - (1 - velo)/(1 + velo)) + dlog(1 + (1 - velo)/(1
     &        + velo))) + (1 + velo**2)* (2*sp1 + sp2 + 2*dlog((1 - velo
     &        )/(1 + velo))* dlog(1 - (1 - velo)/(1 + velo)) + dlog((1 -
     &        velo)/(1 + velo))* dlog(1 + (1 - velo)/(1 + velo)))))/2.d0
     &        ))/3.d0 + ((-3 + velo**2)*(1 + velo**2)* (tr1 - 3*tr2 - 2
     &        *tr3 - 5*tr4 - 4*tr5 + dlog(velo)*dlog((1 - velo)/(1 +
     &        velo))* (-18*dlog(2.d0) - 8*dlog(velo) - 11* dlog((1 -
     &        velo)/(1 + velo)) + 11*dlog(1 - velo**2)) + sp1*dlog((4*(1
     &        - velo**2))/velo**4) + 2*sp2*dlog((1 - velo**2)/(2.d0*velo
     &        **2)) + (Pi**2*(dlog((1 - velo)/(1 + velo)) - dlog((1 -
     &        velo**2)/(4.d0*velo))))/3.d0 + (36*dlog(2.d0)*dlog((1 -
     &        velo)/(1 + velo))**2 - 50*dlog((1 - velo)/(1 + velo))**3 -
     &        36*dlog(2.d0)*dlog((1 - velo)/(1 + velo))* dlog((1 - velo
     &        **2)/4.d0) + 99*dlog((1 - velo)/(1 + velo))**2*dlog((1 -
     &        velo**2)/4.d0) - 48*dlog((1 - velo)/(1 + velo))*dlog((1 -
     &        velo**2)/4.d0)**2 - dlog((1 - velo**2)/4.d0)**3 + 6
     &        *dlog((1 - velo)/(1 + velo))**2* dlog(1 - (1 - velo)**2/(1
     &        + velo)**2))/12.d0 + (11*1.20205690315959428539973816151d0
     &        )/2.d0))/6.d0
      else
         rtmp = -11.d0/8.d0 + zeta3 + 1/4.d0*dlsmu
      endif

      rvnl = cf*tr * rtmp

      return
      end

C-}}}
C-{{{ function rvct:

c..   ------------------------------------------------------------

      function rvct(scms,mq)
c..
c..   two-loop
c..   vector case
c..   double bubble with inner quark = outer quark
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      dlsmu = dlog(scms/mu/mu)

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         if (lrvctexp) then
            rtmp = rvctexp(scms,mq)
         else
            if (mq**2/scms.gt..01d0) then
               rtmp = rvctrv(scms,mq)
     &              + 1.d0/4.d0 * dlog(mq**2/mu**2) * rv1(scms,mq)
            else
c..   use high-energy expansion if mq^2/s is too small 
c..   (rvctreal becomes unstable):
               rtmp = rvctexp(scms,mq)
            endif
         endif
      else
         rtmp = -11.d0/8.d0 + zeta3 + 1/4.d0*dlsmu
      endif

      rvct = cf*tr * rtmp

      return
      end

C-}}}
C-{{{ function rv2ctms:

c..   ------------------------------------------------------------

      function rv2ctms(scms,mq)
c..
c..   two-loop counterterm contribution; OS->MS-bar
c..   vector case
c..
c..   scms = center of mass energy squared
c..   mq   = quark mass
c..
c..   Note, 'rtmp1' and 'rtmp2' contain already the colour factors
c..   but not the overall factor nc
c..   
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      complex*16 spence,csp1,csp2
      real*8 arg
      external spence
      
      include 'common.f'

      velo = dsqrt( 1.d0 - 4.d0*mq**2/scms )

      if ((.not.lmassless).and.(velo.ne.1.d0)) then
         arg = (1.d0 - velo)/(1.d0 + velo)
         csp1 = spence(arg)
         csp2 = spence(arg**2)
         sp1 = dreal(csp1)
         sp2 = dreal(csp2)

c..   expression induced from Born term
         rtmp1 =
     -    (0.0625*(-21.333333333333332 
     -        - 230.25*velo**2 + 524.5*velo**4 - 
     -       272.9166666666667*velo**6 - 48.*velo**2*zeta2 - 
     -       16.*dlog(2.d0)*velo**2*zeta2 + 96.*velo**4*zeta2 + 
     -       32.*dlog(2.d0)*velo**4*zeta2 - 48.*velo**6*zeta2 - 
     -       16.*dlog(2.d0)*velo**6*zeta2 + 4.*velo**2*zeta3 
     -        - 8.*velo**4*zeta3 + 
     -       4.*velo**6*zeta3 - 32.*dlog(mu**2/mq**2) - 
     -       164.33333333333334*velo**2*dlog(mu**2/mq**2) + 
     -       424.6666666666667*velo**4*dlog(mu**2/mq**2) - 
     -       228.33333333333334*velo**6*dlog(mu**2/mq**2) - 
     -       12.*dlog(mu**2/mq**2)**2 
     -        - 43.*velo**2*dlog(mu**2/mq**2)**2 + 
     -       122.*velo**4*dlog(mu**2/mq**2)**2 - 
     -       67.*velo**6*dlog(mu**2/mq**2)**2 + 
     -       (inffin-1)*(11.833333333333334*velo**2 
     -        - 23.666666666666668*velo**4 + 
     -          11.833333333333334*velo**6 + 8.*velo**2*zeta2 - 
     -          16.*velo**4*zeta2 + 8.*velo**6*zeta2 + 
     -          8.666666666666666*velo**2*dlog(mu**2/mq**2) - 
     -          17.333333333333332*velo**4*dlog(mu**2/mq**2) + 
     -          8.666666666666666*velo**6*dlog(mu**2/mq**2) + 
     -          2.*velo**2*dlog(mu**2/mq**2)**2 - 
     -          4.*velo**4*dlog(mu**2/mq**2)**2 + 
     -          2.*velo**6*dlog(mu**2/mq**2)**2)))/velo**3

c..   expression induced from O(as) corrections
         rtmp2 =
     -        (6.222222222222222*velo**3 + 6.222222222222222*velo**4 - 
     -     18.666666666666668*velo**5 - 18.666666666666668*velo**6 + 
     -     18.666666666666668*velo**7 + 18.666666666666668*velo**8 - 
     -     6.222222222222222*velo**9 - 6.222222222222222*velo**10 + 
     -     sp1*(-7.111111111111111*velo**2 - 7.111111111111111*velo**3 + 
     -        28.444444444444443*velo**4 + 28.444444444444443*velo**5 - 
     -        42.666666666666664*velo**6 - 42.666666666666664*velo**7 + 
     -        28.444444444444443*velo**8 + 28.444444444444443*velo**9 - 
     -        7.111111111111111*velo**10 - 7.111111111111111*velo**11 + 
     -        (-5.333333333333333*velo**2 - 5.333333333333333*velo**3 + 
     -           21.333333333333332*velo**4 
     -        + 21.333333333333332*velo**5 - 
     -           32.*velo**6 - 32.*velo**7 
     -        + 21.333333333333332*velo**8 + 
     -           21.333333333333332*velo**9 
     -        - 5.333333333333333*velo**10 - 
     -           5.333333333333333*velo**11)*dlog(mu**2/mq**2)) + 
     -     sp2*(-3.5555555555555554*velo**2 
     -        - 3.5555555555555554*velo**3 + 
     -        14.222222222222221*velo**4 + 14.222222222222221*velo**5 - 
     -        21.333333333333332*velo**6 - 21.333333333333332*velo**7 + 
     -        14.222222222222221*velo**8 + 14.222222222222221*velo**9 - 
     -        3.5555555555555554*velo**10 
     -        - 3.5555555555555554*velo**11 + 
     -        (-2.6666666666666665*velo**2 
     -        - 2.6666666666666665*velo**3 + 
     -           10.666666666666666*velo**4 
     -        + 10.666666666666666*velo**5 - 
     -           16.*velo**6 - 16.*velo**7 
     -        + 10.666666666666666*velo**8 + 
     -           10.666666666666666*velo**9 
     -        - 2.6666666666666665*velo**10 - 
     -           2.6666666666666665*velo**11)*dlog(mu**2/mq**2)) + 
     -     (-10.666666666666666*velo - 10.666666666666666*velo**2 + 
     -        14.222222222222221*velo**3 + 14.222222222222221*velo**4 + 
     -        7.111111111111111*velo**5 + 7.111111111111111*velo**6 - 
     -        14.222222222222221*velo**7 - 14.222222222222221*velo**8 + 
     -        3.5555555555555554*velo**9 + 3.5555555555555554*velo**10)*
     -      dlog(1. - (1.*(-1. + velo)**2)/(1. + velo)**2) + 
     -     (10.666666666666666*velo + 10.666666666666666*velo**2 - 
     -        42.666666666666664*velo**3 - 42.666666666666664*velo**4 + 
     -        64.*velo**5 + 64.*velo**6 - 42.666666666666664*velo**7 - 
     -        42.666666666666664*velo**8 + 10.666666666666666*velo**9 + 
     -        10.666666666666666*velo**10)*
     -      dlog(1. - (1.*(1. - 1.*velo))/(1. + velo)) + 
     -     (10.666666666666666*velo + 10.666666666666666*velo**2 - 
     -        28.444444444444443*velo**3 - 28.444444444444443*velo**4 + 
     -        28.444444444444443*velo**5 + 28.444444444444443*velo**6 - 
     -        14.222222222222221*velo**7 - 14.222222222222221*velo**8 + 
     -        3.5555555555555554*velo**9 + 3.5555555555555554*velo**10)*
     -      dlog(1. + (1. - 1.*velo)/(1. + velo)) + 
     -     dlog((1. - 1.*velo)/(1. + velo))*
     -      (-5.333333333333333 - 5.333333333333333*velo + 
     -        9.333333333333334*velo**2 + 30.666666666666668*velo**3 + 
     -        10.666666666666666*velo**4 - 64.*velo**5 - 
     -        34.666666666666664*velo**6 + 61.333333333333336*velo**7 + 
     -        26.666666666666668*velo**8 - 26.666666666666668*velo**9 - 
     -        6.666666666666667*velo**10 + 4.*velo**11 + 
     -        (-7.111111111111111*velo**2 - 7.111111111111111*velo**3 + 
     -           28.444444444444443*velo**4 
     -        + 28.444444444444443*velo**5 - 
     -           42.666666666666664*velo**6 
     -        - 42.666666666666664*velo**7 + 
     -           28.444444444444443*velo**8 
     -        + 28.444444444444443*velo**9 - 
     -         7.111111111111111*velo**10 - 7.111111111111111*velo**11)*
     -         dlog(1. - (1.*(1. - 1.*velo))/(1. + velo)) + 
     -       (-3.5555555555555554*velo**2 - 3.5555555555555554*velo**3 + 
     -         14.222222222222221*velo**4 + 14.222222222222221*velo**5 - 
     -         21.333333333333332*velo**6 - 21.333333333333332*velo**7 + 
     -         14.222222222222221*velo**8 + 14.222222222222221*velo**9 - 
     -       3.5555555555555554*velo**10 - 3.5555555555555554*velo**11)*
     -         dlog(1. + (1. - 1.*velo)/(1. + velo))) + 
     -     dlog(mu**2/mq**2)*(4.666666666666667*velo**3 + 
     -        4.666666666666667*velo**4 - 14.*velo**5 - 14.*velo**6 + 
     -        14.*velo**7 + 14.*velo**8 - 4.666666666666667*velo**9 - 
     -        4.666666666666667*velo**10 + 
     -        (-8.*velo - 8.*velo**2 + 10.666666666666666*velo**3 + 
     -          10.666666666666666*velo**4 + 5.333333333333333*velo**5 + 
     -          5.333333333333333*velo**6 - 10.666666666666666*velo**7 - 
     -         10.666666666666666*velo**8 + 2.6666666666666665*velo**9 + 
     -           2.6666666666666665*velo**10)*
     -         dlog(1. - (1.*(-1. + velo)**2)/(1. + velo)**2) + 
     -        (8.*velo + 8.*velo**2 - 32.*velo**3 - 32.*velo**4 + 
     -           48.*velo**5 + 48.*velo**6 - 32.*velo**7 - 32.*velo**8 + 
     -           8.*velo**9 + 8.*velo**10)*
     -         dlog(1. - (1.*(1. - 1.*velo))/(1. + velo)) + 
     -        (8.*velo + 8.*velo**2 - 21.333333333333332*velo**3 - 
     -         21.333333333333332*velo**4 + 21.333333333333332*velo**5 + 
     -         21.333333333333332*velo**6 - 10.666666666666666*velo**7 - 
     -         10.666666666666666*velo**8 + 2.6666666666666665*velo**9 + 
     -           2.6666666666666665*velo**10)*
     -         dlog(1. + (1. - 1.*velo)/(1. + velo)) + 
     -        dlog((1. - 1.*velo)/(1. + velo))*
     -         (-4. - 4.*velo + 7.*velo**2 + 23.*velo**3 + 8.*velo**4 - 
     -           48.*velo**5 - 26.*velo**6 + 46.*velo**7 + 20.*velo**8 - 
     -           20.*velo**9 - 5.*velo**10 + 3.*velo**11 + 
     -         (-5.333333333333333*velo**2 - 5.333333333333333*velo**3 + 
     -         21.333333333333332*velo**4 + 21.333333333333332*velo**5 - 
     -         32.*velo**6 - 32.*velo**7 + 21.333333333333332*velo**8 + 
     -         21.333333333333332*velo**9 - 5.333333333333333*velo**10 - 
     -         5.333333333333333*velo**11)*
     -            dlog(1. - (1.*(1. - 1.*velo))/(1. + velo)) + 
     -         (-2.6666666666666665*velo**2 
     -        - 2.6666666666666665*velo**3 + 
     -         10.666666666666666*velo**4 + 10.666666666666666*velo**5 - 
     -         16.*velo**6 - 16.*velo**7 + 10.666666666666666*velo**8 + 
     -              10.666666666666666*velo**9 - 
     -        2.6666666666666665*velo**10 - 2.6666666666666665*velo**11)
     -             *dlog(1. + (1. - 1.*velo)/(1. + velo)))))/
     -   ((-1. + velo)**2*velo**2*(1. + velo)**3)
      else
         rtmp1 = 0.d0
         rtmp2 = 0.d0
      endif

      rv2ctms = rtmp1 + rtmp2

      return
      end

C-}}}
C-{{{ function r0db:

c..   ------------------------------------------------------------

      function r0db(scms,min,label)
c..
c..   two-loop
c..   vector case
c..   double bubble
c..   outer loop massless
c..
c..   scms: cms energy squared
c..   min : mass in inner quark loop
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common.f'

      if (min.eq.0.d0) then    
         r0db = r0nl(scms)
      else
         if (label) then
            if (lmassless) then
               r0db = r0nl(scms)
            else
               r0db = rv2cbl(scms,0.d0,min)
            endif
         else
            if (lpsup) then  ! include power suppressed terms?
               r0db = r0dbh(scms,min)
            else
               r0db = 0.d0
            endif
         endif
      endif

      return
      end

C-}}}
C-{{{ function rvdb:

c..   ------------------------------------------------------------

      function rvdb(scms,mout,min,label)
c..
c..   two-loop
c..   vector case
c..   double bubble
c..
c..   scms: cms energy squared
c..   min : mass of quark coupling to gluons only
c..   mout: mass of quark coupling to external current
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common.f'

      if (lmassless) then
         rvdb = r0db(scms,min,label)
         return
      endif

      if (min.le.mout) then         !  (I)  min <= mout
         if (min.eq.0.d0) then      !  (Ia) min == 0
            rvdb = rvnl(scms,mout)
         elseif (min.eq.mout) then  !  (Ib) min == mout
            rvdb = rvct(scms,mout)
         else                       !  (Ic) min < mout
c..   no mass corrections due to inner loop are included:
            rvdb = rvnl(scms,mout)
         endif
c..         
      else                          !  (II) min > mout
         if (label) then 
c..   (IIa) above threshold
            rvdb = rv2cbl(scms,mout,min)
         else                       !  (IIb) below threshold
            if (lpsup) then
               rvdb = rv2cbh(scms,mout,min)
            else
               rvdb = 0.d0
            endif
         endif
      endif
      
      return
      end

C-}}}
C-{{{ function rv2cbl:

c..   ------------------------------------------------------------

      function rv2cbl(scms,mout,min)
c..
c..   two-loop
c..   vector case
c..   double bubble
c..   light quark in outer loop
c..   valid above threshold, i.e.,
c..
c..   ****    scms > 4*thrq^2    ****
c..
c..   analytical expression
c..
c..   scms: cms energy squared
c..   min : mass of quark coupling to gluons only
c..   mout: mass of quark coupling to external current
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      complex*16 spence,csp1,csp2,csp3,csp4,csp5,csp6,csp7
      complex*16 trilog,ctr1,ctr2,ctr3,ctr4,ctr5,ctr6,ctr7
      real*8 arg1,arg2,arg3,arg4
      external spence,trilog
      
      include 'common.f'

      if (mout.ge.min) then
         write(6,*) '<function rv2cbl>: invalid input'
         write(6,*) '     sqrts = ',dsqrt(scms)
         write(6,*) '     mout  = ',mout
         write(6,*) '     min   = ',min
         stop
      endif

      velo = dsqrt(1.d0-4.d0*min**2/scms)

      dlsmu = dlog(scms/mu/mu)

      arg1 = (1.d0 - velo)/2.d0
      arg2 = (1.d0 + velo)/2.d0
      arg3 = (1.d0 - dsqrt(2.d0-velo**2))/2.d0
      arg4 = (1.d0 + dsqrt(2.d0-velo**2))/2.d0
      arg7 = ( -dsqrt(1.d0/(1.d0-velo**2))+
     &     dsqrt(1.d0+1.d0/(1.d0-velo**2)) )**2
      csp1 = spence(arg1)
      csp2 = spence(arg2)
      csp3 = spence(arg1/arg3)
      csp4 = spence(arg1/arg4)
      csp5 = spence(arg2/arg4)
      csp6 = spence(arg3/arg2)
      csp7 = spence(arg7)
      ctr1 = trilog(arg1)
      ctr2 = trilog(arg2)
      ctr3 = trilog(arg1/arg3)
      ctr4 = trilog(arg1/arg4)
      ctr5 = trilog(arg2/arg4)
      ctr6 = trilog(arg3/arg2)
      ctr7 = trilog(arg7)
      sp1 = dreal(csp1)
      sp2 = dreal(csp2)
      sp3 = dreal(csp3)
      sp4 = dreal(csp4)
      sp5 = dreal(csp5)
      sp6 = dreal(csp6)
      sp7 = dreal(csp7)
      tr1 = dreal(ctr1)
      tr2 = dreal(ctr2)
      tr3 = dreal(ctr3)
      tr4 = dreal(ctr4)
      tr5 = dreal(ctr5)
      tr6 = dreal(ctr6)
      tr7 = dreal(ctr7)

      rtmp = 7.982167648374861d0 - 10.393518518518519d0*velo - 3
     &     .9065840071353524d0*velo**2 + 3.8410493827160495d0*velo**3 +
     &     0 .3005142257898985d0*velo**4 + tr5*(0.8333333333333334d0 +
     &     velo**2 - 0.5d0*velo**4) + tr6*(0.8333333333333334d0 + velo
     &     **2 - 0.5d0*velo **4) + tr1*(0.4166666666666667d0 + 0.5d0
     &     *velo**2 - 0.25d0*velo**4) + tr7*(0.4166666666666667d0 + 0
     &     .5d0*velo**2 - 0.25d0*velo**4) + tr2 *(-0.4166666666666667d0
     &     - 0.5d0*velo**2 + 0.25d0*velo**4) + tr3*(-0
     &     .8333333333333334d0 - 1*velo**2 + 0.5d0*velo**4) + tr4*(-0
     &     .8333333333333334d0 - 1*velo**2 + 0.5d0*velo**4) + 0.25d0
     &     *dlog(scms/mu**2) + (-0.034722222222222224d0 - 0
     &     .041666666666666664d0*velo**2 + 0.020833333333333332d0*velo
     &     **4) * dlog((1 + velo)/(1 - 1*velo))**3 + sp2*(-2
     &     .3055555555555554d0 + 1.5d0*velo**2 - 0.25d0*velo**4 + (0
     &     .4166666666666667d0 + 0.5d0*velo**2 - 0.25d0*velo**4)* dlog(1
     &      + velo) + (0.20833333333333334d0 + 0.25d0*velo**2 - 0.125d0
     &     *velo**4) *dlog(1/(1 - 1*velo**2))) + sp1*(2
     &     .3055555555555554d0 - 1.5d0 *velo**2 + 0.25d0*velo**4 + (-0
     &     .4166666666666667d0 - 0.5d0*velo**2 + 0.25d0*velo**4)* dlog(1
     &      - 1*velo) + (-0.20833333333333334d0 - 0.25d0*velo**2 + 0
     &     .125d0*velo**4)* dlog(1/(1 - 1*velo**2))) + 0 .25d0*dlog(0
     &     .25d0*(1 - 1*velo**2)) + (-7.962962962962963d0 + 3
     &     .0555555555555554d0*velo**2)* dlog(2*dsqrt(1/(1 - 1*velo
     &     **2))) + (1.3707783890401886d0 + 1.6449340668482262d0*velo**2
     &     - 0 .8224670334241131d0*velo**4)* dlog((-1 + dsqrt(2 - 1
     &     *velo**2))/(1 + velo)) + (0.1388888888888889d0 + 0
     &     .16666666666666666d0 *velo**2 - 0.08333333333333333d0*velo**4
     &     )* dlog((-1 + dsqrt(2 - 1*velo**2))/(1 + velo))**3 +
     &     dlog(1 + velo)*(1 .3707783890401886d0 + 1.6449340668482262d0
     &     *velo**2 - 0 .8224670334241131d0*velo**4 + (0
     &     .4166666666666667d0 + 0.5d0*velo**2 - 0.25d0*velo**4)* dlog((
     &     -1 + dsqrt(2 - 1*velo**2))/(1 + velo))**2) + dlog(1/(1 -
     &     1*velo**2))* (0.6853891945200943d0 + 0.8224670334241131d0
     &     *velo**2 - 0.41123351671205655d0*velo**4 + (0
     &     .20833333333333334d0 + 0.25d0*velo**2 - 0.125d0*velo**4)
     &     * dlog((-1 + dsqrt(2 - 1*velo**2))/(1 + velo))**2) + (1
     &     .3707783890401886d0 + 1.6449340668482262d0*velo**2 - 0
     &     .8224670334241131d0*velo**4)* dlog(dsqrt(1/(1 - 1*velo**2))
     &     + dsqrt(1 + 1/(1 - 1*velo**2))) + (-0.2777777777777778d0 -
     &     0 .3333333333333333d0*velo**2 + 0.16666666666666666d0*velo**4
     &     ) * dlog(dsqrt(1/(1 - 1*velo**2)) + dsqrt(1 + 1/(1 - 1
     &     *velo**2)))** 3 + sp7*(3.388888888888889d0*dsqrt(2 - 1*velo
     &     **2) - 1.2777777777777777d0*velo**2*dsqrt(2 - 1*velo**2)) +
     &     dlog(dsqrt(1/(1 - 1*velo**2)) + dsqrt(1 + 1/(1 - 1*velo
     &     **2)))** 2*(3.388888888888889d0*dsqrt(2 - 1*velo**2) - 1
     &     .2777777777777777d0*velo**2*dsqrt(2 - 1*velo**2)) + sp3*((0
     &     .8333333333333334d0 + velo**2 - 0.5d0*velo**4)* dlog(1 - 1
     &     *velo) + (0.4166666666666667d0 + 0.5d0*velo**2 - 0.25d0*velo
     &     **4)* dlog(1 /(1 - 1*velo**2)) + 3.388888888888889d0
     &     *dsqrt(2 - 1*velo**2) - 1.2777777777777777d0*velo**2
     &     *dsqrt(2 - 1*velo**2)) + sp6 *((0.8333333333333334d0 + velo
     &     **2 - 0.5d0*velo**4)*dlog(1 + velo) + (0.4166666666666667d0
     &     + 0.5d0*velo**2 - 0.25d0*velo**4)* dlog(1 /(1 - 1*velo**2))
     &     + 3.388888888888889d0*dsqrt(2 - 1*velo**2) - 1
     &     .2777777777777777d0*velo**2*dsqrt(2 - 1*velo**2)) + sp5 *((
     &     -0.8333333333333334d0 - 1*velo**2 + 0.5d0*velo**4)* dlog(1
     &     + velo) + (-0.4166666666666667d0 - 0.5d0*velo**2 + 0.25d0
     &     *velo**4) *dlog(1/(1 - 1*velo**2)) + 3.388888888888889d0
     &     *dsqrt(2 - 1 *velo**2) - 1.2777777777777777d0*velo**2
     &     *dsqrt(2 - 1*velo**2)) + dlog((-1 + dsqrt(2 - 1*velo**2)
     &     )/(1 + velo))**2* (1 .6944444444444444d0*dsqrt(2 - 1*velo
     &     **2) - 0.6388888888888888d0 *velo**2*dsqrt(2 - 1*velo**2))
     &     + sp4*((0.8333333333333334d0 + velo**2 - 0.5d0*velo**4)
     &     * dlog(1 - 1*velo) + (0 .4166666666666667d0 + 0.5d0*velo**2
     &     - 0.25d0*velo**4)* dlog(1/(1 - 1*velo**2)) - 3
     &     .388888888888889d0*dsqrt(2 - 1*velo**2) + 1
     &     .2777777777777777d0*velo**2*dsqrt(2 - 1*velo**2)) + dlog((1
     &      + velo)/(1 - 1*velo))* (5.1967592592592595d0 - 2
     &     .9305555555555554d0*velo**2 + 0.4375d0*velo**4 + (1
     &     .1527777777777777d0 - 0.75d0*velo**2 + 0.125d0*velo**4)
     &     * dlog(4/(1  - 1*velo**2)) + dlog((-1 + dsqrt(2 - 1
     &     *velo**2))/ (1 + dsqrt(2 - 1*velo**2)))* (3
     &     .388888888888889d0*dsqrt(2 - 1 *velo**2) - 1
     &     .2777777777777777d0*velo**2*dsqrt(2 - 1*velo**2)))


c..
c..   mass terms from outer (light quark) loop:
c..   (m^2-terms are the same of on-shell and ms-bar masses)
c..
      rtmp = rtmp + mout**2/scms * ( -13/2.d0 + 3.d0*dlsmu )

      rv2cbl = cf*tr * rtmp

      return
      end

C-}}}
C-{{{ function rv2cbh:

c..   ------------------------------------------------------------

      function rv2cbh(scms,mout,min)
c..
c..   two-loop
c..   vector case
c..   double bubble contribution with inner and outer loop massive
c..   in the configuration
c..
c..   ****     4*mout^2 <  scms  <<  min^2     ****
c..   ****     scms < 4*thrq^2     ****
c..
c..   scms: cms energy squared
c..   mout: mass of quark coupling to external current
c..   min : mass of quark coupling to gluons only
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      if (.not.lpsup) then
         rv2cbh = 0.d0
         return
      endif

      if (mout.ge.min) then
         write(6,*) '<function rvdbh>: only valid for mout < min'
         write(6,*) '     mout = ',mout
         write(6,*) '     min  = ',min
         stop
      endif

      if (mout.eq.0.d0) then
         rv2cbh = r0dbh(scms,min)
         return
      endif

      if (lmsbar) then
         rv2cbh = rvdbtsms(scms,mout,min)
      else
         rv2cbh = rvdbts(scms,mout,min)
      endif

      return
      end

C-}}}
C-{{{ function rvdbts:

c..   ------------------------------------------------------------

      function rvdbts(scms,mout,min)
c..
c..   two-loop
c..   vector case
c..   double bubble contribution with inner and outer loop massive
c..   in the configuration
c..
c..   ****     4*mout^2 <  scms  <<  min^2     ****
c..
c..   exact coefficient of leading term in the limit of heavy 'min'
c..
c..   scms: cms energy squared
c..   mout: mass of quark coupling to external current
c..   min : mass of quark coupling to gluons only
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

c..   Note: it is assumed that 'min' is decoupled from \alpha_s

c..   analytical expression for mout = 0 expanded for scms/min^2 -> 0
c..   gives good approximation almost up to scms = 4*min^2
c..

c..   Expression from T.Seidensticker, private communication: leading term in 
c..   1/min^2 with exact dependence on mout
      rtmp = ((0.35555555555555557d0*mout**6)/(min**2*scms**2) + (0
     &     .26666666666666666d0*mout**4)/(min**2*scms) - (0
     &     .022222222222222223d0*scms)/min**2)* dlog((1 + dsqrt(1 - (4
     &     *mout**2)/scms))/ (1 - 1*dsqrt(1 - (4*mout**2)/scms))) +((0
     &     .3511111111111111d0*mout**2)/min**2 + (0 .17777777777777778d0
     &     *mout**4)/(min**2*scms) + (0 .09777777777777778d0*scms)/min
     &     **2) * dsqrt(1 - (4*mout**2) /scms) + ((0
     &     .044444444444444446d0 *mout**2)/min**2 + (0
     &     .022222222222222223d0*scms)/min**2) *dlog(min**2/mout**2)
     &     * dsqrt(1 - (4*mout**2)/scms)

      rvdbts = cf*tr * rtmp

      return
      end

C-}}}
C-{{{ function rvdbtsms:

c..   ------------------------------------------------------------

      function rvdbtsms(scms,mout,min)
c..
c..   Note: this is the MS-bar version of 'rvdbts'
c..
c..   two-loop
c..   vector case
c..   double bubble contribution with inner and outer loop massive
c..   in the configuration
c..
c..   ****     4*mout^2 <  scms  <<  min^2     ****
c..
c..   exact coefficient of leading term in the limit of heavy 'min'
c..
c..   scms: cms energy squared
c..   mout: MS-bar mass of quark coupling to external current
c..   min : mass of quark coupling to gluons only
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

c..   Note: it is assumed that 'min' is decoupled from \alpha_s

c..   Expression from T.Seidensticker, diploma thesis (after bug fixing): 
c..   leading term in 1/min^2 with exact dependence on mout
      rtmp = (3.04*mout**6)/(min**2*
     -     dsqrt(1. - (4.*mout**2)/scms)*scms**2) - 
     -   (3.2*mout**6*dlog(1/min))/
     -    (min**2*dsqrt(1. - (4.*mout**2)/scms)*scms**2) + 
     -   (1.6*mout**6*dlog(mout**(-2)))/
     -    (min**2*dsqrt(1. - (4.*mout**2)/scms)*scms**2) + 
     -   ((0.35555555555555557*mout**6)/(min**2*scms**2) + 
     -      (0.26666666666666666*mout**4)/(min**2*scms) - 
     -      (0.022222222222222223*scms)/min**2)*
     -    dlog((1. + dsqrt(1. - (4.*mout**2)/scms))/
     -      (1. - 1.*dsqrt(1. - (4.*mout**2)/scms))) + 
     -   ((0.3511111111111111*mout**2)/min**2 + 
     -      (0.17777777777777778*mout**4)/(min**2*scms) + 
     -      (0.09777777777777778*scms)/min**2)*
     -     dsqrt(1. - (4.*mout**2)/scms)
     -    + ((0.044444444444444446*mout**2)/min**2 + 
     -      (0.022222222222222223*scms)/min**2)*dlog(min**2/mout**2)*
     -    dsqrt(1. - (4.*mout**2)/scms)

      rvdbtsms = cf*tr * rtmp

      return
      end

C-}}}
C-{{{ function r0dbh:

c..   ------------------------------------------------------------

      function r0dbh(scms,min)
c..
c..   two-loop
c..   vector case
c..   outer loop massless
c..   double bubble power suppressed terms, i.e.   
c..
c..   ****     scms < 4*thrq^2     ****
c..
c..   analytical expression
c..
c..   NOTE: WE FIND HUGE CANCELLATIONS AMONG INDIVIDUAL TERMS
c..   IN RTMP! SEE THE VALUES FOR RTMP1, RTMP2, ...
c..   
c..   scms: cms energy squared
c..   min : mass of quark coupling to gluons only
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      complex*16 spence,csp7
      complex*16 trilog,ctr7
      
      include 'common.f'

c..   analytical expression for mout = 0
      arg7 = 
     &     ( -dsqrt(1.d0/(1.d0-(1.d0-4.d0*min**2/scms)))+
     &     dsqrt(1.d0+1.d0/(1.d0-(1.d0-4.d0*min**2/scms)))
     &     )**2
      csp7 = spence(arg7)
      ctr7 = trilog(arg7)
      sp7 = dreal(csp7)
      tr7 = dreal(ctr7)

      rtmp = 7.982167648374861d0- 3.9065840071353524d0*(1 - (4*min **2)
     &     /scms) + 0.3005142257898985d0*(1 - (4*min**2)/scms)**2 + (0
     &     .4166666666666667d0 + 0.5d0*(1 - (4*min**2)/scms) - 0 .25d0
     &     *(1 -(4*min**2)/scms)**2)*tr7 + (-7.962962962962963d0 + 3
     &     .0555555555555554d0*(1 - (4*min**2)/scms))* dlog(2 *dsqrt((0
     &     .25d0*scms)/min**2)) + (1.3707783890401886d0 + 1
     &     .6449340668482262d0*(1 - (4*min**2)/scms) - 0
     &     .8224670334241131d0*(1 - (4*min**2)/scms)**2) * dlog(dsqrt((0
     &     .25d0*scms)/min**2) + dsqrt(1 + (0.25d0*scms) /min**2)) + (-0
     &     .2777777777777778d0 + 0.16666666666666666d0*(1 - (4*min**2)
     &     /scms)**2 + 0.3333333333333333d0*(-1 + (4 *min **2)/scms))
     &     * dlog(dsqrt((0.25d0*scms)/min**2) + dsqrt(1 + (0.25d0*scms)
     &     /min**2))**3 - 5.574498782096766d0*dsqrt(1 + (4 *min**2)/scms
     &     ) +2.101860196528289d0*(1 - (4*min**2) /scms) * dsqrt(1 + (4
     &     *min**2)/scms) + sp7*(3 .388888888888889d0 *dsqrt(1 + (4*min
     &     **2)/scms) - 1 .2777777777777777d0*(1 - (4*min**2)/scms)
     &     * dsqrt(1 + (4 *min**2)/scms)) + dlog(dsqrt((0.25d0*scms)/min
     &     **2) + dsqrt(1 + (0.25d0*scms) /min**2))**2* (3
     &     .388888888888889d0*dsqrt(1 + (4*min**2) /scms) - 1
     &     .2777777777777777d0*(1 - (4*min**2) /scms) * dsqrt(1 + (4*min
     &     **2)/scms))

      r0dbh = cf*tr * rtmp

      return
      end

C-}}}
C-{{{ rvct real+virtual:

C-{{{ function rvctrv:

      function rvctrv(scms,mq)
c..
c..   Adds real and virtual contributions of 2-loop 
c..   double bubble with equal masses.
c..
c..   scms: cms energy squared
c..   mq:   quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      rvctrv = rvctvirt(scms,mq) + rvctreal(scms,mq)

      return
      end

C-}}}
C-{{{ function rvctvirt:

      function rvctvirt(scms,mq)
c..
c..   Virtual contributions to double bubble with equal masses.
c..   
c..   scms:  cms energy squared
c..   mq:    quark mass
c..   
      implicit real*8 (a-z)

      if (scms.le.4*mq*mq) then
         rvctvirt = 0.d0
         return
      endif

      pi= 4.d0*datan(1.d0)
      zeta2 = pi**2/6.d0

      zz = 2*mq/dsqrt(scms)
      zz2 = zz*zz

      vq = dsqrt(1-zz2)
      dlogpp = dlog((1-vq)/(1.d0+vq))

      rtmp = 1/6.d0*( (3 + 10*vq**2 - 5*vq**4)/24.d0*dlogpp**3
     &     + (-3+40*vq**2+16*vq**4-15*vq**6)/12.d0/vq**3*dlogpp**2
     &     + ( (-18+234*vq**2 + 167*vq**4 - 118*vq**6)/18.d0/vq**2
     &     + ( -3 - 10*vq**2 + 5*vq**4)/2.d0*zeta2)*dlogpp
     &     + ( -9+510*vq**2 - 118*vq**4)/9.d0/vq
     &     + vq*(-27 + 5*vq**2)*zeta2 )

      rvctvirt = rtmp

      end

C-}}}
C-{{{ function rvctreal:

      function rvctreal(scms,mq)
c..
c..   Real contributions to double-bubble with equal masses.
c..
c..   scms:  cms energy squared
c..   mq:    quark mass
c..
      implicit real*8 (a-z)
      integer itmx1,itmx2,itmx,ncall1,ncall2,ncall,ndim,nprn
      real avgt,errt,avgi,erri
      external fyzint
      common /dbpar/ zz,zz2
      common/vgres/avgt,errt,avgi,erri

      zz = 2*mq/dsqrt(scms)
      zz2 = zz*zz

      eps = 1.d-8
c..   gcc version 2.95  gives 'false' for the following, if eps=0
      if (2*zz+eps.ge.1.d0) then
         rvctreal = 0.d0
         return
      endif

      ndim = 2
      acc = 1.d-8       !  accuracy for vegas
      itmx1=5         !  interations for vegas run 1
      ncall1=2000      !  calls for vegas run 1
      itmx2=2         !  interations for vegas run 2
      ncall2=5000     !  calls for vegas run 2
      nprn=0          !  =1 -- verbose mode
      igraph = 0
      
      ncall = ncall1
      itmx = itmx1
      
      call vegas(fyzint,acc,ndim,ncall,itmx,nprn,igraph)

      ncall = ncall2
      itmx = itmx2

      call vegas1(fyzint,acc,ndim,ncall,itmx,nprn,igraph)

      rvctreal = avgt

      end

C-}}}
C-{{{ function fyzint:

      real*8 function fyzint(xx,wgt)
c..
c..   Integrand for real contributions to massive double bubble.
c..   
      implicit real*8 (a-z)
      real*8 xx(10)
      external fyz
      common /dbpar/ zz,zz2

      yp = xx(1)
      zp = xx(2)

c..   transformation of variables, so that integration goes from 0 to 1:
      y = (1-2*zz)*yp + zz2
      z = ( (1-dsqrt(y))**2 - zz2 )*zp + zz2

      if (((1-dsqrt(y))**2.le.zz2).or.(z.ge.(1-dsqrt(y))**2)) then
         fyzint = 0.d0
      else
c..   meas1: integration measure as given in Teubner's thesis.
c..   meas2: Jacobian from change of variable.
         meas1 = 1/z*(1 + zz2/2.d0/z) * dsqrt(1 - zz2/z)
         meas2 = (1 - 2*zz)*( (1-dsqrt(y))**2 - zz2 )
         fyzint = 1/3.d0 * meas1 * meas2 * fyz(y,z)
      endif
         
      end

C-}}}
C-{{{ function fyz:

      real*8 function fyz(y,z)
c..
c..   Used by fyzint.
c..
      implicit real*8 (a-z)
      common /dbpar/ zz,zz2

      dslam = dsqrt(1 + y**2 + z**2 - 2*(y + z + y*z))
      wurz = dsqrt(1-zz2/y)

      fyz = (zz2**2/2.d0 + zz2*(1-y+z) - (1-y+z)**2 - 2*(1+z)*y)/
     &     (1-y+z)*dlog(( 1-y+z-wurz*dslam )/(1-y+z + wurz*dslam ))
     &     - wurz*dslam*( 1 + (zz2**2 + 2*zz2 + 4*(1+zz2/2.d0)*z)/
     &     ( (1-y+z)**2 - wurz**2*dslam**2 ) )
      
      end

C-}}}
C-{{{ function rvctexp:

c..   ------------------------------------------------------------

      function rvctexp(scms,mq)
c..
c..   Double bubble with two equal masses.
c..   Uses expansion m^0, m^2, ..., m^12  above (4*mq)**2,
c..   and virtual terms below.
c..   This is faster than rvdbrv, where the real radiation
c..   is integrated numerically.
c..
c..   scms : cms energy squared
c..   mq   : on-shell quark mass 
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      m2s = mq*mq/scms

      dlmus = dlog(mu*mu/scms)
      dlms = dlog(m2s)

      if (scms.lt.(4*mq)**2) then
         rtmp = rvctvirt(scms,mq)
     &        + 1.d0/4.d0 * dlog(mq**2/mu**2) * rv1(scms,mq)
      else
         rtmp =  -1.375 - dlmus/4. + (-6.5 - 3*dlmus)*m2s + m2s**6*(
     &        -1422.0211419753086 + (346981*dlms)/540. + (18937*dlms**2)
     &        /18. + (313.86296296296297 + (3064*dlms)/9.)*dlmus - 1105
     &        *zeta2) + m2s**3*(1.3333333333333333 + (350*dlms)/9. - (76
     &        *dlms**2)/9. + (6.962962962962963 + (116*dlms)/9.)*dlmus +
     &        (352*zeta2)/9.) + m2s**5*(-8.082814814814816 - (1157*dlms)
     &        /75. + (4328*dlms**2)/45. + (91.40592592592593 + (4676
     &        *dlms)/45.)*dlmus + (3592*zeta2)/45.) + m2s**4*(70
     &        .25347222222223 + (1673*dlms)/72. - (25*dlms**2)/4. + (27
     &        .305555555555557 + (203*dlms)/6.)*dlmus + (533*zeta2)/6.)
     &        + m2s**2*(0.6666666666666666 + 12*dlms - 3*dlms**2 + (-2.5
     &        + 6*dlms)*dlmus + 18*zeta2 - 10*zeta3) + zeta3
      endif

      rvctexp = rtmp

      return
      end

C-}}}

C-}}}




