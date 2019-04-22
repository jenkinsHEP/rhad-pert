c..
c..   r34.f
c..
c..   This file is part of `rhad'.
c..
c..   3- and 4-loop results for R_had
c..
C-{{{ function rv3:

      function rv3(scms,mq)

      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      rv3 = rv3m0(scms) + rv3m2(scms,mq) + rv3m4(scms,mq)

      end


C-}}}
C-{{{ function delr03:

      function delr03(scms,mq)

      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      delr03 = delr03m2(scms,mq) + delr03m4(scms,mq)

      end


C-}}}
C-{{{ function rv3m0:

c..   ------------------------------------------------------------

      function rv3m0(scms)
c..
c..   three-loop massless
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlmus = dlog(mu*mu/scms)

      rv3m0 = 302.18402777777777d0 - (7847*inffin)/216.d0 + (151*inffin
     &     **2)/162.d0 +dlmus**2*(7.5625d0 - (11*inffin)/12.d0 + inffin
     &     **2/36.d0) - (121*zeta2)/8.d0+ (11*inffin*zeta2)/6.d0 -
     &     (inffin**2*zeta2)/18.d0 -(1103*zeta3)/4.d0 +(262*inffin*zeta3
     &     )/9.d0 - (19*inffin**2*zeta3)/27.d0 + dlmus*(90
     &     .02083333333333d0 - (785*inffin)/72.d0 + (11*inffin**2)/36.d0
     &     - (121*zeta3)/2.d0 + (22*inffin*zeta3)/3.d0 - (2*inffin**2
     &     *zeta3)/9.d0) + (275*zeta5)/6.d0- (25*inffin*zeta5)/9.d0


      return
      end

C-}}}
C-{{{ function rv3m2:

c..   ------------------------------------------------------------

      function rv3m2(scms,mq)
c..
c..   three-loop m^2 terms
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      if (.not.la3m2) then
         rv3m2 = 0.d0
         return
      endif

      m2s = mq*mq/scms

      dlmus = dlog(mu*mu/scms)
      dlms = dlog(m2s)

      if (lmsbar) then
         rv3m2 = m2s*(2442 - (4846*inffin)/27.d0 + (125*inffin**2)/54.d0
     &        + dlmus**2*(213 .75d0 - 17*inffin + inffin**2/3.d0) +
     &        dlmus*(1126.25d0 -(175*inffin)/2.d0+ (13 *inffin**2)/9.d0)
     &        - (855*zeta2)/2.d0 + 34*inffin *zeta2 - (2*inffin**2*zeta2
     &        ) /3.d0 + (490*zeta3)/3.d0 -(466*inffin *zeta3)/27.d0 -
     &        (5225*zeta5)/6.d0 + (1045*inffin*zeta5)/27.d0)
      else
         rv3m2 = m2s*(1862.5833333333333d0 + 378*dlms - 9*dlms**2 + (596
     &        .25d0 + 132*dlms)*dlmus + (363*dlmus**2)/4.d0 + inffin**2
     &        *(2 .314814814814815d0 + (13*dlmus)/9.d0 + dlmus**2/3.d0 -
     &        (2*zeta2) /3.d0) - (967*zeta2)/2.d0 - 16*dln2*zeta2 + (502
     &        *zeta3)/3.d0 - (5225*zeta5)/6.d0 + inffin*(-156
     &        .09259259259258d0 - (52*dlms)/3.d0 + 2*dlms**2 + (-64
     &        .83333333333333d0 - 8*dlms)*dlmus - 11*dlmus **2 + 42
     &        *zeta2 - (466*zeta3)/27.d0 + (1045*zeta5)/27.d0))
      endif

      return
      end

C-}}}
C-{{{ function rv3m4:

c..   ------------------------------------------------------------

      function rv3m4(scms,mq)
c..
c..   three-loop m^4 terms
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      if (.not.la3m4) then
         rv3m4 = 0.d0
         return
      endif

      m2s = mq*mq/scms

      dlmus = dlog(mu*mu/scms)
      dlms = dlog(m2s)

      if (lmsbar) then
         rv3m4 = m2s**2*(-2926.141203703704d0 + dlms**2*(13 - (2*inffin
     &        )/3.d0) +(130009 *inffin)/648.d0 - (463*inffin**2)/972.d0
     &        + dlmus**2*(-1463.625d0 +91*inffin - (7*inffin**2)/6.d0) +
     &        dlmus**3*(-256.5d0 + (46*inffin)/3.d0 - (2*inffin**2)/9.d0
     &        ) + (12099*zeta2)/4.d0- (574*inffin*zeta2)/3.d0 + (23
     &        *inffin**2*zeta2)/9.d0 + (64123*zeta3)/18.d0 - (1672
     &        *inffin*zeta3)/9.d0 + (28*inffin**2*zeta3)/27.d0 + dlms*(
     &        -218.16666666666666d0 + (199*inffin)/12.d0 -(5*inffin**2)
     &        /27.d0 - 22*zeta3 + (4*inffin*zeta3)/3.d0) + dlmus*(-3348
     &        .75d0 + (877 *inffin)/4.d0 - 2*inffin**2 + dlms*(-61.75d0
     &        +(16*inffin)/3.d0 - inffin**2/9.d0) +1539*zeta2 - 92
     &        *inffin*zeta2 + (4*inffin**2*zeta2)/3.d0 +1026*zeta3 -(124
     &        *inffin*zeta3)/3.d0 + (8*inffin**2*zeta3)/9.d0) + 10
     &        *inffin*zeta4 -(13285*zeta5)/9.d0 +(440*inffin*zeta5)/9.d0
     &        )
      else
         rv3m4 = m2s**2*(842.7314814814815d0 - (608*a4)/3.d0 - (591*dlms
     &        **2)/4.d0 +(15*dlms**3)/2.d0 + (75.625d0 - (363*dlms)/2.d0
     &        )*dlmus**2 -(76 *dln2**4)/9.d0 + (2564287*zeta2)/540.d0 -
     &        (4568*dln2*zeta2)/9.d0 -(128*dln2**2*zeta2)/3.d0 + (56257
     &        *zeta3)/18.d0 -(1439*zeta2 *zeta3)/3.d0 +inffin**2*(1
     &        .9444444444444444d0 + (13*dlms**2)/9.d0 - (2*dlms**3)/9.d0
     &        +(0.2777777777777778d0 - (2*dlms)/3.d0)*dlmus **2 +dlms*(
     &        -3.4814814814814814d0 - (8*zeta2)/3.d0) + (25*zeta2)/3.d0
     &        +dlmus*(1.2962962962962963d0 - 3*dlms + (2*dlms**2)/3.d0 +
     &        4*zeta2 +(8*zeta3)/9.d0) + (112*zeta3)/27.d0) +dlms*(-1845
     &        .3333333333333d0 + 564*zeta2 - 24*dln2*zeta2 +416*zeta3) +
     &        dlmus*(441.4166666666667d0 - (4033*dlms)/4.d0 - (165*dlms
     &        **2) /2.d0 +1199*zeta2 + 88*dln2*zeta2 + 572*zeta3) -
     &        (1565*zeta4)/6.d0 -(3770*zeta5)/3.d0 + inffin*(-97
     &        .27314814814815d0 + (64*a4)/9 .d0 - (157*dlms**2)/6.d0 -(2
     &        *dlms**3)/3.d0 + (-9.166666666666666d0 + 22*dlms)*dlmus**2
     &        +(8*dln2**4)/27.d0 - (3544*zeta2)/9.d0 -(176*dln2*zeta2)/9
     &        .d0 +(32*dln2**2*zeta2)/9.d0 +dlmus*(-52.19444444444444d0
     &        + (361*dlms)/3.d0 - 6*dlms**2 -(416*zeta2)/3.d0 - (16*dln2
     &        *zeta2)/3.d0 - (148*zeta3)/3.d0) -(2323*zeta3)/9.d0 +dlms
     &        *(201.58333333333334d0 + (44*zeta2)/3.d0 + (16*dln2*zeta2)
     &        /3.d0 +(28*zeta3)/3.d0) + (700*zeta4)/9.d0 + (440*zeta5)/9
     &        .d0))
      endif

      return
      end

C-}}}
C-{{{ function delr03m2:

c..   ------------------------------------------------------------

      function delr03m2(scms,mq)
c..
c..   three-loop m^2 terms, outer quark massless, on-shell scheme
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      if (.not.la3m2) then
         delr03m2 = 0.d0
         return
      endif

      m2s = mq*mq/scms

      dlmus = dlog(mu*mu/scms)
      dlms = dlog(m2s)
      
      delr03m2 = m2s*(-80 + (32*inffin)/9.d0 + 60*zeta3 - (8*inffin
     &     *zeta3)/3.d0)

      return
      end

C-}}}
C-{{{ function delr03m4:

c..   ------------------------------------------------------------

      function delr03m4(scms,mq)
c..
c..   three-loop m^4 terms, outer quark massless, on-shell scheme
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      if (.not.la3m4) then
         delr03m4 = 0.d0
         return
      endif

      m2s = mq*mq/scms

      dlmus = dlog(mu*mu/scms)
      dlms = dlog(m2s)
      
      if (lmsbar) then
         delr03m4 = m2s**2*(-67.40972222222223d0 + 2*dlms**2 + (457
     &        *inffin)/108.d0+ 15 *zeta2-(2*inffin*zeta2)/3.d0 + 25
     &        *zeta3 - (22*inffin*zeta3)/9.d0+dlmus*(39
     &        .166666666666664d0 + dlms*(-9.5d0 + inffin/3.d0) - (13
     &        *inffin)/9.d0 -38*zeta3 + (4*inffin*zeta3)/3.d0) +dlms*(3
     &        .5833333333333335d0- (13*inffin)/18.d0 - 22*zeta3 +(4
     &        *inffin*zeta3)/3.d0) + (50*zeta5)/3.d0)
      else
         delr03m4 = m2s**2*(-87.85416666666667d0 - 2*dlms**2 + (457
     &        *inffin)/108.d0 + 15 *zeta2 -(2*inffin*zeta2)/3.d0 + (139
     &        *zeta3)/3.d0- (22*inffin*zeta3)/9.d0 +dlms*(24.25d0 - (13
     &        *inffin)/18.d0 -38*zeta3 +(4*inffin*zeta3)/3.d0) +dlmus
     &        *(23.833333333333332d0+ dlms*(-5.5d0 + inffin/3.d0) - (13
     &        *inffin)/9.d0 -22*zeta3 +(4*inffin*zeta3)/3.d0)+ (50*zeta5
     &        )/3.d0)
      endif

      return
      end

C-}}}
C-{{{ function rv3sing:

c..   ------------------------------------------------------------

      function rv3sing(scms,m1,m2)
c..
c..   three-loop singlet
c..   
c..   vector case
c..
c..   scms : cms energy squared
c..   
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      m12s = m1*m1/scms
      m22s = m2*m2/scms

      rv3sing = 55/216.d0 - 5/9.d0*zeta3

      if (la3m4) then
         rv3sing = rv3sing 
     &        + m12s*m12s * ( -10/9.d0 + 25/3.d0*zeta3 )
     &        + m22s*m22s * ( -10/9.d0 + 25/3.d0*zeta3 )
      endif
         
c..   there are no ml^2 corrections

      rv3sing = rv3sing

      return
      end

C-}}}

C-{{{ function rv4:

      function rv4(scms,mq)

      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      rv4 = rv4m0(scms) + rv4m2(scms,mq)

      end


C-}}}
C-{{{ function rv4m0:

c..   ------------------------------------------------------------

      function rv4m0(scms)
c..
c..   four-loop massless
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      rv4m0 = 
     &        + r04(scms)
     &        + inffin * r04nl(scms)
     &        + inffin**2 * r04nl2(scms)
     &        + inffin**3 * r04nl3(scms)

      return
      end

C-}}}
C-{{{ function rv4m2:

c..   ------------------------------------------------------------

      function rv4m2(scms,mq)
c..
c..   four-loop m^2 terms
c..   only QCD!
c..
c..   scms : cms energy squared
c..   mq   : quark mass
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      if (.not.la4m2) then
         rv4m2 = 0.d0
         return
      endif

      m2s = mq*mq/scms

c..   fit result for nf^0 and nf^1 term
c      r24nf01 = 7114.45d0 - 1432.21d0 * inffin
c      r24nf01 = 7.11d3 - 1.43d3 * inffin

c..   exact result
      r24nf01 =   6225.83637573497986647983825872 - 
     -     1364.03887327264406555898166085*inffin
      
      dlmus = dlog(mu*mu/scms)
      dlms = dlog(m2s)

      if (lmsbar) then
         rv4m2 = 
     -        m2s*((204009*dlmus**2)/32. + (11685*dlmus**3)/16. 
     -        + r24nf01 + 
     -     inffin**3*(-1.3914609053497942 - (13*dlmus**2)/36. - 
     -        dlmus**3/18. + dlmus*(-1.1574074074074074 + zeta2/3.) + 
     -        (13*zeta2)/18.) + 
     -     dlmus*(28444.302083333332 
     -        - (35055*zeta2)/8. + (10045*zeta3)/6. - 
     -        (214225*zeta5)/24.) + 
     -     inffin**2*(189.18312757201647 + (1121*dlmus**2)/36. + 
     -        (143*dlmus**3)/36. - (1121*zeta2)/18. + (13*zeta3)/324. + 
     -        (53*zeta3**2)/9. + 
     -        dlmus*(122.01273148148148 - (143*zeta2)/6. 
     -        + (233*zeta3)/27. - 
     -           (1045*zeta5)/54.) - (3610*zeta5)/81.) + 
     -     inffin*((-19301*dlmus**2)/24. - (2249*dlmus**3)/24. + 
     -        dlmus*(-3471.042824074074 + (2249*zeta2)/4. - 
     -           (15043*zeta3)/54. + (44935*zeta5)/54.)))
      else
         rv4m2 =          
     -        m2s*(-11299.926126742352 
     -        + (178849*dlms)/24. - (2313*dlms**2)/8. + 
     -     (21*dlms**3)/2. + (2669.90625 + (1089*dlms)/2.)*dlmus**2 + 
     -     (3993*dlmus**3)/16. + r24nf01 + 
     -     inffin**3*(-1.3914609053497942 - (13*dlmus**2)/36. - 
     -        dlmus**3/18. + dlmus*(-1.1574074074074074 + zeta2/3.) + 
     -        (13*zeta2)/18.) - 1715.7734016237346*zeta2 - 
     -     620.1837578886431*dlms*zeta2 - (9655*zeta3)/18. + 
     -     (938*dlms*zeta3)/3. + (1439*zeta2*zeta3)/3. 
     -        + (1565*zeta4)/6. + 
     -     dlmus*(16839.03125 + (6849*dlms)/2. - (297*dlms**2)/4. - 
     -      4080.370427833913*zeta2 
     -        + (2761*zeta3)/2. - (57475*zeta5)/8.) + 
     -     inffin**2*(176.31635802469137 
     -        + (199*dlms)/18. - (13*dlms**2)/6. + 
     -        (2*dlms**3)/9. 
     -        + (23.75 + 2*dlms)*dlmus**2 + (11*dlmus**3)/4. - 
     -        (415*zeta2)/6. + (4*dlms*zeta2)/3. - (995*zeta3)/324. + 
     -        (53*zeta3**2)/9. + 
     -        dlmus*(105.13310185185185 + (26*dlms)/3. - dlms**2 - 
     -           (53*zeta2)/2. + (233*zeta3)/27. - (1045*zeta5)/54.) - 
     -        (3610*zeta5)/81.) 
     -        + (18925*zeta5)/9. - (5225*dlms*zeta5)/3. + 
     -     inffin*(929.842802692714 
     -        - (35899*dlms)/54. + (133*dlms**2)/2. - 
     -        (10*dlms**3)/3. + (-455.375 - 66*dlms)*dlmus**2 - 
     -        (363*dlmus**3)/8. + 245.51811785803167*zeta2 + 
     -        17.636548370346958*dlms*zeta2 + (10622*zeta3)/81. - 
     -       (1436*dlms*zeta3)/27. - (610*zeta4)/9. - (8360*zeta5)/81. + 
     -        (2090*dlms*zeta5)/27. + 
     -        dlmus*(-2476.3576388888887 - 370*dlms + 21*dlms**2 + 
     -           593.7951774444796*zeta2 - (4069*zeta3)/18. + 
     -           (13585*zeta5)/18.)))
      endif

      return
      end

C-}}}
C-{{{ function r04:

c..   ------------------------------------------------------------

      function r04(scms)
c..
c..   four-loop, nl=0
c..
c..   scms : cms energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlnmus = dlog(mu**2/scms)

      pi2 = pi*pi

c      r04 = -186.d0
c..   ms: 22Mar09 (exact result in numerical form):
      r04 = -156.608110914430590910388569942d0

c..   logarithms:
      r04 = r04 + (1331*dlnmus**3)/64. + dlnmus**2*(388.8671875 - (3993
     &     *zeta3)/16.) +dlnmus*(2709.2447916666665 - (1331*pi2)/64. -
     &     (38643*zeta3)/16. +(3025*zeta5)/8.)
      

      return
      end

C-}}}
C-{{{ function r04nl:

c..   ------------------------------------------------------------

      function r04nl(scms)
c..
c..   four-loop nl
c..
c..   scms : cms energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlnmus = dlog(mu**2/scms)

      pi2 = pi*pi

c      r04nl = 21.3d0
c..   ms: 22Mar09 (exact result in numerical form):
      r04nl = 18.7747650298249587374020688647d0

c..   logarithms:
      r04nl = r04nl + (-121*dlnmus**3)/32. + dlnmus**2*(-70 .71875 +
     &     (363*zeta3)/8.) +dlnmus*(-490.9401041666667 + (121*pi2) /32.
     &     + (9695*zeta3)/24. -(275*zeta5)/6.)

      return
      end

C-}}}
C-{{{ function r04nl2:

c..   ------------------------------------------------------------

      function r04nl2(scms)
c..
c..   four-loop nl^2
c..
c..   scms : cms energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlnmus = dlog(mu**2/scms)

      pi2 = pi*pi

      r04nl2 = 1045381/15552.d0 - 593/432.d0*pi2 - 40655/864.d0*zeta3
     &     + 11/12.d0*pi2*zeta3 + 5/6.d0*zeta3*zeta3
     &     - 260/27.d0*zeta5

c..   logarithms:
      r04nl2 = r04nl2 + (11*dlnmus**3)/48. + dlnmus**2*(4
     &     .118055555555555 - (11*zeta3)/4.) +dlnmus*(27.39959490740741
     &     - (11*pi2)/48. - (257*zeta3)/12. +(25*zeta5)/18.)

      return
      end

C-}}}
C-{{{ function r04nl3:

c..   ------------------------------------------------------------

      function r04nl3(scms)
c..
c..   four-loop nl^3
c..
c..   scms : cms energy squared
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)
      
      include 'common.f'

      dlnmus = dlog(mu**2/scms)

      pi2 = pi*pi

      r04nl3 = - 6131/5832.d0 + 11/432.d0*pi2 + 203/324.d0*zeta3
     &     - 1/54.d0*pi2*zeta3 + 5/18.d0*zeta5

c..   logarithms:
      r04nl3 = r04nl3 -dlnmus**3/216. + dlnmus**2*(-0.0763888888888889 +
     &     zeta3/18.) +dlnmus*(-0.4660493827160494 + pi2/216. + (19
     &     *zeta3)/54.)


      return
      end

C-}}}






