c..
c..   runal.f
c..
c..   This file is part of `rhad'.
c..
c..   routines for running and decoupling of \alpha_s
c..
C-{{{ function mms2mos:

c..     ------------------------------------------------------------

      subroutine mms2mos(mms,inf,inloop,mos)
c..
c..   converts MS-bar mass, mms(mms), to on-shell mass, mos
c..
c..   mms:    MS-bar m(m)
c..   inf:    number of active quark flavours
c..   inloop: loop-order in \alpha_s
c..   mos:    on-shell mass (output)
c..
c..   Warning: there might be problems if not as(Mz) is used as input
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if(inloop.eq.4)then
         write(6,*) '<function mms2mos>: ',inloop,
     &        ' loop relation not implemented. Use 3-loop relation'
      elseif(inloop.ge.5)then
         write(6,*) '<function mms2mos>: ',inloop,
     &        ' loop relation not implemented.'
         write(6,*) 'Stopped.'
         stop
      endif

c..   Compute \alphas^(inf)(mms)
c..   The quark masses mb and mt are decoupled for mu=mb and mu=mt, resp.

      mubmatch = massb
      mutmatch = masst
      infend = inf
      infact = infini
      alsact = alphasmz
      mu0act = mz
      if(infend.lt.5)then
         call runalpha(alsact/pi,mu0act,mubmatch,infact,iord,0,apiout)
         alsact = pi*apiout
         call decalphams(alsact,massb,mubmatch,infact-1,iord-1,-1,
     &        alsout)
         alsact = alsout
         infact = infact-1
         mu0act = mubmatch
         if(infend.lt.4)then
            write(6,*) '<function mms2mos>: case not implemented'
            stop
         endif
      endif
      if(infend.gt.5)then
         call runalpha(alsact/pi,mu0act,mutmatch,infact,iord,0,apiout)
         alsact = pi*apiout
         call decalphams(alsact,masst,mutmatch,infact-1,iord-1,+1,
     &        alsout)
         alsact = alsout
         infact = infact+1
         mu0act = mutmatch
      endif
      if(infend.gt.6)then
         write(6,*) '<function mms2mos>: case not implemented'
         stop
      endif

      call runalpha(alsact/pi,mu0act,mms,infact,iord,0,apiout)

      apilo = apiout

      mos = mms

      if(inloop.ge.1)then
         mos = mos + mms*(1.3333333333333333333d0*apilo)
      endif
      if(inloop.ge.2)then
         mos = mos +  mms*(apilo**2*(13.443396256931755664d0
     &        - 1.0413669111716310344d0*(-1. + inf)))
      endif
      if(inloop.ge.3)then
         mos = mos +  mms*(apilo**3*(190.59495524378951376d0
     &        + 0.652690749081543758d0*(-1. + inf)**2
     &        - 26.655131692105261133d0*(-1. + inf)))
      endif

      return
      end

C-}}}
C-{{{ function decalphams:

c..     ------------------------------------------------------------

      subroutine decalphams(als,massth,muth,inl,inloop,idir,alsout)
c..
c..   decoupling of heavy quarks from \alpha_s
c..   where the quark mass is defined in the MS-bar scheme
c..
c..   als:    \alpha_s (input)
c..   massth: mass of heavy quark
c..   muth:   scale at which quark should be decoupled
c..   inl:    number of massless quarks (only needed for inloop>=3)
c..   inloop: loop-order for decoupling
c..   idir:   direction: -1: \alphas^(nl+1) -> \alphas^(nl)
c..                      +1: \alphas^(nl) -> \alphas^(nl+1)
c..   alsout: \alpha_s (output)
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if(inloop.ge.4)then
            write(6,*) '<function decalphams>: ',inloop,
     &           ' decoupling not implemented.'
            write(6,*) 'Stopped.'
            stop
      endif

      alsout = 1.d0

c..   \alphas^(nl+1) -> \alphas^(nl)
      if(idir.eq.-1)then
         if(inloop.ge.1)then
            alsout = alsout+als/pi*
     &           (-1.d0/6.d0)*dlog(muth*muth/massth/massth)
         endif
         if(inloop.ge.2)then
            alsout = alsout+als/pi*als/pi*
     &           (11.d0/72.d0
     &           -11.d0/24.d0*dlog(muth*muth/massth/massth)
     &           +1.d0/36.d0*(dlog(muth*muth/massth/massth))**2)
         endif
         if(inloop.ge.3)then
            xlmm = dlog(muth*muth/massth/massth)
            alsout = alsout+als/pi*als/pi*als/pi*
     &           (564731d0/124416d0
     &           - (955d0*xlmm)/576d0
     &           + (53d0*xlmm**2)/576d0
     &           - xlmm**3/216d0
     &           + (-2633d0/31104d0
     &           +   (67d0*xlmm)/576d0 - xlmm**2/36d0)*inl
     &           - (82043d0*zeta3)/27648d0)
         endif
      endif

c..   \alphas^(nl) -> \alphas^(nl+1)
      if(idir.eq.1)then
         if(inloop.ge.1)then
            alsout = alsout+als/pi*
     &           (1.d0/6.d0)*dlog(muth*muth/massth/massth)
         endif
         if(inloop.ge.2)then
            alsout = alsout+als/pi*als/pi*
     &           (-11.d0/72.d0
     &           +11.d0/24.d0*dlog(muth*muth/massth/massth)
     &           +1.d0/36.d0*(dlog(muth*muth/massth/massth))**2)
         endif
         if(inloop.ge.3)then
            xlmm = dlog(muth*muth/massth/massth)
            alsout = alsout+als/pi*als/pi*als/pi*
     &           (-564731d0/124416d0
     &           + (2645d0*xlmm)/1728d0 + (167d0*xlmm**2)/576d0
     &           + xlmm**3/216d0
     &           + (2633d0/31104d0 - (67d0*xlmm)/576d0
     &           + xlmm**2/36d0)*inl
     &           + (82043d0*zeta3)/27648d0)
         endif
      endif

      alsout = als*alsout

      return
      end

C-}}}
C-{{{ function decalpha:

c..     ------------------------------------------------------------

      subroutine decalpha(als,massth,muth,inl,inloop,idir,alsout)
c..
c..   decoupling of heavy quarks from \alpha_s
c..
c..   als:    \alpha_s (input)
c..   massth: mass of heavy quark
c..   muth:   scale at which quark should be decoupled
c..   inl:    number of massless quarks (only needed for inloop>=3)
c..   inloop: loop-order for decoupling
c..   idir:   direction: -1: \alphas^(nl+1) -> \alphas^(nl)
c..                      +1: \alphas^(nl) -> \alphas^(nl+1)
c..   alsout: \alpha_s (output)
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if(inloop.ge.4)then
            write(6,*) '<function decalpha>: ',inloop,'-loop',
     &           ' decoupling not implemented.'
            write(6,*) 'Stopped.'
            stop
      endif

      alsout = 1.d0

c..   \alphas^(nl+1) -> \alphas^(nl)
      if(idir.eq.-1)then
         if(inloop.ge.1)then
            alsout = alsout+als/pi*
     &           (-1.d0/6.d0)*dlog(muth*muth/massth/massth)
         endif
         if(inloop.ge.2)then
            alsout = alsout+als/pi*als/pi*
     &           (-7.d0/24.d0
     &           -19.d0/24.d0*dlog(muth*muth/massth/massth)
     &           +1.d0/36.d0*(dlog(muth*muth/massth/massth))**2)
         endif
         if(inloop.ge.3)then
            xlmm = dlog(muth*muth/massth/massth)
            alsout = alsout+als/pi*als/pi*als/pi*
     &           (-58933d0/124416d0 - (8521d0*xlmm)/1728d0
     &           - (131d0*xlmm**2)/576d0 - xlmm**3/216d0
     &           + inl*(2479/31104d0 + (409d0*xlmm)/1728d0
     &           + zeta2/9d0) - (2d0*zeta2)/3d0
     &           - (2d0*zeta2*dlog(2.d0))/9d0
     &           - (80507d0*zeta3)/27648d0)
         endif
      endif

c..   \alphas^(nl) -> \alphas^(nl+1)
      if(idir.eq.1)then
         if(inloop.ge.1)then
            alsout = alsout+als/pi*
     &           (1.d0/6.d0)*dlog(muth*muth/massth/massth)
         endif
         if(inloop.ge.2)then
            alsout = alsout+als/pi*als/pi*
     &           (7.d0/24.d0
     &           +19.d0/24.d0*dlog(muth*muth/massth/massth)
     &           +1.d0/36.d0*(dlog(muth*muth/massth/massth))**2)
         endif
         if(inloop.ge.3)then
            xlmm = dlog(muth*muth/massth/massth)
            alsout = alsout+als/pi*als/pi*als/pi*
     &           (58933d0/124416d0 + (8941d0*xlmm)/1728d0
     &           + (511d0*xlmm**2)/576d0 + xlmm**3/216d0
     &           - (2479d0*inl)/31104d0 - (409d0*xlmm*inl)/1728d0
     &           + (2d0*zeta2)/3d0 - (inl*zeta2)/9d0
     &           + (2d0*zeta2*dlog(2.d0))/9d0
     &           + (80507d0*zeta3)/27648d0)
         endif
      endif

      alsout = als*alsout

      return
      end

C-}}}
C-{{{ function rundecalpha:

c..     ------------------------------------------------------------

      subroutine rundecalpha(als0,mu0,mu1,alsout)
c..
c..   running and decoupling for \alpha_s
c..
c..   als0:   \alpha_s (input)
c..   mu0:    scale of \alpha_s (input)
c..   mu1:    scale of \alpha_s (output)
c..   alsout: \alpha_s (output)
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if(infini.ne.5)then
            write(6,*) '<function rundecalpha>: '
            write(6,*) '       Initial nf has to be nf=5.'
            write(6,*) '       Stopped.'
            stop
      endif

      infact = infini
      alsact = als0
      mu0act = mu0
      if(inffin.lt.5)then
         call runalpha(alsact/pi,mu0act,mub,infact,iord,0,apiout)
         alsact = pi*apiout
         if (lmsbar) then
            call runalpha(alsact/pi,mub,massb,infact,iord,0,api5mb)
            call runmass(massb,api5mb,alsact/pi,infact,iord,mbmub)
            call decalphams(alsact,mbmub,mub,infact-1,iord-1,-1,alsout)
         else
            call decalpha(alsact,massb,mub,infact-1,iord-1,-1,alsout)
         endif
         alsact = alsout
         infact = infact-1
         mu0act = mub
         if(inffin.lt.4)then
            call runalpha(alsact/pi,mu0act,muc,infact,iord,0,apiout)
            alsact = pi*apiout
         if (lmsbar) then
            call runalpha(alsact/pi,muc,massc,infact,iord,0,api4mc)
            call runmass(massc,api4mc,alsact/pi,infact,iord,mcmuc)
            call decalphams(alsact,mcmuc,muc,infact-1,iord-1,-1,alsout)
         else
            call decalpha(alsact,massc,muc,infact-1,iord-1,-1,alsout)
         endif
            alsact = alsout
            infact = infact-1
            mu0act = muc
         endif
      endif
      if(inffin.gt.5)then
         call runalpha(alsact/pi,mu0act,mut,infact,iord,0,apiout)
         alsact = pi*apiout
         if (lmsbar) then
            call runalpha(alsact/pi,mut,masst,infact,iord,0,api6mt)
            call runmass(masst,api6mt,alsact/pi,infact,iord,mtmut)
            call decalphams(alsact,mtmut,mut,infact+1,iord-1,+1,alsout)
         else
            call decalpha(alsact,masst,mut,infact+1,iord-1,+1,alsout)
         endif
         alsact = alsout
         infact = infact+1
         mu0act = mut
      endif

      call runalpha(alsact/pi,mu0act,mu1,infact,iord,0,apiout)

      alsout = pi*apiout

      return
      end

C-}}}
C-{{{ function decmassms:

c..     ------------------------------------------------------------

      subroutine decmassms(mq,als,massth,muth,inl,inloop,idir,mqout)
c..
c..   decoupling of heavy quarks from m_q
c..
c..   mq:     m_q^(inl)(massth) (input)
c..   als:    \alpha_s
c..   massth: MS-bar mass of heavy quark at scale muth
c..   muth:   scale at which quark should be decoupled
c..   inl:    number of massless quarks (only needed for inloop>=3)
c..   inloop: loop-order for decoupling
c..   idir:   direction: -1: \alphas^(nl+1) -> \alphas^(nl)
c..                      +1: \alphas^(nl) -> \alphas^(nl+1)
c..   mqout:  m_q^(inl+1)(massth) (output)
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if(inloop.ge.4)then
         write(6,*) '<function decmass>: ',inloop,'-loop',
     &        ' decoupling not implemented.'
         write(6,*) 'Stopped.'
         stop
      endif

      mqout = 1.d0

c..   m_q^(nl+1) -> m_q^(nl)
      if(idir.eq.-1)then
         write(6,*) '<function decmass>: ',
     &        ' this direction not implemented.'
         write(6,*) 'Stopped.'
         stop
      endif

c..   m_q^(nl) -> m_q^(nl+1)
      if(idir.eq.1)then
         if(inloop.ge.1)then
            mqout = mqout
         endif
         if(inloop.ge.2)then
            xlmm = dlog(muth*muth/massth/massth)
            mqout = mqout+als/pi*als/pi*
     &           (-89.d0/432.d0 + (5.d0*xlmm)/36.d0
     &           - xlmm**2/12.d0)
         endif
         if(inloop.ge.3)then
            mqout = mqout+als/pi*als/pi*als/pi*
     &           (-2951.d0/2916 + b4/36.d0
     &           - (155.d0*xlmm**2)/432.d0
     &           - (35.d0*xlmm**3)/216.d0
     &           + inl*(-1327.d0/11664.d0 + (53.d0*xlmm)/432.d0
     &           + xlmm**3/108.d0 + (2.d0*zeta3)/27.d0)
     &           + xlmm*(133.d0/2592.d0 + (5.d0*zeta3)/6.d0)
     &           + (407.d0*zeta3)/864.d0 - (5.d0*zeta4)/4.d0)
         endif
      endif

      mqout = mq*mqout

      return
      end

C-}}}
C-{{{ function rundecmass:

c..   ------------------------------------------------------------

      subroutine rundecmass(mq0,inf,mu0,mu1,mqout)
c..
c..   running and decoupling for m_q
c..
c..   Note: this routine is used in order to compute mc(s) with
c..   nf=5 and nf=6 and mb(s) with nf=6 from their respective input
c..   values mc(mc) and mb(mb).
c..
c..   mq0:   m_q^(inf) (input)
c..   inf:   number of active flavours (initially)
c..   mu0:   scale of m_q (initial)
c..   mu1:   scale of m_q (final)
c..   mqout: m_q (output)
c..
      implicit	real*8(a-h,m-z)
      implicit	integer(i,j)
      implicit	character*60(k)
      implicit  logical(l)

      include 'common/common.f'

      if(infini.ne.5)then
            write(6,*) '<function rundecmass>: '
            write(6,*) '       Initial nf has to be nf=5.'
            write(6,*) '       Stopped.'
            stop
      endif

c..   implement only the cases needed for rhad.f
      if(inffin.lt.inf)then
         write(6,*) '<function rundecmass>: nf=',inffin,
     &        ' <',inf,': not allowd.'
         write(6,*) 'Stopped.'
         stop
      endif
      if(.not.((inf.eq.4).or.(inf.eq.5).or.(inf.eq.6))) then
         write(6,*) '<function rundecmass>: ',
     &        'case not implemented.'
         write(6,*) 'Stopped.'
         stop
      endif
      if(.not.((inffin.eq.4).or.(inffin.eq.5).or.(inffin.eq.6))) then
         write(6,*) '<function rundecmass>: ',
     &        'case not implemented.'
         write(6,*) 'Stopped.'
         stop
      endif

c..   initialize alphas
      alsini = alphasmz
      muini  = mz
      inf0   = 5

      infact = inf
      mqact  = mq0
      mu0act = mu0
c..   in case no flavour threshold is crossed:
      call rundecalpha(alsini,muini,mu0act,alsmu0act)
      call rundecalpha(alsini,muini,mu1,alsmu1)
c..   nf=4 -> nf=5
      if((infact.eq.4).and.(inffin.ge.5))then
c..   determine as^(4)(mu0act), as^(4)(mub), as(5)(mub) from as^(5)(Mz)
         call runalpha(alsini/pi,muini,mub,inf0,iord,0,apiout)
         alsact = pi*apiout
         if (lmsbar) then
            call decalphams(alsact,massb,mub,inf0-1,iord-1,-1,alsout)
         else
            call decalpha(alsact,massb,mub,inf0-1,iord-1,-1,alsout)
         endif
         als5mub = alsact
         als4mub = alsout
         call runalpha(als4mub/pi,mub,mu0act,inf0-1,iord,0,apiout)
         als4mu0act = pi*apiout
c     print*,als5mub,als4mub,als4mu0act
c..   determine mq^(4)(mub)
         call runmass(mqact,als4mu0act/pi,als4mub/pi,infact,iord,mqout)
         mqact = mqout
c..   "undo" decoupling of bottom quark
         if (lmsbar) then
            call runalpha(als5mub/pi,mub,massb,inf0,iord,0,api5mb)
            call runmass(massb,api5mb,als5mub/pi,infact,iord,mbmub)
            call decmassms(mqact,als4mu0act,mbmub,mub,infact,iord-1,1,
     &           mqout)
         else
            write(6,*) '<function rundecmass>: ',
     &           'case not implemented.'
            write(6,*) 'Stopped.'
            stop
         endif
         infact    = infact+1
         mqact     = mqout
         mu0act    = mub
         alsmu0act = als5mub
c..   determine as^(infact)(mu1)
         call runalpha(alsmu0act/pi,mu0act,mu1,infact,iord,0,apiout)
         alsmu1 = pi*apiout
      endif

c..   nf=5 -> nf=6
      if((infact.eq.5).and.(inffin.eq.6))then
c..   determine as^(5)(mu0act), as^(5)(mut), as(6)(mut) from as^(5)(Mz)
         call runalpha(alsini/pi,muini,mut,inf0,iord,0,apiout)
         alsact = pi*apiout
         if (lmsbar) then
            call decalphams(alsact,masst,mut,inf0,iord-1,+1,alsout)
         else
            call decalpha(alsact,masst,mut,inf0,iord-1,+1,alsout)
         endif
         als5mut = alsact
         als6mut = alsout
         call runalpha(alsini/pi,muini,mu0act,inf0,iord,0,apiout)
         als5mu0act = pi*apiout
c     print*,als6mut,als5mut,als5mu0act
c..   determine mq^(5)(mut)
         call runmass(mqact,als5mu0act/pi,als5mut/pi,infact,iord,mqout)
         mqact = mqout
c..   "undo" decoupling of top quark
         if (lmsbar) then
            call runalpha(als6mut/pi,mut,masst,inf0,iord,0,api6mt)
            call runmass(masst,api6mt,als6mut/pi,infact,iord,mtmut)
            call decmassms(mqact,als5mu0act,mtmut,mut,infact,iord-1,1,
     &           mqout)
         else
            write(6,*) '<function rundecmass>: ',
     &           'case not implemented.'
            write(6,*) 'Stopped.'
            stop
         endif
         mqact     = mqout
         infact    = infact+1
         mu0act    = mut
         alsmu0act = als6mut
c..   determine as^(infact)(mu1)
         call runalpha(alsmu0act/pi,mu0act,mu1,infact,iord,0,apiout)
         alsmu1 = pi*apiout
      endif
c..   determine mq^(infact)(mu1)
      call runmass(mqact,alsmu0act/pi,alsmu1/pi,infact,iord,mqout)

      return
      end

C-}}}
C-{{{ routines for running of alphas:

C-{{{ subroutine odeint:

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
c..(C) Copr. 1986-92 Numerical Recipes Software 5,".
c..   transscribed to real*8 by R. Harlander, Feb.2002
      implicit real*8 (a-z)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(dabs(x-xsav).gt.dabs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(dabs(hnext).lt.hmin) then
           write(6,*) 'stepsize smaller than minimum in odeint'
           read(5,*)
        endif
        h=hnext
16    continue
      write(6,*) 'too many steps in odeint'
      read(5,*)
      return
      END

C-}}}
C-{{{ subroutine rkck:

      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
c..(C) Copr. 1986-92 Numerical Recipes Software 5,".
c..   transscribed to real*8 by R. Harlander, Feb.2002
      implicit real*8 (a-z)
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END

C-}}}
C-{{{ subroutine rkqs:

      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c..(C) Copr. 1986-92 Numerical Recipes Software 5,".
c..   transscribed to real*8 by R. Harlander, Feb.2002
      implicit real*8 (a-z)
      INTEGER n,NMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,dabs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(dabs(htemp),0.1*dabs(h)),h)
        xnew=x+h
        if(xnew.eq.x) then
           write(6,*) 'stepsize underflow in rkqs'
           read(5,*)
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END

C-}}}
C-{{{ subroutine runalpha:

      subroutine runalpha(api0,mu0,mu,nf,nloop,verb,apiout)
C..
c..   NEEDS:  rkck.f rkqs.f odeint.f  (from Numerical Recipes)
c..
c..   Note:  api = {\alpha_s \over \pi}
C..
c..   purpose : computes the value of api(mu) from api(mu0)
c..   method  : solving RG-equation by adaptive Runge-Kutta method
c..   uses    : odeint.for  from Numerical Recipes
C..
c..   api0  :  api(mu0)
c..   nf    :  number of flavors
c..   nloop :  number of loops
c..   verb  :  0=quiet,  1=verbose
c..   apiout:  api(mu)
C..
      implicit real*8 (a-h,o-z)
      INTEGER KMAXX,NMAX,NVAR
      PARAMETER (KMAXX=200,NMAX=50,NVAR=1)
      INTEGER kmax,kount,nbad,nok,nrhs,nf,verb
      REAL*8 dxsav,eps,h1,hmin,x,y,apif(NVAR),api0,apiout,pi
      real*8 mu,mu0,l0,lf
c..   /path/  is for odeint.for:
      COMMON /path/ kmax,kount,dxsav,x(KMAXX),y(NMAX,KMAXX)
      common /bfunc/ beta0,beta1,beta2,beta3
      COMMON nrhs
      data pi/3.14159265358979323846264338328d0/
      EXTERNAL rhs,rkqs

      if (nloop.eq.0) then
         apiout = api0
         return
      endif

      nrhs=0

c..   integration bounds (note that log(mu^2) is the integration variable)
      l0 = 0.d0
      lf = 2.*dlog(mu/mu0)
      apif(1)=api0

c..   see documentation for odeint.for:
      eps=1.0d-8
      h1=dabs(lf-l0)/10.d0
      hmin=0.0d0
      kmax=100
      dxsav=dabs(lf-l0)/20.d0

c..   initialize beta-function (common block /bfunc/):
      call inibeta(nf,nloop)

c..   check if input values are reasonable
      dlam = mu0*dexp(-1.d0/(2.d0*beta0*api0))
      if (mu.le.dlam) then
         write(6,2001) dlam,mu,mu0,api0*pi
      endif

c..   integrate RG-equation:

      call odeint(apif,NVAR,l0,lf,eps,h1,hmin,nok,nbad,rhs,rkqs)

      if (verb.eq.1) then
         write(6,'(/1x,a,t30,i3)') 'Successful steps:',nok
         write(6,'(1x,a,t30,i3)') 'Bad steps:',nbad
         write(6,'(1x,a,t30,i3)') 'Function evaluations:',nrhs
         write(6,'(1x,a,t30,i3)') 'Stored intermediate values:',kount
      endif

c..   api(mu):
      apiout = apif(1)

 2001 format(' -> <subroutine runalpha>',/,
     &     ' - WARNING: mu-value too low.',/,
     &     ' -     Should be significantly larger than  ',1f8.3,'.',/,
     &     ' -             mu = ',1f8.3,' GeV',/,
     &     ' -            mu0 = ',1f8.3,' GeV',/,
     &     ' -        api0*pi = ',1f8.3,/,
     &     ' -     Integration might break down.',/,
     &     '<- <subroutine runalpha>'
     &     )

      END

C-}}}
C-{{{ subroutine rhs:

      subroutine rhs(logmumu0,api,ainteg)
C..
c..   RG-equation:   (d api)/(d log(mu^2)) = api*beta(api)
C..
      implicit real*8 (a-h,o-z)
      integer nrhs
      real*8 api(*),ainteg(*),logmumu0
      common /bfunc/ beta0,beta1,beta2,beta3
      COMMON nrhs
      nrhs=nrhs+1
      ainteg(1) = api(1)*(- beta0*api(1) - beta1*api(1)**2 - beta2
     &     *api(1)**3 - beta3*api(1)**4)
      end

C-}}}
C-{{{ subroutine inibeta:

      subroutine inibeta(nf,nloopin)
C..
c..   initialize beta function
C..
      implicit real*8 (a-h,o-z)
      data z2/1.6449340668482264364724/,
     &     z3/1.2020569031595942853997/
      common /bfunc/ beta0,beta1,beta2,beta3

      beta0 = (33 - 2*nf)/12.d0
      beta1 = (102 - (38*nf)/3.d0)/16.d0
      beta2 = (2857/2.d0 - (5033*nf)/18.d0 + (325*nf**2)/54.d0)/64.d0
      beta3 = (149753/6.d0 + (1093*nf**3)/729.d0 + 3564*z3 + nf**2
     &     *(50065/162.d0 + (6472*z3)/81.d0) - nf*(1078361/162.d0 +
     &     (6508*z3)/27.d0))/256.d0

      nloop=nloopin

      if (nloop.gt.4) then
         write(6,*) '-> <subroutine inibeta>:'
         write(6,*)
     &        ' - 5-loop beta function unknown. Using 4-loop instead.'
         write(6,*) '<- <subroutine inibeta>'
         nloop=4
      endif
      if (nloop.lt.4) then
         beta3 = 0d0
         if (nloop.lt.3) then
            beta2 = 0d0
            if (nloop.lt.2) then
               beta1 = 0d0
               if (nloop.lt.1) then
                  beta0=0d0
               endif
            endif
         endif
      endif
      end

C-}}}

C-}}}
C-{{{ subroutine runmass:

      subroutine runmass(mass0,api0,apif,nf,nloop,massout)
c..
c..   evaluates the running of the MS-bar quark mass
c..   by expanding the equation
c..
c..   m(mu) = m(mu0) * exp( \int_a0^af dx gammam(x)/x/beta(x) )
c..
c..   in terms of alpha_s. The results agree with RunDec.m.
c..
c..
c..   Input:
c..   ------
c..   mass0  :  m(mu0)
c..   api0   :  alpha_s(mu0)/pi
c..   apif   :  alpha_s(muf)/pi
c..   nf     :  number of flavors
c..   nloop  :  order of calculation (nloop=1..4)
c..
c..   Output:
c..   -------
c..   massout:  m(muf)
c..
      implicit real*8 (a-h,o-z)
      real*8 mass0,massout,massfun
      external massfun
      parameter(accmass=1.d-6)
      common /bfunc/ beta0,beta1,beta2,beta3
      common /gfunc/ gamma0,gamma1,gamma2,gamma3

      if (nloop.eq.0) then
         massout = mass0
         return
      endif

      call inigamma(nf,nloop)
      call inibeta(nf,nloop)

      bb1 = beta1/beta0
      bb2 = beta2/beta0
      bb3 = beta3/beta0

      cc0 = gamma0/beta0
      cc1 = gamma1/beta0
      cc2 = gamma2/beta0
      cc3 = gamma3/beta0

      cfunc1 = 1.d0
      cfunc2 = cc1 - bb1*cc0
      cfunc3 = 1/2.d0*((cc1-bb1*cc0)**2 + cc2 - bb1*cc1 + bb1**2*cc0 -
     &     bb2*cc0)
      cfunc4 = (1/6*(cc1 - bb1*cc0)**3 + 1/2*(cc1 - bb1*cc0)*(cc2 - bb1
     &     *cc1 + bb1**2*cc0 - bb2*cc0) + 1/3*(cc3 - bb1*cc2 + bb1**2
     &     *cc1 - bb2*cc1 - bb1**3*cc0 + 2*bb1*bb2*cc0 - bb3*cc0))

      if (nloop.lt.4) then
         cfunc4 = 0.d0
         if (nloop.lt.3) then
            cfunc3 = 0.d0
            if (nloop.lt.2) then
               cfunc2 = 0.d0
               if (nloop.lt.1) then
                  cfunc1 = 0.d0
               endif
            endif
         endif
      endif

      cfuncmu0 = cfunc1 + cfunc2*api0 + cfunc3*api0**2 + cfunc4*api0**3
      cfuncmuf = cfunc1 + cfunc2*apif + cfunc3*apif**2 + cfunc4*apif**3

      massout = mass0*(apif/api0)**cc0*cfuncmuf/cfuncmu0

      return
      end

C-}}}
C-{{{ subroutine inigamma:

      subroutine inigamma(nfin,nloopin)
C
C     initialize beta function
C
      implicit real*8 (a-h,o-z)
      data z2/1.6449340668482264364724/,
     &     z3/1.2020569031595942853997/,
     &     z5/1.0369277551433699263/,
     &     pi/3.1415926535897932381/
      common /gfunc/ gamma0,gamma1,gamma2,gamma3

      nf = nfin

      gamma0 = 1.d0
      gamma1 = (67.33333333333333d0 - (20*nf)/9.d0)/16.d0
      gamma2 = (1249.d0 - (140*nf**2)/81.d0 + 2*nf*(-20.59259259259259d0
     &     - 48*z3) +(8*nf*(-46 + 48*z3))/9.d0)/64.d0
      gamma3 = (28413.91975308642d0 + (135680*z3)/27.d0 + nf**3*(-1
     &     .3662551440329218d0 + (64*z3)/27.d0) + nf**2*(21
     &     .57201646090535d0 - (16*Pi**4)/27.d0 + (800*z3)/9.d0) - 8800
     &     *z5 + nf*(-3397.1481481481483d0 + (88*Pi**4)/9.d0 - (34192
     &     *z3)/9.d0 + (18400*z5)/9.d0))/256.d0

      nloop=nloopin

      if (nloop.gt.4) then
         write(6,*) '-> <subroutine inigamma>:'
         write(6,*)
     &        ' - 5-loop gamma function unknown. Using 4-loop instead.'
         write(6,*) '<- <subroutine inigamma>'
         nloop=4
      endif
      if (nloop.lt.4) then
         gamma3 = 0d0
         if (nloop.lt.3) then
            gamma2 = 0d0
            if (nloop.lt.2) then
               gamma1 = 0d0
               if (nloop.lt.1) then
                  gamma0 = 0d0
               endif
            endif
         endif
      endif
      end

C-}}}
