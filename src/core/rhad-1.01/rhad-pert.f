






      function Ru_pert(scms, mu_msbar, alpha_s, m_c, m_b, m_t)

       implicit real*8(a-h,m-z)
       implicit integer(i,j)
       implicit character*60(k)
       implicit logical(l)
       include 'common/common.f'

       real*8 :: Ru_pert

       call parameters(scms, mu_msbar, alpha_s, m_c, m_b, m_t)
       call init(scms)

       Ru_pert = ruqrk(scms)

       return
      end

      function Rd_pert(scms, mu_msbar, alpha_s, m_c, m_b, m_t)

       implicit real*8(a-h,m-z)
       implicit integer(i,j)
       implicit character*60(k)
       implicit logical(l)
       include 'common/common.f'

       real*8 :: Rd_pert

       call parameters(scms, mu_msbar, alpha_s, m_c, m_b, m_t)
       call init(scms)

       Rd_pert = rdqrk(scms)

       return
      end

      function Rs_pert(scms, mu_msbar, alpha_s, m_c, m_b, m_t)

       implicit real*8(a-h,m-z)
       implicit integer(i,j)
       implicit character*60(k)
       implicit logical(l)
       include 'common/common.f'

       real*8 :: Rs_pert

       call parameters(scms, mu_msbar, alpha_s, m_c, m_b, m_t)
       call init(scms)

       Rs_pert = rsqrk(scms)

       return
      end

      function Rc_pert(scms, mu_msbar, alpha_s, m_c, m_b, m_t)

       implicit real*8(a-h,m-z)
       implicit integer(i,j)
       implicit character*60(k)
       implicit logical(l)
       include 'common/common.f'

       real*8 :: Rc_pert

       call parameters(scms, mu_msbar, alpha_s, m_c, m_b, m_t)
       call init(scms)

       Rc_pert = rcqrk(scms)

       return
      end
