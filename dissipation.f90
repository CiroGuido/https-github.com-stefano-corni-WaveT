module dissipation 
  use readio
  use random

  use, intrinsic :: iso_c_binding


  implicit none
  real(8)                :: norm, dtot, dsp, dnr, dde 
  complex(16), parameter  :: zeroc=(zero,zero)
  save 
  private
  public norm, dtot, dsp, dnr, dde, add_dis_m, add_dis_nm, loss_norm
  public quan_jump, add_h_rnd, define_h_dis, rnd_noise, add_h_rnd2 
!
  contains

  subroutine add_dis_m(h_dis,nci)
!------------------------------------------------------------------------
! Markovian SSE (eq 25 J. Phys: Condens. Matter vol. 24 (2012) 273201)
! Add the dissipative contribution to H 
! -i/2 \sum_alpha S^dag_alpha S_alpha
! Representation in the CIS basis
!
! Created   : E. Coccia 20 Dec 2016
! Modified  :
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)    :: nci
   real(8), intent(inout) :: h_dis(nci,nci)    
   integer                :: i
   real(8)                :: rate

! Matrix elements of S^+_alpha S_alpha in the system eigenstates basis

   do i=2, nci
! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
      h_dis(i,i) = h_dis(i,i) + sp_gam(i-1)*tmom2_0i(i-1)
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha|
      if (nr_typ.eq.0) then
         rate = nr_gam(i-1)*tmom2_0i(i-1)
      elseif (nr_typ.eq.1) then
         rate = nr_gam(i-1)
      endif
      h_dis(i,i) = h_dis(i,i) + rate 
! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
      h_dis(i,i) = h_dis(i,i) + de_gam(i)
   enddo
   if (idep.eq.0) then
      h_dis(1,1) = de_gam(1) !+ 1.d0
   elseif (idep.eq.1) then
      do i=2,nci
         h_dis(1,1) = h_dis(1,1) + de_gam(i-1)
      enddo
   endif
   h_dis=0.5d0*h_dis

   return

  end subroutine add_dis_m


  subroutine add_dis_nm(h_dis,nci)
!------------------------------------------------------------------------
! Non-Markovian SSE (eq 25 J. Phys: Condens. Matter vol. 24 (2012) 273201)
! Add the dissipative contribution to H 
! 
!
! Created   : E. Coccia 20 Dec 2016
! Modified  :
!------------------------------------------------------------------------

    implicit none
    integer, intent(in)    :: nci
    real(8), intent(inout) :: h_dis(nci,nci)

    write(*,*)
    write(*,*) ' NONMARKOVIAN SSE'
    write(*,*) 'NOT IMPLEMENTED YET'
    write(*,*)

    stop

  end subroutine add_dis_nm

  subroutine loss_norm(c,nci,pjump)
!------------------------------------------------------------------------
! Contributions to the loss of norm 
! Quantum jump from J. Opt. Soc. Am. B. vol. 10 (1993) 524 
! norm = 1 - dtot 
! dtot = dsp + dnr + dde
! dsp -> spontaneous relaxation to the ground state
! dnr -> nonradiative relaxation to the ground state
! dde -> pure dephasing
! 
! Created   : E. Coccia 20 Dec 2016
! Modified  :
!------------------------------------------------------------------------

   implicit none
   complex(16), intent(in)   :: c(nci)
   integer,     intent(in)   :: nci
   real(8),     intent(inout):: pjump(3*nexc+1)
   integer                   :: i
   real(8)                   :: weight, tmp

   norm = dot_product(c,c)

   dsp=0.d0
   dnr=0.d0
   dde=0.d0

   do i=1,nexc
      tmp=abs(c(i+1))
      weight=tmom2_0i(i)*tmp**2
      dsp = dsp + sp_gam(i)*weight
      if (nr_typ.eq.0) then
         dnr = dnr + nr_gam(i)*weight
      elseif (nr_typ.eq.1) then
         dnr = dnr + nr_gam(i)*tmp**2
      endif   
      !dde = dde + de_gam(i+1)*tmp**2
      pjump(i) = sp_gam(i)*weight
      if (nr_typ.eq.0) then
         pjump(i+nexc) = nr_gam(i)*weight
      elseif (nr_typ.eq.1) then
         pjump(i+nexc) = nr_gam(i)*tmp**2
      endif
      !pjump(i+1+2*nexc) = de_gam(i+1)*tmp**2
   enddo
   if (idep.eq.0) then
      pjump(1+2*nexc) = de_gam(1)*abs(c(1))**2 
      dde = dde + pjump(1+2*nexc)
      do i=1,nexc
        tmp=abs(c(i+1))
        pjump(i+1+2*nexc) = de_gam(i+1)*tmp**2
        dde = dde + pjump(i+1+2*nexc)
      enddo
   elseif (idep.eq.1) then
      do i=1,nexc
         tmp=abs(c(i+1))
         pjump(i+2*nexc) = de_gam(i)*(tmp**2 + abs(c(1))**2)
         dde = dde + pjump(i+2*nexc)
      enddo
      
   endif
   dsp = dsp*dt
   dnr = dnr*dt
   dde = (dde+pjump(1+2*nexc))*dt
   dtot = dsp + dnr + dde

   pjump(1:nexc)=pjump(1:nexc)*dt/dsp
   pjump(nexc+1:2*nexc)=pjump(nexc+1:2*nexc)*dt/dnr
   if (idep.eq.0) then
      pjump(2*nexc+1:3*nexc+1)=pjump(2*nexc+1:3*nexc+1)*dt/dde
   elseif (idep.eq.1) then
      pjump(2*nexc+1:3*nexc)=pjump(2*nexc+1:3*nexc)*dt/dde
   endif

   return

  end subroutine loss_norm 

  subroutine quan_jump(c,nci,pjump)
!------------------------------------------------------------------------
! Quantum jump from J. Opt. Soc. Am. B. vol. 10 (1993) 524 
! Random events: dissipation, nonradiative and dephasing
!
! Created   : E. Coccia 20 Dec 2016
! Modified  :
!------------------------------------------------------------------------

   implicit none
   complex(16), intent(inout)   :: c(nci)
   integer(4),  intent(in)      :: nci
   real(8),     intent(in)      :: pjump(3*nexc+1)
   integer(4)                   :: i,istate
   real(8)                      :: eta, eta1, tmp1, tmp2, tmp3, modc, creal, ireal
   real(8)                      :: left, right, ph 

! (0)|------dsp/dtot-----|---dnr/dtot---|--dde/dtot--|(1) 
! Select the type of event according to
! a kinetic Monte Carlo strategy (JCP vol. 111 (1999) 10126)   

   call random_number(eta)

   tmp1 = dsp/dtot
   tmp2 = (dsp+dnr)/dtot
   tmp3 = (dsp+dnr+dde)/dtot

! Select the relaxation channel
! j = n + FLOOR((m+1-n)*rnd), rnd [0,1) -> [n,m] 
! In our case [1,nexc]


! Spontaneous occurring 
   if (eta.ge.0.and.eta.lt.tmp1) then
      call random_number(eta1)
      left=0.d0
      right=pjump(1)
      do i=1,nexc-1
         if (eta1.ge.left.and.eta1.lt.right) then
            istate=i
            exit
         endif 
         left  = right
         right = left + pjump(i+1)
      enddo
      creal = real(c(istate+1))*sqrt(sp_gam(istate)*tmom2_0i(istate))
      ireal = aimag(c(istate+1))*sqrt(sp_gam(istate)*tmom2_0i(istate))
      c(1)  = cmplx(creal,ireal) 
      c(2:nci) = zeroc
      c(1)=c(1)/sqrt(pjump(istate)/dt)
      i_sp=i_sp+1 
! Nonradiative occurring
   elseif (eta.ge.tmp1.and.eta.lt.tmp2) then
      call random_number(eta1)
      left=0.d0
      right=pjump(nexc+1)
      do i=nexc+1,2*nexc-1
         if (eta1.ge.left.and.eta1.lt.right) then
            istate=i-nexc 
            exit
         endif
         left  = right
         right = left + pjump(i+1)
      enddo
      creal = real(c(istate+1))*sqrt(nr_gam(istate)*tmom2_0i(istate))
      ireal = aimag(c(istate+1))*sqrt(nr_gam(istate)*tmom2_0i(istate))
      c(1)  = cmplx(creal,ireal)
      c(2:nci) = zeroc
      c(1)=c(1)/sqrt(pjump(istate+nexc)/dt) 
      i_nr = i_nr +1
! Pure dephasing occurring 
   elseif (eta.ge.tmp2.and.eta.lt.tmp3) then
      call random_number(eta1)
      left=0.d0
      right=pjump(2*nexc+1)
      if (idep.eq.0) then
         do i=2*nexc+1,3*nexc
            if (eta1.ge.left.and.eta1.lt.right) then
               istate=i-2*nexc
               exit  
            endif
            left  = right 
            right = left + pjump(i+1)
         enddo
         ph = atan2(aimag(c(istate)),real(c(istate)))
         creal = abs(c(istate))*cos(ph+delta(istate))*sqrt(de_gam(istate))
         ireal = abs(c(istate))*sin(ph+delta(istate))*sqrt(de_gam(istate))  
         c(istate)  = cmplx(creal,ireal)
         c(1:istate-1) = zeroc 
         c(istate+1:nci) = zeroc 
         c(istate)=c(istate)/sqrt(pjump(istate+2*nexc)/dt)
      elseif (idep.eq.1) then
         do i=2*nexc+1,3*nexc-1
            if (eta1.ge.left.and.eta1.lt.right) then
               istate=i-2*nexc
               exit
            endif
            left  = right
            right = left + pjump(i+1)
         enddo
         creal = abs(c(istate+1))*sqrt(de_gam(istate))
         ireal = abs(c(istate+1))*sqrt(de_gam(istate))
         c(istate+1)  = cmplx(creal,ireal)
         c(1) = - c(1)*sqrt(de_gam(istate))
         c(2:istate) = zeroc
         c(istate+2:nci) = zeroc
      endif
      write(*,*), 'CIAO', istate     
      i_de = i_de + 1 
   endif

   return

  end subroutine quan_jump

  subroutine add_h_rnd(h_rnd,nci,w,w_prev,tdis)
!------------------------------------------------------------------------
! Random term in the Hamiltonian for the stochastic propagation 
! Random events: dissipation, nonradiative and dephasing
!
! Created   : E. Coccia 19 Jan 2017
! Modified  :
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)        :: nci, tdis
   real(8), intent(in)        :: w(3*nci), w_prev(3*nci)
   complex(16), intent(inout) :: h_rnd(nci,nci)
   integer                    :: i
   real(8)                    :: rate, rtmp, itmp 
   real(8)                    :: wrnd(3*nci)


   if (tdis.eq.0) then
      wrnd=w
   elseif (tdis.eq.1) then
      wrnd=w+w_prev
   endif 

! Matrix elements of S_alpha in the basis of the system eigenstates 
   h_rnd=zeroc

   itmp=0.d0
   do i=2, nci
      rtmp=0.d0 
! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
      rtmp = rtmp + sqrt(sp_gam(i-1)*tmom2_0i(i-1))*wrnd(i)
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha|
      if (nr_typ.eq.0) then
         rate = sqrt(nr_gam(i-1)*tmom2_0i(i-1))
      elseif (nr_typ.eq.1) then
         rate = sqrt(nr_gam(i-1))
      endif
      rtmp = rtmp + rate*wrnd(i+nci)
      h_rnd(1,i) = cmplx(rtmp,itmp)
   enddo

   do i=1, nci
! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
      rtmp = sqrt(de_gam(i))*cos(delta(i))*wrnd(i+2*nci)
      itmp = sqrt(de_gam(i))*sin(delta(i))*wrnd(i+2*nci) 
      h_rnd(i,i) = cmplx(rtmp,itmp)
   enddo

   !h_rnd(1,1) = h_rnd(1,1) + wrnd(1) !+ wrnd(1+nci)

   return

  end subroutine add_h_rnd

  subroutine define_h_dis(h_dis,nci,tdis)
!------------------------------------------------------------------------    
! Define the Markovian (tdis=0) or non-Markovian (tdis=1)
! dissipative term in the system Hamiltonian
!
! Created   : E. Coccia 20 Jan 2017
! Modified  :
!------------------------------------------------------------------------

   implicit none  
   integer, intent(in)    :: nci, tdis
   real(8), intent(inout) :: h_dis(nci,nci)
   integer                :: i

   h_dis=zero

   if (tdis.eq.0) then 
      call add_dis_m(h_dis,n_ci)
   elseif (tdis.eq.1) then 
      call add_dis_nm(h_dis,n_ci) 
   endif 

   return
 
  end subroutine define_h_dis

  subroutine rnd_noise(w,w_prev,nci,first,tdis) 
!------------------------------------------------------------------------
! Define the random fluctuating term in the
! stochastic propagator
!
! Created   : E. Coccia 20 Jan 2017
! Modified  :
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)    :: nci, tdis
   real(8), intent(inout) :: w(3*nci), w_prev(3*nci)
   logical, intent(in)    :: first
   integer                :: i,j

   if (tdis.eq.1) then
      if (first) then
         w=0.d0
         w_prev=0.d0
         do i=1,3*nci 
            do j=1,nrnd
               w(i) = w(i) + random_normal()
               w_prev(i) = w_prev(i) + random_normal()
            enddo
         enddo
      else 
         do i=1,3*nci
            w_prev(i) =  w(i)
            w(i) = 0.d0             
            do j=1, nrnd
               w(i) = w(i) + random_normal()
            enddo
         enddo
      endif
   elseif (tdis.eq.0) then
      w=0.d0
      do i=1,3*nci
         do j=1,nrnd
            w(i) = w(i) + random_normal()
         enddo
      enddo
   endif

   return

 end subroutine rnd_noise

  subroutine add_h_rnd2(h_rnd2,nci)
!------------------------------------------------------------------------
! Define the square of the dissipation/dephasing operator 
! Random events: dissipation, nonradiative and dephasing
!
! Created   : E. Coccia 10 Feb 2017
! Modified  :
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)        :: nci
   complex(16), intent(inout) :: h_rnd2(nci,nci)
   integer                    :: i
   real(8)                    :: rtmp, itmp 

! Matrix elements of S^2_alpha in the basis of the system eigenstates 
   h_rnd2=zeroc

! Sp and nr dissipation
! S^2_alpha = 0 (alpha.ne.0, by construction) 
! S^2_alpha = 1 (alpha.eq.0) FALSE
  !h_rnd2(1,1) = 1.d0  

! Pure dephasing (de)
! S^2_alpha = de_gam_alpha exp(i 2*delta_alpha) |Phi_alpha> <Phi_alpha|
   do i=1, nci
      rtmp = de_gam(i)*cos(2.d0*delta(i))
      itmp = de_gam(i)*sin(2.d0*delta(i)) 
      h_rnd2(i,i) = h_rnd2(i,i) + cmplx(rtmp,itmp)
   enddo
   h_rnd2 = 0.5d0*h_rnd2

   return

  end subroutine add_h_rnd2

end module 
