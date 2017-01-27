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
  public quan_jump, add_h_rnd, define_h_dis, rnd_noise 
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

! Matrix elements of S^+_alpha S_alpha in the system eigenstates basis

   do i=2, nci
! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
      h_dis(i,i) = h_dis(i,i) + sp_gam(i-1)*tmom2_0i(i-1)
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha|
      h_dis(i,i) = h_dis(i,i) + nr_gam(i-1)*tmom2_0i(i-1)
! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
      h_dis(i,i) = h_dis(i,i) + de_gam(i-1)
   enddo

   h_dis=0.5d0*h_dis

   !h_dis(1,1)=0.5d0

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

  subroutine loss_norm(c,nci)
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
      dnr = dnr + nr_gam(i)*weight
      dde = dde + de_gam(i)*tmp**2
   enddo
   dsp = dsp*dt
   dnr = dnr*dt
   dde = dde*dt
   dtot = dsp + dnr + dde

   return

  end subroutine loss_norm 

  subroutine quan_jump(c,nci)
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
   integer(4)                   :: istate
   real(8)                      :: eta, eta1, tmp1, tmp2, tmp3, modc, creal, ireal
   logical                      :: state=.false.

! (0)|------dsp/dtot-----|---dnr/dtot---|--dde/dtot--|(1) 
! Select the type of event according to
! a kinetic Monte Carlo strategy (JCP vol. 111 (1999) 10126)   

   call random_number(eta)

   tmp1 = dsp/dtot
   tmp2 = (dsp+dnr)/dtot
   tmp3 = (dsp+dnr+dde)/dtot

   state=.false.

! Select the relaxation channel
! j = n + FLOOR((m+1-n)*rnd), rnd [0,1) -> [n,m] 
! In our case [1,nexc]
   do while (.not.state)
      call random_number(eta1)
      istate = 1 + floor(nexc*eta1)
      modc = sqrt(real(c(istate+1))**2 + aimag(c(istate+1))**2)
      if (modc.ne.0.d0) state=.true.
   enddo

! Spontaneous occurring 
   if (eta.ge.0.and.eta.lt.tmp1) then
      if (istate.eq.1) istate=istate+1
      creal = real(c(istate+1))*sqrt(sp_gam(istate)*tmom2_0i(istate))
      ireal = aimag(c(istate+1))
      c(1)  = cmplx(creal,ireal) 
      c(2:nci) = zeroc
      c=c/abs(c(1))
      i_sp=i_sp+1 
! Nonradiative occurring
   elseif (eta.ge.tmp1.and.eta.lt.tmp2) then
   ! Select the relaxation channel
      if (istate.eq.1) istate=istate+1
      creal = real(c(istate+1))*sqrt(nr_gam(istate)*tmom2_0i(istate))
      ireal = aimag(c(istate+1))
      c(1)  = cmplx(creal,ireal)
      c(2:nci) = zeroc
      c=c/abs(c(1)) 
      i_nr = i_nr +1
! Pure dephasing occurring 
   elseif (eta.ge.tmp2.and.eta.lt.tmp3) then
      creal = real(c(istate+1))*cos(delta(istate))
      ireal = aimag(c(istate+1))*sin(delta(istate))  
      c(istate+1)  = cmplx(creal,ireal)
      c(1:istate) = zeroc 
      c(istate+2:nci) = 0.d0
      c=c/abs(c(istate+1))
      i_de = i_de + 1 
   endif

   return

  end subroutine quan_jump

  subroutine add_h_rnd(h_rnd,nci)
!------------------------------------------------------------------------
! Random term in the Hamiltonian for the stochastic propagation 
! Random events: dissipation, nonradiative and dephasing
!
! Created   : E. Coccia 19 Jan 2017
! Modified  :
!------------------------------------------------------------------------

   implicit none
   integer, intent(in)    :: nci
   real(8), intent(inout) :: h_rnd(nci,nci)
   integer                :: i

! Matrix elements of S_alpha in the basis of the system eigenstates 
   h_rnd=zero

   do i=2, nci
! Relaxation via spontaneous emission (sp)
! S_alpha = sqrt(sp_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha| 
      h_rnd(1,i) = h_rnd(1,i) + sqrt(sp_gam(i-1)*tmom2_0i(i-1))
! Relaxation via nonradiative processes (nr)
! S_alpha = sqrt(nr_gam_alpha) d_(alpha,0)  |Phi_0> <Phi_alpha|
      h_rnd(1,i) = h_rnd(1,i) + sqrt(nr_gam(i-1)*tmom2_0i(i-1))
! Pure dephasing (de)
! S_alpha = sqrt(de_gam_alpha) |Phi_alpha> <Phi_alpha|
      h_rnd(i,i) = h_rnd(i,i) + sqrt(de_gam(i-1))
   enddo

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
   real(8), intent(inout) :: w(nci), w_prev(nci)
   logical, intent(in)    :: first
   integer                :: i

   if (tdis.eq.1) then
      if (first) then
         do i=1,nci 
            w(i) = random_normal()
            w_prev(i) = random_normal()
         enddo
      else 
         do i=1,nci
            w(i) =  w_prev(i)
            w_prev(i) = random_normal()
         enddo
      endif
   elseif (tdis.eq.0) then
      do i=1,nci
         w(i) = random_normal()
      enddo
   endif

   return

 end subroutine rnd_noise

end module dissipation

