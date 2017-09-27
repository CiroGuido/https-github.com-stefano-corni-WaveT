module vib
      use constants

      implicit none
      save

      integer(i4b)              :: nstates,nvib,nmodes,ntot,ncomb,nfc
      real(dbl),    allocatable :: w(:,:),q(:,:)
      real(dbl),    allocatable :: e(:),dip(:,:,:)
      real(dbl),    allocatable :: ef(:),dipf(:,:,:)
      integer(i4b), allocatable :: iv(:,:)
      logical                   :: coupling,mix

      public nstates,nvib,nmodes,w,q,e,dip,ef,dipf,fact,dfact,bin_coef, &
             mix,coupling,iv,ncomb,nfc
       
      contains        

      subroutine read_e_dip()
!------------------------------------------------------------------------
! @brief Read ci_energy.inp and ci_mut.inp from Gamess/Gaussian 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
       integer(i4b) :: i,j
       character(4) :: junk

!    read ci energy
       open(7,file="ci_energy.inp",status="old")
       ! add ground state
       allocate (e(nstates))
       e(1)=0.d0
       do i=2,nstates
         read(7,*) junk,junk,junk,e(i)
         e(i)=e(i)*ev_to_au
         write(6,*) i-1,e(i)
       enddo
       write (6,*)
       close(7)

!    read transition dipoles (also for pcm: useful for analysis)
       open(7,file="ci_mut.inp",status="old")
       allocate (dip(3,nstates,nstates))
       do i=1,nstates
        if (i.le.nstates) then
         read(7,*)junk,junk,junk,junk,dip(1,1,i),dip(2,1,i),dip(3,1,i)
         dip(:,i,1)=dip(:,1,i)
        else
         read(7,*)
        endif
       enddo

       do i=2,nstates
         do j=2,i
          if (i.le.nstates.and.j.le.nstates) then
           read(7,*)junk,junk,junk,junk,dip(1,i,j),dip(2,i,j),dip(3,i,j)
           dip(:,j,i)=dip(:,i,j)
          else
           read(7,*)
          endif
         enddo
       enddo
       close(7)

       return

      end subroutine read_e_dip
    
      subroutine read_input_vib()
!------------------------------------------------------------------------
! @brief Read input for vibrational corrections 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------
 
        integer(i4b)             :: idum,i,j 
 
        namelist /vibrations/ nstates,nvib,nmodes,coupling,mix,nfc

        mix=.false.
        coupling=.false.
        nmodes=1
        nvib=10
        nstates=1
        nfc=1

        ! Read w and q for any vib level
        ! Excited state N
        ! w q for mode 1
        ! w q for mode 2
        ! ...
        ! Excited state N+1
        ! w q for mode 1
        ! w q for mode 2
        ! ...

        read(*,nml=vibrations)

        ncomb=nvib**nmodes

        allocate(w(nstates,nmodes),q(nstates,nmodes))
        allocate(iv(ncomb,nmodes))

        write(*,*)  
        write(*,*) '***************************************'
        write(*,*) '*                                     *'
        write(*,*) '*      Add vibrational levels         *' 
        write(*,*) '*     in harmonic approximation       *'
        write(*,*) '*      to any normal model of         *'
        write(*,*) '*     a given electronic state        *'
        write(*,*) '*                                     *'
        write(*,*) '***************************************'


        ntot=nstates*ncomb
        write(*,*) ''
        write(*,*) 'Total number of states', ntot
        write(*,*) 'Number of electronic states', nstates
        write(*,*) 'Number of normal modes per electronic state', nmodes
        write(*,*) 'Number of vibrational states per normal mode', nvib
        if (mix) then
           write(*,*) 'Duschinsky rotation for normal coordinates' 
        else
           write(*,*) 'Only displacements between normal coordinates'
        endif 
        write(*,*) ''

        w(:,:) = 1000.d0
        q(:,:) = 2.d0

        open(60,file='vib.dat')
        do i=1,nstates
           read(60,*) idum, idum 
           do j=1,nmodes
              read(60,*) w(i,j), q(i,j)
           enddo
        enddo
        close(60)

        w(:,:) = w(:,:)*cm_to_au
 
        allocate(ef(ntot),dipf(3,ntot,ntot)) 

        return

      end subroutine read_input_vib 

      subroutine compute_fc(v,ve,w,we,d,n,fc,mn)
!------------------------------------------------------------------------
! @brief Compute Franck-Condon factors between the vibrational 
! eigenstates (harmonic oscillator) of any electronic ground-excited
! pair 
! Formula from J. Mol. Spectroscopy, vol. 232, 102 (2005)
! 
! @date Created   : E. Coccia 7 Sep 2017
! Modified  :
!------------------------------------------------------------------------

        implicit none

        integer(i4b),  intent(in)  :: v,ve,n
        real(dbl),     intent(in)  :: w,we,d
        real(dbl),     intent(out) :: fc,mn
        integer(i4b)               :: k,ke,kk,k2,vt,vte
        real(dbl)                  :: s,a,b,be,ik,r,nf
        real(dbl)                  :: t1,t2,ikk
        real(dbl)                  :: hv,hve,ww 

! Harmonic oscillator eigenfunction for the electronic ground state
! |v> = Nv*Hv(sqrt(alpha)x)*exp(-1/2*alpha*x^2)
! Nv = sqrt(sqrt(alpha)/(2^v*v!*sqrt(pi)))
! alpha = omega/hbar
! x -> normal coordinate
! Hv -> Hermite polynomial for v

! Harmonic oscillator eigenfunction for the electronic excited state
! |v'> = Nv'*Hv'(sqrt(alpha')x')*exp(-1/2*alpha'*x'^2)
! Nv' = sqrt(sqrt(alpha')/(2^v'*v'!*sqrt(pi)))
! alpha' = omega'/hbar
! x'= x+d -> normal coordinate
! Hv' -> Hermite polynomial for v'

! s = alpha*alpha'*d^2/(alpha+alpha')
! a = 2*sqrt(alpha*alpha')/(alpha+alpha')
! b = -alpha'*sqrt(alpha)*d/(alpha+alpha')
! b'= alpha*sqrt(alpha')*d/(alpha+alpha')
! kk = (k+k')/2
! I(kk) = 0 if k+k' is odd
! I(kk) = (2*kk-1)!!/(alpha+alpha')^kk

! <v|v'> = sqrt(A*exp(-S)/(2^(v+v')*v!*v'!))*sum_k=0^v * sum_k'^v' *
! (v)*(v')*Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*I(kk)
! (k) (k')

! Franck-Condon factor
! |<v|v'>|^2 = A*exp(-S)/(2^(v+v')*v!*v'!)*sum_k=0^v * sum_k'^v' *
! (v)*(v')*Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*I(kk)
! (k) (k')

        !vt=vt-1
        !vte=vte-1 

        ww=1.d0/(w+we) 

        a  = 2.d0*sqrt(w*we)*ww
        s  = w*we*d**2*w
        nf = a*exp(-s)/(2**(vt+vte)*fact(vt)*fact(vte))
        t1 = 2.d0*sqrt(w)
        t2 = 2.d0*sqrt(we)
        r  = -we*d*ww
        b  = -we*sqrt(w)*d*ww
        be = w*sqrt(we)*d*ww 

        fc=0.d0
        mn=0.d0
        do k=0,v
           call hermite(v-k,b,hv)
           do ke=0,ve
              call hermite(ve-ke,be,hve)
              if (mod(k+ke,2).ne.0) then
                 ik=0.d0
              else
                 ik = dfact(k+ke-1)/(w+we)**(0.5d0*(k+ke))
              endif
              fc = fc + bin_coef(v,k)*bin_coef(ve,ke)*hv*hve*t1**k*t2**ke*ik
              do k2=0,n
                 if (mod(k+ke+k2,2).ne.0) then
                     ikk=0.d0
                 else
                     ikk = dfact(k+ke+k2-1)/(w+we)**(0.5d0*(k+ke+k2))
                 endif
                 mn = mn + bin_coef(v,k)*bin_coef(ve,ke)*bin_coef(n,k2)*hv*hve*t1**k*t2**ke*r**(n-k2)*ikk
              enddo
           enddo
        enddo

        fc = nf*fc**2

! <v|x^n|v'> = sqrt(A*exp(-S)/(2^(v+v')*v!*v'!))*sum_k=0^v * sum_k'^v' *
! sum_k''=0^n*
! (v)*(v')*(n)
! *Hv-k(b)*Hv'-k'(b')*(2*sqrt(alpha))^k*(2*sqrt(alpha'))^k'*r^(n-k'')*I(kk)
! (k) (k') (k'')
! r = - alpha'*d/(alpha+alpha')

        mn = sqrt(nf)*mn

        return

      end subroutine compute_fc

      subroutine hermite(nn,x,y)                                   
!------------------------------------------------------------------------
!   @brief Computes the value of the Hermite polynomial of degree nn             
!   at a given point               
!   nn = degree of the polynomial                                   
!   x  = point at which the calculation is done 
!   y  = value of the polynomial in x                                
!
!   @date Created  :  E. Coccia 8 Sep 2017 
!   Modified   :
!------------------------------------------------------------------------
        implicit none                       
         
        integer(i4b), intent(in)  :: nn
        real(dbl),    intent(in)  :: x
        real(dbl),    intent(out) :: y
        integer(i4b)              :: k        
        real(dbl)                 :: yp,dk,ym 
                                                        
        y = 1.d0                                                     
        if (nn.eq.0) return                                              
                                                                        
        y = 2.d0*x                                                   
        if (nn.eq.1) return                                              
                                                                        
        yp=1.d0                                                        
        do 1 k=2,nn                                                        
           dk = dble(k-1)                                               
           ym = y                                                        
           y  = 2.d0*x*y-2.d0*dk*yp                                      
           yp = ym 
1       continue                                                          
                                                                        
        return
                                                            
      end subroutine hermite             

      subroutine deallocate_vib()
!------------------------------------------------------------------------
! @brief Deallocate arrays in module vib 
! 
! @date Created   : E. Coccia 8 Sep 2017
! Modified  :
!------------------------------------------------------------------------

       deallocate(w,q)
       deallocate(e,dip)
       deallocate(ef,dipf)
       deallocate(iv)

       return
      
     end subroutine deallocate_vib 

     real(dbl) function fact(n)
!------------------------------------------------------------------------
! @brief Function computing the factorial n!
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------

       integer(i4b), intent(in) :: n
       integer(i4b)             :: k 
       real(dbl)                :: s

       s=1.d0

       do k=1,n
          s=s*k
       enddo
   
       fact=s

       if (n.eq.0) fact=1 
 
       return 

     end function


     real(dbl) function dfact(n)
!------------------------------------------------------------------------
! @brief Function computing the double factorial n!!
! if n even, n!! = prod_k=1^(n/2) (2k)
! if n odd,  n!! = prod_k=1^(n+1/2) (2k-1)
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------

       integer(i4b), intent(in) :: n 
       integer(i4b)             :: k 
       real(dbl)                :: s

       s=1.d0

       if (mod(n,2).eq.0) then
          do k=1,n/2
             s=s*2.d0*k
          enddo
       else
          do k=1,(n+1)/2
             s=s*(2.d0*k-1)
          enddo
       endif

       if (n.le.0) s=1.d0

       dfact=s

       return
   
     end function 

     real(dbl) function bin_coef(n,k)
!------------------------------------------------------------------------
! @brief Function computing the binomial coefficient for n and k 
! (n) = n! / (k!*(n-k)!) 
! (k)
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------
  
        integer(i4b), intent(in) ::n,k      
        

        if (k.gt.n) then
           write(*,*) ''
           write(*,*) 'ERROR: k must be less than'
           write(*,*) '(or equal to) n: k=', k,'n=',n
           write(*,*) ''
        endif

        bin_coef=fact(n)/(fact(k)*fact(n-k))

        return

     end function

     subroutine compute_e_dip()
!------------------------------------------------------------------------
! @brief Correct electronic energies and dipoles with
! vibrational energies (harmonic approximation) and
! Franck-Condon factors 
! 
! @date Created   : E. Coccia 11 Sep 2017
! Modified  :
!------------------------------------------------------------------------

        integer(i4b)                :: i,j,k,v,v1,kk,kk1,ii,jj

        call gen_map(nmodes,nvib,ncomb,iv)

        ! Only electronic transition i -> j (i>j) are taken into account
        ! No constraint on v and v1 values 

        dipf=0.d0
        ! Neglecting normal mode mixing
        if (.not.mix) then
           kk=0
           do i=1,nstates
              do j=1,ncomb
                 kk=kk+1 
                 kk1=0
                 do ii=1,nstates
                    do jj=1,ncomb
                       kk1=kk1+1
                       call modify_dip(dip(:,i,ii),j,jj,i,ii,nmodes,dipf(:,kk,kk1))                  
                    enddo
                 enddo
              enddo 
           enddo 
        ! Duschinsky rotation
        else
           write(*,*) 'Duschinsky rotation not implemented yet'
           stop
        endif

! The updated energy file has the following structure:
!
! Ground state
!
! E_0 + w_1/2 + w_2/2 + ... + w_k/2  k=1,...,nmodes
!
! nvib**nmodes combinations for the ground state
! E_0 + w_1/2 + w_2*(1+1/2) + ... w_k/2 
!
! ...
!
! E_1 + w_1/2 + w_2/2 + ... + w_k/2
! E_1 + w_1/2 + w_2*(1+1/2) + ... w_k/2 
!
! ...
!
! ...
!
! E_(nstates-1) + w_1/2 + w_2/2 + ... + w_k/2
! E_(nstates-1) + w_1/2 + w_2*(1+1/2) + ... w_k/2 
!
! ...

! Total numbers of states = nstates*nvib**nmodes

        kk=0
        do i=1,nstates 
           do j=1,ncomb
              kk=kk+1
              call add_vibe(e(i),i,nmodes,iv(j,:),ef(kk))
           enddo
        enddo

        ! Print energies and dipoles for WaveT
        open(70,file='ci_energy_new.inp')
        open(71,file='ci_mut_new.inp')        
        kk=0
        do i=1,nstates
           do j=1,ncomb
              kk=kk+1
              if (kk.gt.1) write(70,*) 'Root', kk-1, ':', (ef(kk)-ef(1))/ev_to_au 
           enddo
        enddo

        kk=0
        do i=1,nstates
           do j=1,ncomb
              kk=kk+1
              kk1=0
              do ii=i,nstates
                 do jj=1,ncomb 
                    kk1=kk1+1
                    write(71,*) 'States', kk, 'and',kk1,dipf(:,kk,kk1)
                 enddo
              enddo
           enddo
        enddo

        close(70)
        close(71)

        return
 
     end subroutine compute_e_dip

     subroutine compute_coupling()  
!------------------------------------------------------------------------
! @brief Non adiabatic correction to 
! nradiative decay and dephasing rates
! 
! @date Created   : E. Coccia 12 Sep 2017
! Modified  :
!------------------------------------------------------------------------     

       implicit none


       write(*,*) 'Not implemented yet'
       stop

     end subroutine compute_coupling
     
     subroutine gen_map(nmodes,nvib,ntot,iv)
!------------------------------------------------------------------------
! @brief Defining mapping array 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------     

       implicit none

       integer(i4b),   intent(in)               :: nmodes,nvib,ntot
       integer(i4b),   intent(out)              :: iv(ntot,nmodes)
       integer(i4b),   allocatable              :: state(:)
       integer(i4b),   allocatable              :: x(:,:)

       integer(i4b)                             :: i,j,e,f 


      ! let E = indicates which variable XI the current (recursive) DO loop is for
      ! (e.g. if E = 4, then X4 is current)
      ! let STATE (E) = indicates the current "state" of XI number E (i.e. it's X1 or
      ! X2 value) (e.g. E = 4, then either it's X41 or X42 value)
      ! let F = indicates which combination number (in B) the current "system" of DO
      ! loops represents (where 1 <= F <= nvib**nmodes)

       allocate(state(nmodes))
       allocate(x(nmodes,nvib))

       do i=1,nvib
          x(:,i) = i-1
       enddo

       iv(:,:)=0

       e = 1
       f = 1

       call rec_map(e,state,f,x,iv,nmodes,nvib,ntot)

       !do i=1,ntot
       !   do j=1,nmodes
       !      write(*,*) i,j,iv(i,j)
       !   enddo
       !   write(*,*)
       !enddo

       deallocate(state,x)

       return 

     end subroutine gen_map

     recursive subroutine rec_map(e,state,f,x,b,nmodes,nvib,ntot)
!------------------------------------------------------------------------
! @brief Defining mapping array with a recursive scheme 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------
       implicit none

       integer,   intent(in)    :: e,nmodes,nvib,ntot
       integer,   intent(inout) :: state(nmodes)
       integer,   intent(inout) :: f
       integer,   intent(in)    :: x(nmodes,nvib)
       integer,   intent(inout) :: b(ntot,nmodes)

       integer :: i,j,k,next_e

       if (e.le.nmodes) then
! loop over each "state" of current XI
          do i=1,nvib
             state(e)=i
             next_e=e+1
             call rec_map(next_e,state,f,x,b,nmodes,nvib,ntot)
          enddo
       else
! assign values to the current combination in B
          do j=1,nmodes
             k=state(j)
             b(f,j)=x(j,k)
          enddo
! update the next combination number in B
          f=f+1
          return
       endif

       return

     end subroutine rec_map 

     subroutine add_vibe(e,i,nmodes,ii,ef)
!------------------------------------------------------------------------
! @brief Adding vibrational levels 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------

       implicit none 

       real(dbl),    intent(in)  :: e
       integer(i4b), intent(in)  :: i,nmodes
       integer(i4b), intent(in)  :: ii(nmodes)
       real(dbl),    intent(out) :: ef

       integer(i4b)  :: k

       ef=e
       do k=1,nmodes
          !write(*,*) e,k,w(i,k)*219474, ii(k), ef 
          ef = ef + w(i,k)*(ii(k)+0.5d0) 
       enddo

       return

     end subroutine add_vibe

     subroutine modify_dip(dip,i,j,l,m,nmodes,dipf)
!------------------------------------------------------------------------
! @brief Correct dipoles with Franck-Condon factors 
!                   
! @date Created   : E. Coccia 27 Sep 2017
! Modified  :          
!------------------------------------------------------------------------

       implicit none

       real(dbl),    intent(in)   :: dip(3)
       integer(i4b), intent(in)   :: i,j,l,m,nmodes
       real(dbl),    intent(out)  :: dipf(3)

       integer(i4b)               :: v,v1,k
       real(dbl)                  :: d,fc,mn,tfc

       tfc=1.d0
       do k=1,nmodes 
          v=iv(i,k)
          v1=iv(j,k) 
          d=q(l,k)-q(m,k)
          call compute_fc(v,v1,w(l,k),w(m,k),d,nfc,fc,mn)
          !write(*,*) k,l,k,fc
          tfc=tfc*fc 
       enddo
       dipf(:)=dip(:)*tfc     
 
       return       
 
     end subroutine modify_dip

end module vib
