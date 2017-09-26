module vib
      use constants

      implicit none
      save

      integer(i4b)           :: nstates,nvib,nmodes,ntot
      real(dbl), allocatable :: w(:,:),q(:,:)
      real(dbl), allocatable :: e(:),dip(:,:,:)
      real(dbl), allocatable :: ef(:),dipf(:,:,:)
      logical                :: coupling

      public nstates,nvib,nmodes,w,q,e,dip,ef,dipf,fact,dfact,bin_coef, &
             mix,coupling
       
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
 
        namelist /vibrations/ nstates,nvib,nmodes,coupling,mix

        mix=.false.
        coupling=.false.
        nmodes=1
        nvib=10
        nstates=1

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

        allocate(w(nstates,nmodes),q(nstates,nmodes))

        write(*,*) '***************************************'
        write(*,*) '*                                     *'
        write(*,*) '*      Add vibrational levels         *' 
        write(*,*) '*     in harmonic approximation       *'
        write(*,*) '*      to any normal model of         *'
        write(*,*) '*     a given electronic state        *'
        write(*,*) '*                                     *'
        write(*,*) '***************************************'


        ntot=nstates*(nmodes*nvib) 
        write(*,*) ''
        write(*,*) 'Total number of states', ntot
        write(*,*) 'Number of electronic states', nstates
        write(*,*) 'Number of normal modes per electronic state', nmodes
        write(*,*) 'Number of vibrational states per normal mode', nvib
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
        real(dbl)                  :: hv,hve 
        !real(dbl)                  :: fact,dfact,bin_coef

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

        vt=vt-1
        vte=vte-1 

        nf = A*exp(-s)/(2**(vt+vte)*fact(vt)*fact(vte))
        t1 = 2.d0*sqrt(w)
        t2 = 2.d0*sqrt(we)
        r  = -we*d/(w+we)
        b  = -we*sqrt(w)*d/(w+we)
        be = w*sqrt(we)*d/(w+we) 

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

      subroutine hermite(n,x,y)                                   
!------------------------------------------------------------------------
!   COMPUTES THE VALUE OF THE HERMITE POLYNOMIAL OF DEGREE N            
!   at a given point               
!   n  = degree of the polynomial                                   
!   x  = point at which the calculation is done 
!   y  = value of the polynomial in x                                  
!------------------------------------------------------------------------
        implicit none                       
         
        integer(i4b), intent(in)  :: n
        real(dbl),    intent(in)  :: x
        real(dbl),    intent(out) :: y
        integer(i4b)              :: k        
        real(dbl)                 :: yp,dk,ym 
                                                        
        y = 1.d0                                                     
        if (n.eq.0) return                                              
                                                                        
        y = 2.d0*x                                                   
        if (n.eq.1) return                                              
                                                                        
        yp=1.d0                                                        
        do 1 k=2,n                                                        
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

        integer(i4b)                :: i,j,k,v,v1,kk,kk1,nvibt,n,l
        integer(i4b), allocatable   :: ii(:,:,:)  
        real(dbl)                   :: d,fc,mn

         
        allocate(ii(nstates,nmodes,nvib))
        kk=0
        do i=1,nstates
           do k=1,nmodes
              do v=1,nvib
                 kk=kk+1
                 ii(i,k,v)=kk
              enddo
           enddo
        enddo     

        do j=2,nstates
           do k=1,nmodes
              do v=1,nvib
   
              enddo
           enddo
        enddo  
 

        ! Only electronic transition i -> j (i>j) are taken into account
        ! No constraint on v and v1 values 
        ! Coupling between normal modes allowed

        ! Neglecting normal mode mixing
        if (.not.mix) then
           ! Electronic state i
           do i=nstates,1,-1
           ! Normal mode for state i 
              do k=1,nmodes
           ! Vibrational state v for (i,k)
                 do v=1,nvib
           ! Electronic state j
                    do j=i-1,1,-1
           ! Vibrational state v1 for (j,k) 
                       do v1=1,nvib
                          d=q(i,k)-q(j,k)
                          kk=ii(i,k,v)
                          kk1=ii(j,k,v1)
                          call compute_fc(v,v1,w(i,k),w(j,k),d,n,fc,mn)
                          fc_sum(i,j) = fc_sum(i,j) + fc
                          dipf(:,kk,kk1)=fc*dip(:,i,j)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        ! Duschinsky rotation
        else


        endif 


        stop 

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

        ! Ground state
        kk=0
        do k=1,nmodes
           kk=kk+1
           ef(1) = ef(1) + w(1,k)*0.5d0
        enddo 

        do v=1,nvib
           do k=1,nmodes
              
           enddo
        enddo

        ! Add vibrational energies
        ! Electronic ground state
        kk=0
        do k=1,nmodes
           do v=1,nvib
              kk=kk+1
              ef(kk) = ef(kk) + w(1,k)*(v-1+0.5d0)    
           enddo
        enddo
        ! Electronic excited states
        nvibt=nmodes*nvib
        do i=2,nstates
           kk=(nvibt)*(i-1)
           do k=1,nmodes
              do v=1,nvib
                 kk=kk+1
                 ef(kk)=e(i)+ w(i,k)*(v-1+0.5d0)
              enddo
           enddo
        enddo

        ! Print energies and dipoles for WaveT
        open(70,file='ci_energy_new.inp')
        open(71,file='ci_mut_new.inp')        
        kk=0
        do i=1,nstates
           do k=1,nmodes
              do v=1,nvib
                 kk=kk+1
                 write(70,*) 'Root', kk, ':', ef(kk)/ev_to_au 
              enddo
           enddo
        enddo

        ! Transition dipole moment between different normal modes
        ! in different electronic states is set to zero

        do k=1,nmodes
           do i=1,nstates
              do j=1,nstates
                 do v=1,nvib
                    do v1=1,nvib
                       kk=ii(i,k,v)
                       kk1=ii(j,k,v1)
                       write(71,*) 'States', kk, 'and',kk1,dipf(:,kk,kk1)
                    enddo
                 enddo
              enddo
           enddo
        enddo 


        close(70)
        close(71)

        deallocate(ii)

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
     
     subroutine gen_map(nmodes,nvib,iv)
!------------------------------------------------------------------------
! @brief Defining mapping array 
! 
! @date Created   : E. Coccia 26 Sep 2017
! Modified  :
!------------------------------------------------------------------------     

       implicit none

       integer(i4b),   intent(in)               :: nmodes,nvib
       integer(i4b),   allocatable, intent(out) :: iv(:,:)
       integer(i4b),   intent(out)              :: ntot
       integer(i4b),   allocatable              :: state(:)
       integer(i4b),   allocatable              :: x(:,:)

       integer(i4b)                             :: i,j,e,f 


      ! let E = indicates which variable XI the current (recursive) DO loop is for
      ! (e.g. if E = 4, then X4 is current)
      ! let STATE (E) = indicates the current "state" of XI number E (i.e. it's X1 or
      ! X2 value) (e.g. E = 4, then either it's X41 or X42 value)
      ! let F = indicates which combination number (in B) the current "system" of DO
      ! loops represents (where 1 <= F <= nvib**nmodes)

       ntot=nvib**nmodes

       allocate(state(nmodes))
       allocate(x(nmodes,nvib),iv(ntot,nmodes))

       do i=1,nvib
          x(:,i) = i-1
       enddo

       b(:,:)=0

       e = 1
       f = 1

       call rec_map(e,state,f,x,b,nmodes,nvib,ntot)

       do i=1,ntot
          do j=1,nmodes
             write(*,*) i,j,b(i,j)
          enddo
          write(*,*)
       enddo

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
 


end module vib
