!! Module containing subroutines to parse .wfx, .fchk, .log, and .bas files.
module read_files
   use wfxinfo   
   use geninfo
   use loginfo
   use locatemod
   use numbers
   implicit none
contains 
   !****************************************************************************
   !! Parses a .wfx file to extract molecular and primitive information.
   subroutine filewfx(wfxfilename)   
   character(len=40), intent(in) :: wfxfilename
   character*80 :: line
   character*80 :: zaborra  
   integer :: i,j,k, kk
   integer :: mon  
   double precision :: maxim, minim
   double precision, parameter :: trsh=1.d-16 
   
    corr=.false.
    uhf=.false.
    rhf=.false.
    udens=.false.
    rdens=.false.  
    opsh=.false.
    clsh=.false.
      
   open(unit=1,file=wfxfilename,status='OLD') 
   
   call locate(1,"Number of Nuclei") 
   read(1,*) natoms 
                           
   rewind 1  
   
   call locate(1,"Number of Occupied Molecular Orbitals")
   read(1,*) noccmo
   
   allocate(Occ(noccmo))
   
   rewind 1
   
   call locate(1,"Molecular Orbital Occupation Numbers") 
   do i=1,noccmo                                         
      read(1,*) Occ(i)                                   
   end do
    
   do i=1,noccmo
      if (occ(i).ne.0.d0) then
        if ((Occ(i).ne.1.d0).and.(Occ(i).ne.2.d0)) then
           write(*,*) "Correlated wavefunction"
           corr=.true.
           exit
        end if
      end if  
   end do
   
   maxim=maxval(Occ)                                         
   minim=minval(Occ)   
   if (corr) then
       if (maxim.ge. 2.d0) then             
           if ((maxim-2.d0 .gt. trsh).or.(minim.lt. 0.d0)) then 
              write(*,*) "RELAXED DENSITY (closed shell)"                  
              rdens=.true.
           else     
              write(*,*) "UNRELAXED DENSITY (closed shell)"
              udens=.true.
           end if                    
       else                                                    
           if ((maxim-1.d0 .gt. trsh).or.(minim.lt. 0.d0)) then 
              write(*,*) "RELAXED DENSITY (open shell)"
              rdens=.true.
           else    
              write(*,*) "UNRELAXED DENSITY (open shell)"
              udens=.true.
           end if                    
       end if 
   else                    
      minim=minval(Occ,MASK=occ.gt.0.d0)     
      if (minim.eq.1.d0) then                 
            if (maxim.eq.2.d0) then           
               rhf=.true.
               uhf=.false.
               opsh=.true.
               write(*,*) "Open shell RHF"
            else                                 
               rhf=.false.
               uhf=.true.  
               opsh=.true.
               write(*,*) "UHF"
            end if
      else if (minim.eq.2.d0) then           
            rhf=.true.
            uhf=.false.  
            clsh=.true.
            write(*,*) "Closed shell RHF"
      else 
            write(*,*)   "Error, correlated wf"  
            STOP
      end if    
   end if  
    
   call locate(1,"Net Charge")           
   read(1,*) netch
   
   rewind 1
   
   call locate(1,"Number of Electrons")
   read(1,*) nelec
   
   call locate(1,"Number of Alpha Electrons")
   read(1,*) nalfae
   
   rewind 1
   
   call locate(1,"Number of Beta Electrons")
   read(1,*) nbetae
   
   rewind 1
   
   call locate(1,"Number of Primitives")
   read(1,*) nprim
   
   write(*,*) "nprim=",nprim
   rewind 1
   
   allocate(cartes(natoms,3))
   call locate(1,"Nuclear Cartesian Coordinates")
   read(1,*) ((cartes(i,j) , j=1,3 ), i=1,natoms)
   
   rewind 1
   
   allocate(Ra(nprim))
   call locate(1,"Primitive Centers")
   read(1,*) (Ra(i), i=1,nprim)
   
   rewind 1
   
   allocate(Ptyp(nprim))
   call locate(1,"Primitive Types")
   read(1,*) (Ptyp(i), i=1,nprim)
   
   rewind 1
   
   allocate(Alpha(nprim))
   call locate(1,"Primitive Exponents")
   read(1,*) (Alpha(i), i=1,nprim)
   
   rewind 1
   
   allocate(an(natoms))
   call locate(1,"Atomic Numbers")
   do i=1,natoms
    read(1,*) an(i)
   end do
   
   allocate(chrg(natoms))
   call locate(1,"Nuclear Charges")
   do i=1,natoms
     read(1,*) chrg(i)
   end do  
   
   allocate(T(noccmo,nprim))                        
   call locate(1,"Molecular Orbital Primitive Coefficients")
   do kk=1,noccmo
     read(1,'(a80)')line
     read(1,*) mon 
    read(1,'(a80)') zaborra  
    read(1,*) (T(mon,i), i=1,nprim)
   end do
   
   rewind 1
   
   allocate(TMN(nprim,3))      
   do i=1,nprim                                          
      if (ptyp(i).eq.1) then                                             
        do j=1,3                   
          TMN(i,j)=0
        end do
      else if (ptyp(i).eq.2) then  
         TMN(i,1)=1
         TMN(i,2)=0
         TMN(i,3)=0
      else if (ptyp(i).eq.3) then  
         TMN(i,1)=0
         TMN(i,2)=1
         TMN(i,3)=0   
      else if (ptyp(i).eq.4) then  
         TMN(i,1)=0
         TMN(i,2)=0
         TMN(i,3)=1   
      else if (ptyp(i).eq.5) then  
         TMN(i,1)=2
         TMN(i,2)=0
         TMN(i,3)=0   
      else if (ptyp(i).eq.6) then  
         TMN(i,1)=0
         TMN(i,2)=2
         TMN(i,3)=0  
      else if (ptyp(i).eq.7) then  
         TMN(i,1)=0
         TMN(i,2)=0
         TMN(i,3)=2   
      else if (ptyp(i).eq.8) then  
         TMN(i,1)=1
         TMN(i,2)=1
         TMN(i,3)=0   
      else if (ptyp(i).eq.9) then  
         TMN(i,1)=1
         TMN(i,2)=0
         TMN(i,3)=1
      else if (ptyp(i).eq.10) then 
         TMN(i,1)=0
         TMN(i,2)=1
         TMN(i,3)=1      
      else if (ptyp(i).eq.11) then 
         TMN(i,1)=3
         TMN(i,2)=0
         TMN(i,3)=0
      else if (ptyp(i).eq.12) then 
         TMN(i,1)=0
         TMN(i,2)=3
         TMN(i,3)=0
      else if (ptyp(i).eq.13) then 
         TMN(i,1)=0
         TMN(i,2)=0
         TMN(i,3)=3
      else if (ptyp(i).eq.14) then 
         TMN(i,1)=2
         TMN(i,2)=1
         TMN(i,3)=0
      else if (ptyp(i).eq.15) then 
         TMN(i,1)=2
         TMN(i,2)=0
         TMN(i,3)=1
      else if (ptyp(i).eq.16) then 
         TMN(i,1)=0
         TMN(i,2)=2
         TMN(i,3)=1
      else if (ptyp(i).eq.17) then 
         TMN(i,1)=1
         TMN(i,2)=2
         TMN(i,3)=0
      else if (ptyp(i).eq.18) then 
         TMN(i,1)=1
         TMN(i,2)=0
         TMN(i,3)=2
      else if (ptyp(i).eq.19) then 
         TMN(i,1)=0
         TMN(i,2)=1
         TMN(i,3)=2
      else if (ptyp(i).eq.20) then 
         TMN(i,1)=1
         TMN(i,2)=1
         TMN(i,3)=1
      else
        write(*,*) "G orbitals are not implemented"
        stop
      end if
   end do
                                
   if (uhf) then                
        allocate(T_a(nalfae,nprim))
        allocate(T_b(nbetae,nprim))
        do i=1,nalfae                
             do j=1,nprim    
                   T_a(i,j)= T(i,j)  
             end do
        end do
        k=nalfae+1                   
        do i=1,nbetae     
             do j=1,nprim
                   T_b(i,j) = T(k,j)     
             end do
             k=k+1                   
        end do
        write(*,*) "reading alpha and beta coefficients" 
   end if 
   do i=1,nprim
      write(*,*) cartes(Ra(i),1), cartes(Ra(i),2), cartes(Ra(i),3), Alpha(i), TMN(i,1), TMN(i,2), TMN(i,3)
   end do
   allocate(Xn(nprim));allocate(Yn(nprim));allocate(Zn(nprim))
   do i=1,nprim
      Xn(i)=cartes(Ra(i),1); Yn(i)=cartes(Ra(i),2); Zn(i)=cartes(Ra(i),3)
   end do
   
    close(1) 
   end subroutine filewfx
   
   !****************************************************************************
   !! Parses a .log file to extract AO basis set information and populate related variables.
   subroutine filelog(logfilename)                                             
   use loginfo                    
   use geninfo
   use wfxinfo
   use locatemod
   implicit none
   integer, parameter :: maxao=100
   character(len=*), intent(in) :: logfilename
   character*80 :: line 
   character*80 :: zaborra
   character*10 :: tao 
   integer :: i,j,k,l,m
   integer :: sm, smm, sm2, sm3, sm4 
   integer :: smatom
   double precision :: kk 
   
   open(unit=4 ,file=logfilename, status='OLD') 
   
   allocate(npao(maxao)) 
   allocate(Flg(nprim))  
   allocate(aotyp(maxao))
   
   Flg=0.d0
   
   smm=0 
   sm=0  
   smatom=1
   call locate(4,"AO basis set in the form of general basis input") 
   read(4,'(a80)') zaborra
   
   do while (smatom.le.natoms)   
         read(4,'(a80)') line
         read(line(2:3),*) tao 
         if (tao.eq."**") then 
              smatom=smatom+1
              read(4,'(a80)') zaborra   
         else if (tao.eq."S") then
              sm=sm+1  
              read(line(5:7),*) npao(sm) 
              do i=1,npao(sm)
                   smm=smm+1 
                   read(4,'(a80)') line
                   read(line(24:39),*) Flg(smm)   
              end do
              aotyp(sm)=0  
          else if (tao.eq."SP") then 
               sm=sm+1 
               aotyp(sm)=0
               read(line(5:7),*) npao(sm)      
               do i=1,npao(sm)
                    smm=smm+1  
                    sm2=sm+1     
                    aotyp(sm2)=1
                    sm3=sm+2     
                    aotyp(sm3)=1
                    sm4=sm+3       
                    aotyp(sm4)=1  
                    k=smm+npao(sm)      
                    read(4,*) kk, Flg(smm), Flg(k) 
                    l=k+npao(sm)                   
                    m=l+npao(sm)     
                    Flg(l)=Flg(k)     
                    Flg(m)=Flg(k)     
               end do   
               smm=smm+3*npao(sm)                     
               do i=1,3
                    npao(sm+i)=npao(sm) 
               end do
               sm=sm+3 
          else if (tao.eq."D") then
               sm=sm+1  
               aotyp=2 
               read(line(5:7),*) npao(sm)   
               do i=1,npao(sm) 
                    smm=smm+1  
                    read(4,*) kk, Flg(smm)      
                    do j=1,5  
                           k=smm+j
                           l=sm+j   
                           if (l.le.3) then 
                                    aotyp(l)=2  
                           else
                                    aotyp(l)=3    
                           end if 
                           Flg(k)=Flg(smm)       
                           npao(sm+j)=npao(sm)    
                    end do
               end do    
               smm=smm+6*npao(sm)  
               sm=sm+5     
          end if
   end do
   nao=sm
   
   allocate(N_prim(nprim))  
   do i=1,nprim
         if (ptyp(i).eq.1) then 
               N_prim(i)=(pi*(2*Alpha(i))**(-1))**(-3.d0/4.d0)
         else if ((ptyp(i).gt.1).and.(ptyp(i).le.4)) then
               N_prim(i)= (pi*0.5d0)**(-3.d0/4.d0) * (2*(Alpha(i)**(1.25d0)))   
         else if ((ptyp(i).gt.4).and.(ptyp(i).le.7)) then  
               N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))
         else if ((ptyp(i).gt.7).and.(ptyp(i).le.10)) then  
               N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))/(dsqrt(3.d0))   
         else
               write(*,*) "f type orbital, insert its normalization equation"
         end if    
   end do
   do i=1, nprim              
         Flg(i)=Flg(i)*N_prim(i)
   end do
   
   call cmatrix() 
   
   close(4)
   end subroutine filelog
   
   !*******************************************************************************
   !! Computes the expansion coefficients of MOs in AOs
   subroutine cmatrix()   
   use wfxinfo
   use loginfo
   use geninfo
   implicit none
   integer :: i, j, k
        
   allocate(Ckalk(noccmo,nao))
   write(*,*) nao
   write(*,*) noccmo
   do i=1,noccmo
    k=1
    do j=1,nao
         Ckalk(i,j)=T(i,k)/(Flg(k)) 
         k=k+npao(j)  
     end do
   end do   
   
   end subroutine cmatrix
   
   !****************************************************************************
   !! Parses a .fchk file to extract molecular and primitive information, similar to filewfx but adapted to the .fchk format.
   subroutine filefchk(fchkfilename)
      use wfxinfo
      use geninfo
      use locatemod
      implicit none
      character(len=*), intent(in) :: fchkfilename
      character(len=500) :: line
      character(len=40)  :: typ, method, basis
      integer :: i, j, k, l, l1 , iprim, ncshl, npshl
      integer, allocatable :: shell_types(:), prim_per_shell(:), shell2atom(:)
      double precision, allocatable :: prim_exp(:)
      
      corr=.false.; uhf=.false.; rhf=.false.
      udens=.false.; rdens=.false.
      opsh=.false.; clsh=.false.
   
      open(unit=1,file=fchkfilename,status='OLD')
   
      call getline_with(1,"SP",line)
      read(line,*) typ, method, basis
      call getline_with(1,"Number of alpha electrons",line) 
      read(line,*) typ, typ, typ, typ, typ, nalfae
      call getline_with(1,"Number of beta electrons",line)
      read(line,*) typ, typ, typ, typ, typ, nbetae
      nelec = nalfae + nbetae
      call getline_with(1,"Total Energy",line)
      read(line,*) typ, typ, typ, toteng
      write(*,*) "Total Energy:", toteng
      write(*,*) "Calculation type:", trim(method), trim(basis)
      if (nalfae.eq.nbetae) then
         clsh=.true.
         write(*,*) "Closed shell"
      else
         opsh=.true.
         write(*,*) "Open shell"
      end if
   
      if (method.eq."RHF") then
         rhf=.true.
      else if (method.eq."HF") then
         if (opsh) uhf=.true.
         if (clsh) rhf=.true.
      else if (method.eq."UHF") then  
         write(*,*) "Open shell UHF"
         uhf=.true.
      else
         corr=.true.
         write(*,*) "Correlated wavefunction"      
      end if
   
      write(*,*) "Reading fchk file"
      call getline_with(1,"Number of atoms",line)
      read(line,*) typ, typ, typ, typ, natoms   
   
      allocate(an(natoms))
      call locate(1,"Atomic numbers")
      read(1,*) (an(i), i=1,natoms)
   
      allocate(chrg(natoms))
      chrg(:) = dble(an(:))
   
      allocate(cartes(natoms,3))
      call locate(1,"Current cartesian coord")
      read(1,*) ((cartes(i,j), j=1,3), i=1,natoms)
   
      call getline_with(1,"Number of contracted shells",line)
      read(line,*) typ, typ, typ, typ, typ, ncshl
      call getline_with(1,"Number of primitive shells",line)
      read(line,*) typ, typ, typ, typ, typ, npshl
   
      allocate(shell_types(ncshl))
      call locate(1,"Shell types")
      read(1,*) (shell_types(i), i=1,ncshl)
   
      allocate(prim_per_shell(ncshl))
      call locate(1,"Number of primitives per shell")
      read(1,*) (prim_per_shell(i), i=1,ncshl)
   
      nprim = 0
      do i=1,ncshl
         select case(shell_types(i))
         case(0)  
            nprim = nprim + prim_per_shell(i)
         case(1)  
            nprim = nprim + 3*prim_per_shell(i)
         case(2)  
            nprim = nprim + 6*prim_per_shell(i)
         case(3)  
            nprim = nprim + 10*prim_per_shell(i)   
         case(-1) 
            nprim = nprim + 4*prim_per_shell(i)
         case(-2) 
            nprim = nprim + 6*prim_per_shell(i)
         case(-3) 
            nprim = nprim + 10*prim_per_shell(i)
         case(4)  
            nprim = nprim + 15*prim_per_shell(i)
         case(-4) 
            nprim = nprim + 15*prim_per_shell(i)   
         case default
            write(*,*) "Shell type not implemented:", shell_types(i)
            stop
         end select
      end do
   
      allocate(shell2atom(ncshl))
      call locate(1,"Shell to atom map")
      read(1,*) (shell2atom(i), i=1,ncshl)
   
      allocate(prim_exp(npshl))
      call locate(1,"Primitive exponents")
      read(1,*) (prim_exp(i), i=1,npshl)
     
      allocate(Ra(nprim))
      allocate(Alpha(nprim))
      allocate(Ptyp(nprim))
      allocate(TMN(nprim,3))
      iprim = 0
      l=0 
      do i=1,ncshl
         select case(shell_types(i))
         case(0) 
            do k=1,prim_per_shell(i)
               iprim = iprim+1
               Ra(iprim) = shell2atom(i)
               l=l+1
               Alpha(iprim) = prim_exp(l)
               Ptyp(iprim) = 1
               TMN(iprim,:) = (/0,0,0/)
            end do
         case(1) 
            do k=1,prim_per_shell(i)
               l=l+1
               do j=1,3
                  iprim = iprim+1
                  Ra(iprim) = shell2atom(i)
                  Alpha(iprim) = prim_exp(l)
                  Ptyp(iprim) = 1+j
                  TMN(iprim,:) = 0
                  TMN(iprim,j) = 1
               end do
            end do
         case(2, -2)
            do k=1,prim_per_shell(i)
               l=l+1
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=5; TMN(iprim,:)=(/2,0,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=6; TMN(iprim,:)=(/0,2,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=7; TMN(iprim,:)=(/0,0,2/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=8; TMN(iprim,:)=(/1,1,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=9; TMN(iprim,:)=(/1,0,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=10; TMN(iprim,:)=(/0,1,1/)  
            end do
         case(3, -3)  
            do k=1,prim_per_shell(i)
               l=l+1
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=11; TMN(iprim,:)=(/3,0,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=12; TMN(iprim,:)=(/0,3,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=13; TMN(iprim,:)=(/0,0,3/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=14; TMN(iprim,:)=(/2,1,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=15; TMN(iprim,:)=(/2,0,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=16; TMN(iprim,:)=(/0,2,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=17; TMN(iprim,:)=(/1,2,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=18; TMN(iprim,:)=(/1,0,2/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=19; TMN(iprim,:)=(/0,1,2/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=20; TMN(iprim,:)=(/1,1,1/)  
            end do
         case(4, -4) 
            do k=1,prim_per_shell(i)
               l=l+1
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=21; TMN(iprim,:)=(/4,0,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=22; TMN(iprim,:)=(/0,4,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=23; TMN(iprim,:)=(/0,0,4/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=24; TMN(iprim,:)=(/3,1,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=25; TMN(iprim,:)=(/3,0,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=26; TMN(iprim,:)=(/1,3,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=27; TMN(iprim,:)=(/0,3,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=28; TMN(iprim,:)=(/1,0,3/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=29; TMN(iprim,:)=(/0,1,3/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=30; TMN(iprim,:)=(/2,2,0/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=31; TMN(iprim,:)=(/2,0,2/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=32; TMN(iprim,:)=(/0,2,2/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=33; TMN(iprim,:)=(/2,1,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=34; TMN(iprim,:)=(/1,2,1/)   
               iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
               Ptyp(iprim)=35; TMN(iprim,:)=(/1,1,2/)   
            end do 
         case(-1) 
            l1=l
            do k=1,prim_per_shell(i)
               l=l+1
               iprim = iprim+1
               Ra(iprim) = shell2atom(i)
               Alpha(iprim) = prim_exp(l)
               Ptyp(iprim) = 1
               TMN(iprim,:) = (/0,0,0/)
            end do
            l=l1
            do k=1,prim_per_shell(i)
               l=l+1   
               iprim = iprim+1
               Ra(iprim) = shell2atom(i)
               Alpha(iprim) = prim_exp(l)
               Ptyp(iprim) = 2
               TMN(iprim,:) = (/1,0,0/)
            end do
            l=l1
            do k=1,prim_per_shell(i)
               l=l+1   
               iprim = iprim+1
               Ra(iprim) = shell2atom(i)
               Alpha(iprim) = prim_exp(l)
               Ptyp(iprim) = 3
               TMN(iprim,:) = (/0,1,0/)
            end do
            l=l1
            do k=1,prim_per_shell(i)
               l=l+1   
               iprim = iprim+1
               Ra(iprim) = shell2atom(i)
               Alpha(iprim) = prim_exp(l)
               Ptyp(iprim) = 4
               TMN(iprim,:) = (/0,0,1/)
            end do
         case default
            write(*,*) "H orbitals are not implemented"
            stop  
         end select
      end do
      close(1)
   end subroutine filefchk
   
   !****************************************************************************
   !! Parses a custom .bas file for basis set primitive data.
   subroutine filebas(basfilename)
      use wfxinfo
      use geninfo
      use locatemod
      implicit none
      character(len=*), intent(in) :: basfilename
      integer :: i, j, k
      double precision, allocatable, dimension(:,:) :: bas
      
      write(*,*) "Reading basis set from .bas file"
      open(unit=1,file=basfilename,status='OLD')
      
      allocate(bas(nprim,3))
      allocate(Xn(nprim));allocate(Yn(nprim));allocate(Zn(nprim))
   
      call locate(1,"#")
      do i=1,nprim
         read(1,*) k, bas(i,1), bas(i,2), bas(i,3), Alpha(i), TMN(i,1), TMN(i,2), TMN(i,3)      
         Xn(i)=bas(i,1); Yn(i)=bas(i,2); Zn(i)=bas(i,3)
         do j=1,natoms
            if (((dabs(bas(i,1)-cartes(j,1))).lt.0.00001).and.(dabs((bas(i,2)-cartes(j,2))).lt.0.00001)&
            .and.(dabs((bas(i,3)-cartes(j,3))).lt.0.00001)) then
               Ra(i) = j
            end if
         end do   
      end do
      close(1)
      
      write(*,*) "Reading basis set from .bas file"
      do i=1,nprim
         write(*,*) cartes(int(Ra(i)),1), cartes(int(Ra(i)),2), cartes(int(Ra(i)),3), Alpha(i), TMN(i,1), TMN(i,2), TMN(i,3)
      end do
   end subroutine
end module read_files   