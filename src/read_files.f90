!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   In this subroutine we read an input file where we specify what type of
!   calculations we want to do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readinput() 
use inputdat
use cubeinfo
use locatemod
   implicit none
   character*80 :: name
   integer :: i
   call getarg(1,name)  !gets the name of the input file (writen in the prompt)
   name=trim(name)      !remove the blank spaces of the string 'name'
   open(unit=3,file=name,status='OLD') 
   !set up logical variables as false
   primcube=.false.; aocube=.false.; MOcube=.false.; denscube=.false.; gradient=.false.
   laplacian=.false.; intracalc=.false.
   !read input file
   if (located(3,"$wfxfile")) then
       readwfx=.true.
       read(3,*) wfxfilename
   else if (located(3,"$fchkfile")) then
       readfchk=.true.
       read(3,*) fchkfilename
   else
       write(*,*) "Warning! You must provide at least a .wfx or .fchk file for intracule calculations"
   end if
   rewind 3
   if (located(3,"$basfile")) then
       readbas=.true.
       read(3,*) basfilename
   end if    
   if (located(3,"$logfile")) then
       readlog=.true.
       read(3,*) logfilename
   end if
   rewind 3    
  
   if (located(3,"$cubefile")) then
       read(3,*) (center(i), i=1,3) !cube centered in (x,y,z)
       read(3,*) (step(i), i=1,3) !distance between points in the axis
       read(3,*) (np(i), i=1,3) !number of points for each axis
       rewind 3
         ! --- Cubefile printing options ---
       if (located(3,"$density")) then
           denscube = .true.
       end if
       if (located(3,"$primitive")) then
            primcube = .true.
            read(3,*) npr
       end if     
       rewind(3)
       if (readlog .and. located(3,"$AO")) then
            aocube = .true.
            read(3,*) cao
       end if 
       rewind(3)
       if (located(3,"$MO")) then
            MOcube = .true.
            read(3,*) mo
       end if     
       rewind(3)     
       if (located(3,"$gradient")) then
            gradient = .true.
       end if 
       rewind(3)     
       if (located(3,"$laplacian")) then
           laplacian = .true.
       end if
   end if
    
   if (located(3,"$intracule")) then
       intracalc=.true.
   end if
    
   if (located(3,"$HFhole")) then
       c1calc=.true.
   end if
   end subroutine readinput 
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine readintra()
   !read information from the intracule function     
   use intrainfo
   use locatemod
   use quadratures
   implicit none
   character*80 :: name
   integer :: i, j
   call getarg(1,name)  !gets the name of the input file (writen in the prompt)
   name=trim(name)      !remove the blank spaces of the string 'name'
   open(unit=3,file=name,status='OLD')  
   !read input for intracule
   if (located(3,'$Integral screening threshold')) then !threshold for integral screenings
      read(3,*) thresh
   else  
      thresh=1.0d-12 !default value
   end if    

   rewind(3)
   call locate(3,'$DM2P')
   read(3,*) dm2name           !name of the dm2p file
   !read(3,*) trsh1, trsh2      !thresholds used in DM2prim
   if (located(3,'$outname')) then
      read(3,*) outname           !name of the output file
   else
      outname=trim(name)//'.out' !default name
   end if

   if (located(3,'$nosym')) then
      nosym=.true.
      write(*,*) 'Deactivating symmetry in intracule calculations'
   else
      nosym=.false.
   end if
   rewind(3)
   if (located(3,'$radial_integral')) then
      radial_integral=.true.
      if (located(3,'$definite')) then
         definite=.true.
         read(3,*) a,b
      end if
      if (located(3,'$Multicenter')) then
         if (definite) then
            write(*,*) 'Error! You cannot use definite integrals with multicenter quadrature'
            stop
         end if
         if (located(3,'$manual_grid')) then !number of quadrature centers(not automatically calculated)
            if (located(3,'$number of quadrature centers')) then
               read(3,*) nquad
            else
               write(*,*) 'Error! In manual mode you must provide the number of quadrature centers'
               stop
            end if
            allocate(cent(3,nquad))
            if (located(3,'$position of quadrature centers')) then               
               do i=1,nquad
                  read(3,*) (cent(j,i), j=1,3) !read position of centers
               end do
            else
               write(*,*) 'Error! In manual mode you must provide the positions of the quadrature centers'
               stop
            end if   
            allocate(nradc(nquad))
            allocate(nangc(nquad))
            allocate(sfalpha(nquad))
            if (located(3,'$Gauss-Legendre')) then
               read(3,*) (nradc(i), i=1,nquad)
            else 
               nradc(:)=50 !default value
            end if
            if (located(3,'$scaling factor')) then
               read(3,*) sfalpha(:)
            else
               sfalpha(:)=1.0d0 !default value
            end if  
            if (located(3,'$Gauss-Lebedev')) then
               read(3,*) (nangc(i), i=1,nquad)
            else  
               nangc(:)=590 !default value
            end if  
             allocate(Ps(nquad)) !allocate weight of each centre
            if (located(3,'$Center weights')) then
               read(3,*) Ps(:) !the weigth of a positive center must be equal to the neg.
            else  
               Ps(:)=1.0d0 !default equal weights
            end if  
            rewind(3)
         else !automatically calculate the parameters (default)
            if (located(3,'$Gauss-Legendre')) then
               read(3,*) nrad
            else
               nrad=50
            end if
            if (located(3,'$Gauss-Lebedev')) then
               read(3,*) nang
            else
               nang=590
            end if
            if (located(3,'$hydro')) then
               nohydro=.false.
            else
               nohydro=.true. !do not include H atoms as centres per default
            end if
            call centercalc() !calculate the number of centers and their positions
         end if
      else !single center quadrature (Default)
         write(*,*) 'Using single center quadrature'
         nquad=1
         allocate(cent(3,1))
         allocate(nradc(1)); allocate(nangc(1)); allocate(sfalpha(1)); allocate(Ps(nquad))
         cent(1:3,1)=0.0d0
         if (located(3,'$Gauss-Legendre')) then
            read(3,*) nradc(1) 
         else !default values
            nradc(1)=100
         end if
         if (located(3,'$scaling factor')) then                  
            read(3,*) sfalpha(1)
         else !default values
            sfalpha(1)=1.0d0
         end if
         if (located(3,'$Gauss-Lebedev')) then
            read(3,*) nangc(1)
         else !default values
            nangc(1)=590
         end if
      end if    
   else if (located(3,'$radial_plot')) then
      radial_plot=.true.
      if (located(3,'$plot_name')) then
         read(3,*) r_plot_name
      else
         r_plot_name=trim(name//'.rad')
      end if
      if (located(3,'$number of blocks')) then         
         read(3,*) nblock
      else
         nblock=1 !default value
      end if   
      allocate(n_an_per_part(nblock))
      allocate(tart(2,nblock))
      allocate(stp(nblock))
      rewind(3)
      if (located(3,'$start-end-step-nang')) then
         do i=1,nblock
            read(3,*) tart(:,i), stp(i), n_an_per_part(i)
         end do
      else
         write(*,*) 'Error! You must provide start, end, step and number of angles per block for radial plots'
         stop
      end if
   else if (located(3,'$Vectorial_plot')) then     
      cubeintra=.true.  
      read(3,*) cubeintraname
      read(3,*) (center_i(i), i=1,3) !cube centered in (x,y,z)
      read(3,*) (step_i(i), i=1,3) !distance between points in the axis
      read(3,*) (np_i(i), i=1,3) !number of points for each axis
   else if (located(3,'$Intracule_at_zero')) then
      intracule_at_zero=.true.
   else if (located(3,'$Intracule_two_points')) then
      intracule_two_points=.true.   
      read(3,*) x_point1, y_point1, z_point1
      read(3,*) x_point2, y_point2, z_point2
   else
      write(*,*) 'Warning! You must provide at least a radial or vectorial plot option for intracule calculations'    
   end if
   rewind(3)
   close(3) 
end subroutine readintra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filewfx(wfxfilename)                                             !
!  1-Reads wfxfile and stores the important values in wfxinfo module        !
!    (nelec, noccmo, nprim, ...)                                            !
!  2-Determines the wavefunction type we have: RHF, UHF, correlated WF      !
!  +relaxed/unrelaxed density                                               !
!  3-Builds TMN matrix (matrix with primitive exponents)                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use wfxinfo
use geninfo
use locatemod
implicit none
character*40, intent(in) :: wfxfilename
character*80 :: line
character*80 :: zaborra  !coses que no volem llegir
integer :: i,j,k, kk
integer :: mon !molecular orbital number   
double precision :: maxim, minim
double precision, parameter :: trsh=1.d-16 !threeshold for orbital occupancies

 corr=.false.  
 uhf=.false.
 rhf=.false.
 udens=.false.
 rdens=.false.  
 opsh=.false.
 clsh=.false.
   
open(unit=1,file=wfxfilename,status='OLD') 

call locate(1,"Number of Nuclei") 
read(1,*) natoms ! del contingut de "line", guarda el que hi ha a les 3 primeres posicions en una variable
                        
rewind 1  

call locate(1,"Number of Occupied Molecular Orbitals")
read(1,*) noccmo

allocate(Occ(noccmo))

rewind 1

call locate(1,"Molecular Orbital Occupation Numbers") !check RHF,UHF, or Correlated
do i=1,noccmo                                         !by     
   read(1,*) Occ(i)                                   !reading the MO occupation values
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

maxim=maxval(Occ)                                         !store the maximum occupation value
minim=minval(Occ)   
if (corr) then
    if (maxim.ge. 2.d0) then             !MOs filled by pairs
        if ((maxim-2.d0 .gt. trsh).or.(minim.lt. 0.d0)) then 
           write(*,*) "RELAXED DENSITY (closed shell)"                  
           rdens=.true.
        else     
           write(*,*) "UNRELAXED DENSITY (closed shell)"!all the occupation values between 0 and 2
           udens=.true.
        end if                    
    else                                       !1 electron MOs (Unrestricted)              
        if ((maxim-1.d0 .gt. trsh).or.(minim.lt. 0.d0)) then 
           write(*,*) "RELAXED DENSITY (open shell)"
           rdens=.true.
        else    
           write(*,*) "UNRELAXED DENSITY (open shell)"
           udens=.true.
        end if                    
    end if 
else                    !HF mwavefunction, check if it's RHF or UHF 
   minim=minval(Occ,MASK=occ.gt.0.d0)     
   if (minim.eq.1.d0) then                 !we have 1e- in 1 orbital
         if (maxim.eq.2.d0) then           !at least 1 doubly occ orbital
            rhf=.true.                      !open shell RHF
            uhf=.false.
            opsh=.true.
            write(*,*) "Open shell RHF"
         else                                 !all singly occ. orbitals (alpha and beta)
            rhf=.false.                        !UHF   
            uhf=.true.  
            opsh=.true.           
            write(*,*) "UHF"
         end if
   else if (minim.eq.2.d0) then           !all doubly occ. orbitals
         rhf=.true.                           !closed shell rhf
         uhf=.false.  
         clsh=.true.
         write(*,*) "Closed shell RHF"
   else 
         write(*,*)   "Error, correlated wf"  
         STOP
   end if    
  !deallocate(Occ)                        !with HF wavefunctions we no longer need the occupancy   
end if  
 
call locate(1,"Net Charge")           !store some variables
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

allocate(T(noccmo,nprim))                  !read T matrix (MOs in prims)       
call locate(1,"Molecular Orbital Primitive Coefficients")
do kk=1,noccmo
  read(1,'(a80)')line
  read(1,*) mon                          !molecular orbital number "mo"
 read(1,'(a80)') zaborra  
 read(1,*) (T(mon,i), i=1,nprim)
end do

rewind 1

allocate(TMN(nprim,3))      !create a matrix with t m and n exponents for each primitive (see prim equation)
do i=1,nprim                                                                 !t->lx
   if (ptyp(i).eq.1) then   !s prim                                          !m->ly
     do j=1,3                   !n->lz
       TMN(i,j)=0
     end do
   else if (ptyp(i).eq.2) then  !px
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.3) then  !py
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=0   
   else if (ptyp(i).eq.4) then  !pz
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=1   
   else if (ptyp(i).eq.5) then  !dx²
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=0   
   else if (ptyp(i).eq.6) then  !dy²
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=0  
   else if (ptyp(i).eq.7) then  !dz²
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=2   
   else if (ptyp(i).eq.8) then  !dxy
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=0   
   else if (ptyp(i).eq.9) then  !dxz
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.10) then !dyz 
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=1      
   else if (ptyp(i).eq.11) then !fxxx 
      TMN(i,1)=3
      TMN(i,2)=0
      TMN(i,3)=0
   else if (ptyp(i).eq.12) then !fyyy 
      TMN(i,1)=0
      TMN(i,2)=3
      TMN(i,3)=0
   else if (ptyp(i).eq.13) then !fzzz 
      TMN(i,1)=0
      TMN(i,2)=0
      TMN(i,3)=3
   else if (ptyp(i).eq.14) then !fxxy 
      TMN(i,1)=2
      TMN(i,2)=1
      TMN(i,3)=0
   else if (ptyp(i).eq.15) then !fxxz 
      TMN(i,1)=2
      TMN(i,2)=0
      TMN(i,3)=1
   else if (ptyp(i).eq.16) then !fyyz 
      TMN(i,1)=0
      TMN(i,2)=2
      TMN(i,3)=1
   else if (ptyp(i).eq.17) then !fxyy 
      TMN(i,1)=1
      TMN(i,2)=2
      TMN(i,3)=0
   else if (ptyp(i).eq.18) then !fxzz 
      TMN(i,1)=1
      TMN(i,2)=0
      TMN(i,3)=2
   else if (ptyp(i).eq.19) then !fyzz 
      TMN(i,1)=0
      TMN(i,2)=1
      TMN(i,3)=2
   else if (ptyp(i).eq.20) then !fxyz 
      TMN(i,1)=1
      TMN(i,2)=1
      TMN(i,3)=1
   else
     write(*,*) "G orbitals are not implemented"
     stop
   end if
end do
                             !if we have an unrestricted wf-------What about correlated???
if (uhf) then                !split T matrix into T_a and T_b (alpha and beta)
     allocate(T_a(nalfae,nprim))
     allocate(T_b(nbetae,nprim))
     do i=1,nalfae                !for alpha molecular orbitals
          do j=1,nprim    
                T_a(i,j)= T(i,j)  !store the prim. coefficients of alpha orb.   
          end do
     end do
     k=nalfae+1                   !first beta orbital
     do i=1,nbetae     
          do j=1,nprim
                T_b(i,j) = T(k,j)     
          end do
          k=k+1                   !next beta orbital
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
subroutine filelog(logfilename)                                             !
!1-reads d (Flg) fixed coeficients of the primitives			       !
!2-Normalizes that coeficient (with respect to the primitives, depending on !
!  primitive type (ptyp)						       !
!3-Builds C matrix (expansion of MOs into AOs)                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use loginfo                    
use geninfo
use wfxinfo
use locatemod
implicit none
integer, parameter :: maxao=100
character*10 :: logfilename
character*80 :: line 
character*80 :: zaborra
character*10 :: tao !atomic orbital type
integer :: i,j,k,l,m
integer :: sm, smm, sm2, sm3, sm4 !sums !smp is the sum of the previous primitives
integer :: smatom
double precision :: kk !alpha value of the primitives (already got it from wfx) 

open(unit=4 ,file=logfilename, status='OLD') !obre arxiu amb nom name a la unitat 4

allocate(npao(maxao)) !number of primitives of each AO basis. 
allocate(Flg(nprim))  !d coeficient for each primitive
allocate(aotyp(maxao))!atomic orbital type

Flg=0.d0

smm=0 !sum primitives
sm=0  !sum hybrid AOs (S, SP, D, ...)
smatom=1
call locate(4,"AO basis set in the form of general basis input") 
read(4,'(a80)') zaborra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!read fixed coeficients of the AOs in the primitives!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do while (smatom.le.natoms)   
      read(4,'(a80)') line
      read(line(2:3),*) tao 
      if (tao.eq."**") then !new atom finished 
           smatom=smatom+1
           read(4,'(a80)') zaborra   
      else if (tao.eq."S") then
           sm=sm+1  !one more AO
           read(line(5:7),*) npao(sm) !number of primitives of S ao   
           do i=1,npao(sm)
                smm=smm+1 !one more primitive
                read(4,'(a80)') line
                read(line(24:39),*) Flg(smm) !read d coeficient of smm primitive  
           end do
           aotyp(sm)=0  !S type ao
       else if (tao.eq."SP") then 
            sm=sm+1 !one more ao, S type
            aotyp(sm)=0
            read(line(5:7),*) npao(sm)      
            do i=1,npao(sm)
                 smm=smm+1  !one more primitive
                 sm2=sm+1     !px AO
                 aotyp(sm2)=1
                 sm3=sm+2     !py AO
                 aotyp(sm3)=1
                 sm4=sm+3     !pz AO  
                 aotyp(sm4)=1  
                 k=smm+npao(sm) !k is the p type primitives, smm is for the S type primitives     
                 read(4,*) kk, Flg(smm), Flg(k) !read d coeficients for S and for P primitives
                 l=k+npao(sm) !py primitive number                  
                 m=l+npao(sm) !pz primitive number    
                 Flg(l)=Flg(k)     !same value(C_px1=C_py1=C_pz1, 
                                   !(C_px2=C_py2=C_pz2, ...)  
                 Flg(m)=Flg(k)     
            end do   
            smm=smm+3*npao(sm)     !sum all the primitives used in this P type orbital                
            do i=1,3
                 npao(sm+i)=npao(sm) !In an SP hybrid orbital, px, py and pz AOs have the same number of primitives than S
            end do
            sm=sm+3 !sum px, py and pz to the total number of AOs
       else if (tao.eq."D") then
            sm=sm+1  !one more ao
            aotyp=2 !dx, dy, dz  !aotyp is used to the normalization of d coeficient (see norm subroutine)
            read(line(5:7),*) npao(sm)   
            do i=1,npao(sm) !do for the number of primitives of this AO
                 smm=smm+1  
                 read(4,*) kk, Flg(smm)      
                 !check if we have cartesian or spherical coordinates    
                 do j=1,5  !six AOs with same npao and F coeficient (if cartesian)
                        k=smm+j
                        l=sm+j  !current AO 
                        if (l.le.3) then 
                                 aotyp(l)=2  !dx,dy or dz coordinate
                        else
                                 aotyp(l)=3  !square (different normalization)  
                        end if 
                        Flg(k)=Flg(smm)      !same coeficient for each primitive of this AO 
                        npao(sm+j)=npao(sm)    
                 end do
            end do    
            smm=smm+6*npao(sm)  !sum all the primitives and AOs used. 
            sm=sm+5     
       end if
end do
nao=sm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!Normalize for each d coeficient of each primitive!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(N_prim(nprim))  
do i=1,nprim
      if (ptyp(i).eq.1) then 
            N_prim(i)=(pi*(2*Alpha(i))**(-1))**(-3.d0/4.d0)
           write(*,*) N_prim(i)
      else if ((ptyp(i).gt.1).and.(ptyp(i).le.4)) then
            N_prim(i)= (pi*0.5d0)**(-3.d0/4.d0) * (2*(Alpha(i)**(1.25d0)))   
      else if ((ptyp(i).gt.4).and.(ptyp(i).le.7)) then  !dxy,dxz.dyz
            N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))
      else if ((ptyp(i).gt.7).and.(ptyp(i).le.10)) then  !dx²,dy².dz²
            N_prim(i)=(pi*0.5d0)**(-3.d0/4.d0) * (2**2 * Alpha(i)**(1.75d0))/(dsqrt(3.d0))   
      else
            write(*,*) "f type orbital, insert its normalization equation"
      end if    
end do
do i=1, nprim              !Multiply the normalization of the primitives
      Flg(i)=Flg(i)*N_prim(i)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!Create expansion coef. of MOs in AOs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cmatrix() 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!102 format(30(F4.2,1X))
!103 format(30(I4,1X)) 
 close(4)
end subroutine filelog

!*******************************************************************************
subroutine cmatrix()   
!!calculates C with F and T matrices. C are the contracted coefficients of AOs in
!MOs. F are coefficients of AOs in primitives and T are for MOs in Primirives.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

! 41 format(10(F7.3,1X))
end subroutine cmatrix

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
   
   ! Reset flags
   corr=.false.; uhf=.false.; rhf=.false.
   udens=.false.; rdens=.false.
   opsh=.false.; clsh=.false.

   open(unit=1,file=fchkfilename,status='OLD')
   ! ============================
   ! Determine calculation type
   ! ============================
   call getline_with(1,"SP",line)
   read(line,*) typ, method, basis
   call getline_with(1,"Number of alpha electrons",line) !just to check if UHF or RHF
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
   ! ============================
   ! Number of atoms
   ! ============================
   write(*,*) "Reading fchk file"
   call getline_with(1,"Number of atoms",line)
  ! write(*,*) line
   read(line,*) typ, typ, typ, typ, natoms   ! label has spaces → just grab last value

   ! ============================
   ! Atomic numbers
   ! ============================
   allocate(an(natoms))
   call locate(1,"Atomic numbers")
   read(1,*) (an(i), i=1,natoms)

   ! ============================
   ! Nuclear charges
   ! ============================
   allocate(chrg(natoms))
   chrg(:) = dble(an(:))

   ! ============================
   ! Cartesian coordinates
   ! ============================
   allocate(cartes(natoms,3))
   call locate(1,"Current cartesian coord")
   read(1,*) ((cartes(i,j), j=1,3), i=1,natoms)

   ! =========================================
   ! Number of contracted and primitive shells
   ! ==========================================
   call getline_with(1,"Number of contracted shells",line)
   read(line,*) typ, typ, typ, typ, typ, ncshl
   call getline_with(1,"Number of primitive shells",line)
   read(line,*) typ, typ, typ, typ, typ, npshl
   ! ============================
   ! Shell types
   ! ============================
   allocate(shell_types(ncshl))
   call locate(1,"Shell types")
   read(1,*) (shell_types(i), i=1,ncshl)

   ! ============================
   ! Primitives per shell
   ! ============================
   allocate(prim_per_shell(ncshl))
   call locate(1,"Number of primitives per shell")
   read(1,*) (prim_per_shell(i), i=1,ncshl)
   !calculate total number of pritimives
   !take into account the shell types for that
   nprim = 0
   do i=1,ncshl
      select case(shell_types(i))
      case(0)  ! s shell
         nprim = nprim + prim_per_shell(i)
      case(1)  ! p shell -> expand px, py, pz
         nprim = nprim + 3*prim_per_shell(i)
      case(2)  ! d shell -> Cartesian: xx, yy, zz, xy, xz, yz
         nprim = nprim + 6*prim_per_shell(i)
      case(3)  ! f shell -> Cartesian (10 functions)
         nprim = nprim + 10*prim_per_shell(i)   
      case(-1) ! SP shell: s + p
         nprim = nprim + 4*prim_per_shell(i)
      case(-2) ! spherical d shell
         nprim = nprim + 6*prim_per_shell(i)
      case(-3) ! spherical f shell
         nprim = nprim + 10*prim_per_shell(i)
      case(4)  ! g shell -> Cartesian (15 functions) 
         nprim = nprim + 15*prim_per_shell(i)
      case(-4) ! spherical g shell 
         nprim = nprim + 15*prim_per_shell(i)   
      case default
         write(*,*) "Shell type not implemented:", shell_types(i)
         stop
      end select
   end do
   !write(*,*) "nprim calculated from shells=", nprim
   ! ============================
   ! Shell-to-atom map
   ! ============================
   allocate(shell2atom(ncshl))
   call locate(1,"Shell to atom map")
   read(1,*) (shell2atom(i), i=1,ncshl)
   ! ============================
   ! Primitive exponents
   ! ============================
   allocate(prim_exp(npshl))
   call locate(1,"Primitive exponents")
   read(1,*) (prim_exp(i), i=1,npshl)
   ! ============================
   ! Expand primitives into basis
   ! ============================   
   allocate(Ra(nprim))
   allocate(Alpha(nprim))
   allocate(Ptyp(nprim))
   allocate(TMN(nprim,3))
   iprim = 0
   l=0 !counter for primitive exponents
   do i=1,ncshl
      !counter for primitive shells
      select case(shell_types(i))
      case(0) ! s shell
         do k=1,prim_per_shell(i)
            iprim = iprim+1
            Ra(iprim) = shell2atom(i)
            l=l+1
            Alpha(iprim) = prim_exp(l)
            Ptyp(iprim) = 1
            TMN(iprim,:) = (/0,0,0/)
            !write(*,*) "s shell", iprim
         end do
      case(1) ! p shell -> expand px, py, pz
         do k=1,prim_per_shell(i)
            l=l+1
            do j=1,3
               iprim = iprim+1
               Ra(iprim) = shell2atom(i)
               Alpha(iprim) = prim_exp(l)
               Ptyp(iprim) = 1+j
               TMN(iprim,:) = 0
               TMN(iprim,j) = 1
              ! write(*,*) "p shell", iprim
            end do
         end do
      case(2, -2)
          !consider also case -2   ! d shell -> Cartesian: xx, yy, zz, xy, xz, yz
         do k=1,prim_per_shell(i)
            l=l+1
            ! xx, yy, zz, xy, xz, yz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=5; TMN(iprim,:)=(/2,0,0/)   ! xx
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=6; TMN(iprim,:)=(/0,2,0/)   ! yy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=7; TMN(iprim,:)=(/0,0,2/)   ! zz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=8; TMN(iprim,:)=(/1,1,0/)   ! xy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=9; TMN(iprim,:)=(/1,0,1/)   ! xz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=10; TMN(iprim,:)=(/0,1,1/)  ! yz
            !write(*,*) "d shell", iprim
         end do
      case(3, -3)  ! f shell -> Cartesian (10 functions)
         do k=1,prim_per_shell(i)
            l=l+1
            ! xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=11; TMN(iprim,:)=(/3,0,0/)   ! xxx
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=12; TMN(iprim,:)=(/0,3,0/)   ! yyy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=13; TMN(iprim,:)=(/0,0,3/)   ! zzz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=14; TMN(iprim,:)=(/2,1,0/)   ! xxy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=15; TMN(iprim,:)=(/2,0,1/)   ! xxz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=16; TMN(iprim,:)=(/0,2,1/)   ! yyz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=17; TMN(iprim,:)=(/1,2,0/)   ! xyy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=18; TMN(iprim,:)=(/1,0,2/)   ! xzz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=19; TMN(iprim,:)=(/0,1,2/)   ! yzz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=20; TMN(iprim,:)=(/1,1,1/)  ! xyz
           ! write(*,*) "f shell", iprim
         end do
      case(4, -4) ! g shell -> Cartesian (15 functions) 
         do k=1,prim_per_shell(i)
            l=l+1
            ! xxxx, yyyy, zzzz, xxxy, xxxz, yyyx, yyyz, zzzx, zzzy,
            ! xxyy, xxzz, yyzz, xxyz, xyzz, xyyz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=21; TMN(iprim,:)=(/4,0,0/)   ! xxxx
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=22; TMN(iprim,:)=(/0,4,0/)   ! yyyy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=23; TMN(iprim,:)=(/0,0,4/)   ! zzzz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=24; TMN(iprim,:)=(/3,1,0/)   ! xxxy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=25; TMN(iprim,:)=(/3,0,1/)   ! xxxz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=26; TMN(iprim,:)=(/1,3,0/)   ! yyyx
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=27; TMN(iprim,:)=(/0,3,1/)   ! yyyz     
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=28; TMN(iprim,:)=(/1,0,3/)   ! zzzx
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=29; TMN(iprim,:)=(/0,1,3/)   ! zzzy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=30; TMN(iprim,:)=(/2,2,0/)   ! xxyy
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=31; TMN(iprim,:)=(/2,0,2/)   ! xxzz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=32; TMN(iprim,:)=(/0,2,2/)   ! yyzz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=33; TMN(iprim,:)=(/2,1,1/)   ! xxyz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=34; TMN(iprim,:)=(/1,2,1/)   ! xyyz
            iprim = iprim+1; Ra(iprim)=shell2atom(i); Alpha(iprim)=prim_exp(l)
            Ptyp(iprim)=35; TMN(iprim,:)=(/1,1,2/)   ! xyzz
            !write(*,*) "g shell", iprim            
         end do 
      case(-1) ! SP shell: s + p
         !write(*,*) "SP shell"
         l1=l
         do k=1,prim_per_shell(i)
            l=l+1
            ! first s
            iprim = iprim+1
            Ra(iprim) = shell2atom(i)
            Alpha(iprim) = prim_exp(l)
            Ptyp(iprim) = 1
            TMN(iprim,:) = (/0,0,0/)
         end do
         l=l1
         do k=1,prim_per_shell(i)
            l=l+1   
            ! now px
            iprim = iprim+1
            Ra(iprim) = shell2atom(i)
            Alpha(iprim) = prim_exp(l)
            Ptyp(iprim) = 2
            TMN(iprim,:) = (/1,0,0/)
         end do
         l=l1
         do k=1,prim_per_shell(i)
            l=l+1   
            !py
            iprim = iprim+1
            Ra(iprim) = shell2atom(i)
            Alpha(iprim) = prim_exp(l)
            Ptyp(iprim) = 3
            TMN(iprim,:) = (/0,1,0/)
         end do
         l=l1
         do k=1,prim_per_shell(i)
            l=l+1   
            !pz
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
   allocate(Xn(nprim));allocate(Yn(nprim));allocate(Zn(nprim))
   do i=1,nprim
      Xn(i)=cartes(Ra(i),1); Yn(i)=cartes(Ra(i),2); Zn(i)=cartes(Ra(i),3)
   end do
       
   !   write(*,*) cartes(Ra(i),1), cartes(Ra(i),2), cartes(Ra(i),3), Alpha(i), TMN(i,1), TMN(i,2), TMN(i,3)
   close(1)
end subroutine filefchk
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
   !read basis set info from .bas
   !overwrite .fchk info
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
   !write(*,*) "Reading basis set from .bas file"
   !do i=1,nprim
   !   write(*,*) cartes(int(Ra(i)),1), cartes(int(Ra(i)),2), cartes(int(Ra(i)),3), Alpha(i), TMN(i,1), TMN(i,2), TMN(i,3)
   !end do
end subroutine  
