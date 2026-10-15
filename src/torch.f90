!timreacttorch
!timrattorch
!init_torch
!tfactors(temp)
!screen_torch(y)
!rhstorch(y,rate,dydt,deriva)
!azbar_simple

      module weaktab
      implicit none


!      nuclide table indices of proton and neutron

      integer :: pnuci, nnuci

!      for efficiency we will statically assume a FFN grid
!      which notably uses an integer logarithmic grid for electron density
!      This has 13 temperature values, in K

      integer,parameter :: ntemp = 13
      double precision :: tempgrid(ntemp) = &
          (/ 0.01e9, 0.1e9, 0.2e9, 0.4e9, 0.7e9, 1.0e9, &
             1.5e9 , 2.0e9, 3.0e9, 5.0e9, 10.0e9, 30.0e9, 100.0e9 /)

!      and 11 values of log10 of e^- density
!      (in units of N_A per cm^3 ~= Ye*mass density in g/cc)

      integer, parameter :: nedens = 11
      double precision :: ledens_min=1.0, ledens_max=11.0


! table storage format in memory
!  separate arrays for e^+ decay (emission), e^- capture, pd-ec neutrino loss rates
!                      e^- decay (emission), e^+ capture, ed-pc neutrino loss rates
!  Entries at each density, temperature are stored by target nuclide in the same
!  order as the nuclide table ( aion(:), zion(:) etc ).
!  Entries at one density and temperature for all nuclides are stored contiguously to improve
!  cache-locality during an evalutation at a given density and temperature.

      integer :: numnuc
      integer, allocatable :: map_nuci_to_lepi(:),map_nuci_to_delepi(:)
      double precision, dimension(:,:,:), allocatable :: &
             wtab_pd, wtab_ec, wtab_pdecnu, &
             wtab_ed, wtab_pc, wtab_edpcnu

      integer, dimension(:), allocatable :: lep_prod, delep_prod

      integer :: tempi_saved

      save

      private nucname2num
      contains

! =========================================================================

!      subroutine to initialize memory tables from file
!      Nuclide indexes are derived to match the indexing of the zion and
!      aion arrays passed in.
!      The in-memory tables have one entry per target, but we don't know
!      how many targets we have rates for until we read in the table and
!      compare all the entries to our list of included nuclides.  So the
!      read is done in two stages.  First the rates are read into sparse
!      arrays that have an entry for each included nuclide.  A map is then
!      constructed from nuclide index to a sequence in which nuclides
!      lacking rates do not appear.  There is a separate sequence for
!      leptonization and deleptonization rates, since a given target may
!      have one and not the other.  Finally the rate tables are re-formatted
!      into a compact form based on this compact sequence.
!
!      Why? -- we want the rate tables to be in the same order as the
!      nuclide table, since that is the order the iteration will take
!      during rate calculation.  This improves cache locality.  Further,
!      the tables should be dense for reasons of both cache locality and
!      memory footprint.  We could scan the file twice instead of using
!      the sparse arrays.

      subroutine weaktab_init(netnumnuc,zion,aion)

      implicit none

      integer, intent(in) :: netnumnuc
      double precision,intent(in):: zion(:), aion(:)

!      sparse arrays for reading in data tables
      double precision, dimension(:,:,:), allocatable :: &
             s_wtab_pd, s_wtab_ec, s_wtab_pdecnu, &
             s_wtab_ed, s_wtab_pc, s_wtab_edpcnu

      integer :: i, j, edensi, tempi, delepi, lepi, weakfile, stat
      integer :: z1, a1, z2, a2, nuci1, nuci2
      integer :: numlep, numdelep

      character(len=5) :: name1, name2
      character(len=1024) :: buf

      double precision :: temp, ledens


!     allocate sparse tables based on number of nuclides
      numnuc=netnumnuc


      if (allocated(lep_prod)) deallocate(lep_prod)
      allocate( lep_prod(numnuc) )

      if (allocated(delep_prod)) deallocate(delep_prod)
      allocate( delep_prod(numnuc) )

      if (allocated(s_wtab_pd)) deallocate(s_wtab_pd)
      allocate( s_wtab_pd(numnuc,ntemp,nedens) )

      if (allocated(s_wtab_ec)) deallocate(s_wtab_ec)
      allocate( s_wtab_ec(numnuc,ntemp,nedens) )

      if (allocated(s_wtab_pdecnu)) deallocate(s_wtab_pdecnu)
      allocate( s_wtab_pdecnu(numnuc,ntemp,nedens) )

      if (allocated(s_wtab_ed)) deallocate(s_wtab_ed)
      allocate( s_wtab_ed(numnuc,ntemp,nedens) )

      if (allocated(s_wtab_pc)) deallocate(s_wtab_pc)
      allocate( s_wtab_pc(numnuc,ntemp,nedens) )

      if (allocated(s_wtab_edpcnu)) deallocate(s_wtab_edpcnu)
      allocate( s_wtab_edpcnu(numnuc,ntemp,nedens) )


!      save neutron and proton indices for later use
       do i=1,numnuc
        if ( int(aion(i)).eq. 1 ) then
         if ( int(zion(i)).eq. 0 ) then
          nnuci=i
         else
          pnuci=i
         endif
        endif
       enddo


!      Now read table entries

!     accumulate which nuclides we have rates for (as "targets")
!     Different for leptonization (ed pc) and deleptonization (pd ec) directions
!     will also use later to get products when iterating through targets

      lep_prod(:)   = -1
      delep_prod(:) = -1

!     open file
      weakfile = 33
      open( unit=weakfile, file='weak_rates.dat', status='old', iostat=stat)
      if (stat.ne.0) then
       print *, 'unable to open weak_rates.dat file'
       stop
      endif

!     read first line then iterate through file
!     read nuclide names of leptonization pair

      read(weakfile,"(a)",iostat=stat) buf

!      write(6,*) ' buf = ',trim(buf)


      do while ( stat .ge. 0 )

!      find end of first name on line
       do i=1,1024
        if (buf(i:i)==' ') exit
       enddo

!      find beginning of second name on line
       do j=i+1,1024
        if (buf(j:j)/=' ') exit
       enddo

       if (i>100.or.j>100) stop 'unable parse leptonization pair names'
       name1 = buf(1:i-1)
       name2 = buf(j:j+4)
       call nucname2num(name1, z1,a1)
       call nucname2num(name2, z2,a2)

!       write(6,*) 'after nucname2num z1 a1 name1: ',z1, a1, trim(name1)
!       write(6,*) 'after nucname2num z2 a2 name2: ',z2, a2, trim(name2)
!       read(5,*)


!      look for nuclide indices
       nuci1 = -1
       nuci2 = -1
       do i=1,numnuc
        if ( int(zion(i)).eq.z1 .and. int(aion(i)).eq.a1 ) then
         nuci1 = i
         exit
        endif
       enddo
       do i=1,numnuc
        if ( int(zion(i)).eq.z2 .and. int(aion(i)).eq.a2 ) then
         nuci2 = i
         exit
        endif
       enddo

!       write(6,*) ' nuci1 and nuci2: ',nuci1,nuci2
!       read(5,*)


!      read in rate table
!      Table file format:
!        One table for each nuclide pair The table is gridded in temperature
!        and electron density.
       if ( nuci1 .ne. -1 .and. nuci2 .ne. -1 ) then
        do edensi=1,nedens
         do tempi=1,ntemp
          read(weakfile,*) temp, ledens, &
                             s_wtab_pd(nuci1,tempi,edensi), &
                             s_wtab_ec(nuci1,tempi,edensi), &
                             s_wtab_pdecnu(nuci1,tempi,edensi), &
                             s_wtab_ed(nuci2,tempi,edensi), &
                             s_wtab_pc(nuci2,tempi,edensi), &
                             s_wtab_edpcnu(nuci2,tempi,edensi)
          delep_prod(nuci1) = nuci2
          lep_prod(nuci2) = nuci1
         enddo
        enddo
       else
        do edensi=1,nedens
         do tempi=1,ntemp
          read(weakfile,"(a)") buf
         enddo
        enddo
       endif

!      eat the "@" after each entry
       read(weakfile,"(a)") buf

!      read next leptonization pair
       read(weakfile,"(a)",iostat=stat) buf
      enddo
      close(weakfile)




!     Now pack the tables and create index maps

!     count lep and delep targets
!     and create map from nuclide index to index in lep/delep tables

      if (allocated(map_nuci_to_lepi)) deallocate(map_nuci_to_lepi)
      allocate(map_nuci_to_lepi(numnuc))
      if (allocated(map_nuci_to_delepi)) deallocate(map_nuci_to_delepi)
      allocate(map_nuci_to_delepi(numnuc))

      numlep   = 0
      numdelep = 0

      do i=1,numnuc
       if ( lep_prod(i) > 0 ) then
        numlep = numlep + 1
        map_nuci_to_lepi(i) = numlep
       else
        map_nuci_to_lepi(i) = -1
       endif
       if ( delep_prod(i) > 0 ) then
        numdelep = numdelep + 1
        map_nuci_to_delepi(i) = numdelep
       else
        map_nuci_to_delepi(i) = -1
       endif
      enddo

!     allocate compact tables

      if (allocated(wtab_pd)) deallocate(wtab_pd)
      allocate( wtab_pd(numdelep,ntemp,nedens) )

      if (allocated(wtab_ec)) deallocate(wtab_ec)
      allocate( wtab_ec(numdelep,ntemp,nedens) )

      if (allocated(wtab_pdecnu)) deallocate(wtab_pdecnu)
      allocate( wtab_pdecnu(numdelep,ntemp,nedens) )

      if (allocated(wtab_ed)) deallocate(wtab_ed)
      allocate( wtab_ed(numlep,ntemp,nedens) )

      if (allocated(wtab_pc)) deallocate(wtab_pc)
      allocate( wtab_pc(numlep,ntemp,nedens) )

      if (allocated(wtab_edpcnu)) deallocate(wtab_edpcnu)
      allocate( wtab_edpcnu(numlep,ntemp,nedens) )


!     fill compact tables
      do edensi=1,nedens
       do tempi=1,ntemp
        do i=1,numnuc
         delepi = map_nuci_to_delepi(i)
         if (delepi .gt. 0 ) then
          wtab_pd(delepi,tempi,edensi) = s_wtab_pd(i,tempi,edensi)
          wtab_ec(delepi,tempi,edensi) = s_wtab_ec(i,tempi,edensi)
          wtab_pdecnu(delepi,tempi,edensi)=s_wtab_pdecnu(i,tempi,edensi)
         endif
         lepi = map_nuci_to_lepi(i)
         if (lepi .gt. 0 ) then
          wtab_ed(lepi,tempi,edensi) = s_wtab_ed(i,tempi,edensi)
          wtab_pc(lepi,tempi,edensi) = s_wtab_pc(i,tempi,edensi)
          wtab_edpcnu(lepi,tempi,edensi) =s_wtab_edpcnu(i,tempi,edensi)
         endif
        enddo
       enddo
      enddo

!     deallocate sparse tables and info arrays
      deallocate( s_wtab_pd )
      deallocate( s_wtab_ec )
      deallocate( s_wtab_pdecnu )
      deallocate( s_wtab_ed )
      deallocate( s_wtab_pc )
      deallocate( s_wtab_edpcnu )


!     test table
!     find ti45
!      do i=1,numnuc
!       if (zion(i)==22 .and. aion(i)==45) exit
!      enddo
!      print *, 'ti45 nuci: ', i
!      i = map_nuci_to_delepi(i)
!      print *, 'ti45 values: ',wtab_pd(i,6,3), wtab_ec(i,6,3)
!      stop

      end subroutine weaktab_init





!===================================================================

      subroutine nucname2num(nname, z, a)
      implicit none
      character(len=5),intent(in) :: nname
      integer,intent(out)         :: z, a

      integer, parameter     :: zsymbnum=85
      character(len=2)       :: zsymb(zsymbnum), element

      data  zsymb/'h ','he','li','be','b ','c ','n ','o ','f ','ne', &
                  'na','mg','al','si','p ','s ','cl','ar','k ','ca', &
                  'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
                  'ga','ge','as','se','br','kr','rb','sr','y ','zr', &
                  'nb','mo','tc','ru','rh','pd','ag','cd','in','sn', &
                  'sb','te','i ','xe','cs','ba','la','ce','pr','nd', &
                  'pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
                  'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
                  'tl','pb','bi','po','at'/

!     special cases
      if ( trim(nname)=='p') then
       z=1
       a=1
       return
      endif
      if ( trim(nname)=='n') then
       z=0
       a=1
       return
      endif

!     split name into element and baryon number
!     element name is either 2 or 1 characters
      if ( (ichar(nname(2:2)).ge.ichar('a')) &
                        .and.(ichar(nname(2:2)).le.ichar('z')) ) then
!      2 character element name
       element = nname(1:2)
       read(nname(3:5),*) a
      else
!      1 character element name
       element(1:1)=nname(1:1)
       element(2:2)=' '
       read(nname(2:4),*) a
      endif

      do z=1,zsymbnum
       if ( element .eq. zsymb(z) ) exit
      enddo

      if ( z .gt. zsymbnum) then
       print *, 'unable to find element name: ', element
       stop
      endif

      end subroutine nucname2num




!===================================================================
!     subroutine to fill the relevant entries in the rate arrays for
!     use in the main torch routines.  The e^+/- decay rates from
!     reaclib are assumed to already be put in place, so if we don't
!     have a rate for a particular target nuclide, we should leave it
!     alone.

!     We pass the full rate arrays to keep the compiler from doing
!     any unnecessary copying.  We will only be setting the entries
!     r_tab(5,:) and r_tab(6,:)
!     r_tab(5,i) is the leptonization rate for target nuclide i per nuclei i
!     r_tab(6,i) is the rate at which nuclide i is produced by deleptonization
!                    per parent target nuclei k


      subroutine weaktab_rates(dens,temp,ye, &
                     r_np, drdt_np, drdd_np, &
                     r_pn, drdt_pn, drdd_pn, &
                     r_tab, drdt_tab, drdd_tab)

      implicit none

! declare the pass
      double precision, intent(in) :: dens, temp, ye
      double precision, intent(inout) :: r_np, drdt_np, drdd_np, &
                                          r_pn, drdt_pn, drdd_pn, &
                           r_tab(:,:), drdt_tab(:,:), drdd_tab(:,:)

! local variables

      integer :: tempi, ledensi, i, lepi, delepi, ihi, ilo, prodi
      double precision :: ledens, tw, dw, dt_inv, lnten
      double precision :: dlrdt_ed, dlrdled_ed, r_ed, drdt_ed, drdd_ed
      double precision :: dlrdt_pc, dlrdled_pc, r_pc, drdt_pc, drdd_pc
      double precision :: dlrdt_pd, dlrdled_pd, r_pd, drdt_pd, drdd_pd
      double precision :: dlrdt_ec, dlrdled_ec, r_ec, drdt_ec, drdd_ec

      lnten  = log(10.0e0)
      ledens = log10(ye*dens)


!     cut out early if we are out of bounds

      if ( temp < tempgrid(1) .or. temp > tempgrid(ntemp) .or. &
           ledens < ledens_min .or. ledens > ledens_max ) return


!     Determine indices and weighting factors

      ledensi = min(nedens-1,int(floor(ledens)-ledens_min)+1)


!     Check if saved tempi is still valid
!     otherwise re-find temp index via bisection

      tempi = tempi_saved
      if ( temp<tempgrid(tempi) .or. temp>=tempgrid(tempi+1) ) then
       ilo = 1
       ihi = ntemp-1
       do while (ihi.ne.ilo)
        tempi = ilo+(ihi-ilo+1)/2
        if ( temp<tempgrid(tempi) ) then
         ihi = tempi-1
        else
         ilo = tempi
        endif
       enddo
      endif
      tempi_saved = tempi


!     Calculate weighting factors for linear interpolation
!     these have units since they will be multiplied by derivatives

      tw = temp-tempgrid(tempi)
      dw = ledens-ledensi
      dt_inv = 1.0e0/(tempgrid(tempi+1)-tempgrid(tempi))


!     Step through all nuclides and fill rate structures
!     rates are left unchanged if we don't have an entry


!      write(6,*) 'before loop nnuci pnuci ',nnuci,pnuci

      do i=1,numnuc

!       write(6,*) 'in loop i nnuci pnuci', i,nnuci,pnuci

!  First do leptonizing direction
       lepi = map_nuci_to_lepi(i)
       if (lepi > 0) then
        dlrdt_ed   = (wtab_ed(lepi,tempi+1,ledensi) - wtab_ed(lepi,tempi,ledensi))*dt_inv
        dlrdled_ed = (wtab_ed(lepi,tempi,ledensi+1) - wtab_ed(lepi,tempi,ledensi))
        r_ed       = 10**(wtab_ed(lepi,tempi,ledensi) + tw*dlrdt_ed+dw*dlrdled_ed)
        drdt_ed    = lnten*r_ed*dlrdt_ed
        drdd_ed    = r_ed/dens*dlrdled_ed

        dlrdt_pc   = (wtab_pc(lepi,tempi+1,ledensi) - wtab_pc(lepi,tempi,ledensi))*dt_inv
        dlrdled_pc = (wtab_pc(lepi,tempi,ledensi+1)  - wtab_pc(lepi,tempi,ledensi))
        r_pc       = 10**(wtab_pc(lepi,tempi,ledensi) + tw*dlrdt_pc+dw*dlrdled_pc)
        drdt_pc    = lnten*r_pc*dlrdt_pc
        drdd_pc    = r_pc/dens*dlrdled_pc

!       neutron rates are stored separately
        if (i .eq. nnuci) then
         r_np    = r_ed + r_pc
         drdt_np = drdt_ed + drdt_pc
         drdd_np = drdd_ed + drdd_pc

!         write(6,*) 'r_np ',r_np
!         read(5,*) 

        else
         r_tab(5,i)    = r_ed + r_pc
         drdt_tab(5,i) = drdt_ed + drdt_pc
         drdd_tab(5,i) = drdd_ed + drdd_pc

!         write(6,*) 'setting r_tab(5,i) ',i

        endif

       endif

!  Now do deleptonizing direction
       delepi = map_nuci_to_delepi(i)
       if (delepi > 0) then
        dlrdt_pd   = (wtab_pd(delepi,tempi+1,ledensi) - wtab_pd(delepi,tempi,ledensi))*dt_inv
        dlrdled_pd = (wtab_pd(delepi,tempi,ledensi+1) - wtab_pd(delepi,tempi,ledensi))
        r_pd       = 10**(wtab_pd(delepi,tempi,ledensi) + tw*dlrdt_pd+dw*dlrdled_pd)
        drdt_pd    = lnten*r_pd*dlrdt_pd
        drdd_pd    = r_pd/dens*dlrdled_pd

        dlrdt_ec   = (wtab_ec(delepi,tempi+1,ledensi) - wtab_ec(delepi,tempi,ledensi))*dt_inv
        dlrdled_ec = (wtab_ec(delepi,tempi,ledensi+1) - wtab_ec(delepi,tempi,ledensi))
        r_ec       = 10**(wtab_ec(delepi,tempi,ledensi) + tw*dlrdt_ec+dw*dlrdled_ec)
        drdt_ec    = lnten*r_ec*dlrdt_ec
        drdd_ec    = r_ec/dens*dlrdled_ec

!       proton rates are stored separately
        if (i .eq. pnuci) then
         r_pn    = r_pd+r_ec
         drdt_pn = drdt_pd+drdt_ec
         drdd_pn = drdt_pd+drdt_ec

!         write(6,*) 'r_pn ', r_pn
!         read(5,*) 


        else
         prodi            = delep_prod(i)
         r_tab(6,prodi)   = r_pd+r_ec
         drdt_tab(6,prodi) = drdt_pd+drdt_ec
         drdd_tab(6,prodi) = drdd_pd+drdd_ec

!         write(6,*) 'setting r_tab(6,i) ',i
        endif

       endif


!      to spot check
!       if (i==154) then
!        print *,'log10 ni56 pd ec ed pc rates:',
!     1     log10(r_pd), log10(r_ec), log10(r_ed), log10(r_pc)
!        print *,lepi,delepi
!       endif

      enddo

      end subroutine

      end module weaktab
      




  !     Get Timmes's reaction rates from MG solution
      subroutine torch_react(nd, doweak, doneut, pa, ydot)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'vector_eos.dek'      
      include 'burn_common.dek'
!     Input
      integer   :: nd,doweak, doneut
      real*8    :: pa(110), ydot(100)
!     Local
      logical   :: deriva 
      integer   :: i      
      real*8    :: y(250), dydt(250),temp, den, pres
      real*8  ::  enuc,taud,taut,y1,r1,y2,r2,y3,r3,snupp, &
                  b1a,b1b,b2a,b2b,b3a,b3b,y4,b4a,b4b,y5,b5a,b5b, &
                  snucno,z,abar,zbar,ye,snuda,snudz
      deriva = .false.
! weak reactions on/off
      weak_on = doweak
!     Set temperature, density and pressure
      temp = pa(ionmax+4+nd)
      den = pa(1)
      pres = pa(nd+2)
      call tfactors(temp)
      bden = den
      btemp = temp
      bpres = pres
!     Set abundances
      do i=1,ionmax
         y(i) = pa(2+nd+i)
      enddo
!     Get rates
      call azbar_simple(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)
      call torchrat(ye)
!     Do screening
      screen_on = 1
      call screen_torch(y)      
      call rhstorch(y,sig,ratdum,dydt)
      do i=1,ionmax
         ydot(i) = dydt(i)
      enddo
      call ener_gener_rate(dydt, enuc)
      if (doneut .eq. 1) then
         call sneut5(btemp,bden,abar,zbar, &
              sneut,dsneutdt,dsneutdd,snuda,snudz)

         ! get the specific pp neutrino losses
         y1 = 0.0d0
         r1 = 0.0d0
         y2 = 0.0d0
         r2 = 0.0d0
         y3 = 0.0d0
         r3 = 0.0d0
         if (ih2 .ne. 0) then
            y1 = y(iprot)
            r1 = 0.5d0*ratdum(irpp)
         end if
         if (ibe7 .ne. 0  .and. ili7 .ne. 0) then
            y2 = y(ibe7)
            r2 = sig(6,ili7)
         end if
         if (ib8 .ne. 0) then
            y3 = y(ib8)
            r3 = ratdum(irb8ep)
         end if
         sneutpp = snupp(y1,r1,y2,r2,y3,r3)


         ! get the cno specific neutrino losses
         y1  = 0.0d0
         b1a = 0.0d0
         b1b = 0.0d0
         if (in13 .ne. 0  .and. ic13 .ne. 0) then
            y1  = y(in13)
            b1a = bion(ic13)
            b1b = bion(in13)
         end if

         y2  = 0.0d0
         b2a = 0.0d0
         b2b = 0.0d0
         if (io14 .ne. 0  .and. in14 .ne. 0) then
            y2  = y(io14)
            b2a = bion(in14)
            b2b = bion(io14)
         end if

         y3  = 0.0d0
         b3a = 0.0d0
         b3b = 0.0d0
         if (io15 .ne. 0  .and. in15 .ne. 0) then
            y3  = y(io15)
            b3a = bion(in15)
            b3b = bion(io15)
         end if

         y4  = 0.0d0
         b4a = 0.0d0
         b4b = 0.0d0
         if (if17 .ne. 0  .and. io17 .ne. 0) then
            y4  = y(if17)
            b4a = bion(io17)
            b4b = bion(if17)
         end if

         y5  = 0.0d0
         b5a = 0.0d0
         b5b = 0.0d0
         if (if18 .ne. 0  .and. io18 .ne. 0) then
            y5  = y(if18)
            b5a = bion(io18)
            b5b = bion(if18)
         end if

         sneutcno = snucno(y1,b1a,b1b,y2,b2a,b2b, &
              y3,b3a,b3b,y4,b4a,b4b, &
              y5,b5a,b5b)


         sneut = sneut + sneutpp + sneutcno
         enuc = enuc - sneut
      endif
      ydot(ionmax+1) = enuc
      return
    end subroutine torch_react

      subroutine ener_gener_rate(dydt,enuc)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! computes the instantaneous energy generation rate

! declare the pass
      real*8  :: dydt(*),enuc

! local variables
       integer :: i
       real*8 :: a1,a2,a3,xx

       integer :: itmp(abignet)
       real*8  :: tmp(abignet)

! conversion factors for the nuclear energy generation rate
! detlap is the mass excess of the proton in MeV
! detlan is the mass excess of the neutron in MeV

      real*8, parameter :: enuc_conv  = ev2erg*1.0d6*avo, &
                           enuc_conv2 = -avo*clight*clight, &
                           deltap     = 7.288969d0, &
                           deltan     = 8.071323d0, &
                           mev2gr = 1.7826617585133994e-27
      
      
      xx = 0.0q0
      if (weak_on .eq. 1) then
         do i=1,ionmax
            a1 = dydt(i)
            a2 = mion(i)
            xx = xx + a1*a2
         enddo
      else
         do i=1,ionmax
            a1 = dydt(i)
            a2 = - bion(i)*mev2gr
            xx = xx + a1*a2
         enddo
      endif
      enuc = enuc_conv2 * xx
!      write(6,'(1x,1pe24.16)') enuc

      return
    end subroutine ener_gener_rate

      subroutine sneut5(temp,den,abar,zbar, &
                        snu,dsnudt,dsnudd,dsnuda,dsnudz)
      include 'implno.dek'
      include 'const.dek'

! this routine computes thermal neutrino losses from the analytic fits of
! itoh et al. apjs 102, 411, 1996 and also returns their derivatives.

! input:
! temp = temperature
! den  = density
! abar = mean atomic weight
! zbar = mean charge

! output:
! snu    = total neutrino loss rate in erg/g/sec
! dsnudt = derivative of snu with temperature
! dsnudd = derivative of snu with density
! dsnuda = derivative of snu with abar
! dsnudz = derivative of snu with zbar


! declare the pass
      double precision temp,den,abar,zbar, &
                       snu,dsnudt,dsnudd,dsnuda,dsnudz

! local variables
      double precision spair,spairdt,spairdd,spairda,spairdz, &
                       splas,splasdt,splasdd,splasda,splasdz, &
                       sphot,sphotdt,sphotdd,sphotda,sphotdz, &
                       sbrem,sbremdt,sbremdd,sbremda,sbremdz, &
                       sreco,srecodt,srecodd,srecoda,srecodz

      double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9, &
                       xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi, &
                       a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,c05,c06, &
                       c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
                       c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
                       dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0, &
                       f1,deni,tempi,abari,zbari,f2,f3,z,xmue,ye, &
                       dum,dumdt,dumdd,dumda,dumdz, &
                       gum,gumdt,gumdd,gumda,gumdz


! pair production
      double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt, &
                       zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3, &
                       xnum,xnumdt,xnumdd,xnumda,xnumdz, &
                       xden,xdendt,xdendd,xdenda,xdendz, &
                       fpair,fpairdt,fpairdd,fpairda,fpairdz, &
                       qpair,qpairdt,qpairdd,qpairda,qpairdz

! plasma
      double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
                       ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
                       fxy,fxydt,fxydd,fxyda,fxydz

! photo
      double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
                       sin3,sin4,sin5,last,xast, &
                       fphot,fphotdt,fphotdd,fphotda,fphotdz, &
                       qphot,qphotdt,qphotdd,qphotda,qphotdz

! brem
      double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5, &
                       t8m6, &
                       eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
                       fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
                       gbrem,gbremdt,gbremdd,gbremda,gbremdz, &
                       u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
                       fliq,fliqdt,fliqdd,fliqda,fliqdz, &
                       gliq,gliqdt,gliqdd,gliqda,gliqdz

! recomb
      double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz, &
                       nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz



! numerical constants
      double precision fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
      parameter        (fac1   = 5.0d0 * pi / 3.0d0, &
                        fac2   = 10.0d0 * pi, &
                        fac3   = pi / 5.0d0, &
                        oneth  = 1.0d0/3.0d0, &
                        twoth  = 2.0d0/3.0d0, &
                        con1   = 1.0d0/5.9302d0, &
                        sixth  = 1.0d0/6.0d0, &
                        iln10  = 4.342944819032518d-1)


! theta is sin**2(theta_weinberg) = 0.2223 (2019)
! xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
! change theta and xnufam if need be, and the changes will automatically
! propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3, &
                       tfac4,tfac5,tfac6
      parameter        (theta  = 0.2223d0, &
                        xnufam = 3.0d0, &
                        cv     = 0.5d0 + 2.0d0 * theta, &
                        cvp    = 1.0d0 - cv, &
                        ca     = 0.5d0, &
                        cap    = 1.0d0 - ca, &
                        tfac1  = cv*cv + ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp+cap*cap), &
                        tfac2  = cv*cv - ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp - cap*cap), &
                        tfac3  = tfac2/tfac1, &
                        tfac4  = 0.5d0 * tfac1, &
                        tfac5  = 0.5d0 * tfac2, &
                        tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)* &
                                 (cvp*cvp + 1.5d0*cap*cap))



! initialize
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0

      if (temp .lt. 1.0e7) return


! to avoid lots of divisions
      deni  = 1.0d0/den
      tempi = 1.0d0/temp
      abari = 1.0d0/abar
      zbari = 1.0d0/zbar


! some composition variables
      ye    = zbar*abari
      xmue  = abar*zbari




! some frequent factors
      t9     = temp * 1.0d-9
      xl     = t9 * con1
      xldt   = 1.0d-9 * con1
      xlp5   = sqrt(xl)
      xl2    = xl*xl
      xl3    = xl2*xl
      xl4    = xl3*xl
      xl5    = xl4*xl
      xl6    = xl5*xl
      xl7    = xl6*xl
      xl8    = xl7*xl
      xl9    = xl8*xl
      xlmp5  = 1.0d0/xlp5
      xlm1   = 1.0d0/xl
      xlm2   = xlm1*xlm1
      xlm3   = xlm1*xlm2
      xlm4   = xlm1*xlm3

      rm     = den*ye
      rmdd   = ye
      rmda   = -rm*abari
      rmdz   = den*abari
      rmi    = 1.0d0/rm

      a0     = rm * 1.0d-9
      a1     = a0**oneth
      zeta   = a1 * xlm1
      zetadt = -a1 * xlm2 * xldt
      a2     = oneth * a1*rmi * xlm1
      zetadd = a2 * rmdd
      zetada = a2 * rmda
      zetadz = a2 * rmdz

      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta




! pair neutrino section
! for reactions like e+ + e- => nu_e + nubar_e

! equation 2.8
      gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
      gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)

! equation 2.7

      a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
      a2     = 2.084d20 + 2.0d0*1.872d21*zeta

      if (t9 .lt. 10.0) then
       b1     = exp(-5.5924d0*zeta)
       b2     = -b1*5.5924d0
      else
       b1     = exp(-4.9924d0*zeta)
       b2     = -b1*4.9924d0
      end if

      xnum   = a1 * b1
      c      = a2*b1 + a1*b2
      xnumdt = c*zetadt
      xnumdd = c*zetadd
      xnumda = c*zetada
      xnumdz = c*zetadz

      if (t9 .lt. 10.0) then
       a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
       a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
      else
       a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2
       a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3
      end if

      b1   = 3.0d0*zeta2

      xden   = zeta3 + a1
      xdendt = b1*zetadt + a2*xldt
      xdendd = b1*zetadd
      xdenda = b1*zetada
      xdendz = b1*zetadz

      a1      = 1.0d0/xden
      fpair   = xnum*a1
      fpairdt = (xnumdt - fpair*xdendt)*a1
      fpairdd = (xnumdd - fpair*xdendd)*a1
      fpairda = (xnumda - fpair*xdenda)*a1
      fpairdz = (xnumdz - fpair*xdendz)*a1


! equation 2.6
      a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
      a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5)
      xnum   = 1.0d0/a1
      xnumdt = -xnum*xnum*a2

      a1     = 7.692d7*xl3 + 9.715d6*xlp5
      a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)

      c      = 1.0d0/a1
      b1     = 1.0d0 + rm*c

      xden   = b1**(-0.3d0)

      d      = -0.3d0*xden/b1
      xdendt = -d*rm*c*c*a2
      xdendd = d*rmdd*c
      xdenda = d*rmda*c
      xdendz = d*rmdz*c

      qpair   = xnum*xden
      qpairdt = xnumdt*xden + xnum*xdendt
      qpairdd = xnum*xdendd
      qpairda = xnum*xdenda
      qpairdz = xnum*xdendz



! equation 2.5
      a1    = exp(-2.0d0*xlm1)
      a2    = a1*2.0d0*xlm2*xldt

      spair   = a1*fpair
      spairdt = a2*fpair + a1*fpairdt
      spairdd = a1*fpairdd
      spairda = a1*fpairda
      spairdz = a1*fpairdz

      a1      = spair
      spair   = gl*a1
      spairdt = gl*spairdt + gldt*a1
      spairdd = gl*spairdd
      spairda = gl*spairda
      spairdz = gl*spairdz

      a1      = tfac4*(1.0d0 + tfac3 * qpair)
      a2      = tfac4*tfac3

      a3      = spair
      spair   = a1*a3
      spairdt = a1*spairdt + a2*qpairdt*a3
      spairdd = a1*spairdd + a2*qpairdd*a3
      spairda = a1*spairda + a2*qpairda*a3
      spairdz = a1*spairdz + a2*qpairdz*a3




! plasma neutrino section
! for collective reactions like gamma_plasmon => nu_e + nubar_e
! equation 4.6

      a1   = 1.019d-6*rm
      a2   = a1**twoth
      a3   = twoth*a2/a1

      b1   =  sqrt(1.0d0 + a2)
      b2   = 1.0d0/b1

      c00  = 1.0d0/(temp*temp*b1)

      gl2   = 1.1095d11 * rm * c00

      gl2dt = -2.0d0*gl2*tempi
      d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
      gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
      gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
      gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)


      gl    = sqrt(gl2)
      gl12  = sqrt(gl)
      gl32  = gl * gl12
      gl72  = gl2 * gl32
      gl6   = gl2 * gl2 * gl2


! equation 4.7
      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
      gum  = 1.0d0/gl2
      a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
      ftdt = a1*gl2dt
      ftdd = a1*gl2dd
      ftda = a1*gl2da
      ftdz = a1*gl2dz


! equation 4.8
      a1   = 8.6d0*gl2 + 1.35d0*gl72
      a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

      b1   = 225.0d0 - 17.0d0*gl + gl2
      b2   = -0.5d0*17.0d0*gl*gum + 1.0d0

      c    = 1.0d0/b1
      fl   = a1*c

      d    = (a2 - fl*b2)*c
      fldt = d*gl2dt
      fldd = d*gl2dd
      flda = d*gl2da
      fldz = d*gl2dz


! equation 4.9 and 4.10
      cc   = log10(2.0d0*rm)
      xlnt = log10(temp)

      xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
      xnumdt = -iln10*0.5d0*tempi
      a2     = iln10*sixth*rmi
      xnumdd = a2*rmdd
      xnumda = a2*rmda
      xnumdz = a2*rmdz

      xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
      xdendt = iln10*0.5d0*tempi
      xdendd = a2*rmdd
      xdenda = a2*rmda
      xdendz = a2*rmdz


! equation 4.11
      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
       fxy   = 1.0d0
       fxydt = 0.0d0
       fxydd = 0.0d0
       fxydz = 0.0d0
       fxyda = 0.0d0

      else

       a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
       a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)

       b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
       b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0

       c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
       if (c .eq. 0.0) then
        dumdt = 0.0d0
        dumdd = 0.0d0
        dumda = 0.0d0
        dumdz = 0.0d0
       else
        dumdt = xdendt + 1.25d0*xnumdt
        dumdd = xdendd + 1.25d0*xnumdd
        dumda = xdenda + 1.25d0*xnumda
        dumdz = xdendz + 1.25d0*xnumdz
       end if

       d   = 0.57d0 - 0.25d0*xnum
       a3  = c/d
       c00 = exp(-1.0d0*a3**2)

       f1  = -c00*2.0d0*a3/d
       c01 = f1*(dumdt + a3*0.25d0*xnumdt)
       c02 = f1*(dumdd + a3*0.25d0*xnumdd)
       c03 = f1*(dumda + a3*0.25d0*xnumda)
       c04 = f1*(dumdz + a3*0.25d0*xnumdz)

       fxy   = 1.05d0 + (a1 - b1)*c00
       fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
       fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
       fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
       fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04

      end if



! equation 4.1 and 4.5
      splas   = (ft + fl) * fxy
      splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
      splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
      splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
      splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

      a2      = exp(-gl)
      a3      = -0.5d0*a2*gl*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1

      a2      = gl6
      a3      = 3.0d0*gl6*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1


      a2      = 0.93153d0 * 3.0d21 * xl9
      a3      = 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*a1
      splasdd = a2*splasdd
      splasda = a2*splasda
      splasdz = a2*splasdz




! photoneutrino process section
! for reactions like e- + gamma => e- + nu_e + nubar_e
!                    e+ + gamma => e+ + nu_e + nubar_e
! equation 3.8 for tau, equation 3.6 for cc,
! and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau   =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10
       c01  = -4.524d8
       c02  = -6.088d6
       c03  =  4.269d7
       c04  =  5.172d7
       c05  =  4.910d7
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9
       c12  = -3.304d9
       c13  = -1.031d9
       c14  = -1.764d9
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9
       c23  = -1.695d9
       c24  = -2.865d9
       c25  = -3.395d9
       c26  = -3.418d9
       dd01 = -1.135d8
       dd02 =  1.256d8
       dd03 =  5.149d7
       dd04 =  3.436d7
       dd05 =  1.005d7
       dd11 =  1.652d9
       dd12 = -3.119d9
       dd13 = -1.839d9
       dd14 = -1.458d9
       dd15 = -8.956d8
       dd21 = -1.548d10
       dd22 = -9.338d9
       dd23 = -5.899d9
       dd24 = -3.035d9
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8
       c03  =  2.236d8
       c04  =  1.580d8
       c05  =  2.165d8
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11
       c13  = -1.765d11
       c14  = -1.867d11
       c15  = -1.983d11
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9
       c23  = -7.967d9
       c24  = -7.932d9
       c25  = -7.987d9
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8
       dd03 =  2.242d8
       dd04 =  7.937d7
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11
       dd14 = -1.273d11
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

      taudt = iln10*tempi


! equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)
      xast = sin(fac2*tau)

      a0 = 0.5d0*c00 &
           + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
           + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
           + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

      f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 &
           + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 &
           - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
           - c05*sin5*5.0d0 + dd05*cos5*5.0d0) &
           - 0.5d0*c06*xast*fac2*taudt

      a1 = 0.5d0*c10 &
           + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
           + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
           + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

      f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 &
           + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 &
           - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 &
           + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

      a2 = 0.5d0*c20 &
           + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
           + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
           + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

      f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 &
           + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 &
           - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 &
           + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

! equation 3.4
      dum   = a0 + a1*zeta + a2*zeta2
      dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
      dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
      dumda = a1*zetada + 2.0d0*a2*zeta*zetada
      dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz

      z      = exp(-cc*zeta)

      xnum   = dum*z
      xnumdt = dumdt*z - dum*z*cc*zetadt
      xnumdd = dumdd*z - dum*z*cc*zetadd
      xnumda = dumda*z - dum*z*cc*zetada
      xnumdz = dumdz*z - dum*z*cc*zetadz

      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3

      dum    = 3.0d0*zeta2
      xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 &
               + 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
      xdendd = dum*zetadd
      xdenda = dum*zetada
      xdendz = dum*zetadz

      dum      = 1.0d0/xden
      fphot   = xnum*dum
      fphotdt = (xnumdt - fphot*xdendt)*dum
      fphotdd = (xnumdd - fphot*xdendd)*dum
      fphotda = (xnumda - fphot*xdenda)*dum
      fphotdz = (xnumdz - fphot*xdendz)*dum


! equation 3.3
      a0     = 1.0d0 + 2.045d0 * xl
      xnum   = 0.666d0*a0**(-2.066d0)
      xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt

      dum    = 1.875d8*xl + 1.653d8*xl2 + 8.499d8*xl3 - 1.604d8*xl4
      dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 &
               - 4.0d0*1.604d8*xl3)

      z      = 1.0d0/dum
      xden   = 1.0d0 + rm*z
      xdendt =  -rm*z*z*dumdt
      xdendd =  rmdd*z
      xdenda =  rmda*z
      xdendz =  rmdz*z

      z      = 1.0d0/xden
      qphot = xnum*z
      qphotdt = (xnumdt - qphot*xdendt)*z
      dum      = -qphot*z
      qphotdd = dum*xdendd
      qphotda = dum*xdenda
      qphotdz = dum*xdendz

! equation 3.2
      sphot   = xl5 * fphot
      sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
      sphotdd = xl5*fphotdd
      sphotda = xl5*fphotda
      sphotdz = xl5*fphotdz

      a1      = sphot
      sphot   = rm*a1
      sphotdt = rm*sphotdt
      sphotdd = rm*sphotdd + rmdd*a1
      sphotda = rm*sphotda + rmda*a1
      sphotdz = rm*sphotdz + rmdz*a1

      a1      = tfac4*(1.0d0 - tfac3 * qphot)
      a2      = -tfac4*tfac3

      a3      = sphot
      sphot   = a1*a3
      sphotdt = a1*sphotdt + a2*qphotdt*a3
      sphotdd = a1*sphotdd + a2*qphotdd*a3
      sphotda = a1*sphotda + a2*qphotda*a3
      sphotdz = a1*sphotdz + a2*qphotdz*a3

      if (sphot .le. 0.0) then
       sphot   = 0.0d0
       sphotdt = 0.0d0
       sphotdd = 0.0d0
       sphotda = 0.0d0
       sphotdz = 0.0d0
      end if





! bremsstrahlung neutrino section
! for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!                    n  + n     => n + n + nu + nubar
!                    n  + p     => n + p + nu + nubar
! equation 4.3

      den6   = den * 1.0d-6
      t8     = temp * 1.0d-8
      t812   = sqrt(t8)
      t832   = t8 * t812
      t82    = t8*t8
      t83    = t82*t8
      t85    = t82*t83
      t86    = t85*t8
      t8m1   = 1.0d0/t8
      t8m2   = t8m1*t8m1
      t8m3   = t8m2*t8m1
      t8m5   = t8m3*t8m2
      t8m6   = t8m5*t8m1


      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)

! "weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then

! equation 5.3
       dum   = 7.05d6 * t832 + 5.12d4 * t83
       dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8

       z     = 1.0d0/dum
       eta   = rm*z
       etadt = -rm*z*z*dumdt
       etadd = rmdd*z
       etada = rmda*z
       etadz = rmdz*z

       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1
       etam3 = etam2 * etam1


! equation 5.2
       a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
       f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
       xnum  = 1.0d0/a0

       dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
       z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
       dumdt = z*etadt
       dumdd = z*etadd
       dumda = z*etada
       dumdz = z*etadz

       c00   = 1.26d0 * (1.0d0+etam1)
       z     = -1.26d0*etam2
       c01   = z*etadt
       c02   = z*etadd
       c03   = z*etada
       c04   = z*etadz

       z      = 1.0d0/dum
       xden   = c00*z
       xdendt = (c01 - xden*dumdt)*z
       xdendd = (c02 - xden*dumdd)*z
       xdenda = (c03 - xden*dumda)*z
       xdendz = (c04 - xden*dumdz)*z

       fbrem   = xnum + xden
       fbremdt = -xnum*xnum*f0 + xdendt
       fbremdd = xdendd
       fbremda = xdenda
       fbremdz = xdendz


! equation 5.9
       a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
       f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8

       z     = 1.0d0 + rm*1.0d-9
       dum   = a0*z
       dumdt = f0*z
       z     = a0*1.0d-9
       dumdd = z*rmdd
       dumda = z*rmda
       dumdz = z*rmdz

       xnum   = 1.0d0/dum
       z      = -xnum*xnum
       xnumdt = z*dumdt
       xnumdd = z*dumdd
       xnumda = z*dumda
       xnumdz = z*dumdz

       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8

       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8

       c02   = 4.59d-5 * t8**(-0.110d0)
       dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8

       z     = den**(0.656d0)
       dum   = c00*rmi  + c01  + c02*z
       dumdt = dd00*rmi + dd01 + dd02*z
       z     = -c00*rmi*rmi
       dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
       dumda = z*rmda
       dumdz = z*rmdz

       xden  = 1.0d0/dum
       z      = -xden*xden
       xdendt = z*dumdt
       xdendd = z*dumdd
       xdenda = z*dumda
       xdendz = z*dumdz

       gbrem   = xnum + xden
       gbremdt = xnumdt + xdendt
       gbremdd = xnumdd + xdendd
       gbremda = xnumda + xdenda
       gbremdz = xnumdz + xdendz


! equation 5.1
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fbrem - tfac5*gbrem
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
       sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
       sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
       sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)




! liquid metal with c12 parameters (not too different for other elements)
! equation 5.18 and 5.16

      else
       u     = fac3 * (log10(den) - 3.0d0)
       a0    = iln10*fac3*deni

! compute the expensive trig functions of equation 5.21 only once
       cos1 = cos(u)
       cos2 = cos(2.0d0*u)
       cos3 = cos(3.0d0*u)
       cos4 = cos(4.0d0*u)
       cos5 = cos(5.0d0*u)

       sin1 = sin(u)
       sin2 = sin(2.0d0*u)
       sin3 = sin(3.0d0*u)
       sin4 = sin(4.0d0*u)
       sin5 = sin(5.0d0*u)

! equation 5.21
       fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0 &
             - 0.05821d0*cos1 - 0.04969d0*sin1 &
             - 0.01089d0*cos2 - 0.01584d0*sin2 &
             - 0.01147d0*cos3 - 0.00504d0*sin3 &
             - 0.00656d0*cos4 - 0.00281d0*sin4 &
             - 0.00519d0*cos5

       c00 =  a0*(0.00945d0 &
             + 0.05821d0*sin1       - 0.04969d0*cos1 &
             + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
             + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
             + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
             + 0.00519d0*sin5*5.0d0)


! equation 5.22
       ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
             - 0.00944d0*cos1 - 0.02213d0*sin1 &
             - 0.01289d0*cos2 - 0.01136d0*sin2 &
             - 0.00589d0*cos3 - 0.00467d0*sin3 &
             - 0.00404d0*cos4 - 0.00131d0*sin4 &
             - 0.00330d0*cos5

       c01 = a0*(-0.02342d0 &
             + 0.00944d0*sin1       - 0.02213d0*cos1 &
             + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
             + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
             + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
             + 0.00330d0*sin5*5.0d0)


! equation 5.23
       gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
             - 0.00710d0*cos1 + 0.02300d0*sin1 &
             - 0.00028d0*cos2 - 0.01078d0*sin2 &
             + 0.00232d0*cos3 + 0.00118d0*sin3 &
             + 0.00044d0*cos4 - 0.00089d0*sin4 &
             + 0.00158d0*cos5

       c02 = a0*(-0.01259d0 &
             + 0.00710d0*sin1       + 0.02300d0*cos1 &
             + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
             - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
             - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
             - 0.00158d0*sin5*5.0d0)


! equation 5.24
       gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
             + 0.00356d0*cos1 + 0.01052d0*sin1 &
             - 0.00184d0*cos2 - 0.00354d0*sin2 &
             + 0.00146d0*cos3 - 0.00014d0*sin3 &
             + 0.00031d0*cos4 - 0.00018d0*sin4 &
             + 0.00069d0*cos5

       c03 = a0*(-0.00829d0 &
             - 0.00356d0*sin1       + 0.01052d0*cos1 &
             + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
             - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
             - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
             - 0.00069d0*sin5*5.0d0)


       dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
       dumdt = -dum*tempi
       dumdd = oneth*dum*deni
       dumda = -oneth*dum*abari
       dumdz = 2.0d0*dum*zbari

       gm1   = 1.0d0/dum
       gm2   = gm1*gm1
       gm13  = gm1**oneth
       gm23  = gm13 * gm13
       gm43  = gm13*gm1
       gm53  = gm23*gm1


! equation 5.25 and 5.26
       v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2

       w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
       a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2


! equation 5.19 and 5.20
       fliq   = v*fb + (1.0d0 - v)*ft
       fliqdt = a0*dumdt*(fb - ft)
       fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01
       fliqda = a0*dumda*(fb - ft)
       fliqdz = a0*dumdz*(fb - ft)

       gliq   = w*gb + (1.0d0 - w)*gt
       gliqdt = a1*dumdt*(gb - gt)
       gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
       gliqda = a1*dumda*(gb - gt)
       gliqdz = a1*dumdz*(gb - gt)


! equation 5.17
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fliq - tfac5*gliq
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
       sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
       sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
       sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)

      end if




! recombination neutrino section
! for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
! equation 6.11 solved for nu
      xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
      xnumdt = -1.50d0*xnum*tempi
      xnumdd = xnum*deni
      xnumda = -xnum*abari
      xnumdz = xnum*zbari

! the chemical potential
      nu   = ifermi12(xnum)

! a0 is d(nu)/d(xnum)
      a0 = 1.0d0/(0.5d0*zfermim12(nu))
      nudt = a0*xnumdt
      nudd = a0*xnumdd
      nuda = a0*xnumda
      nudz = a0*xnumdz

      nu2  = nu * nu
      nu3  = nu2 * nu

! table 12
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if


! equation 6.7, 6.13 and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then

       zeta   = 1.579d5*zbar*zbar*tempi
       zetadt = -zeta*tempi
       zetadd = 0.0d0
       zetada = 0.0d0
       zetadz = 2.0d0*zeta*zbari

       c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
       c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
       dum    = zeta*c00
       dumdt  = zetadt*c00 + zeta*c01*nudt
       dumdd  = zeta*c01*nudd
       dumda  = zeta*c01*nuda
       dumdz  = zetadz*c00 + zeta*c01*nudz


       z      = 1.0d0/dum
       dd00   = dum**(-2.25)
       dd01   = dum**(-4.55)
       c00    = a1*z + a2*dd00 + a3*dd01
       c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z


       z      = exp(c*nu)
       dd00   = b*z*(1.0d0 + d*dum)
       gum    = 1.0d0 + dd00
       gumdt  = dd00*c*nudt + b*z*d*dumdt
       gumdd  = dd00*c*nudd + b*z*d*dumdd
       gumda  = dd00*c*nuda + b*z*d*dumda
       gumdz  = dd00*c*nudz + b*z*d*dumdz


       z   = exp(nu)
       a1  = 1.0d0/gum

       bigj   = c00 * z * a1
       bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
       bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
       bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
       bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz


! equation 6.5
       z     = exp(zeta + nu)
       dum   = 1.0d0 + z
       a1    = 1.0d0/dum
       a2    = 1.0d0/bigj

       sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
       srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
       srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
       srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
       srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

      end if


! convert from erg/cm^3/s to erg/g/s
! comment these out to duplicate the itoh et al plots

      spair   = spair*deni
      spairdt = spairdt*deni
      spairdd = spairdd*deni - spair*deni
      spairda = spairda*deni
      spairdz = spairdz*deni

      splas   = splas*deni
      splasdt = splasdt*deni
      splasdd = splasdd*deni - splas*deni
      splasda = splasda*deni
      splasdz = splasdz*deni

      sphot   = sphot*deni
      sphotdt = sphotdt*deni
      sphotdd = sphotdd*deni - sphot*deni
      sphotda = sphotda*deni
      sphotdz = sphotdz*deni

      sbrem   = sbrem*deni
      sbremdt = sbremdt*deni
      sbremdd = sbremdd*deni - sbrem*deni
      sbremda = sbremda*deni
      sbremdz = sbremdz*deni

      sreco   = sreco*deni
      srecodt = srecodt*deni
      srecodd = srecodd*deni - sreco*deni
      srecoda = srecoda*deni
      srecodz = srecodz*deni


! the total neutrino loss rate
      snu    =  splas + spair + sphot + sbrem + sreco
      dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt
      dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd
      dsnuda =  splasda + spairda + sphotda + sbremda + srecoda
      dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz

      return
      end subroutine sneut5

      double precision function snupp(yp,ratepp,ybe7,ratebeec, &
                                      yb8,rateb8epnu)
      include 'implno.dek'
      include 'const.dek'

! computes approximate neutrino losses from pp chain reactions
! see page 142 of astro 289j notes for these loss formulas

! input:
! yp         = proton molar abbundance
! ratepp     = pp reaction rate
! ybe7       = be7 molar abundance
! ratebeec   = be7 electron capture reaction rate
! yb8        = b8 molar abundance
! rateb8epnu = b8 decay reaction rate


! declare the pass
      double precision yp,ratepp,ybe7,ratebeec,yb8,rateb8epnu


! local variables
      double precision pp1nu,pp2nu,pp3nu,conv
      parameter        (conv = ev2erg*1.0d6*avo)


! nu losses from p(p,e-nu)h2
      pp1nu  = yp*yp*ratepp * 0.5d0 * 0.263d0


! nu losses from be7(n=>p)li7
      pp2nu  = ybe7 * ratebeec * 0.81d0


! nu losses from b8(p=>n)be8=>2a
      pp3nu  = yb8 * rateb8epnu * 7.73d0

! sum the pp-chain neutrino losses and convert to erg/g/s
      snupp  = (pp1nu + pp2nu + pp3nu) * conv

      return
      end function snupp
      
      double precision function snucno(yn13,bc13,bn13,yo14,bn14,bo14, &
                                       yo15,bn15,bo15,yf17,bo17,bf17, &
                                       yf18,bo18,bf18)
      include 'implno.dek'
      include 'const.dek'

! computes approximate neutrino losses from cno cycle  reactions
! see page 142 of astro 289j notes for these loss formulas

! input:
! yn13 = n13 molar abundance
! bc13 = c13 binding energy in mev
! bn13 = n13 binding energy in mev
! yo14 = o14 molar abundance
! bn14 = n14 binding energy in mev
! bo14 = o14 binding energy in mev
! yo15 = o15 molar abundance
! bn15 = n15 binding energy in mev
! bo15 = o15 binding energy in mev
! yf17 = f17 molar abundance
! bo17 = o17 binding energy in mev
! bf17 = f17 binding energy in mev
! yf18 = f18 molar abundance
! bo18 = o18 binding energy in mev
! bf18 = f18 binding energy in mev


! declare the pass
      double precision yn13,bc13,bn13,yo14,bn14,bo14, &
                       yo15,bn15,bo15,yf17,bo17,bf17, &
                       yf18,bo18,bf18

! local variables
      double precision sum,sum2,enu13n,enu14o,enu15o,enu17f,enu18f, &
                       conv,lntwo,tm1,tm2,tm3,tm4,tm5
      parameter        (conv  = ev2erg*1.0d6*avo, &
                        lntwo = 0.693147181d0, &
                        tm1   = lntwo/597.9d0, &
                        tm2   = lntwo/70.606d0, &
                        tm3   = lntwo/124.0, &
                        tm4   = lntwo/64.49, &
                        tm5   = lntwo/6586.2)


! 13n(e+nu)13c
      sum    = bc13 - bn13 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu13n = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu13n = yn13 * enu13n * tm1



! hot cno cycle 14o(e+nu)14n
      sum    = bn14 - bo14 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu14o = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu14o = yo14 * enu14o * tm2


! 15o(e+nu)15n
      sum    = bn15 - bo15 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu15o = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu15o = yo15 * enu15o * tm3


! 17f(e+nu)17o
      sum    = bo17 - bf17 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu17f = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu17f  = yf17 * enu17f * tm4


! 18f(e+nu)18o
      sum    = bo18 - bf18 - 0.782d0 - 1.022d0
      sum    = 1.0d0 + sum/0.511d0
      sum2   = sum*sum
      enu18f = 0.5d0 * sum * 0.511d0 * (1.0d0 - 1.0d0/sum2) &
               * (1.0d0 - 1.0d0/(4.0d0*sum) - 1.0d0/(9.0d0*sum2))
      enu18f = yf18 * enu18f * tm5


! sum the cno cycle losses and convert to erg/g/s
      snucno = (enu13n + enu14o + enu15o + enu17f + enu18f) * conv

      return
      end function snucno

      real*8 function ifermi12(f)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the inverse
! fermi function of order 1/2 when it is equal to f.
! maximum error is 4.19d-9.   reference: antia apjs 84,101 1993
!
! declare
      integer  :: i,m1,k1,m2,k2
      real*8   :: f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
                           6.610132843877d2,   3.818838129486d1, &
                           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
                           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
                           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
    end function ifermi12

      real*8 function zfermi12(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi
! function of order 1/2 evaluated at x. maximum error is 5.47d-13.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx


! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 7, 7, 10, 11/
      data  (a1(i),i=1,8)/5.75834152995465d6,   1.30964880355883d7, &
                          1.07608632249013d7,   3.93536421893014d6, &
                          6.42493233715640d5,   4.16031909245777d4, &
                          7.77238678539648d2,   1.0d0/
      data  (b1(i),i=1,8)/6.49759261942269d6,   1.70750501625775d7, &
                          1.69288134856160d7,   7.95192647756086d6, &
                          1.83167424554505d6,   1.95155948326832d5, &
                          8.17922106644547d3,   9.02129136642157d1/
      data (a2(i),i=1,11)/4.85378381173415d-14, 1.64429113030738d-11, &
                          3.76794942277806d-9,  4.69233883900644d-7, &
                          3.40679845803144d-5,  1.32212995937796d-3, &
                          2.60768398973913d-2,  2.48653216266227d-1, &
                          1.08037861921488d0,   1.91247528779676d0, &
                          1.0d0/
      data (b2(i),i=1,12)/7.28067571760518d-14, 2.45745452167585d-11, &
                          5.62152894375277d-9,  6.96888634549649d-7, &
                          5.02360015186394d-5,  1.92040136756592d-3, &
                          3.66887808002874d-2,  3.24095226486468d-1, &
                          1.16434871200131d0,   1.34981244060549d0, &
                          2.01311836975930d-1, -2.14562434782759d-2/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermi12 = xx * rn/den

      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermi12 = x*sqrt(x)*rn/den
      end if
      return
      end function zfermi12

      real*8 function zfermim12(x)
      include 'implno.dek'
!
! this routine applies a rational function expansion to get the fermi function
! of order -1/2 evaluated at x. maximum error is 1.23d-12.
! reference: antia apjs 84,101 1993
!
! declare
      integer :: i,m1,k1,m2,k2
      real*8  :: x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
!
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end function zfermim12

      subroutine torch_rat(nd, pa, RR)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'burn_common.dek'
!     Input
      integer   ::  nd            
      real*8    ::  pa(*), RR(*)
!     Local
      logical   :: deriva
      integer   :: i      
      real*8    :: y(100), temp, den, pres, abar, zbar, ye
      deriva = .false.
      !     Set temperature, density and pressure
      temp = pa(ionmax+4+nd)      
      den = pa(1)
      pres = pa(nd+2)      
      call tfactors(temp)
      bden = den
      btemp = temp
      bpres = pres
!     Set abundances
      do i=1,46
         y(i) = pa(2+nd+i)
      enddo
      !Get rates
      call azbar_simple(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)
      call torchrat(ye)      
!     Do screening
      screen_on = 1
      call screen_torch(y)
      do i=1,290
         RR(i) = ratdum(i)
      enddo      
      return
    end subroutine torch_rat

      subroutine torch_init(niso, tnrat, eaion, ezion, ebion, mgnetname)
      use weaktab
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
! this routine initializes stuff for the torch network
!     Input
      character*200 :: mgnetname
!     Output
      integer :: niso, tnrat
      real*8 :: eaion(100), ezion(100), ebion(100)
! declare      
      integer, parameter :: fil14=14
      real*8, parameter  :: mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/clight**2
      logical            :: ibhere
      character*132      :: string,word
      character*10       :: aname
      integer            :: i,j,k,l,ii,jj,ll,llmin,llmax,mm,j1,nz,na,nn,kk, &
                            inta,intz,izpf,iapf,inpf,nppf
      logical, parameter :: icoords = .true.
      real*8             :: dum,qful,xx
      netname = mgnetname

! popular format statements
01    format(a,i4)
06    format(2i5,f10.4)
07    format(f10.3)
08    format(6f10.3)
09    format(1x,i4,i4,i4,'    ',a5)

! decide on the arrow orientation
! downarrow true puts neutron, protons, alfa to the end
! downarrow false (i.e uparrow) puts neutron, protons, alfa at the beginning
! in general downarrowtrue  is faster for dense linear algebra, and either
! orientation for sparse linear algebra. gift routines, however, do
! much better with uparrow (downarrow false).

      downarrow = .true.
      if (downarrow) then
       ionbeg = 1
      else
       ionbeg = 4
      end if

! set which_bdat
      which_bdat = 1
! imput MESA file      
      call zet_read_file


! open the nuclear reaction rate data file

      call read_partition_functions
      open(unit=fil14,file='BDATNEW',status='old')
      nppf  = 6

! put neutrons, protons and alfa first for up-arrow orientation
      if (.not.downarrow) then
       ineut        = 1
       aion(ineut)  = 1.0d0
       nion(ineut)  = 1.0d0
       zion(ineut)  = 0.0d0
       bion(ineut)  = 0.0d0
!       mion(ineut)  = nion(ineut)*mn + zion(ineut)*(mp+me)  - bion(ineut)*mev2gr
       mion(ineut)  = nion(ineut)*mn + zion(ineut)*mp  - bion(ineut)*mev2gr       
       wion(ineut)  = avo * mion(ineut)
       ionam(ineut) = 'neut'

       iprot        = 2
       aion(iprot)  = 1.0d0
       nion(iprot)  = 0.0d0
       zion(iprot)  = 1.0d0
       bion(iprot)  = 0.0d0
!       mion(iprot)  = nion(iprot)*mn + zion(iprot)*(mp+me) - bion(iprot)*mev2gr
       mion(iprot)  = nion(iprot)*mn + zion(iprot)*mp - bion(iprot)*mev2gr
       wion(iprot)  = avo * mion(iprot)
       ionam(iprot) = 'prot'

       ihe4         = 3
       aion(ihe4)   = 4.0d0
       nion(ihe4)   = 2.0d0
       zion(ihe4)   = 2.0d0
       bion(ihe4)   = 28.29603d0
!       mion(ihe4)   = nion(ihe4)*mn + zion(ihe4)*(mp+me) - bion(ihe4)*mev2gr
       mion(ihe4)   = nion(ihe4)*mn + zion(ihe4)*mp - bion(ihe4)*mev2gr
       wion(ihe4)  = avo * mion(ihe4)
       ionam(ihe4)  = 'he4'
      end if


! now start reading the nuclear reaction rate data
! i is the code number of element z(i),n(i).
! forward reaction rates are 
! j = 1 = ng                 j = 6 = an
! j = 2 = pn                 j = 7  = ag
! j = 3 = ground state b-    j = 8 for semi-empirical electron capture
! j = 4 = pg                 j = 9 for semi-empirical positron decay
! j = 5 = ap                 j = 10 for semi-empirical beta decay
! ic1(j,i) = type formula to be used to calculate rate
! ic2(j,i) = number of constants in fitting reaction j on species i
! cx(k) holds the fitting coefficients
! ic3(j,i) = where to start counting cx from



! initialize counters for the ensuing goto loop construction
      k        = 1
      i        = ionbeg - 1


! we keep returning here from various goto and loop constructions
60    i = i+1
      if (i .gt. abignet) stop 'abignet too small in init_torch'
!      write(6,*) 'reading', i


! read in z and a and the two fitting information arrays
      read(fil14,02) nz,na,(ic1(j,i),ic2(j,i), j=1,10)
02    format(2i6,20i3)
!      write(6,*) 'read line 1', nz,na
      zion(i) = nz
      aion(i) = na
      if (nz .eq. 99 .and. na .eq. 99) go to 120


! read temperature dependent partition function coefficients
! not used now as the rauscher & thieleman table is interpolated

      llmin = nppf*(i-1) + 1
      llmax = llmin + (nppf - 1)
      if (llmax .gt. 6*abignet) stop 'array as past bounds'
      if (nppf .eq. 5) read(fil14,03) nz,nn,bion(i),(as(ll),ll=llmin,llmax),jj,aname
      if (nppf .eq. 6) read(fil14,04) nz,nn,bion(i),(as(ll),ll=llmin,llmax),jj,aname
03    format(2i3,f11.4,f5.1,4e12.3,i2,a10)
04    format(2i3,f11.4,f5.1,5e12.3,i2,a10)
      nion(i) = nn
      call sqeeze(aname)
!      write(6,*) 'read line 2 ',aname,' ',nn


! decide if this isotope is in the network and branch accordingly

! no he4 from jina reaclib/bdat
      if (nz .eq. 2 .and. na .eq. 4) goto 80

! no be8 from jina reaclib/bdat
      if (nz .eq. 4 .and. na .eq. 8) goto 80


      do jj=1,inetin
       if (int(zion(i)) .eq. izzz(jj)  .and. &
           int(nion(i)) .ge. inmin(jj) .and. &
           int(nion(i)) .le. inmax(jj)) goto 90
      enddo


! not using this isotope, but
! do the read, backup i by one, and go back to 60 for another isotope

80    continue
!      write(6,117) 'not using ',i,int(zion(i)),int(aion(i)),aname

      do jj=1,10
       if (ic1(jj,i) .gt. 0) then
        read(fil14,240) (dum,j1=1,ic2(jj,i))
240     format (7e13.6)
       end if
      enddo
      i = i - 1
      go to 60


! using this isotope, read coefficients of fitting formula into cx array
90    continue


!      write(6,117) 'using ',i,int(zion(i)),int(aion(i)),int(nion(i)),aname
117   format(1x,a,4i4,' ',a)


      do j=1,8
       ic3(j,i) = k
       if (ic1(j,i) .gt. 0) then
        kmax = k + ic2(j,i) - 1

        if (kmax .gt. cxdim) then
         write(6,*) 'kmax =', kmax,' cxdim =',cxdim
         stop 'kmax > cxdim in routine init_torch'
        end if

        read(fil14,240)  (cx(j1), j1=k,kmax)

!        write(6,*) j,ic1(j,i),ic2(j,i)
!        write(6,*) k,kmax,kmax-k+1
!        write(6,118) (cx(j1), j1=k,kmax)
!118     format(1x,1p7e16.8)

        k = kmax + 1
       end if
      enddo



! and go back for another isotope
      go to 60


! all done with the this part of the loading
120   continue
      close(unit=fil14)
!      write(6,*) 'done reading bdat'


! append neutrons, protons and alfa if down arrow
! set the ending isotope index ionend
      if (downarrow) then
       ionmax       = i
       ineut        = ionmax
       aion(ineut)  = 1.0d0
       nion(ineut)  = 1.0d0
       zion(ineut)  = 0.0d0
       bion(ineut)  = 0.0d0
!       mion(ineut)  = nion(ineut)*mn + zion(ineut)*(mp+me) - bion(ineut)*mev2gr
       mion(ineut)  = nion(ineut)*mn + zion(ineut)*mp - bion(ineut)*mev2gr       
       wion(ineut)  = avo * mion(ineut)
       ionam(ineut) = 'neut'
!       write(6,117) 'using ',ineut,int(zion(ineut)),int(aion(ineut)),int(nion(ineut)),ionam(ineut)

       ionmax       = ionmax + 1
       iprot        = ionmax
       aion(iprot)  = 1.0d0
       nion(iprot)  = 0.0d0
       zion(iprot)  = 1.0d0
       bion(iprot)  = 0.0d0
!       mion(iprot)  = nion(iprot)*mn + zion(iprot)*(mp+me) - bion(iprot)*mev2gr
       mion(iprot)  = nion(iprot)*mn + zion(iprot)*mp - bion(iprot)*mev2gr
       wion(iprot)  = avo * mion(iprot)
       ionam(iprot) = 'prot'
!       write(6,117) 'using ',iprot,int(zion(iprot)),int(aion(iprot)),int(nion(iprot)),ionam(iprot)

       ionmax       = ionmax + 1
       ihe4         = ionmax
       aion(ihe4)   = 4.0d0
       nion(ihe4)   = 2.0d0
       zion(ihe4)   = 2.0d0
       bion(ihe4)   = 28.29603d0
!       mion(ihe4)   = nion(ihe4)*mn + zion(ihe4)*(mp+me) - bion(ihe4)*mev2gr
       mion(ihe4)   = nion(ihe4)*mn + zion(ihe4)*mp - bion(ihe4)*mev2gr
       wion(ihe4)   = avo * mion(ihe4)
       ionam(ihe4)  = 'he4'
!       write(6,117) 'using ',ihe4,int(zion(ihe4)),int(aion(ihe4)),int(nion(ihe4)),ionam(ihe4)

       ionend       = ionmax - 3

! for up-arrow configurations
      else
       ionmax = i - 1
       ionend = ionmax
      end if


! for either orientation, append energy, temperature, and density pointers
       iener = ionmax + 1
       itemp = ionmax + 2
       iden  = ionmax + 3
       ivelx = ionmax + 4
       iposx = ionmax + 5
       neqs  = iposx

       ionam(iener) = 'ener '
       ionam(itemp) = 'temp '
       ionam(iden)  = 'den  '
       ionam(ivelx) = 'velx '
       ionam(iposx) = 'posx '


! number of neutrons
      nion(ionbeg:ionend) = aion(ionbeg:ionend) - zion(ionbeg:ionend)

! mass of each isotope assuming fully ionized
!      mion(ionbeg:ionend) = nion(ionbeg:ionend)*mn + zion(ionbeg:ionend)*(mp+me) - bion(ionbeg:ionend)*mev2gr
      mion(ionbeg:ionend) = nion(ionbeg:ionend)*mn + zion(ionbeg:ionend)*mp - bion(ionbeg:ionend)*mev2gr

! molar mass of each isotope
      wion(ionbeg:ionend) = avo * mion(ionbeg:ionend)

! here is a common approximation for the molar mass
      wion(ionbeg:ionmax) = aion(ionbeg:ionmax)

      ! set niso, eaion, ezion, ebion
      niso = ionmax
      eaion(1:ionmax) = aion(1:ionmax)
      ezion(1:ionmax) = zion(1:ionmax)
      ebion(1:ionmax) = bion(1:ionmax)
   
! build the links between the isotopes in the network

      call naray



! load the ffn, oda and langanke & martin-pinda weak rates

      call weaktab_init(ionmax,zion,aion)



! set the isotope names and pointers
!       write(6,*) ' '
!       write(6,*) ' using isotopes:'
!       write(6,*) '   i   z   a      name'

       if (.not.downarrow) then
        i    = ineut
        inta = aion(i)
        intz = zion(i)
!        write(6,09) ineut,intz,inta,ineut,ionam(i)

        i    = iprot
        inta = aion(i)
        intz = zion(i)
!        write(6,09) iprot,intz,inta,ionam(i)

        i    = ihe4
        inta = aion(i)
        intz = zion(i)
!        write(6,09) ihe4,intz,inta,ionam(i)
       endif


       do i=ionbeg,ionend
        inta = aion(i)
        intz = zion(i)
        if (intz .ge. 1  .and.  intz .le. zmax) then
         if (inta .ge. aidmin(intz)  .and.  inta .le. aidmax(intz)) then
          do ii = aidmin(intz),aidmax(intz)
           if (ii .eq. inta) then

            write(string,01) zsymb(intz),inta
            call sqeeze(string)
            ionam(i) = string

! using this isotope
!            write(6,09) i,intz,inta,ionam(i)

           end if
          enddo
         else
          write(6,*) ' bad aion',inta,' in routine init_torch'
          write(6,*) ' zion=',intz
          write(6,*) ' amin=',aidmin(intz),' amax=',aidmax(intz)
          stop 'error: bad inta in routine init_torch'
         end if
        else
         write(6,*) 'bad zion',intz,' in routine init_torch'
         write(6,*) 'inta =',inta,' zmax=', zmax
         stop 'error: bad intz in routine init_torch'
        end if
       enddo

       if (downarrow) then
        i    = ineut
        inta = aion(i)
        intz = zion(i)
!        write(6,09) ionmax-2,intz,inta,ionam(i)

        i    = iprot
        inta = aion(i)
        intz = zion(i)
!        write(6,09) ionmax-1,intz,inta,ionam(i)

        i    = ihe4
        inta = aion(i)
        intz = zion(i)
!        write(6,09) ionmax,intz,inta,ionam(i)
       endif



! write a isotope file - useful for plotting the network
       if (icoords) then
        open(unit=44,file='isotope_listing.dat',status='unknown')
        do i=1,ionmax
         write(44,'(1x,3i4,a,a)') int(zion(i)), int(aion(i)), int(nion(i)),' ',ionam(i)
        enddo
        close(unit=44)
!        write(6,*) 'wrote file isotope_listing.dat'
!        write(6,*)
       end if

!      do i=1,ionmax
!       write(6,888) ionam(i), &
!                   int(zion(i)),int(nion(i)),int(aion(i)), &
!                   mion(i)*avo,(mion(i)*avo-aion(i)), &
!                   (mion(i)*avo-aion(i))/(mev2gr*avo)
! 888   format(1x,a,3i4,1p5e18.10)
!      enddo
!      read(5,*)



! set the id numbers of certain key isotopes
       do i=ionbeg,ionend
        if (ionam(i) .eq. 'h2   ')   then
         ih2 = i
        else if (ionam(i) .eq. 'h3   ') then
         ih3 = i
        else if (ionam(i) .eq. 'he3  ') then
         ihe3 = i
        else if (ionam(i) .eq. 'li6  ') then
         ili6 = i
        else if (ionam(i) .eq. 'li7  ') then
         ili7 = i
        else if (ionam(i) .eq. 'li8  ') then
         ili8 = i
        else if (ionam(i) .eq. 'be7  ') then
         ibe7 = i
        else if (ionam(i) .eq. 'be9  ') then
         ibe9 = i
        else if (ionam(i) .eq. 'b8   ') then
         ib8 = i
        else if (ionam(i) .eq. 'b9   ') then
         ib9 = i
        else if (ionam(i) .eq. 'b10  ') then
         ib10 = i
        else if (ionam(i) .eq. 'b11  ') then
         ib11 = i
        else if (ionam(i) .eq. 'c11  ') then
         ic11 = i
        else if (ionam(i) .eq. 'c12  ') then
         ic12 = i
        else if (ionam(i) .eq. 'c13  ') then
         ic13 = i
        else if (ionam(i) .eq. 'c14  ') then
         ic14 = i
        else if (ionam(i) .eq. 'n12  ') then
         in12 = i
        else if (ionam(i) .eq. 'n13  ') then
         in13 = i
        else if (ionam(i) .eq. 'n14  ') then
         in14 = i
        else if (ionam(i) .eq. 'n15  ') then
         in15 = i
        else if (ionam(i) .eq. 'o14  ') then
         io14 = i
        else if (ionam(i) .eq. 'o15  ') then
         io15 = i
        else if (ionam(i) .eq. 'o16  ') then
         io16 = i
        else if (ionam(i) .eq. 'o17  ') then
         io17 = i
        else if (ionam(i) .eq. 'o18  ') then
         io18 = i
        else if (ionam(i) .eq. 'f17  ') then
         if17 = i
        else if (ionam(i) .eq. 'f18  ') then
         if18 = i
        else if (ionam(i) .eq. 'f19  ') then
         if19 = i
        else if (ionam(i) .eq. 'ne18 ') then
         ine18 = i
        else if (ionam(i) .eq. 'ne19 ') then
         ine19 = i
        else if (ionam(i) .eq. 'ne20 ') then
         ine20 = i
        else if (ionam(i) .eq. 'ne21 ') then
         ine21 = i
        else if (ionam(i) .eq. 'ne22 ') then
         ine22 = i
        else if (ionam(i) .eq. 'ne23 ') then
         ine23 = i
        else if (ionam(i) .eq. 'na20 ') then
         ina20 = i
        else if (ionam(i) .eq. 'na21 ') then
         ina21 = i
        else if (ionam(i) .eq. 'na22 ') then
         ina22 = i
        else if (ionam(i) .eq. 'na23 ') then
         ina23 = i
        else if (ionam(i) .eq. 'mg22 ') then
         img22 = i
        else if (ionam(i) .eq. 'mg23 ') then
         img23 = i
        else if (ionam(i) .eq. 'mg24 ') then
         img24 = i
        else if (ionam(i) .eq. 'mg25 ') then
         img25 = i
        else if (ionam(i) .eq. 'mg26 ') then
         img26 = i
        else if (ionam(i) .eq. 'al25 ') then
         ial25 = i
        else if (ionam(i) .eq. 'al26 ') then
         ial26 = i
        else if (ionam(i) .eq. 'al27 ') then
         ial27 = i
        else if (ionam(i) .eq. 'si27 ') then
         isi27 = i
        else if (ionam(i) .eq. 'si28 ') then
         isi28 = i
        else if (ionam(i) .eq. 'si29 ') then
         isi29 = i
        else if (ionam(i) .eq. 'si30 ') then
         isi30 = i
        else if (ionam(i) .eq. 'p30  ') then
         ip30  = i
        else if (ionam(i) .eq. 'p31  ') then
         ip31  = i
        else if (ionam(i) .eq. 's30  ') then
         is30  = i
        else if (ionam(i) .eq. 's31  ') then
         is31  = i
        else if (ionam(i) .eq. 's32  ') then
         is32  = i
        else if (ionam(i) .eq. 'cl35 ') then
         icl35  = i
        else if (ionam(i) .eq. 'ar36 ') then
         iar36  = i
        else if (ionam(i) .eq. 'ar37 ') then
         iar37  = i
        else if (ionam(i) .eq. 'ar38 ') then
         iar38  = i
        else if (ionam(i) .eq. 'k39  ') then
         ik39  = i
        else if (ionam(i) .eq. 'ca40 ') then
         ica40  = i
        else if (ionam(i) .eq. 'sc43 ') then
         isc43  = i
        else if (ionam(i) .eq. 'sc45 ') then
         isc45  = i
        else if (ionam(i) .eq. 'ti44 ') then
         iti44  = i
        else if (ionam(i) .eq. 'ti46 ') then
         iti46  = i
        else if (ionam(i) .eq. 'ti48 ') then
         iti48  = i
        else if (ionam(i) .eq. 'ti50 ') then
         iti50  = i
        else if (ionam(i) .eq. 'v46  ') then
         iv46  = i
        else if (ionam(i) .eq. 'v47  ') then
         iv47  = i
        else if (ionam(i) .eq. 'v48  ') then
         iv48  = i
        else if (ionam(i) .eq. 'v51  ') then
         iv51  = i
        else if (ionam(i) .eq. 'cr47 ') then
         icr47  = i
        else if (ionam(i) .eq. 'cr48 ') then
         icr48  = i
        else if (ionam(i) .eq. 'cr49 ') then
         icr49  = i
        else if (ionam(i) .eq. 'cr50 ') then
         icr50  = i
        else if (ionam(i) .eq. 'cr51 ') then
         icr51  = i
        else if (ionam(i) .eq. 'cr52 ') then
         icr52  = i
        else if (ionam(i) .eq. 'cr53 ') then
         icr53  = i
        else if (ionam(i) .eq. 'cr54 ') then
         icr54  = i
        else if (ionam(i) .eq. 'mn50 ') then
         imn50  = i
        else if (ionam(i) .eq. 'mn51 ') then
         imn51  = i
        else if (ionam(i) .eq. 'mn52 ') then
         imn52  = i
        else if (ionam(i) .eq. 'mn55 ') then
         imn55  = i
        else if (ionam(i) .eq. 'fe52 ') then
         ife52 = i
        else if (ionam(i) .eq. 'fe54 ') then
         ife54 = i
        else if (ionam(i) .eq. 'fe55 ') then
         ife55 = i
        else if (ionam(i) .eq. 'fe56 ') then
         ife56 = i
        else if (ionam(i) .eq. 'fe57 ') then
         ife57 = i
        else if (ionam(i) .eq. 'fe58 ') then
         ife58 = i
        else if (ionam(i) .eq. 'fe59 ') then
         ife59 = i
        else if (ionam(i) .eq. 'co54 ') then
         ico54 = i
        else if (ionam(i) .eq. 'co55 ') then
         ico55 = i
        else if (ionam(i) .eq. 'co56 ') then
         ico56 = i
        else if (ionam(i) .eq. 'co59 ') then
         ico59 = i
        else if (ionam(i) .eq. 'ni56 ') then
         ini56 = i
        else if (ionam(i) .eq. 'ni57 ') then
         ini57 = i
        else if (ionam(i) .eq. 'ni58 ') then
         ini58 = i
        else if (ionam(i) .eq. 'ni59 ') then
         ini59 = i
        else if (ionam(i) .eq. 'ni64 ') then
         ini64 = i
        else if (ionam(i) .eq. 'ni66 ') then
         ini66 = i
        else if (ionam(i) .eq. 'cu63 ') then
         icu63 = i
        else if (ionam(i) .eq. 'zn60 ') then
         izn60 = i
        else if (ionam(i) .eq. 'zn64 ') then
         izn64 = i
        end if
       enddo


! now set reaction rate names and pointers

       if (icoords) open(unit=44,file='isotope_linkages.dat',status='unknown')

       nrat = 0
       ratnam(1:abigrat) = ' ' 
       do j=ionbeg,ionend

        if (nrat+14 .gt. abigrat) stop 'nrat > abigrat in init_torch'

! set up the (n,g) and (g,n) names and a few reaction rate pointers
        k = nrr(1,j)
        if (k .gt. 0) then

         if (icoords) call coords(j,k)

         string = 'r_'//ionam(j)//'_ng_'//ionam(k)
!         string = ionam(j)//'(n,g)'//ionam(k)
!         string = 'n '//ionam(j)//' g '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         string = 'r_'//ionam(k)//'_gn_'//ionam(j)
!         string = ionam(k)//'(g,n)'//ionam(j)
!         string = 'g '//ionam(k)//' n '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ih2   .and. k .eq. ih3)  then
          irdng = nrat-1
          irtgn = nrat
         else if (j .eq. ili6  .and. k .eq. ili7) then
          irli6ng = nrat-1
          irli7gn = nrat
         else if (j .eq. ili7  .and. k .eq. ili8) then
          irli7ng = nrat-1
          irli8gn = nrat
         else if (j .eq. ine20  .and. k .eq. ine21) then
          irne20ng = nrat-1
          irne21gn = nrat
         end if
        end if


! set up the (p,n) (n,p) beta- beta+ decay names and a few reaction rate pointers
        k = nrr(2,j)
        if (k .gt. 0) then

         if (icoords) call coords(j,k)

         string = 'r_'//ionam(j)//'_pn_'//ionam(k)
!         string = ionam(j)//'(p,n)'//ionam(k)
!         string = 'p '//ionam(j)//' n '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string


         string = 'r_'//ionam(k)//'_np_'//ionam(j)
!         string = ionam(k)//'(n,p)'//ionam(j)
!         string = 'n '//ionam(k)//' p '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ih3   .and. k .eq. ihe3) then
          irtpn = nrat - 1
          irhe3np = nrat
         else if (j .eq. ili7  .and. k .eq. ibe7) then
          irli7pn = nrat - 1
          irbe7np = nrat
         else if (j .eq. ibe9  .and. k .eq. ib9)  then
          irbe9pn = nrat - 1
          irb9np = nrat
         else if (j .eq. ib11  .and. k .eq. ic11) then
          irb11pn = nrat - 1
          irc11np = nrat
         else if (j .eq. ic13  .and. k .eq. in13) then
          irc13pn = nrat - 1
          irn13np = nrat
         else if (j .eq. ic14  .and. k .eq. in14) then
          irc14pn = nrat - 1
          irn13np = nrat
         else if (j .eq. in14  .and. k .eq. io14) then
          irn14pn = nrat - 1
          iro14np = nrat
         else if (j .eq. in15  .and. k .eq. io15) then
          irn15pn = nrat - 1
          iro15np = nrat
         else if (j .eq. if19  .and. k .eq. ine19) then
          irf19pn = nrat - 1
          irne19np = nrat
         else if (j .eq. ine22  .and. k .eq. ina22) then
          irne22pn = nrat - 1
          irna22np = nrat
         else if (j .eq. ina23  .and. k .eq. img23) then
          irna23pn = nrat - 1
          irmg23np = nrat
         end if

         string = 'r_'//ionam(j)//'_wk-minus_'//ionam(k)
!         string = ionam(j)//'(n=>p)'//ionam(k)
!         string = ionam(j)//' e+ + nub '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ibe7   .and. k .eq. ili7) irbeec = nrat

         string = 'r_'//ionam(k)//'_wk_'//ionam(j)
!         string = ionam(k)//'(p=>n)'//ionam(j)
!         string = ionam(k)//' e- + nu '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

        end if


! set up the (p,g) and (g,p) names and a few reaction rate pointers
        k = nrr(3,j)
        if (k .gt. 0) then

         if (icoords) call coords(j,k)

         string = 'r_'//ionam(j)//'_pg_'//ionam(k)
!         string = ionam(j)//'(p,g)'//ionam(k)
!         string = 'p '//ionam(j)//' g '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         string = 'r_'//ionam(k)//'_gp_'//ionam(j)
!         string = ionam(k)//'(g,p)'//ionam(j)
!         string = 'g '//ionam(k)//' p '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ih2   .and. k .eq. ihe3) then
          irdpg = nrat - 1
          irhe3gp = nrat
         else if (j .eq. ili6  .and. k .eq. ibe7) then
          irli6pg = nrat - 1
          irbe7gp = nrat
         else if (j .eq. ibe7  .and. k .eq. ib8)  then
          irbe7pg = nrat - 1
          irb8gp = nrat
         else if (j .eq. ibe9  .and. k .eq. ib10) then
          irbe9pg = nrat - 1
          irb10gp = nrat
         else if (j .eq. ib10  .and. k .eq. ic11) then
          irb10pg = nrat - 1
          irc11gp = nrat
         else if (j .eq. ib11  .and. k .eq. ic12) then
          irb11pg = nrat - 1
          irc12gp = nrat
         else if (j .eq. ic11  .and. k .eq. in12) then
          irc11pg = nrat - 1
          irn12gp = nrat
         else if (j .eq. ic14  .and. k .eq. in15) then
          irc14pg = nrat - 1
          irn15gp = nrat
         else if (j .eq. in13  .and. k .eq. io14) then
          irn13pg = nrat - 1
          iro14gp = nrat
         else if (j .eq. in14  .and. k .eq. io15) then
          irn14pg = nrat - 1
          iro15gp = nrat
         else if (j .eq. in15  .and. k .eq. io16) then
          irn15pg = nrat - 1
          iro16gp = nrat
         else if (j .eq. io16  .and. k .eq. if17) then
          iro16pg = nrat - 1
          irf17gp = nrat
         else if (j .eq. io17  .and. k .eq. if18) then
          iro17pg = nrat - 1
          irf18gp = nrat
         else if (j .eq. io18  .and. k .eq. if19) then
          iro18pg = nrat - 1
          irf19gp = nrat
         else if (j .eq. if17  .and. k .eq. ine18) then
          irf17pg  = nrat - 1
          irne18gp = nrat
         else if (j .eq. if18  .and. k .eq. ine19) then
          irf18pg  = nrat - 1
          irne19gp = nrat
         else if (j .eq. if19  .and. k .eq. ine20) then
          irf19pg  = nrat - 1
          irne20gp = nrat
         else if (j .eq. ine19  .and. k .eq. ina20) then
          irne19pg = nrat - 1
          irna20gp = nrat
         else if (j .eq. ine20  .and. k .eq. ina21) then
          irne20pg = nrat - 1
          irna21gp = nrat
         else if (j .eq. ine21  .and. k .eq. ina22) then
          irne21pg = nrat - 1
          irna22gp = nrat
         else if (j .eq. ine22  .and. k .eq. ina23) then
          irne22pg = nrat - 1
          irna23gp = nrat
         else if (j .eq. ina21  .and. k .eq. img22) then
          irna21pg = nrat - 1
          irmg22gp = nrat
         else if (j .eq. ina22  .and. k .eq. img23) then
          irna22pg = nrat - 1
          irmg23gp = nrat
         else if (j .eq. ina23  .and. k .eq. img24) then
          irna23pg = nrat - 1
          irmg24gp = nrat
         else if (j .eq. img24  .and. k .eq. ial25) then
          irmg24pg = nrat - 1
          iral25gp = nrat
         else if (j .eq. img25  .and. k .eq. ial26) then
          irmg25pg = nrat - 1
          iral26gp = nrat
         else if (j .eq. img26  .and. k .eq. ial27) then
          irmg26pg = nrat - 1
          iral27gp = nrat
         else if (j .eq. ial25  .and. k .eq. isi26) then
          iral25pg = nrat - 1
          irsi26gp = nrat
         else if (j .eq. ial26  .and. k .eq. isi27) then
          iral26pg = nrat - 1
          irsi27gp = nrat
         else if (j .eq. ial27  .and. k .eq. isi28) then
          iral27pg = nrat - 1
          irsi28gp = nrat
         else if (j .eq. isi27  .and. k .eq. ip28) then
          irsi27pg = nrat - 1
          irp28gp = nrat
         else if (j .eq. isi28  .and. k .eq. ip29) then
          irsi28pg = nrat - 1
          irp29gp = nrat
         else if (j .eq. isi29  .and. k .eq. ip30) then
          irsi29pg = nrat - 1
          irp30gp = nrat
         else if (j .eq. isi30  .and. k .eq. ip31) then
          irsi30pg = nrat - 1
          irp31gp = nrat
         end if
        end if


! set up the (a,p) and (p,a) names and a few reaction rate pointers
        k = nrr(4,j)
        if (k .gt. 0) then

         if (icoords) call coords(j,k)

         string = 'r_'//ionam(j)//'_ap_'//ionam(k)
!         string = ionam(j)//'(a,p)'//ionam(k)
!         string = 'he4 '//ionam(j)//' p '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         string = 'r_'//ionam(k)//'_pa_'//ionam(j)
!         string = ionam(k)//'(p,a)'//ionam(j)
!         string = 'p '//ionam(k)//' he4 '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ihe3  .and. k .eq. ili6) then
          irhe3ap = nrat - 1
          irli6pa = nrat
         else if (j .eq. ili6  .and. k .eq. ibe9) then
          irli6ap = nrat - 1
          irbe9pa = nrat
         else if (j .eq. ibe7  .and. k .eq. ib10) then
          irbe7ap = nrat - 1
          irb10pa = nrat
         else if (j .eq. ib8   .and. k .eq. ic11) then
          irb8ap = nrat - 1
          irc11pa = nrat
         else if (j .eq. ib11  .and. k .eq. ic14) then
          irb11ap = nrat - 1
          irc14pa = nrat
         else if (j .eq. ic11  .and. k .eq. in14) then
          irc11ap = nrat - 1
          irn14pa = nrat
         else if (j .eq. in13  .and. k .eq. io16) then
          irn13ap = nrat - 1
          iro16pa = nrat
         else if (j .eq. in14  .and. k .eq. io17) then
          irn14ap = nrat - 1
          iro17pa = nrat
         else if (j .eq. in15  .and. k .eq. io18) then
          irn15ap = nrat - 1
          iro18pa = nrat
         else if (j .eq. io14  .and. k .eq. if17) then
          iro14ap = nrat - 1
          irf17pa = nrat
         else if (j .eq. io15  .and. k .eq. if18) then
          iro15ap = nrat - 1
          irf18pa = nrat
         else if (j .eq. io16  .and. k .eq. if19) then
          iro16ap = nrat - 1
          irf19pa = nrat
         else if (j .eq. if17  .and. k .eq. ine20) then
          irf17ap  = nrat - 1
          irne20pa = nrat
         else if (j .eq. if19  .and. k .eq. ine22) then
          irf19ap  = nrat - 1
          irne22pa = nrat
         else if (j .eq. ine20  .and. k .eq. ina23) then
          irne20ap  = nrat - 1
          irna23pa  = nrat
         else if (j .eq. ina21  .and. k .eq. img24) then
          irna21ap  = nrat - 1
          irmg24pa  = nrat
         else if (j .eq. img24  .and. k .eq. ial27) then
          irmg24ap  = nrat - 1
          iral27pa  = nrat
         else if (j .eq. img25  .and. k .eq. ial28) then
          irmg25ap  = nrat - 1
          iral28pa  = nrat
         end if
        end if


! set up the (a,n) and (n,a) names and a few reaction rate pointers
        k = nrr(5,j)
        if (k .gt. 0) then

         if (icoords) call coords(j,k)

         string = 'r_'//ionam(j)//'_an_'//ionam(k)
!         string = ionam(j)//'(a,n)'//ionam(k)
!         string = 'he4 '//ionam(j)//' n '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         string = 'r_'//ionam(k)//'_na_'//ionam(j)
!         string = ionam(k)//'(n,a)'//ionam(j)
!         string = 'n '//ionam(k)//' he4 '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ih3   .and. k .eq. ili6) then
          irtan = nrat - 1
          irli6na = nrat
         else if (j .eq. ili7  .and. k .eq. ib10) then
          irli7an = nrat - 1
          irb10na = nrat
         else if (j .eq. ibe9  .and. k .eq. ic12) then
          irbe9an = nrat - 1
          irc12na = nrat
         else if (j .eq. ib10  .and. k .eq. in13) then
          irb10an = nrat - 1
          irn13na = nrat
         else if (j .eq. ib11  .and. k .eq. in14) then
          irb11an = nrat - 1
          irn14na = nrat
         else if (j .eq. ic12  .and. k .eq. io15) then
          irc12an = nrat - 1
          iro15na = nrat
         else if (j .eq. ic13  .and. k .eq. io16) then
          irc13an = nrat - 1
          iro16na = nrat
         else if (j .eq. in14  .and. k .eq. if17) then
          irn14an = nrat - 1
          irf17na = nrat
         else if (j .eq. in15  .and. k .eq. if18) then
          irn15an = nrat - 1
          irf18na = nrat
         else if (j .eq. io17  .and. k .eq. ine20) then
          iro17an = nrat - 1
          irne20na = nrat
         else if (j .eq. io18  .and. k .eq. ine21) then
          iro18an = nrat - 1
          irne21na = nrat
         else if (j .eq. if19  .and. k .eq. ina22) then
          irf19an  = nrat - 1
          irna22na = nrat
         else if (j .eq. ine21  .and. k .eq. img24) then
          irne21an  = nrat - 1
          irmg24na  = nrat
         else if (j .eq. ine22  .and. k .eq. img25) then
          irne22an  = nrat - 1
          irmg25na  = nrat
         else if (j .eq. img25  .and. k .eq. isi28) then
          irmg25an  = nrat - 1
          irsi28na  = nrat
         else if (j .eq. img26  .and. k .eq. isi29) then
          irmg26an  = nrat - 1
          irsi29na  = nrat
         else if (j .eq. ial27  .and. k .eq. ip30) then
          iral27an  = nrat - 1
          irp30na  = nrat
         end if
        end if


! and the (a,g) and (g,a) names and a few reaction rate pointers
        k = nrr(6,j)
        if (k .gt. 0) then

         if (icoords) call coords(j,k)

         string = 'r_'//ionam(j)//'_ag_'//ionam(k)
!         string = ionam(j)//'(a,g)'//ionam(k)
!         string = 'he4 '//ionam(j)//' g '//ionam(k)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         string = 'r_'//ionam(k)//'_ga_'//ionam(j)
!         string = ionam(k)//'(g,a)'//ionam(j)
!         string = 'g '//ionam(k)//' he4 '//ionam(j)

         call sqeeze(string)
         nrat = nrat + 1
         ratnam(nrat) = string

         if (j .eq. ih2   .and. k .eq. ili6) then
          irdag = nrat - 1
          irli6ga = nrat
         else if (j .eq. ih3   .and. k .eq. ili7) then
          irtag = nrat - 1
          irli7ga = nrat
         else if (j .eq. ihe3  .and. k .eq. ibe7) then
          irhe3ag = nrat - 1
          irbe7ga = nrat
         else if (j .eq. ili6  .and. k .eq. ib10) then
          irli6ag = nrat - 1
          irb10ga = nrat
         else if (j .eq. ili7  .and. k .eq. ib11) then
          irli7ag = nrat - 1
          irb11ga = nrat
         else if (j .eq. ibe7  .and. k .eq. ic11) then
          irbe7ag = nrat - 1
          irc11ga = nrat
         else if (j .eq. ic12  .and. k .eq. io16) then
          ircag = nrat - 1
          iroga = nrat
         else if (j .eq. ic14  .and. k .eq. io18) then
          irc14ag = nrat - 1
          iro18ga = nrat
         else if (j .eq. in14  .and. k .eq. if18) then
          irn14ag = nrat - 1
          irf18ga = nrat
         else if (j .eq. in15  .and. k .eq. if19) then
          irn15ag = nrat - 1
          irf19ga = nrat
         else if (j .eq. io14  .and. k .eq. ine18) then
          iro14ag = nrat - 1
          irne18ga = nrat
         else if (j .eq. io15  .and. k .eq. ine19) then
          iro15ag = nrat - 1
          irne19ga = nrat
         else if (j .eq. io16  .and. k .eq. ine20) then
          iroag = nrat - 1
          irnega = nrat
         else if (j .eq. io17  .and. k .eq. ine21) then
          iro17ag = nrat - 1
          irne21ga = nrat
         else if (j .eq. io18  .and. k .eq. ine22) then
          iro18ag = nrat - 1
          irne22ga = nrat
         else if (j .eq. ine20  .and. k .eq. img24) then
          irneag = nrat - 1
          irmgga = nrat
         else if (j .eq. ine21  .and. k .eq. img25) then
          irne21ag = nrat - 1
          irmg25ga = nrat
         else if (j .eq. ine22  .and. k .eq. img26) then
          irne22ag = nrat - 1
          irmg26ga = nrat
         else if (j .eq. img24  .and. k .eq. isi28) then
          irmgag = nrat - 1
          irsiga = nrat
         else if (j .eq. img25  .and. k .eq. isi29) then
          irmg25ag = nrat - 1
          irsi29ga = nrat
         else if (j .eq. img26  .and. k .eq. isi30) then
          irmg26ag = nrat - 1
          irsi30ga = nrat
         else if (j .eq. isi28  .and. k .eq. is32) then
          irsiag = nrat - 1
          irsga  = nrat
         end if
        end if
       enddo


! make sure we have the space before adding in the hardcoded rates

       if (nrat + 67 .gt. abigrat) stop 'nrat > abigrat in init_torch'


! for p(e-,nu)n and n(e+,nub)p reactions
       nrat  = nrat + 1
       irpen = nrat
       ratnam(irpen) = 'r_h1_wk_neut'

       nrat = nrat + 1
       irnep = nrat
       ratnam(irnep) = 'r_neut_wk-minus_h1'

       if (icoords) call coords(ineut,iprot)



! c12 reactions, first half of triple alpha
       if (ic12 .ne. 0) then
        nrat = nrat + 1
        ir3a = nrat
        ratnam(ir3a) = 'r_he4_he4_he4_to_c12'

        nrat = nrat + 1
        irg3a = nrat
        ratnam(irg3a) = 'r_c12_to_he4_he4_he4'

        if (icoords) call coords(ihe4,ic12)


! c12+c12 reactions; must have ne20, na23, mg23 in the network
        if (ine20 .ne. 0  .and. ina23 .ne. 0  .and.  img23 .ne. 0) then

         nrat = nrat + 1
         ir1212n = nrat
         ratnam(ir1212n) = 'r_c12_c12_to_neut_mg23'

         nrat = nrat + 1
         irmg23nc = nrat
         ratnam(irmg23nc) = 'r_neut_mg23_to_c12_c12'

         nrat = nrat + 1
         ir1212p = nrat
         ratnam(ir1212p) = 'r_c12_c12_to_h1_na23'

         nrat = nrat + 1
         irna23pc = nrat
         ratnam(irna23pc) = 'r_h1_na23_to_c12_c12'

         nrat = nrat + 1
         ir1212a = nrat
         ratnam(ir1212a) = 'r_c12_c12_to_he4_ne20'

         nrat = nrat + 1
         irne20ac = nrat
         ratnam(irne20ac) = 'r_he4_ne20_to_c12_c12'

         if (icoords) then
          call coords(ic12,ine20)
          call coords(ic12,ina23)
          call coords(ic12,img23)
         end if

        end if
       end if


! o16+o16 reactions; must have si28, p30, p31 and s31 in the network
       if (io16 .ne. 0 .and. isi28 .ne. 0 .and. ip30 .ne. 0 .and. ip31 .ne. 0 .and. is31 .ne. 0) then

        nrat = nrat + 1
        ir1616n = nrat
        ratnam(nrat) = 'r_o16_o16_to_neut_s31'

        nrat = nrat + 1
        irs31no = nrat
        ratnam(irs31no) = 'r_neut_s31_to_o16_o16'

        nrat = nrat + 1
        ir1616p = nrat
        ratnam(ir1616p) = 'r_o16_o16_to_h1_p31'

        nrat = nrat + 1
        irp31po = nrat
        ratnam(irp31po) = 'r_h1_p31_to_o16_o16'

        nrat = nrat + 1
        ir1616a = nrat
        ratnam(ir1616a) = 'r_o16_o16_to_he4_si28'

        nrat = nrat + 1
        irsi28ao = nrat
        ratnam(irsi28ao) = 'r_he4_si28_to_o16_o16'

        nrat = nrat + 1
        ir1616d = nrat
        ratnam(ir1616d) = 'r_o16_o16_to_h2_p30'

        nrat = nrat + 1
        irp30do = nrat
        ratnam(irp30do) = 'r_h2_p30_to_o16_o16'

         if (icoords) then
          call coords(io16,isi28)
          call coords(io16,ip30)
          call coords(io16,ip31)
          call coords(io16,is31)
         end if
       end if


! c12+o16 reactions; must have mg24, al27, si27 in the network
       if (ic12 .ne. 0 .and. io16 .ne. 0 .and. img24 .ne. 0  .and. ial27 .ne. 0  .and.  isi27 .ne. 0) then

        nrat = nrat + 1
        ir1216n = nrat
        ratnam(ir1216n) = 'r_c12_o16_to_neut_si27'

        nrat = nrat + 1
        irsi27no = nrat
        ratnam(irsi27no) = 'r_neut_si27_to_c12_o16'

        nrat = nrat + 1
        ir1216p = nrat
        ratnam(ir1216p) = 'r_c12_o16_to_h1_al27'

        nrat = nrat + 1
        iral27po = nrat
        ratnam(iral27po) = 'r_h1_al27_to_c12_o16'

        nrat = nrat + 1
        ir1216a = nrat
        ratnam(ir1216a) = 'r_c12_o16_to_he4_mg24'

        nrat = nrat + 1
        irmg24ao = nrat
        ratnam(irmg24ao) = 'r_he4_mg24_to_c12_o16'

        if (icoords) then
         call coords(ic12,io16)
         call coords(ic12,img24)
         call coords(ic12,ial27)
         call coords(ic12,isi27)
         call coords(io16,img24)
         call coords(io16,ial27)
         call coords(io16,isi27)
        end if
       end if


! if we have deuterium
       if (ih2 .ne. 0) then

! pp
        nrat = nrat + 1
        irpp = nrat
        ratnam(irpp) = 'r_h1_h1_wk_h2'

        nrat = nrat + 1
        irpep = nrat
        ratnam(irpep) = 'r_h1_h1_ec_h2'

! p(n,g)d
        nrat = nrat + 1
        irpng = nrat
        ratnam(irpng) = 'r_h1_neut_to_h2'

        nrat = nrat + 1
        irdgn = nrat
        ratnam(irdgn) = 'r_h2_to_h1_neut'

        if (icoords) then
         call coords(iprot,ih2)
         call coords(ineut,ih2)
        end if


! d(p,n)2p
        nrat = nrat + 1
        irdpn = nrat
        ratnam(irdpn) = 'r_h2_h1_to_neut_h1_h1'

        nrat = nrat + 1
        ir2pnp = nrat
        ratnam(ir2pnp) = 'r_h1_h1_neut_to_h1_h2'

! d(d,g)he4
        nrat = nrat + 1
        irddg = nrat
        ratnam(irddg) = 'r_h2_h2_to_he4'

        nrat = nrat + 1
        irhe4gd = nrat
        ratnam(irhe4gd) = 'r_he4_to_h2_h2'

        if (icoords) call coords(ih2,ihe4)

       end if


! if we have deuterium and tritium
       if (ih2 .ne. 0  .and. ih3 .ne. 0) then

! t(d,n)he4 dt reaction
        nrat = nrat + 1
        irtdn = nrat
        ratnam(irtdn) = 'r_h2_h3_to_neut_he4'

        nrat = nrat + 1
        irhe4nd = nrat
        ratnam(irhe4nd) = 'r_neut_he4_to_h2_h3'

        if (icoords) then
         call coords(ineut,ih2)
         call coords(ineut,ih3)
         call coords(ineut,ihe4)
         call coords(ih2,ih3)
         call coords(ih2,ihe4)
         call coords(ih3,ihe4)
        end if

! d(d,p)t
        nrat = nrat + 1
        irddp = nrat
        ratnam(irddp) = 'r_h2_h2_to_h1_h3'

        nrat = nrat + 1
        irtpd = nrat
        ratnam(irtpd) = 'r_h1_h3_to_h2_h2'

        if (icoords) then
         call coords(iprot,ih2)
         call coords(iprot,ih3)
        end if

       end if


! if we have tritium
       if (ih3 .ne. 0) then

! t(p,g)he4
        nrat = nrat + 1
        irh3pg = nrat
        ratnam(irh3pg) = 'r_h1_h3_to_he4'

        nrat = nrat + 1
        irhe4gp = nrat
        ratnam(irhe4gp) = 'r_he4_to_h1_h3'

        if (icoords) then
         call coords(iprot,ih3)
         call coords(iprot,ihe4)
         call coords(ih3,ihe4)
        end if

! t(t,2n)he4
        nrat = nrat + 1
        irtt2n = nrat
        ratnam(irtt2n) = 'r_h3_h3_to_neut_neut_he4'

        nrat = nrat + 1
        irhe42nt = nrat
        ratnam(irhe42nt) = 'r_neut_neut_he4__to_h2_h3'

        if (icoords) then
         call coords(ineut,ih3)
         call coords(ineut,ihe4)
         call coords(ih3,ihe4)
        end if
       end if





! if we have he3
       if (ihe3 .ne. 0) then

! he3(he3,2p)he4
        nrat = nrat + 1
        ir33 = nrat
        ratnam(ir33) = 'r_he3_he3_to_h1_h1_he4'

        nrat = nrat + 1
        ir33inv = nrat
        ratnam(ir33inv) = 'r_h1_h1_he4_to_he3_he3'

        if (icoords) then
         call coords(iprot,ihe3)
         call coords(iprot,ihe4)
         call coords(ihe3,ihe4)
        end if


! he3(p,e+nu)he4
        nrat = nrat + 1
        irhep = nrat
        ratnam(irhep) = 'r_h1_he3_wk_he4'


! he3(n,g)he4
        nrat = nrat + 1
        irhe3ng = nrat
        ratnam(irhe3ng) = 'r_he3_ng_he4'

        nrat = nrat + 1
        irhe4gn = nrat
        ratnam(irhe4gn) = 'r_he4_gn_he3'

        if (icoords) then
         call coords(ineut,ihe3)
         call coords(ineut,ihe4)
         call coords(ineut,ihe4)
        end if
       end if



! if we have he3 and deuterium

      if (ihe3 .ne. 0 .and. ih2 .ne. 0) then

! he3(d,p)he4
        nrat = nrat + 1
        irhe3dp = nrat
        ratnam(irhe3dp) = 'r_h2_he3_to_h1_he4'

        nrat = nrat + 1
        irhe4pd = nrat
        ratnam(irhe4pd) = 'r_h1_he4_to_h2_he3'

        if (icoords) then
         call coords(iprot,ih2)
         call coords(iprot,ihe3)
         call coords(iprot,ihe4)
         call coords(ih2,ihe3)
         call coords(ih2,ihe4)
         call coords(ihe3,ihe4)
        end if

! d(d,n)he3
        nrat = nrat + 1
        irddn = nrat
        ratnam(irddn) = 'r_h2_h2_to_neut_he3'

        nrat = nrat + 1
        irhe3nd = nrat
        ratnam(irhe3nd) = 'r_neut_he3_to_h2_h2'

        if (icoords) then
         call coords(ineut,ih2)
         call coords(ineut,ihe3)
         call coords(ih2,ihe4)
        end if

      end if



! if we have he3, deuterium, tritium

      if (ihe3 .ne. 0  .and. ih2 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,d)he4
        nrat = nrat + 1
        irhe3td = nrat
        ratnam(irhe3td) = 'r_he3_h3_to_d_he4'

        nrat = nrat + 1
        irhe4dt = nrat
        ratnam(irhe4dt) = 'r_h2_he4_to_h3_he3'

        if (icoords) then
         call coords(ih2,ih3)
         call coords(ih2,ihe3)
         call coords(ih2,ihe4)
         call coords(ih3,ihe3)
         call coords(ih3,ihe4)
         call coords(ihe3,ihe4)
        end if

      end if


! if we have he3 and tritium

      if (ihe3 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,np)he4
        nrat = nrat + 1
        irhe3tnp = nrat
        ratnam(irhe3tnp) = 'r_he3_h3_to_neut_h1_he4'

        nrat = nrat + 1
        irhe4pnt = nrat
        ratnam(irhe4pnt) = 'r_neut_h1_he4_to_h3_he3'

        if (icoords) then
         call coords(iprot,ineut)
         call coords(iprot,ih3)
         call coords(iprot,ihe3)
         call coords(iprot,ihe4)
         call coords(ineut,ih3)
         call coords(ineut,ihe3)
         call coords(ineut,ihe4)
         call coords(ih3,ihe3)
         call coords(ih3,ihe4)
         call coords(ihe3,ihe4)
        end if

      end if



! if we have li7

       if (ili7 .ne. 0) then

! li7(p,g)be8 and li7(p,a)he4
        nrat = nrat + 1
        irli7pag = nrat
        ratnam(irli7pag) = 'r_h1_li7_to_he4_he4'

        nrat = nrat + 1
        ir2he4ga = nrat
        ratnam(ir2he4ga) = 'r_he4_he4_to_h1_li7'

        if (icoords) then
         call coords(iprot,ili7)
         call coords(iprot,ihe4)
         call coords(ili7,ihe4)
        end if

       end if

! if we have li7 and deuterium

      if (ili7 .ne. 0 .and. ih2 .ne. 0) then

! li7(d,n)2a
        nrat = nrat + 1
        irli7dn = nrat
        ratnam(irli7dn) = 'r_h2_li7_to_neut_he4_he4'

        nrat = nrat + 1
        ir2and = nrat
        ratnam(ir2and) = 'r_neut_he4_he4_to_h2_li7'

        if (icoords) then
         call coords(ineut,ih2)
         call coords(ineut,ihe4)
         call coords(ineut,ili7)
         call coords(ih2,ihe4)
         call coords(ih2,ili7)
         call coords(ihe4,ili7)
        end if

      end if

! if we have li7 and tritium

      if (ili7 .ne. 0 .and. ih3 .ne. 0) then

! li7(t,2n)2a
        nrat = nrat + 1
        irli7t2n = nrat
        ratnam(irli7t2n) = 'r_h3_li7_to_neut_neut_he4_he4'

        nrat = nrat + 1
        ir2a2nt = nrat
        ratnam(ir2a2nt) = 'r_neut_neut_he4_he4_to_h3_li7'

        if (icoords) then
         call coords(ineut,ih3)
         call coords(ineut,ihe4)
         call coords(ineut,ili7)
         call coords(ih3,ihe4)
         call coords(ih3,ili7)
         call coords(ihe4,ili7)
        end if

      end if


! if we have li7 and he3

      if (ili7 .ne. 0 .and. ihe3 .ne. 0) then

! li7(he3,np)2a
        nrat = nrat + 1
        irli7he3np = nrat
        ratnam(irli7he3np) = 'r_he3_li7_to_neut_h1_he4_he4'

        nrat = nrat + 1
        ir2apnhe3 = nrat
        ratnam(ir2apnhe3) = 'r_neut_h1_he4_he4_to_he3_li7'

        if (icoords) then
         call coords(ineut,iprot)
         call coords(ineut,ihe3)
         call coords(ineut,ihe4)
         call coords(ineut,ili7)
         call coords(iprot,ihe3)
         call coords(iprot,ihe4)
         call coords(iprot,ili7)
         call coords(ihe3,ihe4)
         call coords(ihe3,ili7)
         call coords(ihe4,ili7)
        end if

      end if


! if we have be7 and deutrium
! be7(d,p)2a
       if (ibe7 .ne. 0 .and. ih2 .ne. 0) then
        nrat = nrat + 1
        irbe7dp = nrat
        ratnam(irbe7dp) = 'r_h2_be7_to_h1_he4_he4'

        nrat = nrat + 1
        ir2apd = nrat
        ratnam(ir2apd) = 'r_h1_he4_he4_to_h2_be7'

        if (icoords) then
         call coords(iprot,ih2)
         call coords(iprot,ihe4)
         call coords(iprot,ibe7)
         call coords(ih2,ihe4)
         call coords(ih2,ibe7)
         call coords(ihe4,ibe7)
        end if

       end if


! if we have be7 and tritium
! be7(t,np)2a
       if (ibe7 .ne. 0 .and. ih3 .ne. 0) then
        nrat = nrat + 1
        irbe7tnp = nrat
        ratnam(irbe7tnp) = 'r_h3_be7_to_neut_h1_he4_he4'

        nrat = nrat + 1
        ir2apnt = nrat
        ratnam(ir2apnt) = 'r_neut_h1_he4_he4_to_h3_be7'

        if (icoords) then
         call coords(ineut,iprot)
         call coords(ineut,ih3)
         call coords(ineut,ihe4)
         call coords(ineut,ibe7)
         call coords(iprot,ih3)
         call coords(iprot,ihe4)
         call coords(iprot,ibe7)
         call coords(ih3,ihe4)
         call coords(ih3,ibe7)
         call coords(ihe4,ibe7)
        end if

       end if


! if we have be7 and he3
! be7(he3,2p)2a
       if (ibe7 .ne. 0 .and. ihe3 .ne. 0) then
        nrat = nrat + 1
        irbe7he32p = nrat
        ratnam(irbe7he32p) = 'r_he3_be7_to_h1_h1_he4_he4'

        nrat = nrat + 1
        ir2a2phe3 = nrat
        ratnam(ir2a2phe3) = 'r_h1_h1_he4_he4_to_he3_be7'

        if (icoords) then
         call coords(iprot,ihe3)
         call coords(iprot,ihe4)
         call coords(iprot,ibe7)
         call coords(ihe3,ihe4)
         call coords(ihe3,ibe7)
         call coords(ihe4,ibe7)
        end if

       end if



! if we have be9
       if (ibe9 .ne. 0) then

! a(an,g)be9
        nrat = nrat + 1
        iraan = nrat
        ratnam(iraan) = 'r_neut_he4_he4_to_be9'

        nrat = nrat + 1
        irgaan = nrat
        ratnam(irgaan) = 'r_be9_to_neut_he4_he4'

        if (icoords) then
         call coords(ineut,ihe4)
         call coords(ineut,ibe9)
         call coords(ihe4,ibe9)
        end if
       end if

! if we have be9 and deuterium
       if (ibe9 .ne. 0 .and. ih2 .ne. 0) then

! be9(p,d)be8 =>2a
        nrat = nrat + 1
        irbe9pd = nrat
        ratnam(irbe9pd) = 'r_h1_be9_to_h2_he4_he4'

        nrat = nrat + 1
        ir2adp = nrat
        ratnam(ir2adp) = 'r_h2_he4_he4_to_h1_be9'

        if (icoords) then
         call coords(iprot,ih2)
         call coords(iprot,ihe4)
         call coords(iprot,ibe9)
         call coords(ih2,ihe4)
         call coords(ih2,ibe9)
         call coords(ihe4,ibe9)
        end if
       end if


! if we have b8
       if (ib8 .ne. 0) then

! b8(p=>n)be8 =>2a reactions
        nrat = nrat + 1
        irb8ep = nrat
        ratnam(irb8ep) = 'r_b8_wk_he4_he4'

        if (icoords) then
         call coords(ib8,ihe4)
        end if
       end if


! if we have b11
       if (ib11 .ne. 0) then

! b11(p,a)be8 => 2a reactions
        nrat = nrat + 1
        irb11pa = nrat
        ratnam(irb11pa) = 'r_h1_b11_to_he4_he4_he4'

        nrat = nrat + 1
        ir3ap = nrat
        ratnam(ir3ap) = 'r_he4_he4_he4_to_h1_b11'

        if (icoords) then
         call coords(iprot,ihe4)
         call coords(iprot,ib11)
         call coords(ihe4,ib11)
        end if
       end if


! if we have c11
       if (ic11 .ne. 0) then

! c11(n,a)be8 => 2a
        nrat = nrat + 1
        irc11na = nrat
        ratnam(irc11na) = 'r_c11_neut_to_he4_he4_he4'

        nrat = nrat + 1
        ir3an = nrat
        ratnam(ir3an) = 'r_he4_he4_he4_to_c11_neut'

        if (icoords) then
         call coords(ineut,ihe4)
         call coords(ineut,ic11)
         call coords(ihe4,ic11)
        end if
       end if

       if (icoords) then
        close(unit=44)
!        write(6,*) 'wrote isotope_linkages.dat'
!        write(6,*)
       end if
!      set no of rates for export
       tnrat = nrat
! Extra initialisation (helm etc)       
       init_screen = 1
       init_table  = 1
       init_pivot  = 1
       init_weak   = 1
       call read_helm_table
       return
       end subroutine torch_init
    
      subroutine tfactors(temp)
      include 'implno.dek'
      include 'tfactors.dek'

! sets various popular temperature factors into common block
! this routine must be called before any of the rates are called

! declare the pass
      double precision temp

! all these are in common block

      t9    = temp * 1.0d-9
      t92   = t9*t9
      t93   = t9*t92
      t94   = t9*t93
      t95   = t9*t94
      t96   = t9*t95

      t912  = sqrt(t9)
      t932  = t9*t912
      t952  = t9*t932
      t972  = t9*t952

      t913  = t9**oneth
      t923  = t913*t913
      t943  = t9*t913
      t953  = t9*t923
      t973  = t953*t923
      t9113 = t973*t943

      t914  = t9**(0.25d0)
      t934  = t914*t914*t914
      t954  = t9*t914
      t974  = t9*t934

      t915  = t9**onefif
      t935  = t915*t915*t915
      t945  = t915 * t935
      t965  = t9 * t915

      t916  = t9**onesix
      t976  = t9 * t916
      t9i76 = 1.0d0/t976

      t917  = t9**onesev
      t927  = t917*t917
      t947  = t927*t927

      t918  = sqrt(t914)
      t938  = t918*t918*t918
      t958  = t938*t918*t918

      t9i   = 1.0d0/t9
      t9i2  = t9i*t9i
      t9i3  = t9i2*t9i

      t9i12 = 1.0d0/t912
      t9i32 = t9i*t9i12
      t9i52 = t9i*t9i32
      t9i72 = t9i*t9i52

      t9i13 = 1.0d0/t913
      t9i23 = t9i13*t9i13
      t9i43 = t9i*t9i13
      t9i53 = t9i*t9i23

      t9i14 = 1.0d0/t914
      t9i34 = t9i14*t9i14*t9i14
      t9i54 = t9i*t9i14

      t9i15 = 1.0d0/t915
      t9i35 = t9i15*t9i15*t9i15
      t9i45 = t9i15 * t9i35
      t9i65 = t9i*t9i15

      t9i17 = 1.0d0/t917
      t9i27 = t9i17*t9i17
      t9i47 = t9i27*t9i27

      t9i18 = 1.0d0/t918
      t9i38 = t9i18*t9i18*t9i18
      t9i58 = t9i38*t9i18*t9i18

      return
    end subroutine tfactors

      subroutine zet_read_file
      include 'implno.dek'
      include 'network.dek'

! this routine reads a flat file of isotope names to setup a torch network

! for reading network file
      character*80      :: string,word,atoken
      logical           :: err
      integer           :: i,j,k,l,ip,ntoken,ain,getnam
      character*2       :: junk
      real*8            :: value2


! popular format statements
01    format(1x,a,a,a)
02    format(a)


! get file name
10    continue
!      write(6,*) 'give network file name =>'
!      read(5,2) string
      open(unit=2, file=netname(1:len_trim(netname)), status='old', err=10)
      idnet   = idtorch
!      netname = string
      inetin = 0


! for every line in the file, read it, remove leading blanks, get the non-blank length
      do i=1,10000
       read(2,2,end=40) string

! replace any tab with a space
       j = len_trim(string)
       do k=1,j
        if (iachar(string(k:k)) .eq. 9) string(k:k) = ' '
       enddo

! strip any mesa comments, which are preceeded by a !
       word = adjustl(string)
       j = len_trim(word)
       do k=1,j
        if (iachar(word(k:k)) .eq. 33) then
         do ip=k,j
          word(ip:ip) = ' '
         end do
        end if
       end do

! get another line if this one is now empty 
       j = len_trim(word)
       if (j .eq. 0) cycle

!       write(6,*) trim(word)


! we'll deal with two cases, chemical name and atomic weight, such as o16,
! or a chemical name and a range of atmic weights, such os o 14 19. this 
! means either 1 or 3 tokens on the line.
! this takes me back to metasoftware days.


! get the number of tokens on the line
! get the next substring; j=0 means that we have found the end of word

       ntoken = 0
       ip = 1
       do k=1,3
        l = getnam(word,atoken,ip)
        if (l .eq. 0) cycle
        ntoken = ntoken + 1 
       end do
!       write(6,*) 'number of tokens: ',ntoken



! single chemical symbol
       if (ntoken .eq. 1) then 


! special cases for n, p, alpha, d & t
        if (word(1:j) .eq. 'n' .or. word(1:j) .eq. 'neut' .or. &
            word(1:j) .eq. 'p' .or. word(1:j) .eq. 'h1'   .or. & 
            word(1:j) .eq. 'he4') goto 30

! special mesa cases
        if (word(1:4) .eq. 'add_' .or. word(1:1) .eq. ')') goto 30

! convert common notation for deuterium and tritium
        if (word(1:j) .eq. 'd') then
         word(1:2) = 'h2'
         j = 2
        else if (word(1:j) .eq. 't') then
         word(1:2) = 'h3'
         j = 2
        end if

! find the location of first integer
        do k=1,j
         if (iachar(word(k:k)) .ge. 49 .and. iachar(word(k:k)) .le. 57) then      
          ip = k
          l  = k - 1
          goto 20
         end if
        end do 
        write(6,01) 'did not find isotope number =>',word(1:j)
        goto 10
20      continue

! special cases of single letter chemical names
        if (l .eq. 1) then
         junk(1:2) = word(1:1)//' '
        else
         junk(1:2) = word(1:2)
        end if

! look for the symbol in the chemical names table, fill in z and n
        do k=1,zmax
         if (junk .eq. zsymb(k)) then
          ain = value2(word(ip:j),err)
          if (.not.err) then
           inetin = inetin + 1
           izzz(inetin)   = k
           inmin(inetin)  = ain - k 
           inmax(inetin)  = ain - k
!           write(6,*) izzz(inetin),ain,ain-k,'  ',word(1:j)
           goto 30
          else
           write(6,01) 'error converting number =>',word(ip:j)
           goto 10
          end if
         end if
        enddo

! did not find chemcal symbol name
        write(6,01) 'did not find isotope =>',word(1:j)
        goto 10



! chemical symbol and range of isotopes
       else if (ntoken .eq. 3) then
     
        ip = 1
        l = getnam(word,atoken,ip)
        do k=1,zmax
         if (atoken(1:2) .eq. zsymb(k)) then
          inetin = inetin + 1
          izzz(inetin)   = k

          l = getnam(word,atoken,ip)
          ain = value2(atoken,err)
          if (.not.err) then
            inmin(inetin)  = ain - k 
!            write(6,*) izzz(inetin),ain,ain-k,'  ',word(1:j) 
          else
            write(6,*) 
           write(6,01) 'error converting number =>',trim(atoken)
           goto 10
          end if

          l = getnam(word,atoken,ip)
          ain = value2(atoken,err)
          if (.not.err) then
            inmax(inetin)  = ain - k 
!            write(6,*) izzz(inetin),ain,ain-k,'  ',word(1:j) 
          else
            write(6,*) 
           write(6,01) 'error converting number =>',trim(atoken)
           goto 10
          end if

          goto 30

         end if
        end do

! did not find chemcal symbol name
        write(6,01) 'did not find isotope =>',word(1:j)
        goto 10



! otherwise there is a wierd number of tokens
       else 
        write(6,01) 'did not find isotopes =>',word(1:j)
        goto 10
       end if
        


! back for another line in the file 
30     continue
      enddo

! close file and return
40    close(unit=2)
      return 
    end subroutine zet_read_file

      real*8 function value2(string,err)
      implicit none

! this routine takes a character string and converts it to a real number.
! on error during the conversion, the err flag is set to true
!

! declare the pass
      character*(*)  :: string
      logical        :: err

! local variables
      logical        :: pflag
      integer        :: noblnk,long,ipoint,power,psign,j,z,i
      real*8         ::  x,sign,factor,temp

      character*1, parameter   :: plus = '+'  , minus = '-' , decmal = '.'   , &
                                  blank = ' ' , se = 'e'    , sd = 'd'       , &
                                  se1 = 'E'   , sd1 = 'D'   
      real*8, parameter        :: rten =  10.0d0
      integer, parameter       :: iten = 10 


! initialize
      err    = .false.
      x      =  0.0d0
      sign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)

! remove any leading blanks and get the sign of the number
      do z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto  30
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0d0
        end if
        goto 10
       end if
      enddo

! main number conversion loop
10    do i = noblnk,long
       ipoint = i + 1

! if a blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x     = x * sign
        value2 = x
        return

! if an exponent character, process the whole exponent, and return
       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or. &
                string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0 .and. ipoint.eq.2)     x = 1.0d0
        if (sign .eq. -1.0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value2 = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) goto 30
          power= (power * iten)  + j
         end if
        enddo

! if an ascii number character, process ie
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) goto 30
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         goto 20
        end if

! must be a decimal point if none of the above
! check that there are not two decimal points!
       else
        if (pflag) goto 30
        pflag = .true.
       end if
20    continue
      end do

! if we got through the do loop ok, then we must be done
      x      = x * sign
      value2 = x
      return

! error processing the number
30    err = .true.
      return
      end function value2

      integer function getnam(string,word,ipos)
      implicit none

! this routine is a basic string parser, only spaces and equal signs are
! used as delimiters.

! input is the string to chop and where to start the parse from ipos.

! output is the parsed word and the updated position on the string ipos.
! the length of the word is returned as getnam. if no word is found,
! meaning that there are no more words in the string, then getnam is
! returned as zero and the word is filled with spaces.

! declare the pass
      character*(*) :: string,word
      integer       :: ipos

! local variables
      integer       :: kpos,k,lend,nbegin


! get the length of the input line, blank word and save ipos
      lend = len(string)
      word = ' '
      kpos = ipos


! find begining of the word; if nothing found return getnam as zero
      do k = kpos,lend
       nbegin = k
       if (string(k:k).ne.' ' .and. string(k:k).ne.'=') goto  25
      enddo
      getnam = 0
      return

! find end of the word
 25   continue
      do k = nbegin,lend
       ipos = k
       if( string(k:k).eq.' ' .or. string(k:k).eq.'=') exit
      enddo

! build the word, set getnam and return
      word   = string(nbegin:ipos-1)
      getnam = ipos - nbegin
      return
      end function getnam

      
      subroutine read_helm_table
      include 'implno.dek'
      include 'helm_table_storage.dek'

! this routine reads the helmholtz eos file for electrons/positrons.
! must be called once before the helmeos routine is invoked.

! declare local variables
      logical          file_exists
      integer          i,j
      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                       dd,dd2,ddi,dd2i,dd3i


! for standard table limits
! temperature
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp

! density
       dlo   = -12.0d0
       dhi   = 15.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp


! store the temperature and density arrays
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
       enddo
       do i=1,imax
        dsav = dlo + (i-1)*dstp
        d(i) = 10.0d0**(dsav)
       enddo


! construct the temperature and density deltas and their inverses
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt3i        = dt2i*dti
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
        dt3i_sav(j) = dt3i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd3i        = dd2i*ddi
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
        dd3i_sav(i) = dd3i
       enddo



! now we read the tables

! see if a binary version exists

       inquire(file='helm_table.bin', exist=file_exists)

! if the binary file exists, read it
       if (file_exists) then
!        write(6,*) 'reading binary helm table file'
        open(unit=18,file='helm_table.bin',form='unformatted',status='old')
        read(18) f
        read(18) fd
        read(18) ft
        read(18) fdd
        read(18) ftt
        read(18) fdt
        read(18) fddt
        read(18) fdtt
        read(18) fddtt
        
        read(18) dpdf
        read(18) dpdfd
        read(18) dpdft
        read(18) dpdfdt

        read(18) ef
        read(18) efd
        read(18) eft
        read(18) efdt

        read(18) xf
        read(18) xfd
        read(18) xft
        read(18) xfdt
        close(unit=18)

! binary file does not exist
       else


! see if the ascii file exists
        inquire(file='torch_helm_table.dat', exist=file_exists)
        if (.not.file_exists) then
          stop 'helm table does not exist'
        else

! ascii file exists. first read it
!         write(6,*) 'reading ascii helm table file'
         open(unit=19,file='torch_helm_table.dat',status='old')

! read the helmholtz free energy and its derivatives
         do j=1,jmax
          do i=1,imax
           read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                    fddt(i,j),fdtt(i,j),fddtt(i,j)
          enddo
         enddo

! read the pressure derivative with density table
         do j=1,jmax
          do i=1,imax
           read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
          enddo
         enddo

! read the electron chemical potential table
         do j=1,jmax
          do i=1,imax
           read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
          enddo
         enddo

! read the number density table
         do j=1,jmax
          do i=1,imax
           read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
          enddo
         enddo

! close the ascii file
         close(unit=19)

! create a binary file for future use
!         write(6,*) 'writing binary file'
         open(unit=18,file='helm_table.bin',form='unformatted')
         write(18) f
         write(18) fd
         write(18) ft
         write(18) fdd
         write(18) ftt
         write(18) fdt
         write(18) fddt
         write(18) fdtt
         write(18) fddtt
        
         write(18) dpdf
         write(18) dpdfd
         write(18) dpdft
         write(18) dpdfdt

         write(18) ef
         write(18) efd
         write(18) eft
         write(18) efdt

         write(18) xf
         write(18) xfd
         write(18) xft
         write(18) xfdt
         close(unit=18)
        end if
       end if


!      write(6,*)
!      write(6,*) 'finished reading eos table'
!      write(6,04) 'imax=',imax,' jmax=',jmax
!04    format(1x,4(a,i4))
!      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
!      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
!03    format(1x,4(a,1pe11.3))
!      write(6,*)

      return
      end subroutine read_helm_table



      subroutine helmeos
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'helm_table_storage.dek'


! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar and zbar, this routine returns most of the other
! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
! their derivatives with respect to temperature, density, abar, and zbar.
! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.
!
! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. interpolation in a table of the helmholtz free energy
! is used to return the electron-positron thermodynamic quantities.
! all other derivatives are analytic.
!
! references: cox & giuli chapter 24 ; timmes & swesty apj 1999


! declare
      integer          i,j
      double precision temp,den,abar,zbar,ytot1,ye, &
                       x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                       dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                       dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                       deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                       dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                       sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                       dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                       gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                       detadt,detadd,xnefer,dxnedt,dxnedd,s

      double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
                       egas,degasdd,degasdt,degasda,degasdz, &
                       sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
                       cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
                       chit_gas,chid_gas,nabad_gas,sound_gas


      double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                        forth   = 4.0d0/3.0d0, &
                        forpi   = 4.0d0 * pi, &
                        kergavo = kerg * avo, &
                        ikavo   = 1.0d0/kergavo, &
                        asoli3  = asol/3.0d0, &
                        light2  = clight * clight)

! for the abar derivatives
      double precision dpradda,deradda,dsradda, &
                       dpionda,deionda,dsionda, &
                       dpepda,deepda,dsepda, &
                       dpresda,denerda,dentrda, &
                       detada,dxneda

! for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz, &
                       dpiondz,deiondz,dsiondz, &
                       dpepdz,deepdz,dsepdz, &
                       dpresdz,denerdz,dentrdz, &
                       detadz,dxnedz

! for the interpolations
      integer          iat,jat
      double precision free,df_d,df_t,df_dd,df_tt,df_dt
      double precision xt,xd,mxt,mxd, &
                       si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                       si0d,si1d,si2d,si0md,si1md,si2md, &
                       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                       ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                       z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                       dpsi2,ddpsi2,din,h5,fi(36), &
                       xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
                       w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                       w0d,w1d,w2d,w0md,w1md,w2md


! for the uniform background coulomb correction
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
                       plasg,plasgdd,plasgdt,plasgda,plasgdz, &
                       ecoul,decouldd,decouldt,decoulda,decouldz, &
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
                       a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
      parameter        (a1    = -0.898004d0, &
                        b1    =  0.96786d0, &
                        c1    =  0.220703d0, &
                        d1    = -0.86097d0, &
                        e1    =  2.5269d0, &
                        a2    =  0.29561d0, &
                        b2    =  1.9885d0, &
                        c2    =  0.288675d0, &
                        third =  1.0d0/3.0d0, &
                        esqu  =  qe * qe)


! quintic hermite polynomial statement functions
! psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


! psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


! psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


! biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
             fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



! cubic hermite polynomial statement functions
! psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)


! psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


! bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
             fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



! popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))



! start of pipeline loop, normal execution starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

!       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
!       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = max(1.0d-16,ytot1 * zbar)



! initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


! radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0

       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0


! ion section:
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnida  = -xni * ytot1

        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = dxnida * kt
        dpiondz = 0.0d0

        eion    = 1.5d0 * pion*deni
        deiondd = (1.5d0 * dpiondd - eion)*deni
        deiondt = 1.5d0 * dpiondt*deni
        deionda = 1.5d0 * dpionda*deni
        deiondz = 0.0d0


! sackur-tetrode equation for the ion entropy of
! a single ideal gas characterized by abar
        x       = abar*abar*sqrt(abar) * deni/avo
        s       = sioncon * temp
        z       = x * s * sqrt(s)
        y       = log(z)

!        y       = 1.0d0/(abar*kt)
!        yy      = y * sqrt(y)
!        z       = xni * sifac * yy
!        etaion  = log(z)


        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                   - kergavo * deni * ytot1
        dsiondt = (dpiondt*deni + deiondt)*tempi - &
                  (pion*deni + eion) * tempi*tempi &
                  + 1.5d0 * kergavo * tempi*ytot1
        x       = avo*kerg/abar
        dsionda = (dpionda*deni + deionda)*tempi &
                  + kergavo*ytot1*ytot1* (2.5d0 - y)
        dsiondz = 0.0d0



! electron-positron section:


! assume complete ionization
        xnem    = xni * zbar


! enter the table with ye*den
        din = ye*den


! bomb proof the input
        if (temp .gt. t(jmax)) then
         write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
         write(6,*) 'temp too hot, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (temp .lt. t(1)) then
         write(6,01) 'temp=',temp,' t(1)=',t(1)
         write(6,*) 'temp too cold, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .gt. d(imax)) then
         write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
         write(6,*) 'ye*den too big, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .lt. d(1)) then
         write(6,01) 'ye*den=',din,' d(1)=',d(1)
         write(6,*) 'ye*den too small, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if

! hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,imax-1))

! various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

! the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt)*dt_sav(jat)
        si2t =   psi2(xt)*dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt)*dt_sav(jat)
        si2mt =  psi2(mxt)*dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd)*dd_sav(iat)
        si2d =   psi2(xd)*dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd)*dd_sav(iat)
        si2md =  psi2(mxd)*dd2_sav(iat)

! derivatives of the weight functions
        dsi0t =   dpsi0(xt)*dti_sav(jat)
        dsi1t =   dpsi1(xt)
        dsi2t =   dpsi2(xt)*dt_sav(jat)

        dsi0mt = -dpsi0(mxt)*dti_sav(jat)
        dsi1mt =  dpsi1(mxt)
        dsi2mt = -dpsi2(mxt)*dt_sav(jat)

        dsi0d =   dpsi0(xd)*ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd)*dd_sav(iat)

        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd)*dd_sav(iat)

! second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)

        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

! now get the pressure derivative with density, chemical potential, and
! electron positron number densities
! get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)


! derivatives of weight functions
        dsi0t  = xdpsi0(xt)*dti_sav(jat)
        dsi1t  = xdpsi1(xt)

        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        dsi0d  = xdpsi0(xd)*ddi_sav(iat)
        dsi1d  = xdpsi1(xd)

        dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
        dsi1md = xdpsi1(mxd)


! look in the electron chemical potential table only once
        fi(1)  = ef(iat,jat)
        fi(2)  = ef(iat+1,jat)
        fi(3)  = ef(iat,jat+1)
        fi(4)  = ef(iat+1,jat+1)
        fi(5)  = eft(iat,jat)
        fi(6)  = eft(iat+1,jat)
        fi(7)  = eft(iat,jat+1)
        fi(8)  = eft(iat+1,jat+1)
        fi(9)  = efd(iat,jat)
        fi(10) = efd(iat+1,jat)
        fi(11) = efd(iat,jat+1)
        fi(12) = efd(iat+1,jat+1)
        fi(13) = efdt(iat,jat)
        fi(14) = efdt(iat+1,jat)
        fi(15) = efdt(iat,jat+1)
        fi(16) = efdt(iat+1,jat+1)


! electron chemical potential etaele
        etaele  = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)


! derivative with respect to density
        x       = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detadd  = ye * x

! derivative with respect to temperature
        detadt  = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd

! end of pipeline loop
      enddo
      return
      end subroutine helmeos
    
      subroutine naray
      include 'implno.dek'
      include 'network.dek'

! this routine builds the nrr(7,i) and nrrneut(7,i) arrays, which specify
! the location of isotopes coupled to i by various reactions. that is,
! it build what becomes the band diagonal part of the jacobian for binary reactions.

! the first index on nrr refers to reactions of the form
!      1=ng   2=pn   3=pg   4=ap   5=an   6=ag   7=b-

! second index on nrr refers to the entrance channel isotope
! the value of nrr hold the exit chanel isotope

! for example i(a,g)k is encoded as nrr(6,i) = k

!
! nrreut if for neutrino nucleosynthesis, which is disabled for now
! the first index on nrrneut refers to reactions of the form
!      1=nu,e-,n  2=nu,e-  3=nu,e-,p   4=nu e+,n   5=nu,e+   6=nu,e+,p


! declare
      integer  ::  i,k,n,kz,kn,jz(7),jn(7)


! initialize for nrr
      jz(1) = 0
      jz(2) = 1
      jz(3) = 1
      jz(4) = 1
      jz(5) = 2
      jz(6) = 2
      jz(7) = -1

      jn(1) = 1
      jn(2) = -1
      jn(3) = 0
      jn(4) = 2
      jn(5) = 1
      jn(6) = 2
      jn(7) = 1

! build nrr
      do i=ionbeg,ionend
       do n=1,7
        nrr(n,i) = 0
        kz = int(zion(i)) + jz(n)
        kn = int(nion(i)) + jn(n)
        do k=ionbeg,ionend
         if (kz.eq.int(zion(k)) .and. kn.eq.int(nion(k))) then
          nrr(n,i)=k

! these are the n-z plane coordinates connecting isotop i with isotope k
! does not pick up the specials though
!          write(6,'(1x,0p4f8.1)') nion(i), zion(i), nion(k), zion(k)
!          write(6,'(1x,4i4)') int(nion(i)), int(zion(i)), int(nion(k)), int(zion(k))
         end if

        enddo
       enddo
      enddo


! initialize for nrrneut for neutrino reactions
!      jz(1) =  1
!      jz(2) =  1
!      jz(3) =  0
!      jz(4) = -1
!      jz(5) = -1
!      jz(6) = -2
!      jz(7) = -2
!      jn(1) = -2
!      jn(2) = -1
!      jn(3) = -1
!      jn(4) = 0
!      jn(5) = 1
!      jn(6) = 1
!      jn(7) = -2

! build nrrneut
!      do i=ionbeg,ionend
!       do n=1,7
!        nrrneut(n,i) = 0
!        kz = int(zion(i)) + jz(n)
!        kn = int(nion(i)) + jn(n)
!        do k=ionbeg,ionend
!         if (kz.eq.int(zion(k)) .and. kn.eq.int(nion(k))) nrrneut(n,i)=k
!        enddo
!       enddo
!      enddo

      return
      end subroutine naray

      subroutine sqeeze(line)
      implicit none

! this routine takes line and removes all blanks, such as
! those from writing to string with fortran format statements
!
! declare the pass
      character*(*) :: line

! local variables
      character*1  :: achar
      integer      :: l,n,k,lend,lsiz


! find the end of the line
      lsiz = len(line)
      lend = len_trim(line)
      n    = 0
      l    = 0

! do the compression in place
10    continue
      l = l + 1
      achar = line(l:l)
      if (achar .eq. ' ') goto 10
      n = n + 1
      line(n:n) = achar
      if (l .lt. lend) goto 10

! blank the rest of the line
      do k=n+1,lsiz
       line(k:k) = ' '
      enddo
      return
      end subroutine sqeeze

      subroutine torchrat(ye)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'vector_eos.dek'
      include 'network.dek'
      include 'tfactors.dek'


! this routine generates nuclear reaction rates for the torch network
!
! declare the pass
      real*8 ::  ye


! local variables
      integer ::i,j,k,jd1,jd2,jd3,jd4,jd5,jn,j1,j2,k0,k1,k2,k3, &
                k4,k5,k6,k7,k8,k9,jxx,inap, &
                nrate,lk0,jk,lk1,lk2,lk3,lk4,lk5,lk6,lk7,lk8

      real*8 :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                frate,dfratedt,dfratedd,rrate,drratedt,drratedd,z,g0,dg0,sp1,sn1

! initialize
      inap = 0
      sigraw(1:14,1:ionmax)   = 0.0d0
      sigrawdt(1:14,1:ionmax) = 0.0d0
      sigrawdd(1:14,1:ionmax) = 0.0d0

      if (btemp .lt. 1.0e6) return


! get the temperature factors
      call tfactors(btemp)


! get the partition functions, arrays wpart, zwork1, and zwork2 filled
! these are used in the evaluate_rate calls below

      if (which_bdat .eq. jina_bdat) then      
       call get_partition_functions(t9)

      else if (which_bdat .eq. pre_jina_bdat) then      
       call bdat_partition_functions

      else
       stop 'unknown bdat option in torchrat'
      end if



! generate strong, electromagnetic and ground state beta decay rates
! j=1=(ng) j=2=(pn) j=3=b-  j=4=(pg) j=5=(ap)  j=6=(an)  j=7=(ag) j=8=a decay


      do i=ionbeg,ionend

! set up the y(j)(n,g)y(k) and y(k)(g,n)y(j) components

       if (nrr(1,i) .gt. 0) then
        call evaluate_rate(1,i,sigraw(1,i),sigrawdt(1,i),sigrawdd(1,i), &
                               sigraw(2,i),sigrawdt(2,i),sigrawdd(2,i))
        inap = inap + 2

        if (shut_down_ng .eq. 1) then
          sigraw(1,i)=0.0d0 ; sigrawdt(1,i)=0.0d0 ; sigrawdd(1,i)=0.0d0
          sigraw(2,i)=0.0d0 ; sigrawdt(2,i)=0.0d0 ; sigrawdd(2,i)=0.0d0
        end if
       end if


! set up the (p,n) (n,p) and beta- beta+ decay components

       if (nrr(2,i) .gt. 0) then
        call evaluate_rate(3,i,sigraw(3,i),sigrawdt(3,i),sigrawdd(3,i), &
                               sigraw(4,i),sigrawdt(4,i),sigrawdd(4,i))
        if (weak_on .eq. 1) then        
         call evaluate_rate(5,i,sigraw(5,i),sigrawdt(5,i),sigrawdd(5,i), &
                               sigraw(6,i),sigrawdt(6,i),sigrawdd(6,i))
        endif
        inap = inap + 4

        if (shut_down_pn .eq. 1) then
          sigraw(3,i)=0.0d0 ; sigrawdt(3,i)=0.0d0 ; sigrawdd(3,i)=0.0d0
          sigraw(4,i)=0.0d0 ; sigrawdt(4,i)=0.0d0 ; sigrawdd(4,i)=0.0d0
        end if

! zero the weak rates
         if (weak_on .eq. 0) then
          sigraw(5,i)   = 0.0d0
          sigrawdt(5,i) = 0.0d0
          sigrawdd(5,i) = 0.0d0
          sigraw(6,i)   = 0.0d0
          sigrawdt(6,i) = 0.0d0
          sigrawdd(6,i) = 0.0d0
         end if
       end if



! set up the (p,g) and (g,p) components

       if (nrr(3,i) .gt. 0) then
        call evaluate_rate(7,i,sigraw(7,i),sigrawdt(7,i),sigrawdd(7,i), &
                               sigraw(8,i),sigrawdt(8,i),sigrawdd(8,i))
        inap = inap + 2

        if (shut_down_pg .eq. 1) then
          sigraw(7,i)=0.0d0 ; sigrawdt(7,i)=0.0d0 ; sigrawdd(7,i)=0.0d0
          sigraw(8,i)=0.0d0 ; sigrawdt(8,i)=0.0d0 ; sigrawdd(8,i)=0.0d0
        end if
       end if


! set up the (a,p) and (p,a) reactions

       if (nrr(4,i) .gt. 0) then
        call evaluate_rate(9,i,sigraw(9,i),sigrawdt(9,i),sigrawdd(9,i), &
                               sigraw(10,i),sigrawdt(10,i),sigrawdd(10,i))
        inap = inap + 2

        if (shut_down_ap .eq. 1) then
          sigraw(9,i)=0.0d0  ; sigrawdt(9,i)=0.0d0  ; sigrawdd(9,i)=0.0d0
          sigraw(10,i)=0.0d0 ; sigrawdt(10,i)=0.0d0 ; sigrawdd(10,i)=0.0d0
        end if
       end if


! set up the (a,n) and (n,a) components

       if (nrr(5,i) .gt. 0) then
        call evaluate_rate(11,i,sigraw(11,i),sigrawdt(11,i),sigrawdd(11,i), &
                                sigraw(12,i),sigrawdt(12,i),sigrawdd(12,i))
        inap = inap + 2

        if (shut_down_an .eq. 1) then
          sigraw(11,i)=0.0d0 ; sigrawdt(11,i)=0.0d0 ; sigrawdd(11,i)=0.0d0
          sigraw(12,i)=0.0d0 ; sigrawdt(12,i)=0.0d0 ; sigrawdd(12,i)=0.0d0
        end if
       end if


! and the (a,g) and (g,a) components

       if (nrr(6,i) .gt. 0) then
        call evaluate_rate(13,i,sigraw(13,i),sigrawdt(13,i),sigrawdd(13,i), &
                        sigraw(14,i),sigrawdt(14,i),sigrawdd(14,i))
        inap = inap + 2

        if (shut_down_ag .eq. 1) then
          sigraw(13,i)=0.0d0 ; sigrawdt(13,i)=0.0d0 ; sigrawdd(13,i)=0.0d0
          sigraw(14,i)=0.0d0 ; sigrawdt(14,i)=0.0d0 ; sigrawdd(14,i)=0.0d0
        end if
       end if

      enddo





! now add some specific reaction rates
! versions of BDAT before jina's requires setting the following specific reactions. 


! for p(e-,nu)n and n(e+,nub)p reactions
! count them, but don't compute them here
      inap = inap + 2


! triple alpha to c12
      if (ic12 .ne. 0) then
       call rate_tripalf(btemp,bden, &
                         ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
                         ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))
       inap = inap + 2


!       ratraw(ir3a)     = 0.0d0
!       dratrawdt(ir3a)  = 0.0d0
!       dratrawdd(ir3a)  = 0.0d0
!       ratraw(irg3a)    = 0.0d0
!       dratrawdt(irg3a) = 0.0d0
!       dratrawdd(irg3a) = 0.0d0




! c12 + c12
       if (ine20.ne.0 .and. ina23.ne.0 .and. img23.ne.0) then

        call rate_c12c12npa(btemp,bden, &
                            ratraw(ir1212n),dratrawdt(ir1212n),dratrawdd(ir1212n),    &
                            ratraw(irmg23nc),dratrawdt(irmg23nc),dratrawdd(irmg23nc), &
                            ratraw(ir1212p),dratrawdt(ir1212p),dratrawdd(ir1212p),    &
                            ratraw(irna23pc),dratrawdt(irna23pc),dratrawdd(irna23pc), &
                            ratraw(ir1212a),dratrawdt(ir1212a),dratrawdd(ir1212a),    &
                            ratraw(irne20ac),dratrawdt(irne20ac),dratrawdd(irne20ac))

        inap = inap + 6


!       ratraw(ir1212n)     = 0.0d0
!       dratrawdt(ir1212n)  = 0.0d0
!       dratrawdd(ir1212n)  = 0.0d0
!       ratraw(irmg23nc)    = 0.0d0
!       dratrawdt(irmg23nc) = 0.0d0
!       dratrawdd(irmg23nc) = 0.0d0
!       ratraw(ir1212p)     = 0.0d0
!       dratrawdt(ir1212p)  = 0.0d0
!       dratrawdd(ir1212p)  = 0.0d0
!       ratraw(irna23pc)    = 0.0d0
!       dratrawdt(irna23pc) = 0.0d0
!       dratrawdd(irna23pc) = 0.0d0
!       ratraw(ir1212a)     = 0.0d0
!       dratrawdt(ir1212a)  = 0.0d0
!       dratrawdd(ir1212a)  = 0.0d0
!       ratraw(irne20ac)    = 0.0d0
!       dratrawdt(irne20ac) = 0.0d0
!       dratrawdd(irne20ac) = 0.0d0

       end if
      end if




! o16 + o16
      if (io16 .ne. 0 .and. isi28 .ne. 0 .and. ip30.ne.0 .and. &
          ip31 .ne. 0 .and. is31 .ne. 0) then

       call rate_o16o16npad(btemp,bden, &
                            ratraw(ir1616n),dratrawdt(ir1616n),dratrawdd(ir1616n), &
                            ratraw(irs31no),dratrawdt(irs31no),dratrawdd(irs31no), &
                            ratraw(ir1616p),dratrawdt(ir1616p),dratrawdd(ir1616p), &
                            ratraw(irp31po),dratrawdt(irp31po),dratrawdd(irp31po), &
                            ratraw(ir1616a),dratrawdt(ir1616a),dratrawdd(ir1616a), &
                            ratraw(irsi28ao),dratrawdt(irsi28ao),dratrawdd(irsi28ao), &
                            ratraw(ir1616d),dratrawdt(ir1616d),dratrawdd(ir1616d), &
                            ratraw(irp30do),dratrawdt(irp30do),dratrawdd(irp30do))

       inap = inap + 8

!       ratraw(ir1616n)     = 0.0d0
!       dratrawdt(ir1616n)  = 0.0d0
!       dratrawdd(ir1616n)  = 0.0d0
!       ratraw(irs31no)     = 0.0d0
!       dratrawdt(irs31no)  = 0.0d0
!       dratrawdd(irs31no)  = 0.0d0
!       ratraw(ir1616p)     = 0.0d0
!       dratrawdt(ir1616p)  = 0.0d0
!       dratrawdd(ir1616p)  = 0.0d0
!       ratraw(irp31po)     = 0.0d0
!       dratrawdt(irp31po)  = 0.0d0
!       dratrawdd(irp31po)  = 0.0d0
!       ratraw(ir1616a)     = 0.0d0
!       dratrawdt(ir1616a)  = 0.0d0
!       dratrawdd(ir1616a)  = 0.0d0
!       ratraw(irsi28ao)    = 0.0d0
!       dratrawdt(irsi28ao) = 0.0d0
!       dratrawdd(irsi28ao) = 0.0d0
!       ratraw(ir1616d)     = 0.0d0
!       dratrawdt(ir1616d)  = 0.0d0
!       dratrawdd(ir1616d)  = 0.0d0
!       ratraw(irp30do)     = 0.0d0
!       dratrawdt(irp30do)  = 0.0d0
!       dratrawdd(irp30do)  = 0.0d0

      end if




! c12 + o16 reaction
      if (ic12 .ne. 0 .and. io16 .ne. 0 .and. &
          img24 .ne. 0  .and. ial27 .ne. 0  .and.  isi27 .ne. 0) then

       call rate_c12o16npa(btemp,bden, &
                           ratraw(ir1216n),dratrawdt(ir1216n),dratrawdd(ir1216n), &
                           ratraw(irsi27no),dratrawdt(irsi27no),dratrawdd(irsi27no), &
                           ratraw(ir1216p),dratrawdt(ir1216p),dratrawdd(ir1216p), &
                           ratraw(iral27po),dratrawdt(iral27po),dratrawdd(iral27po), &
                           ratraw(ir1216a),dratrawdt(ir1216a),dratrawdd(ir1216a), &
                           ratraw(irmg24ao),dratrawdt(irmg24ao),dratrawdd(irmg24ao))

       inap = inap + 6

!       ratraw(ir1216n)     = 0.0d0
!       dratrawdt(ir1216n)  = 0.0d0
!       dratrawdd(ir1216n)  = 0.0d0
!       ratraw(irsi27no)    = 0.0d0
!       dratrawdt(irsi27no) = 0.0d0
!       dratrawdd(irsi27no) = 0.0d0
!       ratraw(ir1216p)     = 0.0d0
!       dratrawdt(ir1216p)  = 0.0d0
!       dratrawdd(ir1216p)  = 0.0d0
!       ratraw(iral27po)    = 0.0d0
!       dratrawdt(iral27po) = 0.0d0
!       dratrawdd(iral27po) = 0.0d0
!       ratraw(ir1216a)     = 0.0d0
!       dratrawdt(ir1216a)  = 0.0d0
!       dratrawdd(ir1216a)  = 0.0d0
!       ratraw(irmg24ao)    = 0.0d0
!       dratrawdt(irmg24ao) = 0.0d0
!       dratrawdd(irmg24ao) = 0.0d0

      end if



! if have deuterium
      if (ih2 .ne. 0) then

! p(p,e+nu)d
      call rate_pp(btemp,bden, &
                   ratraw(irpp),dratrawdt(irpp),dratrawdd(irpp), &
                   rrate,drratedt,drratedd)

       inap = inap + 1

!       if (weak_on .eq. 0) then
!        ratraw(irpp)    = 0.0d0
!        dratrawdt(irpp) = 0.0d0
!        dratrawdd(irpp) = 0.0d0
!       end if


! p(e-p,nu)d
      call rate_pep(btemp,bden,ye, &
                   ratraw(irpep),dratrawdt(irpep),dratrawdd(irpep), &
                   rrate,drratedt,drratedd)
       inap = inap + 1

!       if (weak_on .eq. 0) then
!        ratraw(irpep)    = 0.0d0
!        dratrawdt(irpep) = 0.0d0
!        dratrawdd(irpep) = 0.0d0
!       end if


! p(n,g)d
      call rate_png(btemp,bden, &
                    ratraw(irpng),dratrawdt(irpng),dratrawdd(irpng), &
                    ratraw(irdgn),dratrawdt(irdgn),dratrawdd(irdgn))
       inap = inap + 2


! d(p,n)2p
      call rate_dpn(btemp,bden, &
                    ratraw(irdpn),dratrawdt(irdpn),dratrawdd(irdpn), &
                    ratraw(ir2pnp),dratrawdt(ir2pnp),dratrawdd(ir2pnp))
       inap = inap + 2


! d(d,g)he4
      call rate_ddg(btemp,bden, &
                    ratraw(irddg),dratrawdt(irddg),dratrawdd(irddg), &
                    ratraw(irhe4gd),dratrawdt(irhe4gd),dratrawdd(irhe4gd))
       inap = inap + 2
      end if



! if we have deuterium and tritium

      if (ih2 .ne. 0 .and. ih3 .ne. 0) then

! t(d,n)he4 dt reaction
      call rate_tdn(btemp,bden, &
                    ratraw(irtdn),dratrawdt(irtdn),dratrawdd(irtdn), &
                    ratraw(irhe4nd),dratrawdt(irhe4nd),dratrawdd(irhe4nd))
       inap = inap + 2

! d(d,p)t
      call rate_ddp(btemp,bden, &
                    ratraw(irddp),dratrawdt(irddp),dratrawdd(irddp), &
                    ratraw(irtpd),dratrawdt(irtpd),dratrawdd(irtpd))
       inap = inap + 2


! d(n,g)t
      if (old_bdat_format)      & 
      call rate_dng(btemp,bden, &
                    sigraw(1,ih2),sigrawdt(1,ih2),sigrawdd(1,ih2), &
                    sigraw(2,ih2),sigrawdt(2,ih2),sigrawdd(2,ih2))

      end if



! if we have tritium
      if (ih3 .ne. 0) then

! t(p,g)he4  
      call rate_tpg(btemp,bden, &
                    ratraw(irh3pg),dratrawdt(irh3pg),dratrawdd(irh3pg), &
                    ratraw(irhe4gp),dratrawdt(irhe4gp),dratrawdd(irhe4gp))
       inap = inap + 2


! t(t,2n)he4
      call rate_tt2n(btemp,bden, &
                     ratraw(irtt2n),dratrawdt(irtt2n),dratrawdd(irtt2n), &
                     ratraw(irhe42nt),dratrawdt(irhe42nt),dratrawdd(irhe42nt))
       inap = inap + 2



! t(p,n)he3
      if (old_bdat_format)      & 
      call rate_tpn(btemp,bden, &
                    sigraw(3,ih3),sigrawdt(3,ih3),sigrawdd(3,ih3), &
                    sigraw(4,ih3),sigrawdt(4,ih3),sigrawdd(4,ih3))
      end if



! if we have he3
      if (ihe3 .ne. 0) then

! he3(he3,2p)he4
       call rate_he3he3(btemp,bden, &
           ratraw(ir33),dratrawdt(ir33),dratrawdd(ir33), &
           ratraw(ir33inv),dratrawdt(ir33inv),dratrawdd(ir33inv))
       inap = inap + 2


! he3(p,e+nu)he4
       call rate_hep(btemp,bden, &
                   ratraw(irhep),dratrawdt(irhep),dratrawdd(irhep), &
                   rrate,drratedt,drratedd)
       inap = inap + 1

!        if (weak_on .eq. 0) then
!         ratraw(irhep)    = 0.0d0
!         dratrawdt(irhep) = 0.0d0
!         dratrawdd(irhep) = 0.0d0
!        end if


! he3(n,g)he4
       call rate_he3ng(btemp,bden, &
                       ratraw(irhe3ng),dratrawdt(irhe3ng),dratrawdd(irhe3ng), &
                       ratraw(irhe4gn),dratrawdt(irhe4gn),dratrawdd(irhe4gn))
       inap = inap + 2



! he3(a,g)be7
      if (old_bdat_format)      & 
       call rate_he3he4(btemp,bden, &
                    sigraw(13,ihe3),sigrawdt(13,ihe3),sigrawdd(13,ihe3), &
                    sigraw(14,ihe3),sigrawdt(14,ihe3),sigrawdd(14,ihe3))
      end if



! if we have he3 and deuterium

      if (ihe3 .ne. 0 .and. ih2 .ne. 0) then

! he3(d,p)he4
       call rate_he3dp(btemp,bden, &
                       ratraw(irhe3dp),dratrawdt(irhe3dp),dratrawdd(irhe3dp), &
                       ratraw(irhe4pd),dratrawdt(irhe4pd),dratrawdd(irhe4pd))
       inap = inap + 2

! d(d,n)he3
       call rate_ddn(btemp,bden, &
                     ratraw(irddn),dratrawdt(irddn),dratrawdd(irddn), &
                     ratraw(irhe3nd),dratrawdt(irhe3nd),dratrawdd(irhe3nd))
       inap = inap + 2


! d(p,g)he3
       if (old_bdat_format)      & 
       call rate_dpg(btemp,bden, &
                     sigraw(7,ih2),sigrawdt(7,ih2),sigrawdd(7,ih2), &
                     sigraw(8,ih2),sigrawdt(8,ih2),sigrawdd(8,ih2))
      end if



! if we have he3, deuterium, tritium

      if (ihe3 .ne. 0  .and. ih2 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,d)he4
       call rate_he3td(btemp,bden, &
           ratraw(irhe3td),dratrawdt(irhe3td),dratrawdd(irhe3td), &
           ratraw(irhe4dt),dratrawdt(irhe4dt),dratrawdd(irhe4dt))
       inap = inap + 2
      end if


! if we have he3 and tritium

      if (ihe3 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,np)he4
       call rate_he3tnp(btemp,bden, &
                        ratraw(irhe3tnp),dratrawdt(irhe3tnp),dratrawdd(irhe3tnp), &
                        ratraw(irhe4pnt),dratrawdt(irhe4pnt),dratrawdd(irhe4pnt))
       inap = inap + 2
      end if




! if we have li6
      if (ili6 .ne. 0) then

! d(a,g)li6
       if (ili6 .ne. 0 .and. ih2 .ne. 0) then
        if (old_bdat_format)      & 
        call rate_he4dg(btemp,bden, &
                        sigraw(13,ih2),sigrawdt(13,ih2),sigrawdd(13,ih2), &
                        sigraw(14,ih2),sigrawdt(14,ih2),sigrawdd(14,ih2))
       end if

! t(a,n)li6
       if (ili6 .ne. 0 .and. ih3 .ne. 0) then
        if (old_bdat_format)      & 
        call rate_he4tn(btemp,bden, &
                        sigraw(11,ih3),sigrawdt(11,ih3),sigrawdd(11,ih3), &
                        sigraw(12,ih3),sigrawdt(12,ih3),sigrawdd(12,ih3))
       end if

! he3(a,p)li6 - forward rate is given as li6 + p, so note the flip of sigraw(9,i) and sigraw(10,i)
       if (ili6 .ne. 0 .and. ihe3 .ne. 0) then
        if (old_bdat_format)      & 
         call rate_li6phe3(btemp,bden, &
                      sigraw(10,ihe3),sigrawdt(10,ihe3),sigrawdd(10,ihe3), &
                      sigraw(9,ihe3),sigrawdt(9,ihe3),sigrawdd(9,ihe3))
       end if

! li6(n,g)li7
       if (ili6 .ne. 0 .and. ili7 .ne. 0) then
        if (old_bdat_format)      & 
        call rate_li6ng(btemp,bden, &
                      sigraw(1,ili6),sigrawdt(1,ili6),sigrawdd(1,ili6), &
                      sigraw(2,ili6),sigrawdt(2,ili6),sigrawdd(2,ili6))
       end if

! li6(p,g)be7
       if (ili6 .ne. 0 .and. ibe7 .ne. 0) then
        if (old_bdat_format)      & 
        call rate_li6pg(btemp,bden, &
                        sigraw(7,ili6),sigrawdt(7,ili6),sigrawdd(7,ili6), &
                        sigraw(8,ili6),sigrawdt(8,ili6),sigrawdd(8,ili6))
       end if

! li6(a,p)be9
       if (ili6 .ne. 0 .and. ibe9 .ne. 0) then
        if (old_bdat_format)      & 
        call rate_be9pa(btemp,bden, &
                      sigraw(10,ili6),sigrawdt(10,ili6),sigrawdd(10,ili6), &
                      sigraw(9,ili6),sigrawdt(9,ili6),sigrawdd(9,ili6))
       end if

! li6(a,g)b10
       if (ili6 .ne. 0 .and. ib10 .ne. 0) then
        if (old_bdat_format)      & 
        call rate_li6ag(btemp,bden, &
                      sigraw(13,ili6),sigrawdt(13,ili6),sigrawdd(13,ili6), &
                      sigraw(14,ili6),sigrawdt(14,ili6),sigrawdd(14,ili6))
       end if
      end if



! if we have li7 
      if (ili7 .ne. 0)then


! li7(p,g)be8 and li7(p,a)he4
      call rate_li7pag(btemp,bden, &
                       ratraw(irli7pag),dratrawdt(irli7pag),dratrawdd(irli7pag), &
                       ratraw(ir2he4ga),dratrawdt(ir2he4ga),dratrawdd(ir2he4ga))
       inap = inap + 2


! li7(n,g)li8
      if (old_bdat_format)      & 
      call rate_li7ng(btemp,bden, &
                      sigraw(1,ili7),sigrawdt(1,ili7),sigrawdd(1,ili7), &
                      sigraw(2,ili7),sigrawdt(2,ili7),sigrawdd(2,ili7))

! li7(p,n)be7
      if (old_bdat_format)      & 
      call rate_li7pn(btemp,bden, &
                      sigraw(3,ili7),sigrawdt(3,ili7),sigrawdd(3,ili7), &
                      sigraw(4,ili7),sigrawdt(4,ili7),sigrawdd(4,ili7))
! li7(a,g)b11
      if (old_bdat_format)      & 
      call rate_li7ag(btemp,bden, &
                      sigraw(13,ili7),sigrawdt(13,ili7),sigrawdd(13,ili7), &
                      sigraw(14,ili7),sigrawdt(14,ili7),sigrawdd(14,ili7))

! li7(a,n)b10   q = -2.790
      if (old_bdat_format)      & 
      call rate_li7an(btemp,bden, &
                      sigraw(11,ili7),sigrawdt(11,ili7),sigrawdd(11,ili7), &
                      sigraw(12,ili7),sigrawdt(12,ili7),sigrawdd(12,ili7))


! be7(e-,nu+g)li7
      if (old_bdat_format)      & 
      call rate_be7em(btemp,bden,ye, &
                      sigraw(6,ili7),sigrawdt(6,ili7),sigrawdd(6,ili7), &
                      rrate,drratedt,drratedd)

      end if



! if we have li7 and deuterium

      if (ili7 .ne. 0 .and. ih2 .ne. 0) then

! li7(d,n)2a
      call rate_li7dn(btemp,bden, &
                      ratraw(irli7dn),dratrawdt(irli7dn),dratrawdd(irli7dn), &
                      ratraw(ir2and),dratrawdt(ir2and),dratrawdd(ir2and))
       inap = inap + 2
      end if


! if we have li7 and tritium

      if (ili7 .ne. 0 .and. ih3 .ne. 0) then

! li7(t,2n)2a
      call rate_li7t2n(btemp,bden, &
                       ratraw(irli7t2n),dratrawdt(irli7t2n),dratrawdd(irli7t2n), &
                       ratraw(ir2a2nt),dratrawdt(ir2a2nt),dratrawdd(ir2a2nt))
       inap = inap + 2

! t(a,g)li7
      if (old_bdat_format)      & 
      call rate_he4tg(btemp,bden, &
                      sigraw(13,ih3),sigrawdt(13,ih3),sigrawdd(13,ih3), &
                      sigraw(14,ih3),sigrawdt(14,ih3),sigrawdd(14,ih3))
      end if


! if we have li7 and he3

      if (ili7 .ne. 0 .and. ihe3 .ne. 0) then

! li7(he3,np)2a
       call rate_li7he3np(btemp,bden, &
                          ratraw(irli7he3np),dratrawdt(irli7he3np),dratrawdd(irli7he3np), &
                          ratraw(ir2apnhe3),dratrawdt(ir2apnhe3),dratrawdd(ir2apnhe3))
       inap = inap + 2

!       ratraw(irli7he3np) = 0.0d0 ; dratrawdt(irli7he3np) = 0.0d0 ; dratrawdd(irli7he3np) = 0.0d0
!       ratraw(ir2apnhe3) = 0.0d0  ; dratrawdt(ir2apnhe3) = 0.0d0  ; dratrawdd(ir2apnhe3) = 0.0d0
      end if



! if we have be7 and deutrium
! be7(d,p)2a
      if (ibe7 .ne. 0  .and. ih2 .ne. 0) then
       call rate_be7dp(btemp,bden, &
                       ratraw(irbe7dp),dratrawdt(irbe7dp),dratrawdd(irbe7dp), &
                       ratraw(ir2apd),dratrawdt(ir2apd),dratrawdd(ir2apd))
       inap = inap + 2
      end if

! if we have be7 and tritium
! be7(t,np)2a
      if (ibe7 .ne. 0  .and. ih3 .ne. 0) then
       call rate_be7tnp(btemp,bden, &
                        ratraw(irbe7tnp),dratrawdt(irbe7tnp),dratrawdd(irbe7tnp), &
                        ratraw(ir2apnt),dratrawdt(ir2apnt),dratrawdd(ir2apnt))
       inap = inap + 2

!       ratraw(irbe7tnp) = 0.0d0 ; dratrawdt(irbe7tnp) = 0.0d0 ; dratrawdd(irbe7tnp) = 0.0d0
!       ratraw(ir2apnt) = 0.0d0  ; dratrawdt(ir2apnt)= 0.0d0   ; dratrawdd(ir2apnt) = 0.0d0
      end if

! if we have be7 and he3
! be7(he3,2p)2a
      if (ibe7 .ne. 0  .and. ihe3 .ne. 0) then
       call rate_be7he32p(btemp,bden, &
                          ratraw(irbe7he32p),dratrawdt(irbe7he32p),dratrawdd(irbe7he32p), &
                          ratraw(ir2a2phe3),dratrawdt(ir2a2phe3),dratrawdd(ir2a2phe3))
       inap = inap + 2

!       ratraw(irbe7he32p) = 0.0d0 ; dratrawdt(irbe7he32p) = 0.0d0 ; dratrawdd(irbe7he32p) = 0.0d0
!       ratraw(ir2a2phe3) = 0.0d0  ; dratrawdt(ir2a2phe3) = 0.0d0  ; dratrawdd(ir2a2phe3) = 0.0d0
      end if

! be7(p,g)b8
      if (ibe7 .ne. 0  .and. ib8 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_be7pg(btemp,bden, &
                       sigraw(7,ibe7),sigrawdt(7,ibe7),sigrawdd(7,ibe7), &
                       sigraw(8,ibe7),sigrawdt(8,ibe7),sigrawdd(8,ibe7))
      end if

! b10(p,a)be7
      if (ibe7 .ne. 0  .and. ib10 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b10pa(btemp,bden, &
                     sigraw(10,ibe7),sigrawdt(10,ibe7),sigrawdd(10,ibe7), &
                     sigraw(9,ibe7),sigrawdt(9,ibe7),sigrawdd(9,ibe7))
      end if

! be7(a,g)c11
      if (ibe7 .ne. 0  .and. ic11 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_be7ag(btemp,bden, &
                     sigraw(13,ibe7),sigrawdt(13,ibe7),sigrawdd(13,ibe7), &
                     sigraw(14,ibe7),sigrawdt(14,ibe7),sigrawdd(14,ibe7))
      end if




! if we have be9
! a(an,g)be9
      if (ibe9 .ne. 0) then
       call rate_aan(btemp,bden, &
                     ratraw(iraan),dratrawdt(iraan),dratrawdd(iraan), &
                     ratraw(irgaan),dratrawdt(irgaan),dratrawdd(irgaan))
       inap = inap + 2
      end if


! if we have be9 and deuterium
! be9(p,d)be8
      if (ibe9 .ne. 0 .and. ih2 .ne. 0) then
       call rate_be9pd(btemp,bden, &
                       ratraw(irbe9pd),dratrawdt(irbe9pd),dratrawdd(irbe9pd), &
                       ratraw(ir2adp),dratrawdt(ir2adp),dratrawdd(ir2adp))
       inap = inap + 2
      end if


! be9(p,n)b9
      if (ibe9 .ne. 0 .and. ib9 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_be9pn(btemp,bden, &
                       sigraw(3,ibe9),sigrawdt(3,ibe9),sigrawdd(3,ibe9), &
                       sigraw(4,ibe9),sigrawdt(4,ibe9),sigrawdd(4,ibe9))
      end if

! be9(p,g)b10
      if (ibe9 .ne. 0 .and. ib10 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_be9pg(btemp,bden, &
                       sigraw(7,ibe9),sigrawdt(7,ibe9),sigrawdd(7,ibe9), &
                       sigraw(8,ibe9),sigrawdt(8,ibe9),sigrawdd(8,ibe9))
      end if


! be9(a,n)c12
      if (ibe9 .ne. 0 .and. ic12 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_be9an(btemp,bden, &
                       sigraw(11,ibe9),sigrawdt(11,ibe9),sigrawdd(11,ibe9), &
                       sigraw(12,ibe9),sigrawdt(12,ibe9),sigrawdd(12,ibe9))
      end if



! if we have b8
! b8(e+,nu)be8 => 2a
      if (ib8 .ne. 0) then
       call rate_b8ep(btemp,bden, &
                      ratraw(irb8ep),dratrawdt(irb8ep),dratrawdd(irb8ep), &
                      rrate,drratedt,drratedd)
       inap = inap + 1

!       if (weak_on .eq. 0) then
!        ratraw(irb8ep)    = 0.0d0
!        dratrawdt(irb8ep) = 0.0d0
!        dratrawdd(irb8ep) = 0.0d0
!       end if
      end if

! b8(a,p)c11
      if (ib8 .ne. 0 .and. ic11 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b8ap(btemp,bden, &
                      sigraw(9,ib8),sigrawdt(9,ib8),sigrawdd(9,ib8), &
                      sigraw(10,ib8),sigrawdt(10,ib8),sigrawdd(10,ib8))
      end if





! if we have b10
! b10(p,g)c11
      if (ib10 .ne. 0 .and. ic11 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b10pg(btemp,bden, &
                     sigraw(7,ib10),sigrawdt(7,ib10),sigrawdd(7,ib10), &
                     sigraw(8,ib10),sigrawdt(8,ib10),sigrawdd(8,ib10))
       end if

! b10(a,n)n13
      if (ib10 .ne. 0 .and. in13 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b10an(btemp,bden, &
                     sigraw(11,ib10),sigrawdt(11,ib10),sigrawdd(11,ib10), &
                     sigraw(12,ib10),sigrawdt(12,ib10),sigrawdd(12,ib10))
      end if



! if we have b11
! b11(p,a)be8=>2a
      if (ib11 .ne. 0) then
      call rate_b11pa(btemp,bden, &
                      ratraw(irb11pa),dratrawdt(irb11pa),dratrawdd(irb11pa), &
                      ratraw(ir3ap),dratrawdt(ir3ap),dratrawdd(ir3ap))
       inap = inap + 2
      end if



! b11(p,n)c11
      if (ib11 .ne. 0 .and. ic11 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b11pn(btemp,bden, &
                       sigraw(3,ib11),sigrawdt(3,ib11),sigrawdd(3,ib11), &
                       sigraw(4,ib11),sigrawdt(4,ib11),sigrawdd(4,ib11))
      end if


! b11(p,g)c12
      if (ib11 .ne. 0 .and. ic12 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b11pg(btemp,bden, &
                       sigraw(7,ib11),sigrawdt(7,ib11),sigrawdd(7,ib11), &
                       sigraw(8,ib11),sigrawdt(8,ib11),sigrawdd(8,ib11))
      end if



! b11(a,p)c14
      if (ib11 .ne. 0 .and. ic14 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b11ap(btemp,bden, &
                       sigraw(9,ib11),sigrawdt(9,ib11),sigrawdd(9,ib11), &
                       sigraw(10,ib11),sigrawdt(10,ib11),sigrawdd(10,ib11))
       end if


! b11(a,n)n14
      if (ib11 .ne. 0 .and. in14 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_b11an(btemp,bden, &
                       sigraw(11,ib11),sigrawdt(11,ib11),sigrawdd(11,ib11), &
                       sigraw(12,ib11),sigrawdt(12,ib11),sigrawdd(12,ib11))
      end if




! if we have c11
! c11(n,a)be8 => 2a
      if (ic11 .ne. 0) then
       call rate_c11na(btemp,bden, &
                       ratraw(irc11na),dratrawdt(irc11na),dratrawdd(irc11na), &
                       ratraw(ir3an),dratrawdt(ir3an),dratrawdd(ir3an))
       inap = inap + 2
      end if

! n14(p,a)c11
      if (ic11 .ne. 0 .and. in14 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_n14pa(btemp,bden, &
                       sigraw(10,ic11),sigrawdt(10,ic11),sigrawdd(10,ic11), &
                       sigraw(9,ic11),sigrawdt(9,ic11),sigrawdd(9,ic11))
      end if

! c11(p,g)n12
      if (ic11 .ne. 0 .and. in12 .ne. 0) then
       if (old_bdat_format)      & 
       call rate_c11pg(btemp,bden, &
                       sigraw(7,ic11),sigrawdt(7,ic11),sigrawdd(7,ic11), &
                       sigraw(8,ic11),sigrawdt(8,ic11),sigrawdd(8,ic11))
      end if






! for c12
      if (ic12 .ne. 0) then

! c12(p,g)n13
      if (old_bdat_format)      & 
      call rate_c12pg(btemp,bden, &
                    sigraw(7,ic12),sigrawdt(7,ic12),sigrawdd(7,ic12), &
                    sigraw(8,ic12),sigrawdt(8,ic12),sigrawdd(8,ic12))

! c12(a,n)o15
      if (old_bdat_format)      & 
      call rate_c12an(btemp,bden, &
                    sigraw(11,ic12),sigrawdt(11,ic12),sigrawdd(11,ic12), &
                    sigraw(12,ic12),sigrawdt(12,ic12),sigrawdd(12,ic12))

! c12(a,g)o16
      if (old_bdat_format)      & 
      call rate_c12ag(btemp,bden, &
                    sigraw(13,ic12),sigrawdt(13,ic12),sigrawdd(13,ic12), &
                    sigraw(14,ic12),sigrawdt(14,ic12),sigrawdd(14,ic12))

! n15(p,a)c12
      if (old_bdat_format)      & 
      call rate_n15pa(btemp,bden, &
                    sigraw(10,ic12),sigrawdt(10,ic12),sigrawdd(10,ic12), &
                    sigraw(9,ic12),sigrawdt(9,ic12),sigrawdd(9,ic12))

      end if





! for c13
      if (ic13 .ne. 0) then

! c13(p,g)n14
      if (old_bdat_format)      & 
      call rate_c13pg(btemp,bden, &
                    sigraw(7,ic13),sigrawdt(7,ic13),sigrawdd(7,ic13), &
                    sigraw(8,ic13),sigrawdt(8,ic13),sigrawdd(8,ic13))

! c13(a,n)o16
      if (old_bdat_format)      & 
      call rate_c13an(btemp,bden, &
                    sigraw(11,ic13),sigrawdt(11,ic13),sigrawdd(11,ic13), &
                    sigraw(12,ic13),sigrawdt(12,ic13),sigrawdd(12,ic13))

! c13(p,n)n13
      if (old_bdat_format)      & 
      call rate_c13pn(btemp,bden, &
                    sigraw(3,ic13),sigrawdt(3,ic13),sigrawdd(3,ic13), &
                    sigraw(4,ic13),sigrawdt(4,ic13),sigrawdd(4,ic13))

      end if





! for c14
      if (ic14 .ne. 0) then

! c14(p,g)n15
      if (old_bdat_format)      & 
      call rate_c14pg(btemp,bden, &
                    sigraw(7,ic14),sigrawdt(7,ic14),sigrawdd(7,ic14), &
                    sigraw(8,ic14),sigrawdt(8,ic14),sigrawdd(8,ic14))

! c14(p,n)n14
      if (old_bdat_format)      & 
      call rate_c14pn(btemp,bden, &
                    sigraw(3,ic14),sigrawdt(3,ic14),sigrawdd(3,ic14), &
                    sigraw(4,ic14),sigrawdt(4,ic14),sigrawdd(4,ic14))

! c14(a,g)o18
      if (old_bdat_format)      & 
      call rate_c14ag(btemp,bden, &
                    sigraw(13,ic14),sigrawdt(13,ic14),sigrawdd(13,ic14), &
                    sigraw(14,ic14),sigrawdt(14,ic14),sigrawdd(14,ic14))

      end if




! for n13
       if (in13 .ne. 0) then

! n13(p,g)o14
      if (old_bdat_format)      & 
      call rate_n13pg(btemp,bden, &
                    sigraw(7,in13),sigrawdt(7,in13),sigrawdd(7,in13), &
                    sigraw(8,in13),sigrawdt(8,in13),sigrawdd(8,in13))

! o16(p,a)n13
      if (old_bdat_format)      & 
      call rate_o16pa(btemp,bden, &
                    sigraw(10,in13),sigrawdt(10,in13),sigrawdd(10,in13), &
                    sigraw(9,in13),sigrawdt(9,in13),sigrawdd(9,in13))

       end if




! for n14
      if (in14 .ne. 0) then

! n14(p,n)o14
      if (old_bdat_format)      & 
      call rate_n14pn(btemp,bden, &
                    sigraw(3,in14),sigrawdt(3,in14),sigrawdd(3,in14), &
                    sigraw(4,in14),sigrawdt(4,in14),sigrawdd(4,in14))

! n14(p,g)o15
      if (old_bdat_format)      & 
      call rate_n14pg(btemp,bden, &
                    sigraw(7,in14),sigrawdt(7,in14),sigrawdd(7,in14), &
                    sigraw(8,in14),sigrawdt(8,in14),sigrawdd(8,in14))

! n14(a,n)f17
      if (old_bdat_format)      & 
      call rate_n14an(btemp,bden, &
                    sigraw(11,in14),sigrawdt(11,in14),sigrawdd(11,in14), &
                    sigraw(12,in14),sigrawdt(12,in14),sigrawdd(12,in14))


! n14(a,g)f18
      if (old_bdat_format)      & 
      call rate_n14ag(btemp,bden, &
                    sigraw(13,in14),sigrawdt(13,in14),sigrawdd(13,in14), &
                    sigraw(14,in14),sigrawdt(14,in14),sigrawdd(14,in14))

! o17(p,a)n14
      if (old_bdat_format)      & 
      call rate_o17pa(btemp,bden, &
                    sigraw(10,in14),sigrawdt(10,in14),sigrawdd(10,in14), &
                    sigraw(9,in14),sigrawdt(9,in14),sigrawdd(9,in14))

      end if





! for n15
      if (in15 .ne. 0)  then

! n15(p,g)o16
      if (old_bdat_format)      & 
      call rate_n15pg(btemp,bden, &
                    sigraw(7,in15),sigrawdt(7,in15),sigrawdd(7,in15), &
                    sigraw(8,in15),sigrawdt(8,in15),sigrawdd(8,in15))

! n15(p,n)o15
      if (old_bdat_format)      & 
      call rate_n15pn(btemp,bden, &
                    sigraw(3,in15),sigrawdt(3,in15),sigrawdd(3,in15), &
                    sigraw(4,in15),sigrawdt(4,in15),sigrawdd(4,in15))

! n15(a,n)f18
      if (old_bdat_format)      & 
      call rate_n15an(btemp,bden, &
                    sigraw(11,in15),sigrawdt(11,in15),sigrawdd(11,in15), &
                    sigraw(12,in15),sigrawdt(12,in15),sigrawdd(12,in15))

! n15(a,g)f19
      if (old_bdat_format)      & 
      call rate_n15ag(btemp,bden, &
                    sigraw(13,in15),sigrawdt(13,in15),sigrawdd(13,in15), &
                    sigraw(14,in15),sigrawdt(14,in15),sigrawdd(14,in15))

! o18(p,a)n15
      if (old_bdat_format)      & 
      call rate_o18pa(btemp,bden, &
                    sigraw(10,in15),sigrawdt(10,in15),sigrawdd(10,in15), &
                    sigraw(9,in15),sigrawdt(9,in15),sigrawdd(9,in15))

      end if




! for o14
      if (io14 .ne. 0) then

! o14(a,g)ne18
      if (old_bdat_format)      & 
      call rate_o14ag(btemp,bden, &
                    sigraw(13,io14),sigrawdt(13,io14),sigrawdd(13,io14), &
                    sigraw(14,io14),sigrawdt(14,io14),sigrawdd(14,io14))

! o14(a,p)f17
      if (old_bdat_format)      & 
      call rate_o14ap(btemp,bden, &
                    sigraw(9,io14),sigrawdt(9,io14),sigrawdd(9,io14), &
                    sigraw(10,io14),sigrawdt(10,io14),sigrawdd(10,io14))

      end if




! for o15
      if (io15 .ne. 0) then

! o15(a,g)ne19
      if (old_bdat_format)      & 
      call rate_o15ag(btemp,bden, &
                    sigraw(13,io15),sigrawdt(13,io15),sigrawdd(13,io15), &
                    sigraw(14,io15),sigrawdt(14,io15),sigrawdd(14,io15))


! f18(p,a)o15
      if (old_bdat_format)      & 
      call rate_f18pa(btemp,bden, &
                    sigraw(10,io15),sigrawdt(10,io15),sigrawdd(10,io15), &
                    sigraw(9,io15),sigrawdt(9,io15),sigrawdd(9,io15))

      end if





! for o16
      if (io16 .ne. 0) then

! o16(p,g)f17
      if (old_bdat_format)      & 
      call rate_o16pg(btemp,bden, &
                    sigraw(7,io16),sigrawdt(7,io16),sigrawdd(7,io16), &
                    sigraw(8,io16),sigrawdt(8,io16),sigrawdd(8,io16))

! o16(a,g)ne20
      if (old_bdat_format)      & 
      call rate_o16ag(btemp,bden, &
                    sigraw(13,io16),sigrawdt(13,io16),sigrawdd(13,io16), &
                    sigraw(14,io16),sigrawdt(14,io16),sigrawdd(14,io16))

! f19(p,a)o16
      if (old_bdat_format)      & 
      call rate_f19pa(btemp,bden, &
                    sigraw(10,io16),sigrawdt(10,io16),sigrawdd(10,io16), &
                    sigraw(9,io16),sigrawdt(9,io16),sigrawdd(9,io16))

      end if





! for o17
      if (io17 .ne. 0) then

! o17(p,g)f18
      if (old_bdat_format)      & 
      call rate_o17pg(btemp,bden, &
                    sigraw(7,io17),sigrawdt(7,io17),sigrawdd(7,io17), &
                    sigraw(8,io17),sigrawdt(8,io17),sigrawdd(8,io17))

! o17(a,g)ne21
      if (old_bdat_format)      & 
      call rate_o17ag(btemp,bden, &
                    sigraw(13,io17),sigrawdt(13,io17),sigrawdd(13,io17), &
                    sigraw(14,io17),sigrawdt(14,io17),sigrawdd(14,io17))

! o17(a,n)ne20
      if (old_bdat_format)      & 
      call rate_o17an(btemp,bden, &
                    sigraw(11,io17),sigrawdt(11,io17),sigrawdd(11,io17), &
                    sigraw(12,io17),sigrawdt(12,io17),sigrawdd(12,io17))

      end if





! for o18
      if (io18 .ne. 0) then

! o18(p,g)f19
      if (old_bdat_format)      & 
      call rate_o18pg(btemp,bden, &
                    sigraw(7,io18),sigrawdt(7,io18),sigrawdd(7,io18), &
                    sigraw(8,io18),sigrawdt(8,io18),sigrawdd(8,io18))

! o18(a,g)ne22
      if (old_bdat_format)      & 
      call rate_o18ag(btemp,bden, &
                    sigraw(13,io18),sigrawdt(13,io18),sigrawdd(13,io18), &
                    sigraw(14,io18),sigrawdt(14,io18),sigrawdd(14,io18))

! o18(a,n)ne21
      if (old_bdat_format)      & 
      call rate_o18an(btemp,bden, &
                    sigraw(11,io18),sigrawdt(11,io18),sigrawdd(11,io18), &
                    sigraw(12,io18),sigrawdt(12,io18),sigrawdd(12,io18))

      end if





! for f17
      if (if17 .ne. 0) then

! f17(p,g)ne18
      if (old_bdat_format)      & 
      call rate_f17pg(btemp,bden, &
                    sigraw(7,if17),sigrawdt(7,if17),sigrawdd(7,if17), &
                    sigraw(8,if17),sigrawdt(8,if17),sigrawdd(8,if17))


! ne20(p,a)f17
      if (old_bdat_format)      & 
      call rate_ne20pa(btemp,bden, &
                    sigraw(10,if17),sigrawdt(10,if17),sigrawdd(10,if17), &
                    sigraw(9,if17),sigrawdt(9,if17),sigrawdd(9,if17))

      end if





! for f18
      if (if18 .ne. 0) then

! f18(p,g)ne19
      if (old_bdat_format)      & 
      call rate_f18pg(btemp,bden, &
                    sigraw(7,if18),sigrawdt(7,if18),sigrawdd(7,if18), &
                    sigraw(8,if18),sigrawdt(8,if18),sigrawdd(8,if18))

      end if





! for f19
      if (if19 .ne. 0) then

! f19(p,g)ne20
      if (old_bdat_format)      & 
      call rate_f19pg(btemp,bden, &
                    sigraw(7,if19),sigrawdt(7,if19),sigrawdd(7,if19), &
                    sigraw(8,if19),sigrawdt(8,if19),sigrawdd(8,if19))

! f19(p,n)ne19
      if (old_bdat_format)      & 
      call rate_f19pn(btemp,bden, &
                    sigraw(3,if19),sigrawdt(3,if19),sigrawdd(3,if19), &
                    sigraw(4,if19),sigrawdt(4,if19),sigrawdd(4,if19))

! f19(a,p)ne22
      if (old_bdat_format)      & 
      call rate_f19ap(btemp,bden, &
                    sigraw(9,if19),sigrawdt(9,if19),sigrawdd(9,if19), &
                    sigraw(10,if19),sigrawdt(10,if19),sigrawdd(10,if19))

! na22(n,a)f19
      if (old_bdat_format)      & 
      call rate_na22na(btemp,bden, &
                    sigraw(12,if19),sigrawdt(12,if19),sigrawdd(12,if19), &
                    sigraw(11,if19),sigrawdt(11,if19),sigrawdd(11,if19))

      end if





! for ne19
      if (ine19 .ne. 0) then

! ne19(p,g)na20
      if (old_bdat_format)      & 
      call rate_ne19pg(btemp,bden, &
                    sigraw(7,ine19),sigrawdt(7,ine19),sigrawdd(7,ine19), &
                    sigraw(8,ine19),sigrawdt(8,ine19),sigrawdd(8,ine19))

      end if





! for ne20
      if (ine20 .ne. 0) then

! ne20(p,g)na21
      if (old_bdat_format)      & 
      call rate_ne20pg(btemp,bden, &
                 sigraw(7,ine20),sigrawdt(7,ine20),sigrawdd(7,ine20), &
                 sigraw(8,ine20),sigrawdt(8,ine20),sigrawdd(8,ine20))

! ne20(a,g)mg24
      if (old_bdat_format)      & 
      call rate_ne20ag(btemp,bden, &
                 sigraw(13,ine20),sigrawdt(13,ine20),sigrawdd(13,ine20), &
                 sigraw(14,ine20),sigrawdt(14,ine20),sigrawdd(14,ine20))

! na23(p,a)ne20
      if (old_bdat_format)      & 
      call rate_na23pa(btemp,bden, &
                 sigraw(10,ine20),sigrawdt(10,ine20),sigrawdd(10,ine20), &
                 sigraw(9,ine20),sigrawdt(9,ine20),sigrawdd(9,ine20))

! ne20(n,g)ne21
      if (old_bdat_format)      & 
      call rate_ne20ng(btemp,bden, &
                 sigraw(1,ine20),sigrawdt(1,ine20),sigrawdd(1,ine20), &
                 sigraw(2,ine20),sigrawdt(2,ine20),sigrawdd(2,ine20))

      end if




! for ne21
      if (ine21 .ne. 0) then

! ne21(p,g)na22
      if (old_bdat_format)      & 
      call rate_ne21pg(btemp,bden, &
                 sigraw(7,ine21),sigrawdt(7,ine21),sigrawdd(7,ine21), &
                 sigraw(8,ine21),sigrawdt(8,ine21),sigrawdd(8,ine21))

! ne21(a,g)mg25
      if (old_bdat_format)      & 
      call rate_ne21ag(btemp,bden, &
                 sigraw(13,ine21),sigrawdt(13,ine21),sigrawdd(13,ine21), &
                 sigraw(14,ine21),sigrawdt(14,ine21),sigrawdd(14,ine21))

! ne21(a,n)mg24
      if (old_bdat_format)      & 
      call rate_ne21an(btemp,bden, &
                 sigraw(11,ine21),sigrawdt(11,ine21),sigrawdd(11,ine21), &
                 sigraw(12,ine21),sigrawdt(12,ine21),sigrawdd(12,ine21))

      end if





! for ne22
      if (ine22 .ne. 0) then

! ne22(p,g)na23
      if (old_bdat_format)      & 
      call rate_ne22pg(btemp,bden, &
                 sigraw(7,ine22),sigrawdt(7,ine22),sigrawdd(7,ine22), &
                 sigraw(8,ine22),sigrawdt(8,ine22),sigrawdd(8,ine22))

! ne22(a,g)mg26
      if (old_bdat_format)      & 
      call rate_ne22ag(btemp,bden, &
                 sigraw(13,ine22),sigrawdt(13,ine22),sigrawdd(13,ine22), &
                 sigraw(14,ine22),sigrawdt(14,ine22),sigrawdd(14,ine22))

! na22(n,p)ne22
      if (old_bdat_format)      & 
      call rate_na22np(btemp,bden, &
                 sigraw(4,ine22),sigrawdt(4,ine22),sigrawdd(4,ine22), &
                 sigraw(3,ine22),sigrawdt(3,ine22),sigrawdd(3,ine22))

! ne22(a,n)mg25
      if (old_bdat_format)      & 
      call rate_ne22an(btemp,bden, &
                 sigraw(11,ine22),sigrawdt(11,ine22),sigrawdd(11,ine22), &
                 sigraw(12,ine22),sigrawdt(12,ine22),sigrawdd(12,ine22))

      end if






! for na21
      if (ina21 .ne. 0) then

! na21(p,g)mg22
      if (old_bdat_format)      & 
      call rate_na21pg(btemp,bden, &
                 sigraw(7,ina21),sigrawdt(7,ina21),sigrawdd(7,ina21), &
                 sigraw(8,ina21),sigrawdt(8,ina21),sigrawdd(8,ina21))

! mg24(p,a)na21
      if (old_bdat_format)      & 
      call rate_mg24pa(btemp,bden, &
                 sigraw(10,ina21),sigrawdt(10,ina21),sigrawdd(10,ina21), &
                 sigraw(9,ina21),sigrawdt(9,ina21),sigrawdd(9,ina21))

      end if



! for na22
      if (ina22 .ne. 0) then

! na22(p,g)mg23
      if (old_bdat_format)      & 
      call rate_na22pg(btemp,bden, &
                 sigraw(7,ina22),sigrawdt(7,ina22),sigrawdd(7,ina22), &
                 sigraw(8,ina22),sigrawdt(8,ina22),sigrawdd(8,ina22))

      end if




! for na23
      if (ina23 .ne. 0) then

! na23(p,g)mg24
      if (old_bdat_format)      & 
      call rate_na23pg(btemp,bden, &
                 sigraw(7,ina23),sigrawdt(7,ina23),sigrawdd(7,ina23), &
                 sigraw(8,ina23),sigrawdt(8,ina23),sigrawdd(8,ina23))

! na23(p,n)mg23
      if (old_bdat_format)      & 
      call rate_na23pn(btemp,bden, &
                 sigraw(3,ina23),sigrawdt(3,ina23),sigrawdd(3,ina23), &
                 sigraw(4,ina23),sigrawdt(4,ina23),sigrawdd(4,ina23))

      end if





! for mg24
      if (img24 .ne. 0) then

! mg24(p,g)al25
      if (old_bdat_format)      & 
      call rate_mg24pg(btemp,bden, &
                 sigraw(7,img24),sigrawdt(7,img24),sigrawdd(7,img24), &
                 sigraw(8,img24),sigrawdt(8,img24),sigrawdd(8,img24))

! mg24(a,g)si28
      if (old_bdat_format)      & 
      call rate_mg24ag(btemp,bden, &
                 sigraw(13,img24),sigrawdt(13,img24),sigrawdd(13,img24), &
                 sigraw(14,img24),sigrawdt(14,img24),sigrawdd(14,img24))

! al27(p,a)mg24
      if (old_bdat_format)      & 
      call rate_al27pa(btemp,bden, &
                 sigraw(10,img24),sigrawdt(10,img24),sigrawdd(10,img24), &
                 sigraw(9,img24),sigrawdt(9,img24),sigrawdd(9,img24))

      end if






! for mg25
      if (img25 .ne. 0) then

! mg25(p,g)al26
      if (old_bdat_format)      & 
      call rate_mg25pg(btemp,bden, &
                 sigraw(7,img25),sigrawdt(7,img25),sigrawdd(7,img25), &
                 sigraw(8,img25),sigrawdt(8,img25),sigrawdd(8,img25))

! mg25(a,p)al28
      if (old_bdat_format)      & 
      call rate_mg25ap(btemp,bden, &
                 sigraw(9,img25),sigrawdt(9,img25),sigrawdd(9,img25), &
                 sigraw(10,img25),sigrawdt(10,img25),sigrawdd(10,img25))

! mg25(a,g)si29
      if (old_bdat_format)      & 
      call rate_mg25ag(btemp,bden, &
                 sigraw(13,img25),sigrawdt(13,img25),sigrawdd(13,img25), &
                 sigraw(14,img25),sigrawdt(14,img25),sigrawdd(14,img25))

! mg25(a,n)si28
      if (old_bdat_format)      & 
      call rate_mg25an(btemp,bden, &
                 sigraw(11,img25),sigrawdt(11,img25),sigrawdd(11,img25), &
                 sigraw(12,img25),sigrawdt(12,img25),sigrawdd(12,img25))

      end if






! for mg26
      if (img26 .ne. 0) then

! mg26(p,g)al27
      if (old_bdat_format)      & 
      call rate_mg26pg(btemp,bden, &
                 sigraw(7,img26),sigrawdt(7,img26),sigrawdd(7,img26), &
                 sigraw(8,img26),sigrawdt(8,img26),sigrawdd(8,img26))

! mg26(a,g)si30
      if (old_bdat_format)      & 
      call rate_mg26ag(btemp,bden, &
                 sigraw(13,img26),sigrawdt(13,img26),sigrawdd(13,img26), &
                 sigraw(14,img26),sigrawdt(14,img26),sigrawdd(14,img26))

! mg26(a,n)si29
      if (old_bdat_format)      & 
      call rate_mg26an(btemp,bden, &
                 sigraw(11,img26),sigrawdt(11,img26),sigrawdd(11,img26), &
                 sigraw(12,img26),sigrawdt(12,img26),sigrawdd(12,img26))
      end if




! for al25
      if (ial25 .ne. 0) then

! al25(p,g)si26
      if (old_bdat_format)      & 
      call rate_al25pg(btemp,bden, &
                 sigraw(7,ial25),sigrawdt(7,ial25),sigrawdd(7,ial25), &
                 sigraw(8,ial25),sigrawdt(8,ial25),sigrawdd(8,ial25))

      endif





! for al26
      if (ial26 .ne. 0) then

! al26(p,g)si27
      if (old_bdat_format)      & 
      call rate_al26pg(btemp,bden, &
                 sigraw(7,ial26),sigrawdt(7,ial26),sigrawdd(7,ial26), &
                 sigraw(8,ial26),sigrawdt(8,ial26),sigrawdd(8,ial26))

      end if





! for al27
      if (ial27 .ne. 0) then

! al27(p,g)si28
      if (old_bdat_format)      & 
      call rate_al27pg(btemp,bden, &
                 sigraw(7,ial27),sigrawdt(7,ial27),sigrawdd(7,ial27), &
                 sigraw(8,ial27),sigrawdt(8,ial27),sigrawdd(8,ial27))

! al27(a,n)p30
      if (old_bdat_format)      & 
      call rate_al27an(btemp,bden, &
                 sigraw(11,ial27),sigrawdt(11,ial27),sigrawdd(11,ial27), &
                 sigraw(12,ial27),sigrawdt(12,ial27),sigrawdd(12,ial27))

      end if





! for si27
      if (isi27 .ne. 0) then

! si27(p,g)p28
      if (old_bdat_format)      & 
      call rate_si27pg(btemp,bden, &
                 sigraw(7,isi27),sigrawdt(7,isi27),sigrawdd(7,isi27), &
                 sigraw(8,isi27),sigrawdt(8,isi27),sigrawdd(8,isi27))

      end if




! for si28
      if (isi28 .ne. 0) then


! si28(a,g)s32
      if (old_bdat_format)      & 
      call rate_si28ag(btemp,bden, &
                 sigraw(13,isi28),sigrawdt(13,isi28),sigrawdd(13,isi28), &
                 sigraw(14,isi28),sigrawdt(14,isi28),sigrawdd(14,isi28))


! si28(p,g)p29
      if (old_bdat_format)      & 
      call rate_si28pg(btemp,bden, &
                 sigraw(7,isi28),sigrawdt(7,isi28),sigrawdd(7,isi28), &
                 sigraw(8,isi28),sigrawdt(8,isi28),sigrawdd(8,isi28))

      end if



! for si29
      if (isi29 .ne. 0) then


! si29(p,g)p30
      if (old_bdat_format)      & 
      call rate_si29pg(btemp,bden, &
                 sigraw(7,isi29),sigrawdt(7,isi29),sigrawdd(7,isi29), &
                 sigraw(8,isi29),sigrawdt(8,isi29),sigrawdd(8,isi29))
      end if



! for si30
      if (isi30 .ne. 0) then

! si30(p,g)p31
      if (old_bdat_format)      & 
      call rate_si30pg(btemp,bden, &
                 sigraw(7,isi30),sigrawdt(7,isi30),sigrawdd(7,isi30), &
                 sigraw(8,isi30),sigrawdt(8,isi30),sigrawdd(8,isi30))
      end if



! bullet check the counting
      if (inap .ne. nrat) then
       write(6,*)
       write(6,*) 'in torchrat inap =',inap
       write(6,*) 'in torchrat nrat =',nrat
       write(6,*) 'inap is not equal to nrat'
       write(6,*) 'fatal counting error'
       write(6,*)
       stop 'fatal counting error in torchrat'
      end if

      return
      end subroutine torchrat

      subroutine evaluate_rate(jat,iat,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'tfactors.dek'


! this routine generates nuclear reaction rates for the torch network
!
! declare the pass
      integer :: jat,iat
      real*8  :: fr,dfrdt,dfrdd,rr,drrdt,drrdd


! local variables
      integer :: j,k,jd1,jd2,jd3,jd4,jd5,jn,j1,j2,k0,k1,k2,k3, &
                 k4,k5,k6,k7,k8,k9,jxx, &
                 nrate,lk0,jk,lk1,lk2,lk3,lk4,lk5,lk6,lk7,lk8
      real*8  :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee,z



! initialize
      fr    = 0.0d0
      dfrdt = 0.0d0
      dfrdd = 0.0d0
      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0
      if (btemp .lt. 1.0e6) return

 
! set the reaction index j based on jat being
! even (reverse rates) or odd (foward rates)

      if (int(jat/2) .eq. int((jat+1)/2)) then
       j = jat/2
      else
       j = max(1,(jat - 1)/2 + 1)
      end if

! set the location index jn
      jn = j
      if (jn - 3 .ge. 0) jn = j-1


! bail if there is no link to the isotope or no stored formula coefficients

      if (nrr(jn,iat) .eq. 0  .or.  ic1(j,iat) .eq. 0) return


! at most 10 fitting constants for strong and electromagnetic rates
      k  = ic3(j,iat)
      k0 = k
      k1 = k + 1
      k2 = k + 2
      k3 = k + 3
      k4 = k + 4
      k5 = k + 5
      k6 = k + 6
      k7 = k + 7
      k8 = k + 8
      k9 = k + 9


! ground state plus first excited state positron decay
! excited state assumed to be in thermal equilibrium
      if (ic1(j,iat) .eq. 6) then

       aa  = cx(k2) * cx(k3) * exp(-11.60485d0*cx(k4)*t9i)
       daa = aa*11.60485d0*cx(k4)*t9i2

       bb  = cx(k0)*cx(k1) + aa
       dbb = daa

       cc  = cx(k2)*exp(-11.60485d0*cx(k4)*t9i)
       dcc = cc*11.60485d0*cx(k4)*t9i2

       dd  = cx(k0) + cc
       ddd = dcc

       z   = 1.0d0/dd
       ee  = bb*z
       dee = (dbb - ee*ddd)*z

       rr   = ee
       drrdt = dee * 1.0d-9
       drrdd = 0.0d0

! ground state beta decay
      else if (ic1(j,iat) .eq. 7) then
       fr    = cx(k0)
       dfrdt = 0.0d0
       dfrdd = 0.0d0


! ground state positron decay or electron capture
      else if (ic1(j,iat) .eq. 8) then
       rr    = cx(k0)
       drrdt = 0.0d0
       drrdd = 0.0d0


! rates from f.k.thielemann reaclib deck
! exoergic frwd (14), rev (15), beta- (16), and beta+ or ec (17)
! do not evaluate if t9.lt.0.01
      else if (ic1(j,iat) .eq. 14  .or. &
               ic1(j,iat) .eq. 15  .or. &
               ic1(j,iat) .eq. 16  .or. &
               ic1(j,iat) .eq. 17 ) then


       if (t9 .ge. 0.01) then
        nrate = (ic2(j,iat)-2)/7
        lk0 = k
        do jk=1,nrate
         lk1 = lk0 + 1
         lk2 = lk1 + 1
         lk3 = lk2 + 1
         lk4 = lk3 + 1
         lk5 = lk4 + 1
         lk6 = lk5 + 1
         lk7 = lk6 + 1
         lk8 = lk7 + 1

         aa   =  cx(lk0) + cx(lk1)*t9i + cx(lk2)*t9i13 &
               + cx(lk3)*t913 + cx(lk4)*t9 + cx(lk5)*t953 &
               + cx(lk6)*log(t9)

         daa  =   -cx(lk1)*t9i2 - oneth*cx(lk2)*t9i43 &
               + oneth*cx(lk3)*t9i23 + cx(lk4) + fiveth*cx(lk5)*t923 &
               + cx(lk6)*t9i

         if (aa .lt. 200.0) then
          term    = exp(aa)
          dtermdt = term*daa*1.0d-9
         else
          term    = exp(200.0d0)
          dtermdt = 0.0d0
         end if

         if (ic1(j,iat).eq.14) then
          fr  = fr + term
          dfrdt = dfrdt + dtermdt
         end if

         if (ic1(j,iat).eq.15) then
          rr = rr + term
          drrdt = drrdt + dtermdt
         end if
         if (ic1(j,iat).eq.16) then
          fr = fr + term
          dfrdt = dfrdt + dtermdt
         end if
         if (ic1(j,iat).eq.17) then
          rr  = rr + term
          drrdt = drrdt + dtermdt
         end if

         lk0 = lk0+7
        enddo


! no rev rate for weak
        if (ic1(j,iat) .ge. 16 .or. ic1(j,iat) .eq. 17) return

        if (j.eq.2 .or. j.eq.5 .or. j.eq.6) then

         if (ic1(j,iat).eq.14) then
          z   = 1.0d0/zwork1(nrr(jn,iat))
          aa  = zwork1(iat)*z
          daa = (zwork2(iat) - aa*zwork2(nrr(jn,iat)))*z
          bb  = cx(lk7) * exp(-cx(lk8)*t9i)
          dbb = bb*cx(lk8)*t9i2
          cc  = aa*bb
          dcc = (daa*bb + aa*dbb) * 1.0d-9

! forward rate
          dfrdd = fr
          fr    = fr*bden
          dfrdt = dfrdt*bden

! reverse rate
          rr    = fr*cc
          drrdt = dfrdt*cc + fr*dcc
          drrdd = dfrdd*cc

         else
          z   = 1.0d0/zwork1(iat)
          aa  = zwork1(nrr(jn,iat))*z
          daa = (zwork2(nrr(jn,iat)) - aa*zwork2(iat))*z
          bb  = cx(lk7) * exp(-cx(lk8)*t9i)
          dbb = bb*cx(lk8)*t9i2
          cc  = aa*bb
          dcc = (daa*bb + aa*dbb) * 1.0d-9

! reverse rate
          drrdd = rr
          rr    = rr *bden
          drrdt = drrdt*bden

! forward rate
          fr    = rr*cc
          dfrdt = drrdt*cc + rr*dcc
          drrdd = drrdd*cc
         endif

        else

         if (ic1(j,iat).eq.14) then
          aa  = zwork1(iat)/zwork1(nrr(jn,iat))
          daa = (zwork2(iat) - aa*zwork2(nrr(jn,iat))) &
                  / zwork1(nrr(jn,iat))
          bb  = cx(lk7) * t932 * exp(-cx(lk8)*t9i)
          dbb = 1.5d0*bb*t9i + bb*cx(lk8)*t9i2
          cc  = aa*bb
          dcc = (daa*bb + aa*dbb) * 1.0d-9

! reverse rate
          rr    = fr*cc
          drrdt = dfrdt*cc + fr*dcc
          drrdd = 0.0d0

! forward rate
          dfrdd = fr
          fr    = fr*bden
          dfrdt = dfrdt*bden
         else
          aa  = zwork1(nrr(jn,iat))/zwork1(iat)
          daa = (zwork2(nrr(jn,iat)) - aa*zwork2(iat))/zwork1(iat)
          bb  = cx(lk7) * t932 * exp(-cx(lk8)*t9i)
          dbb = 1.5d0*bb*t9i + bb*cx(lk8)*t9i2
          cc  = aa*bb
          dcc = (daa*bb + aa*dbb) * 1.0d-9

! forward rate
          fr    = rr *cc
          dfrdt = drrdt*cc + rr*dcc
          drrdd = 0.0d0

! reverse rate
          drrdd = rr
          rr    = rr*bden
          drrdt = drrdt*bden
         end if
        end if
       endif


! particle capture rates (n,g), (p,g), (a,g) from rath 2000,
! these are treated specially because they have negative Q values

      else if (ic1(j,iat) .eq. 18) then
       if (t9 .ge. 0.01) then
        lk0 = k
        lk1 = lk0+1
        lk2 = lk1+1
        lk3 = lk2+1
        lk4 = lk3+1
        lk5 = lk4+1
        lk6 = lk5+1
        lk7 = lk6+1
        lk8 = lk7+1

! forward (x,g) rate
        aa  =   cx(lk0) + cx(lk1)*t9i + cx(lk2)*t9i13 &
              + cx(lk3)*t913 + cx(lk4)*t9 + cx(lk5)*t953 &
              + cx(lk6)*log(t9)

        daa =   -cx(lk1)*t9i2 - oneth*cx(lk2)*t9i43 &
              + oneth*cx(lk3)*t9i23 + cx(lk4) + fiveth*cx(lk5)*t923 &
              + cx(lk6)*t9i

        bb  = exp(aa)
        dbb = bb*daa

        fr    = bden * bb
        dfrdt = bden * dbb * 1.0d-9
        dfrdd = bb

! reverse (g,x) rate
        z   = 1.0d0/zwork1(nrr(jn,iat))
        aa  = zwork1(iat)*z
        daa = (zwork2(iat) - aa*zwork2(nrr(jn,iat)))*z

        bb  =   cx(lk7) + cx(lk8)*t9i + cx(lk2)*t9i13 &
              + cx(lk3)*t913 + cx(lk4)*t9 + cx(lk5)*t953 &
              +(cx(lk6)+1.5)*log(t9)

        dbb = -cx(lk8)*t9i2 - oneth*cx(lk2)*t9i43 &
             + oneth*cx(lk3)*t9i23 + cx(lk4) + fiveth*cx(lk5)*t923 &
             + (cx(lk6)+1.5)*t9i

        cc = exp(bb)
        dcc = cc*dbb

        rr    = cc*aa
        drrdt = (dcc*aa + cc*daa)*1.0d-9
        drrdd = 0.0d0
       endif


! g.s. 1 parameter alpha decay from nwc
! apply by adding to (g,a) reaction channel.
! The rate for alpha-decay is stored in nucleus
! card deck (z-2,a-4) in bdat. Ex: te106ad found
! in sn102 card deck just after sn102(a,g)

! just adding a constant, no changes to the derivatives

      else if (ic1(j,iat) .eq. 19) then
       fr = fr + cx(k0)


! ground state plus first excited state beta- decay.
! excited state assumed in thermal equilibrium.

      else if (ic1(j,iat) .eq. 20) then

       aa  = cx(k2) * cx(k3) * exp(-11.60485d0*cx(k4)*t9i)
       daa = aa*11.60485d0*cx(k4)*t9i2

       bb  = cx(k0)*cx(k1) + aa
       dbb = daa

       cc  = cx(k2)*exp(-11.60485d0*cx(k4)*t9i)
       dcc = cc*11.60485d0*cx(k4)*t9i2

       dd  = cx(k0) + cc
       ddd = dcc

       z   = 1.0d0/dd
       ee  = bb*z
       dee = (dbb - ee*ddd)*z

       fr   = ee
       dfrdt = dee * 1.0d-9
       dfrdd = 0.0d0

! or an unknown reaction
      else
       write(6,*) 'unknown ic1(j,iat)', j, iat, ic1(j,iat)
       stop 'unknown ic1(j,iat) in routine torchrat'
      end if

      return
      end subroutine evaluate_rate

      subroutine screen_torch(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine computes the screening factors and applies them to the raw reaction rates,
! producing the final reaction rates used by the right hand sides and jacobian matrix elements.

! notionally one looks up the raw forward reaction rate, screens the forward rate 
! and then computes the (screened) reverse rate using detailed balance.

! we enter this routine with the raw forward and reverse rates already calculated
! and in detail balance. two options are avaliable through the switch screen_on:
! 1) asymmetric screening - this is the traditional route, where things like (g,p) reactiosn are not screened.
! 2) symmetric screening - a screening factor is applied to forward and reverse rates.

! screening is expensive computationally - real powers, exponentials and divisions
! so we are going to go through some gyrations to minimize the number of screening
! factors that must be calculated. its worth the time saving at the cost of bookkeeping.

! declare the pass
      real*8 :: y(*)


! local variables
      integer :: i,j,k,jscr,jscr_check,maxrate_loc(1),minrate_loc(1)
      real*8  :: abar,zbar,ye,z2bar, &
                 scfacp,scfacpdt,scfacpdd,scfaca,scfacadt,scfacadd

      integer, parameter :: nscreen_max = 2*abignet + 40

      real*8  :: zsc1(nscreen_max),asc1(nscreen_max), &
                 zsc2(nscreen_max),asc2(nscreen_max), &
                 scvec(nscreen_max),scvecdt(nscreen_max),scvecdd(nscreen_max)

      real*8 :: scpp,scppdt,scppdd,sc1a,sc1adt,sc1add, &
                scdp,scdpdt,scdpdd,scdd,scdddt,scdddd, & 
                sc2a,sc2adt,sc2add,sc3a,sc3adt,sc3add, &
                sc1212,sc1212dt,sc1212dd, &
                sc1216,sc1216dt,sc1216dd, &
                sc1616,sc1616dt,sc1616dd, &
                scdt,scdtdt,scdtdd,sctp,sctpdt,sctpdd, &
                sctt,scttdt,scttdd,sc33,sc33dt,sc33dd, &
                sche3d,sche3ddt,sche3ddd, &
                sche3t,sche3tdt,sche3tdd, &
                sche4p,sche4pdt,sche4pdd, &
                scli7d,scli7ddt,scli7ddd, &
                scli7t,scli7tdt,scli7tdd, &
                scli7he3,scli7he3dt,scli7he3dd, &
                scbe7d,scbe7ddt,scbe7ddd, &
                scbe7t,scbe7tdt,scbe7tdd, &
                scbe7he3,scbe7he3dt,scbe7he3dd, &
                sche4d,sche4ddt,sche4ddd,sche3p,sche3pdt,sche3pdd


! initialize the screening factors with no-screening values

      scvec(ionbeg:ionend)   = 1.0d0
      scvecdt(ionbeg:ionend) = 0.0d0
      scvecdd(ionbeg:ionend) = 0.0d0

      scfacp   = 1.0d0  ; scfacpdt   = 0.0d0   ; scfacpdd = 0.0d0
      scfaca   = 1.0d0  ; scfacadt   = 0.0d0   ; scfacadd = 0.0d0
      scpp     = 1.0d0  ; scppdt     = 0.0d0   ; scppdd   = 0.0d0
      scdp     = 1.0d0  ; scdpdt     = 0.0d0   ; scdpdd   = 0.0d0
      scdd     = 1.0d0  ; scdddt     = 0.0d0   ; scdddd   = 0.0d0
      sc1a     = 1.0d0  ; sc1adt     = 0.0d0   ; sc1add   = 0.0d0
      sc2a     = 1.0d0  ; sc2adt     = 0.0d0   ; sc2add   = 0.0d0
      sc3a     = 1.0d0  ; sc3adt     = 0.0d0   ; sc3add   = 0.0d0
      sc1212   = 1.0d0  ; sc1212dt   = 0.0d0   ; sc1212dd = 0.0d0
      sc1216   = 1.0d0  ; sc1216dt   = 0.0d0   ; sc1216dd = 0.0d0
      sc1616   = 1.0d0  ; sc1616dt   = 0.0d0   ; sc1616dd = 0.0d0
      scdt     = 1.0d0  ; scdtdt     = 0.0d0   ; scdtdd   = 0.0d0
      sctp     = 1.0d0  ; sctpdt     = 0.0d0   ; sctpdd   = 0.0d0
      sctt     = 1.0d0  ; scttdt     = 0.0d0   ; scttdd   = 0.0d0
      sc33     = 1.0d0  ; sc33dt     = 0.0d0   ; sc33dd   = 0.0d0
      sche3p   = 1.0d0  ; sche3pdt   = 0.0d0   ; sche3pdd = 0.0d0
      sche3d   = 1.0d0  ; sche3ddt   = 0.0d0   ; sche3ddd = 0.0d0
      sche3t   = 1.0d0  ; sche3tdt   = 0.0d0   ; sche3tdd = 0.0d0
      sche4d   = 1.0d0  ; sche4ddt   = 0.0d0   ; sche4ddd = 0.0d0
      sche4p   = 1.0d0  ; sche4pdt   = 0.0d0   ; sche4pdd = 0.0d0
      scli7d   = 1.0d0  ; scli7ddt   = 0.0d0   ; scli7ddd = 0.0d0
      scli7t   = 1.0d0  ; scli7tdt   = 0.0d0   ; scli7tdd = 0.0d0
      scli7he3 = 1.0d0  ; scli7he3dt = 0.0d0 ; scli7he3dd = 0.0d0
      scbe7d   = 1.0d0  ; scbe7ddt   = 0.0d0   ; scbe7ddd = 0.0d0
      scbe7t   = 1.0d0  ; scbe7tdt   = 0.0d0   ; scbe7tdd = 0.0d0
      scbe7he3 = 1.0d0  ; scbe7he3dt = 0.0d0 ; scbe7he3dd = 0.0d0



! if screening corrections are on

      if (screen_on .eq. asymmetric_screening .or. screen_on .eq. symmetric_screening) then


! load the screening pipeline. needs to be done once for a given network,
! but do not reset the init_screen flag as that will be done in routine screen6 below

       if (init_screen .eq. 1) then
        jscr = 0

        do j=ionbeg,ionend
         jscr = jscr + 1
         zsc1(jscr) = zion(j)  ; asc1(jscr) = aion(j)  ; zsc2(jscr) = zion(iprot) ; asc2(jscr) = aion(iprot)

         jscr = jscr + 1
         zsc1(jscr) = zion(j)  ; asc1(jscr) = aion(j)  ; zsc2(jscr) = zion(ihe4)  ; asc2(jscr) = aion(ihe4)
        enddo


! now the specials

! triple alpha
        if (ic12 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe4) ; asc1(jscr) = aion(ihe4) ; zsc2(jscr) = zion(ihe4) ; asc2(jscr) = aion(ihe4)
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe4) ; asc1(jscr) = aion(ihe4) ; zsc2(jscr) = 4.0d0 ; asc2(jscr) = 8.0d0
        

! c12 + c12
         if (ine20 .ne. 0  .and. ina23 .ne. 0  .and.  img23 .ne. 0) then
          jscr = jscr + 1
          zsc1(jscr) = zion(ic12)  ; asc1(jscr) = aion(ic12)  ; zsc2(jscr) = zion(ic12)  ; asc2(jscr) = aion(ic12)
         end if
        end if


! o16 + o16
      if (io16 .ne. 0 .and. isi28 .ne. 0 .and. ip30 .ne. 0 .and. ip31 .ne. 0 .and. is31 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(io16) ; asc1(jscr) = aion(io16) ; zsc2(jscr) = zion(io16) ; asc2(jscr) = aion(io16)
        end if

! c12 + o16
        if (ic12 .ne. 0 .and. io16 .ne. 0 .and. img24 .ne. 0  .and. ial27 .ne. 0  .and. isi27 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ic12)  ; asc1(jscr) = aion(ic12) ; zsc2(jscr) = zion(io16) ; asc2(jscr) = aion(io16)
        end if


! deuterium
        if (ih2 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(iprot) ; asc1(jscr) = aion(iprot) ; zsc2(jscr) = zion(iprot) ; asc2(jscr) = aion(iprot)

         jscr = jscr + 1
         zsc1(jscr) = zion(ih2) ; asc1(jscr) = aion(ih2) ; zsc2(jscr) = zion(iprot) ; asc2(jscr) = aion(iprot)

         jscr = jscr + 1
         zsc1(jscr) = zion(ih2) ; asc1(jscr) = aion(ih2) ; zsc2(jscr) = zion(ih2) ; asc2(jscr) = aion(ih2)
        end if


! d + t
        if (ih2 .ne. 0  .and. ih3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ih2)  ; asc1(jscr) = aion(ih2) ; zsc2(jscr) = zion(ih3) ; asc2(jscr) = aion(ih3)
        end if


! h3
        if (ih3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ih3) ; asc1(jscr) = aion(ih3) ; zsc2(jscr) = zion(iprot) ; asc2(jscr) = aion(iprot)

         jscr = jscr + 1
         zsc1(jscr) = zion(ih3) ; asc1(jscr) = aion(ih3) ; zsc2(jscr) = zion(ih3) ; asc2(jscr) = aion(ih3)
        end if

! he3
        if (ihe3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe3) ; asc1(jscr) = aion(ihe3) ; zsc2(jscr) = zion(ihe3) ; asc2(jscr) = aion(ihe3)

         jscr = jscr + 1
         zsc1(jscr) = zion(ihe3)  ; asc1(jscr) = aion(ihe3)  ; zsc2(jscr) = zion(iprot)  ; asc2(jscr) = aion(iprot)
        end if

! he3 + d
        if (ihe3 .ne. 0 .and. ih2 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe3) ; asc1(jscr) = aion(ihe3) ; zsc2(jscr) = zion(ih2) ; asc2(jscr) = aion(ih2)
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe4) ; asc1(jscr) = aion(ihe4) ; zsc2(jscr) = zion(iprot) ; asc2(jscr) = aion(iprot)
        end if

! he3 + t
        if (ihe3 .ne. 0 .and. ih3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe3) ; asc1(jscr) = aion(ihe3) ; zsc2(jscr) = zion(ih3) ; asc2(jscr) = aion(ih3)
         jscr = jscr + 1
         zsc1(jscr) = zion(ihe4) ; asc1(jscr) = aion(ihe4) ; zsc2(jscr) = zion(ih2) ; asc2(jscr) = aion(ih2)
        end if

! li7 + d
        if (ili7 .ne. 0 .and. ih2 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ili7) ; asc1(jscr) = aion(ili7) ; zsc2(jscr) = zion(ih2) ; asc2(jscr) = aion(ih2)
        end if

! li7 + t
        if (ili7 .ne. 0 .and. ih3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ili7) ; asc1(jscr) = aion(ili7) ; zsc2(jscr) = zion(ih3) ; asc2(jscr) = aion(ih3)
        end if


! li7 + he3
        if (ili7 .ne. 0 .and. ihe3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ili7) ; asc1(jscr) = aion(ili7) ; zsc2(jscr) = zion(ihe3) ; asc2(jscr) = aion(ihe3)
        end if

! be7 + d
        if (ibe7 .ne. 0 .and. ih2 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ibe7) ; asc1(jscr) = aion(ibe7) ; zsc2(jscr) = zion(ih2) ; asc2(jscr) = aion(ih2)
        end if

! be7 + t
        if (ibe7 .ne. 0 .and. ih3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ibe7) ; asc1(jscr) = aion(ibe7) ; zsc2(jscr) = zion(ih3) ; asc2(jscr) = aion(ih3)
        end if

! be7 + he3
        if (ibe7 .ne. 0 .and. ihe3 .ne. 0) then
         jscr = jscr + 1
         zsc1(jscr) = zion(ibe7) ; asc1(jscr) = aion(ibe7) ; zsc2(jscr) = zion(ihe3) ; asc2(jscr) = aion(ihe3)
        end if

! make sure we are within the bounds allocated
        if (jscr .gt. nscreen_max) then
         write(6,*) 
         write(6,*) 'jscr=',jscr,' nscreen_max=',nscreen_max
         write(6,*) 
         stop 'jscr > nscreen_max in routine screen_torch'
        end if

        jscr_check = jscr

! end of init_screen block
       end if




! with the passed composition compute abar, zbar and z**2 bar

       call azbar_simple(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)
       z2bar  = abar * sum( zion(1:ionmax) * zion(1:ionmax) * y(1:ionmax) )


! get all the screening factors in the pipeline

       call screen6(jscr, &
                    btemp,bden,zbar,abar,z2bar, &
                    zsc1,asc1,zsc2,asc2, &
                    scvec,scvecdt,scvecdd)


! end of screen_on if
      end if




! apply the screening factors to the sigarray 
      i    = 0
      jscr = 0
      do j=ionbeg,ionend


       if (screen_on .eq. asymmetric_screening .or. screen_on .eq. symmetric_screening) then
        jscr     = jscr + 1
        scfacp   = scvec(jscr)  ; scfacpdt = scvecdt(jscr)  ; scfacpdd = scvecdd(jscr)

        jscr     = jscr + 1
        scfaca   = scvec(jscr)  ; scfacadt = scvecdt(jscr)  ; scfacadd = scvecdd(jscr)
       end if



! (n,g) and (g,n) reactions
! no screening here
       k = nrr(1,j)
       if (k .gt. 0) then
        i            = i + 1
        sig(1,j)     = sigraw(1,j)  ; sigdt(1,j)   = sigrawdt(1,j)  ; sigdd(1,j)   = sigrawdd(1,j)
        ratraw(i)    = sigraw(1,j)  ; dratrawdt(i) = sigrawdt(1,j)  ; dratrawdd(i) = sigrawdd(1,j)
        ratdum(i)    = sig(1,j)     ; dratdumdt(i) = sigdt(1,j)     ; dratdumdd(i) = sigdd(1,j)
        scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdd(i)  = 0.0d0

        i            = i + 1
        sig(2,j)     = sigraw(2,j)  ; sigdt(2,j)   = sigrawdt(2,j)  ; sigdd(2,j)   = sigrawdd(2,j)
        ratraw(i)    = sigraw(2,j)  ; dratrawdt(i) = sigrawdt(2,j)  ; dratrawdd(i) = sigrawdd(2,j)
        ratdum(i)    = sig(2,j)     ; dratdumdt(i) = sigdt(2,j)     ; dratdumdd(i) = sigdd(2,j)
        scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0
       end if



! set up the (p,n) beta- beta+ decay components

       k = nrr(2,j)
       if (k .gt. 0) then

        i          = i + 1
        sig(3,j)   = sigraw(3,j) * scfacp
        sigdt(3,j) = sigrawdt(3,j)*scfacp + sigraw(3,j)*scfacpdt
        sigdd(3,j) = sigrawdd(3,j)*scfacp + sigraw(3,j)*scfacpdd

        ratraw(i)   = sigraw(3,j)  ; dratrawdt(i) = sigrawdt(3,j)  ; dratrawdd(i) = sigrawdd(3,j)
        ratdum(i)   = sig(3,j)     ; dratdumdt(i) = sigdt(3,j)     ; dratdumdd(i) = sigdd(3,j)
        scfac(i)    = scfacp       ; dscfacdt(i)  = scfacpdt       ; dscfacdd(i)  = scfacpdd

        i            = i + 1
        if (screen_on .eq. asymmetric_screening .or. screen_on .eq. no_screening) then
         sig(4,j)     = sigraw(4,j) 
         sigdt(4,j)   = sigrawdt(4,j)
         sigdd(4,j)   = sigrawdd(4,j)
         scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0
        else if (screen_on .eq. symmetric_screening) then
         sig(4,j)     = sigraw(4,j) * scfacp
         sigdt(4,j)   = sigrawdt(4,j)*scfacp + sigraw(4,j)*scfacpdt
         sigdd(4,j)   = sigrawdd(4,j)*scfacp + sigraw(4,j)*scfacpdd
         scfac(i)     = scfacp       ; dscfacdt(i)  = scfacpdt       ; dscfacdd(i)  = scfacpdd
        end if
        ratraw(i)    = sigraw(4,j)  ; dratrawdt(i) = sigrawdt(4,j)  ; dratrawdd(i) = sigrawdd(4,j)
        ratdum(i)    = sig(4,j)     ; dratdumdt(i) = sigdt(4,j)     ; dratdumdd(i) = sigdd(4,j)


        i            = i + 1
        sig(5,j)     = sigraw(5,j)  ; sigdt(5,j)   = sigrawdt(5,j)  ; sigdd(5,j)   = sigrawdd(5,j)
        ratraw(i)    = sigraw(5,j)  ; dratrawdt(i) = sigrawdt(5,j)  ; dratrawdd(i) = sigrawdd(5,j)
        ratdum(i)    = sig(5,j)     ; dratdumdt(i) = sigdt(5,j)     ; dratdumdd(i) = sigdd(5,j)
        scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0

        i            = i + 1
        sig(6,j)     = sigraw(6,j)  ; sigdt(6,j)   = sigrawdt(6,j)  ; sigdd(6,j)   = sigrawdd(6,j)
        ratraw(i)    = sigraw(6,j)  ; dratrawdt(i) = sigrawdt(6,j)  ; dratrawdd(i) = sigrawdd(6,j)
        ratdum(i)    = sig(6,j)     ; dratdumdt(i) = sigdt(6,j)     ; dratdumdd(i) = sigdd(6,j)
        scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0
       end if



! (p,g) and (g,p) reactions

       k = nrr(3,j)
       if (k .gt. 0) then

        i          = i + 1
        sig(7,j)   = sigraw(7,j) * scfacp
        sigdt(7,j) = sigrawdt(7,j)*scfacp + sigraw(7,j)*scfacpdt
        sigdd(7,j) = sigrawdd(7,j)*scfacp + sigraw(7,j)*scfacpdd

        ratraw(i)   = sigraw(7,j)  ; dratrawdt(i) = sigrawdt(7,j)  ; dratrawdd(i) = sigrawdd(7,j)
        ratdum(i)   = sig(7,j)     ; dratdumdt(i) = sigdt(7,j)     ; dratdumdd(i) = sigdd(7,j)
        scfac(i)    = scfacp       ; dscfacdt(i)  = scfacpdt       ; dscfacdd(i)  = scfacpdd

        i            = i + 1
        if (screen_on .eq. asymmetric_screening .or. screen_on .eq. no_screening) then
         sig(8,j)     = sigraw(8,j) 
         sigdt(8,j)   = sigrawdt(8,j)
         sigdd(8,j)   = sigrawdd(8,j)
         scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0
        else if (screen_on .eq. symmetric_screening) then
         sig(8,j)     = sigraw(8,j) * scfacp
         sigdt(8,j)   = sigrawdt(8,j)*scfacp + sigraw(8,j)*scfacpdt
         sigdd(8,j)   = sigrawdd(8,j)*scfacp + sigraw(8,j)*scfacpdd
         scfac(i)     = scfacp       ; dscfacdt(i)  = scfacpdt       ; dscfacdd(i)  = scfacpdd
        end if
        ratraw(i)    = sigraw(8,j)  ; dratrawdt(i) = sigrawdt(8,j)  ; dratrawdd(i) = sigrawdd(8,j)
        ratdum(i)    = sig(8,j)     ; dratdumdt(i) = sigdt(8,j)     ; dratdumdd(i) = sigdd(8,j)
       end if



! (a,p) and (p,a) reactions
! this case is ambiguous. for now we will take (p,a) as the forward direction 

       k = nrr(4,j)
       if (k .gt. 0) then

        i           = i + 1
        if (screen_on .eq. asymmetric_screening .or. screen_on .eq. no_screening) then
         sig(9,j)   = sigraw(9,j) * scfaca
         sigdt(9,j) = sigrawdt(9,j)*scfaca + sigraw(9,j)*scfacadt
         sigdd(9,j) = sigrawdd(9,j)*scfaca + sigraw(9,j)*scfacadd
         scfac(i)    = scfaca       ; dscfacdt(i)  = scfacadt       ; dscfacdd(i)  = scfacadd
        else if (screen_on .eq. symmetric_screening) then
         sig(9,j)   = sigraw(9,j) * scfacp
         sigdt(9,j) = sigrawdt(9,j)*scfacp + sigraw(9,j)*scfacpdt
         sigdd(9,j) = sigrawdd(9,j)*scfacp + sigraw(9,j)*scfacpdd
         scfac(i)   = scfacp       ; dscfacdt(i)  = scfacpdt       ; dscfacdd(i)  = scfacpdd
        end if
        ratraw(i)    = sigraw(9,j)  ; dratrawdt(i) = sigrawdt(9,j)  ; dratrawdd(i) = sigrawdd(9,j)
        ratdum(i)    = sig(9,j)     ; dratdumdt(i) = sigdt(9,j)     ; dratdumdd(i) = sigdd(9,j)

        i           = i + 1
        sig(10,j)   = sigraw(10,j) * scfacp
        sigdt(10,j) = sigrawdt(10,j)*scfacp + sigraw(10,j)*scfacpdt
        sigdd(10,j) = sigrawdd(10,j)*scfacp + sigraw(10,j)*scfacpdd

        ratraw(i)    = sigraw(10,j)  ; dratrawdt(i) = sigrawdt(10,j)  ; dratrawdd(i) = sigrawdd(10,j)
        ratdum(i)    = sig(10,j)     ; dratdumdt(i) = sigdt(10,j)     ; dratdumdd(i) = sigdd(10,j)
        scfac(i)     = scfacp        ; dscfacdt(i)  = scfacpdt        ; dscfacdd(i)  = scfacpdd
       end if


! (a,n) and (n,a) reactions

       k = nrr(5,j)
       if (k .gt. 0) then

        i           = i + 1
        sig(11,j)   = sigraw(11,j) * scfaca
        sigdt(11,j) = sigrawdt(11,j)*scfaca + sigraw(11,j)*scfacadt
        sigdd(11,j) = sigrawdd(11,j)*scfaca + sigraw(11,j)*scfacadd

        ratraw(i)   = sigraw(11,j)  ; dratrawdt(i) = sigrawdt(11,j)  ; dratrawdd(i) = sigrawdd(11,j)
        ratdum(i)   = sig(11,j)     ; dratdumdt(i) = sigdt(11,j)     ; dratdumdd(i) = sigdd(11,j)
        scfac(i)    = scfaca        ; dscfacdt(i)  = scfacadt        ; dscfacdd(i)  = scfacadd


        i           = i + 1
        if (screen_on .eq. asymmetric_screening .or. screen_on .eq. no_screening) then
         sig(12,j)   = sigraw(12,j) 
         sigdt(12,j) = sigrawdt(12,j)
         sigdd(12,j) = sigrawdd(12,j)
         scfac(i)     = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0
        else if (screen_on .eq. symmetric_screening) then
         sig(12,j)   = sigraw(12,j) * scfaca
         sigdt(12,j) = sigrawdt(12,j)*scfaca + sigraw(12,j)*scfacadt
         sigdd(12,j) = sigrawdd(12,j)*scfaca + sigraw(12,j)*scfacadd
         scfac(i)     = scfaca        ; dscfacdt(i)  = scfacadt        ; dscfacdd(i)  = scfacadd
        end if
        ratraw(i)    = sigraw(12,j)  ; dratrawdt(i) = sigrawdt(12,j)  ; dratrawdd(i) = sigrawdd(12,j)
        ratdum(i)    = sig(12,j)     ; dratdumdt(i) = sigdt(12,j)     ; dratdumdd(i) = sigdd(12,j)
       end if


! (a,g) and (g,a) reactions

       k = nrr(6,j)
       if (k .gt. 0) then

        i           = i + 1
        sig(13,j)   = sigraw(13,j) * scfaca
        sigdt(13,j) = sigrawdt(13,j)*scfaca + sigraw(13,j)*scfacadt
        sigdd(13,j) = sigrawdd(13,j)*scfaca + sigraw(13,j)*scfacadd

        ratraw(i)   = sigraw(13,j)  ; dratrawdt(i) = sigrawdt(13,j)  ; dratrawdd(i) = sigrawdd(13,j)
        ratdum(i)   = sig(13,j)     ; dratdumdt(i) = sigdt(13,j)     ; dratdumdd(i) = sigdd(13,j)
        scfac(i)    = scfaca        ; dscfacdt(i)  = scfacadt        ; dscfacdd(i)  = scfacadd

        i           = i + 1
        if (screen_on .eq. asymmetric_screening .or. screen_on .eq. no_screening) then
         sig(14,j)   = sigraw(14,j) 
         sigdt(14,j) = sigrawdt(14,j)
         sigdd(14,j) = sigrawdd(14,j)
         scfac(i)    = 1.0d0        ; dscfacdt(i)  = 0.0d0          ; dscfacdt(i)  = 0.0d0
        else if (screen_on .eq. symmetric_screening) then
         sig(14,j)   = sigraw(14,j) * scfaca
         sigdt(14,j) = sigrawdt(14,j)*scfaca + sigraw(14,j)*scfacadt
         sigdd(14,j) = sigrawdd(14,j)*scfaca + sigraw(14,j)*scfacadd
         scfac(i)    = scfaca        ; dscfacdt(i)  = scfacadt        ; dscfacdd(i)  = scfacadd
        end if
        ratraw(i)    = sigraw(14,j)  ; dratrawdt(i) = sigrawdt(14,j)  ; dratrawdd(i) = sigrawdd(14,j)
        ratdum(i)    = sig(14,j)     ; dratdumdt(i) = sigdt(14,j)     ; dratdumdd(i) = sigdd(14,j)
       end if

! end of reaction type loop
      enddo




! if screening corrections are on, load the rest of the screening vector 
! the into the appropriate scalars

      if (screen_on .eq. asymmetric_screening .or. screen_on .eq. symmetric_screening) then

! triple alpha
       if (ic12 .ne. 0) then
        jscr   = jscr + 1  ; sc1a   = scvec(jscr)  ; sc1adt = scvecdt(jscr)  ; sc1add = scvecdd(jscr)
        jscr   = jscr + 1  ; sc2a   = scvec(jscr)  ; sc2adt = scvecdt(jscr)  ; sc2add = scvecdd(jscr)  ;
        sc3a   = sc1a * sc2a  ; sc3adt = sc1adt*sc2a + sc1a*sc2adt  ; sc3add = sc1add*sc2a + sc1a*sc2add

! c12 + c12
        if (ine20 .ne. 0  .and. ina23 .ne. 0  .and.  img23 .ne. 0) then
         jscr  = jscr + 1  ; sc1212   = scvec(jscr)  ; sc1212dt = scvecdt(jscr)  ; sc1212dd = scvecdd(jscr)
        end if
       end if

! o16 + o16
      if (io16 .ne. 0 .and. isi28 .ne. 0 .and. ip30 .ne. 0 .and. ip31 .ne. 0 .and. is31 .ne. 0) then
        jscr  = jscr + 1  ; sc1616   = scvec(jscr)  ; sc1616dt = scvecdt(jscr)  ; sc1616dd = scvecdd(jscr)
       end if

! c12 + o16
       if (ic12 .ne. 0 .and. io16 .ne. 0 .and. img24 .ne. 0  .and. ial27 .ne. 0  .and. isi27 .ne. 0) then
        jscr  = jscr + 1  ; sc1216   = scvec(jscr)  ; sc1216dt = scvecdt(jscr)  ; sc1216dd = scvecdd(jscr)
       end if

! deuterium
       if (ih2 .ne. 0) then
        jscr   = jscr + 1  ; scpp   = scvec(jscr)  ; scppdt = scvecdt(jscr)  ; scppdd = scvecdd(jscr)
        jscr   = jscr + 1  ; scdp   = scvec(jscr)  ; scdpdt = scvecdt(jscr)  ; scdpdd = scvecdd(jscr)
        jscr   = jscr + 1  ; scdd   = scvec(jscr)  ; scdddt = scvecdt(jscr)  ; scdddd = scvecdd(jscr)
       end if

! d + t
       if (ih2 .ne. 0  .and. ih3 .ne. 0) then
        jscr  = jscr + 1  ;  scdt   = scvec(jscr)  ; scdtdt = scvecdt(jscr)  ; scdtdd = scvecdd(jscr)
       end if

! h3 
       if (ih3 .ne. 0) then
        jscr  = jscr + 1  ; sctp   = scvec(jscr)  ; sctpdt = scvecdt(jscr)  ; sctpdd = scvecdd(jscr)
        jscr  = jscr + 1 ;  sctt    = scvec(jscr)  ; scttdt  = scvecdt(jscr)  ; scttdd  = scvecdd(jscr)
       end if

! he3
       if (ihe3 .ne. 0) then
        jscr  = jscr + 1  ; sc33   = scvec(jscr)  ; sc33dt = scvecdt(jscr)  ; sc33dd = scvecdd(jscr)
        jscr  = jscr + 1 ; sche3p   = scvec(jscr)  ; sche3pdt = scvecdt(jscr)  ; sche3pdd = scvecdd(jscr)
       end if

! he3 + d
       if (ihe3 .ne. 0 .and. ih2 .ne. 0) then
        jscr  = jscr + 1  ; sche3d   = scvec(jscr)  ; sche3ddt = scvecdt(jscr)  ; sche3ddd = scvecdd(jscr)
        jscr  = jscr + 1  ; sche4p   = scvec(jscr)  ; sche4pdt = scvecdt(jscr)  ; sche4pdd = scvecdd(jscr)
       end if

! he3 + t
       if (ihe3 .ne. 0 .and. ih3 .ne. 0) then
        jscr  = jscr + 1  ; sche3t   = scvec(jscr)  ; sche3tdt = scvecdt(jscr)  ; sche3tdd = scvecdd(jscr)
        jscr  = jscr + 1  ; sche4d   = scvec(jscr)  ; sche4ddt = scvecdt(jscr)  ; sche4ddd = scvecdd(jscr)
       end if

! li + d
       if (ili7 .ne. 0 .and. ih2 .ne. 0) then
        jscr  = jscr + 1  ; scli7d   = scvec(jscr)  ; scli7ddt = scvecdt(jscr)  ; scli7ddd = scvecdd(jscr)
       end if

! li7 + t
       if (ili7 .ne. 0 .and. ih3 .ne. 0) then
        jscr  = jscr + 1  ; scli7t   = scvec(jscr)  ; scli7tdt = scvecdt(jscr)  ; scli7tdd = scvecdd(jscr)
       end if

! li7 + he3
       if (ili7 .ne. 0 .and. ihe3 .ne. 0) then
        jscr  = jscr + 1  ; scli7he3   = scvec(jscr)  ; scli7he3dt = scvecdt(jscr)  ; scli7he3dd = scvecdd(jscr)
       end if

! be7 + d
       if (ibe7 .ne. 0 .and. ih2 .ne. 0) then
        jscr  = jscr + 1  ; scbe7d   = scvec(jscr)  ; scbe7ddt = scvecdt(jscr)  ; scbe7ddd = scvecdd(jscr)
       end if

! be7 + t
       if (ibe7 .ne. 0 .and. ih3 .ne. 0) then
        jscr  = jscr + 1  ; scbe7t   = scvec(jscr)  ; scbe7tdt = scvecdt(jscr)  ; scbe7tdd = scvecdd(jscr)
       end if

! be7 + he3
       if (ibe7 .ne. 0 .and. ihe3 .ne. 0) then
        jscr  = jscr + 1  ; scbe7he3   = scvec(jscr)  ; scbe7he3dt = scvecdt(jscr)  ; scbe7he3dd = scvecdd(jscr)
       end if

! make sure we are counted correctedly
       if (jscr .ne. jscr_check) then
        write(6,*) 'jscr=',jscr,' jscr_check=',jscr_check
        stop 'fatal counting error in routine screen_torch'
       end if
      end if



! finish recording the rates and multipliers
! for p(e-,nu)n and n(e+,nub)p reactions

      i  = i + 1
      ratdum(irpen)    = ratraw(irpen)
      dratdumdt(irpen) = dratrawdt(irpen)
      dratdumdd(irpen) = dratrawdd(irpen)
      scfac(irpen)     = 1.0d0
      dscfacdt(irpen)  = 0.0d0
      dscfacdd(irpen)  = 0.0d0

      i  = i + 1
      ratdum(irnep)    = ratraw(irnep)
      dratdumdt(irnep) = dratrawdt(irnep)
      dratdumdd(irnep) = dratrawdd(irnep)
      scfac(irnep)     = 1.0d0
      dscfacdt(irnep)  = 0.0d0
      dscfacdd(irnep)  = 0.0d0


! c12 reactions, first triple alpha
      if (ic12 .ne. 0) then
       i = i + 1
       ratdum(ir3a)    = ratraw(ir3a) * sc3a
       dratdumdt(ir3a) = dratrawdt(ir3a) * sc3a + ratraw(ir3a) * sc3adt
       dratdumdd(ir3a) = dratrawdd(ir3a) * sc3a + ratraw(ir3a) * sc3add
       scfac(ir3a)     = sc3a
       dscfacdt(ir3a)  = sc3adt
       dscfacdd(ir3a)  = sc3add

       i = i + 1
       ratdum(irg3a)    = ratraw(irg3a) * sc3a
       dratdumdt(irg3a) = dratrawdt(irg3a)*sc3a + ratraw(irg3a)*sc3adt
       dratdumdd(irg3a) = dratrawdd(irg3a)*sc3a + ratraw(irg3a)*sc3add
       scfac(irg3a)     = sc3a
       dscfacdt(irg3a)  = sc3adt
       dscfacdd(irg3a)  = sc3add



! c12+c12 reactions; must have ne20, na23, mg23 in the network
       if (ine20 .ne. 0  .and. ina23 .ne. 0  .and.  img23 .ne. 0) then

        i = i + 1
        ratdum(ir1212n)    = ratraw(ir1212n) * sc1212
        dratdumdt(ir1212n) = dratrawdt(ir1212n) * sc1212 + ratraw(ir1212n) * sc1212dt
        dratdumdd(ir1212n) = dratrawdd(ir1212n) * sc1212 + ratraw(ir1212n) * sc1212dd
        scfac(ir1212n)     = sc1212
        dscfacdt(ir1212n)  = sc1212dt
        dscfacdd(ir1212n)  = sc1212dd

        i = i + 1
        ratdum(irmg23nc)    = ratraw(irmg23nc) * sc1212
        dratdumdt(irmg23nc) = dratrawdt(irmg23nc) * sc1212 + ratraw(irmg23nc) * sc1212dt
        dratdumdd(irmg23nc) = dratrawdd(irmg23nc) * sc1212 &
                            + ratraw(irmg23nc) * sc1212dd
        scfac(irmg23nc)     = sc1212
        dscfacdt(irmg23nc)  = sc1212dt
        dscfacdd(irmg23nc)  = sc1212dd

        i = i + 1
        ratdum(ir1212p)    = ratraw(ir1212p) * sc1212
        dratdumdt(ir1212p) = dratrawdt(ir1212p) * sc1212 + ratraw(ir1212p) * sc1212dt
        dratdumdd(ir1212p) = dratrawdd(ir1212p) * sc1212 + ratraw(ir1212p) * sc1212dd
        scfac(ir1212p)     = sc1212
        dscfacdt(ir1212p)  = sc1212dt
        dscfacdd(ir1212p)  = sc1212dd

        i = i + 1
        ratdum(irna23pc)    = ratraw(irna23pc) * sc1212
        dratdumdt(irna23pc) = dratrawdt(irna23pc) * sc1212 + ratraw(irna23pc) * sc1212dt
        dratdumdd(irna23pc) = dratrawdd(irna23pc) * sc1212 + ratraw(irna23pc) * sc1212dd
        scfac(irna23pc)     = sc1212
        dscfacdt(irna23pc)  = sc1212dt
        dscfacdd(irna23pc)  = sc1212dd

        i = i + 1
        ratdum(ir1212a)    = ratraw(ir1212a) * sc1212
        dratdumdt(ir1212a) = dratrawdt(ir1212a) * sc1212 + ratraw(ir1212a) * sc1212dt
        dratdumdd(ir1212a) = dratrawdd(ir1212a) * sc1212 + ratraw(ir1212a) * sc1212dd
        scfac(ir1212a)     = sc1212
        dscfacdt(ir1212a)  = sc1212dt
        dscfacdd(ir1212a)  = sc1212dd

        i = i + 1
        ratdum(irne20ac)    = ratraw(irne20ac) * sc1212
        dratdumdt(irne20ac) = dratrawdt(irne20ac) * sc1212 + ratraw(irne20ac) * sc1212dt
        dratdumdd(irne20ac) = dratrawdd(irne20ac) * sc1212 + ratraw(irne20ac) * sc1212dd
        scfac(irne20ac)     = sc1212
        dscfacdt(irne20ac)  = sc1212dt
        dscfacdd(irne20ac)  = sc1212dd

       end if
      end if



! o16+o16 reactions; must have si28, p 30, p31 and s31 in the network
      if (io16 .ne. 0 .and. isi28 .ne. 0 .and. ip30 .ne. 0 .and. ip31 .ne. 0 .and. is31 .ne. 0) then

       i = i + 1
       ratdum(ir1616n)    = ratraw(ir1616n) * sc1616
       dratdumdt(ir1616n) = dratrawdt(ir1616n) * sc1616 + ratraw(ir1616n) * sc1616dt
       dratdumdd(ir1616n) = dratrawdd(ir1616n) * sc1616 + ratraw(ir1616n) * sc1616dd
       scfac(ir1616n)     = sc1616
       dscfacdt(ir1616n)  = sc1616dt
       dscfacdd(ir1616n)  = sc1616dd

       i = i + 1
       ratdum(irs31no)    = ratraw(irs31no) * sc1616
       dratdumdt(irs31no) = dratrawdt(irs31no) * sc1616 + ratraw(irs31no) * sc1616dt
       dratdumdd(irs31no) = dratrawdd(irs31no) * sc1616 + ratraw(irs31no) * sc1616dd
       scfac(irs31no)     = sc1616
       dscfacdt(irs31no)  = sc1616dt
       dscfacdd(irs31no)  = sc1616dd

       i = i + 1
       ratdum(ir1616p)    = ratraw(ir1616p) * sc1616
       dratdumdt(ir1616p) = dratrawdt(ir1616p) * sc1616 + ratraw(ir1616p) * sc1616dt
       dratdumdd(ir1616p) = dratrawdd(ir1616p) * sc1616 + ratraw(ir1616p) * sc1616dd
       scfac(ir1616p)     = sc1616
       dscfacdt(ir1616p)  = sc1616dt
       dscfacdd(ir1616p)  = sc1616dd

       i = i + 1
       ratdum(irp31po)    = ratraw(irp31po) * sc1616
       dratdumdt(irp31po) = dratrawdt(irp31po) * sc1616 + ratraw(irp31po) * sc1616dt
       dratdumdd(irp31po) = dratrawdd(irp31po) * sc1616 + ratraw(irp31po) * sc1616dd
       scfac(irp31po)     = sc1616
       dscfacdt(irp31po)  = sc1616dt
       dscfacdd(irp31po)  = sc1616dd

       i = i + 1
       ratdum(ir1616a)    = ratraw(ir1616a) * sc1616
       dratdumdt(ir1616a) = dratrawdt(ir1616a) * sc1616 + ratraw(ir1616a) * sc1616dt
       dratdumdd(ir1616a) = dratrawdd(ir1616a) * sc1616 + ratraw(ir1616a) * sc1616dd
       scfac(ir1616a)     = sc1616
       dscfacdt(ir1616a)  = sc1616dt
       dscfacdd(ir1616a)  = sc1616dd

       i = i + 1
       ratdum(irsi28ao)    = ratraw(irsi28ao) * sc1616
       dratdumdt(irsi28ao) = dratrawdt(irsi28ao) * sc1616 + ratraw(irsi28ao) * sc1616dt
       dratdumdd(irsi28ao) = dratrawdd(irsi28ao) * sc1616 + ratraw(irsi28ao) * sc1616dd
       scfac(irsi28ao)     = sc1616
       dscfacdt(irsi28ao)  = sc1616dt
       dscfacdd(irsi28ao)  = sc1616dd

       i = i + 1
       ratdum(ir1616d)    = ratraw(ir1616d) * sc1616
       dratdumdt(ir1616d) = dratrawdt(ir1616d) * sc1616 + ratraw(ir1616d) * sc1616dt
       dratdumdd(ir1616d) = dratrawdd(ir1616d) * sc1616 + ratraw(ir1616d) * sc1616dd
       scfac(ir1616d)     = sc1616
       dscfacdt(ir1616d)  = sc1616dt
       dscfacdd(ir1616d)  = sc1616dd

       i = i + 1
       ratdum(irp30do)    = ratraw(irp30do) * sc1616
       dratdumdt(irp30do) = dratrawdt(irp30do) * sc1616 + ratraw(irp30do) * sc1616dt
       dratdumdd(irp30do) = dratrawdd(irp30do) * sc1616 + ratraw(irp30do) * sc1616dd
       scfac(irp30do)     = sc1616
       dscfacdt(irp30do)  = sc1616dt
       dscfacdd(irp30do)  = sc1616dd

      end if



! c12+o16 reactions; must have mg24, al27, si27 in the network
      if (ic12 .ne. 0 .and. io16 .ne. 0 .and. img24 .ne. 0  .and. ial27 .ne. 0  .and. isi27 .ne. 0) then

       i = i + 1
       ratdum(ir1216n)    = ratraw(ir1216n) * sc1216
       dratdumdt(ir1216n) = dratrawdt(ir1216n) * sc1216 + ratraw(ir1216n) * sc1216dt
       dratdumdd(ir1216n) = dratrawdd(ir1216n) * sc1216 + ratraw(ir1216n) * sc1216dd
       scfac(ir1216n)     = sc1216
       dscfacdt(ir1216n)  = sc1216dt
       dscfacdd(ir1216n)  = sc1216dd

       i = i + 1
       ratdum(irsi27no)     = ratraw(irsi27no) * sc1216
       dratdumdt(irsi27no)  = dratrawdt(irsi27no) * sc1216 + ratraw(irsi27no) * sc1216dt
       dratdumdd(irsi27no)  = dratrawdd(irsi27no) * sc1216 + ratraw(irsi27no) * sc1216dd
       scfac(irsi27no)      = sc1216
       dscfacdt(irsi27no)   = sc1216dt
       dscfacdd(irsi27no)   = sc1216dd

       i = i + 1
       ratdum(ir1216p)     = ratraw(ir1216p) * sc1216 
       dratdumdt(ir1216p)  = dratrawdt(ir1216p) * sc1216 + ratraw(ir1216p) * sc1216dt
       dratdumdd(ir1216p)  = dratrawdd(ir1216p) * sc1216 + ratraw(ir1216p) * sc1216dd
       scfac(ir1216p)      = sc1216
       dscfacdt(ir1216p)   = sc1216dt
       dscfacdd(ir1216p)   = sc1216dd

       i = i + 1
       ratdum(iral27po)    = ratraw(iral27po) * sc1216
       dratdumdt(iral27po) = dratrawdt(iral27po) * sc1216 + ratraw(iral27po) * sc1216dt
       dratdumdd(iral27po) = dratrawdd(iral27po) * sc1216 + ratraw(iral27po) * sc1216dd
       scfac(iral27po)     = sc1216
       dscfacdt(iral27po)  = sc1216dt
       dscfacdd(iral27po)  = sc1216dd

       i = i + 1
       ratdum(ir1216a)     = ratraw(ir1216a) * sc1216
       dratdumdt(ir1216a)  = dratrawdt(ir1216a) * sc1216 + ratraw(ir1216a) * sc1216dt
       dratdumdd(ir1216a)  = dratrawdd(ir1216a) * sc1216 + ratraw(ir1216a) * sc1216dd
       scfac(ir1216a)      = sc1216
       dscfacdt(ir1216a)   = sc1216dt
       dscfacdd(ir1216a)   = sc1216dd

       i = i + 1
       ratdum(irmg24ao)    = ratraw(irmg24ao) * sc1216
       dratdumdt(irmg24ao) = dratrawdt(irmg24ao) * sc1216 + ratraw(irmg24ao) * sc1216dt
       dratdumdd(irmg24ao) = dratrawdd(irmg24ao) * sc1216 + ratraw(irmg24ao) * sc1216dd
       scfac(irmg24ao)     = sc1216
       dscfacdt(irmg24ao)  = sc1216dt
       dscfacdd(irmg24ao)  = sc1216dd

      end if


! if we have deuterium

      if (ih2 .ne. 0) then

! pp
       i = i + 1
       ratdum(irpp)    = ratraw(irpp) * scpp
       dratdumdt(irpp) = dratrawdt(irpp)*scpp + ratraw(irpp)*scppdt
       dratdumdd(irpp) = dratrawdd(irpp)*scpp + ratraw(irpp)*scppdd
       scfac(irpp)     = scpp
       dscfacdt(irpp)  = scppdt
       dscfacdd(irpp)  = scppdd

       i = i + 1
       ratdum(irpep)    = ratraw(irpep) * scpp
       dratdumdt(irpep) = dratrawdt(irpep)*scpp + ratraw(irpep)*scppdt
       dratdumdd(irpep) = dratrawdd(irpep)*scpp + ratraw(irpep)*scppdd
       scfac(irpep)     = scpp
       dscfacdt(irpep)  = scppdt
       dscfacdd(irpep)  = scppdd

! p(n,g)d
       i = i + 1
       ratdum(irpng)     = ratraw(irpng)
       dratdumdt(irpng)  = dratrawdt(irpng)
       dratdumdd(irpng)  = dratrawdd(irpng)
       scfac(irpng)      = 1.0d0
       dscfacdt(irpng)   = 0.0d0
       dscfacdd(irpng)   = 0.0d0

       i = i + 1
       ratdum(irdgn)     = ratraw(irdgn)
       dratdumdt(irdgn)  = dratrawdt(irdgn)
       dratdumdd(irdgn)  = dratrawdd(irdgn)
       scfac(irdgn)      = 1.0d0
       dscfacdt(irdgn)   = 0.0d0
       dscfacdd(irdgn)   = 0.0d0

! d(p,n)2p
       i = i + 1
       ratdum(irdpn)     = ratraw(irdpn) * scdp
       dratdumdt(irdpn)  = dratrawdt(irdpn)*scdp + ratraw(irdpn)*scdpdt
       dratdumdd(irdpn)  = dratrawdd(irdpn)*scdp + ratraw(irdpn)*scdpdd
       scfac(irdpn)      = scdp
       dscfacdt(irdpn)   = scdpdt
       dscfacdd(irdpn)   = scdpdd

       i = i + 1
       ratdum(ir2pnp)     = ratraw(ir2pnp) * scdp
       dratdumdt(ir2pnp)  = dratrawdt(ir2pnp)*scdp + ratraw(ir2pnp)*scdpdt
       dratdumdd(ir2pnp)  = dratrawdd(ir2pnp)*scdp + ratraw(ir2pnp)*scdpdd
       scfac(ir2pnp)      = scdp
       dscfacdt(ir2pnp)   = scdpdt
       dscfacdd(ir2pnp)   = scdpdd

! d(d,g)he4
       i = i + 1
       ratdum(irddg)     = ratraw(irddg) * scdd
       dratdumdt(irddg)  = dratrawdt(irddg)*scdd + ratraw(irddg)*scdddt
       dratdumdd(irddg)  = dratrawdd(irddg)*scdd + ratraw(irddg)*scdddd
       scfac(irddg)      = scdd
       dscfacdt(irddg)   = scdddt
       dscfacdd(irddg)   = scdddd

       i = i + 1
       ratdum(irhe4gd)     = ratraw(irhe4gd) * scdd
       dratdumdt(irhe4gd)  = dratrawdt(irhe4gd)*scdd + ratraw(irhe4gd)*scdddt
       dratdumdd(irhe4gd)  = dratrawdd(irhe4gd)*scdd + ratraw(irhe4gd)*scdddd
       scfac(irhe4gd)      = scdd
       dscfacdt(irhe4gd)   = scdddt
       dscfacdd(irhe4gd)   = scdddd
      end if


! if we have deutrium and tritium

      if (ih2 .ne. 0  .and. ih3 .ne. 0) then

! t(d,n)he4  dt reaction
       i = i + 1
       ratdum(irtdn)    = ratraw(irtdn) * scdt
       dratdumdt(irtdn) = dratrawdt(irtdn)*scdt + ratraw(irtdn)*scdtdt
       dratdumdd(irtdn) = dratrawdd(irtdn)*scdt + ratraw(irtdn)*scdtdd
       scfac(irtdn)     = scdt
       dscfacdt(irtdn)  = scdtdt
       dscfacdd(irtdn)  = scdtdd

       i = i + 1
       ratdum(irhe4nd)    = ratraw(irhe4nd) * scdt
       dratdumdt(irhe4nd) = dratrawdt(irhe4nd)*scdt
       dratdumdd(irhe4nd) = dratrawdd(irhe4nd)*scdt
       scfac(irhe4nd)     = scdt
       dscfacdt(irhe4nd)  = scdtdt
       dscfacdd(irhe4nd)  = scdtdd


! d(d,p)t
       i = i + 1
       ratdum(irddp)    = ratraw(irddp) * scdd
       dratdumdt(irddp) = dratrawdt(irddp)*scdd + ratraw(irddp)*scdddt
       dratdumdd(irddp) = dratrawdd(irddp)*scdd + ratraw(irddp)*scdddd
       scfac(irddp)     = scdd
       dscfacdt(irddp)  = scdddt
       dscfacdd(irddp)  = scdddd

       i = i + 1
       ratdum(irtpd)    = ratraw(irtpd) * scdd
       dratdumdt(irtpd) = dratrawdt(irtpd)*scdd + ratraw(irtpd)*scdddt
       dratdumdd(irtpd) = dratrawdd(irtpd)*scdd + ratraw(irtpd)*scdddd
       scfac(irtpd)     = scdd
       dscfacdt(irtpd)  = scdddt
       dscfacdd(irtpd)  = scdddd
      end if


! if we have tritium

      if (ih3 .ne. 0) then

! t(p,g)he4
       i = i + 1
       ratdum(irh3pg)    = ratraw(irh3pg) * sctp
       dratdumdt(irh3pg) = dratrawdt(irh3pg)*sctp + ratraw(irh3pg)*sctpdt
       dratdumdd(irh3pg) = dratrawdd(irh3pg)*sctp + ratraw(irh3pg)*sctpdd
       scfac(irh3pg)     = sctp
       dscfacdt(irh3pg)  = sctpdt
       dscfacdd(irh3pg)  = sctpdd

       i = i + 1
       ratdum(irhe4gp)    = ratraw(irhe4gp) * sctp
       dratdumdt(irhe4gp) = dratrawdt(irhe4gp)*sctp + ratraw(irhe4gp)*sctpdt
       dratdumdd(irhe4gp) = dratrawdd(irhe4gp)*sctp + ratraw(irhe4gp)*sctpdd
       scfac(irhe4gp)     = sctp
       dscfacdt(irhe4gp)  = sctpdt
       dscfacdd(irhe4gp)  = sctpdd

! t(t,2n)he4
       i = i + 1
       ratdum(irtt2n)    = ratraw(irtt2n) * sctt
       dratdumdt(irtt2n) = dratrawdt(irtt2n)*sctt + ratraw(irtt2n)*scttdt
       dratdumdd(irtt2n) = dratrawdd(irtt2n)*sctt + ratraw(irtt2n)*scttdd
       scfac(irtt2n)     = sctt
       dscfacdt(irtt2n)  = scttdt
       dscfacdd(irtt2n)  = scttdd

       i = i + 1
       ratdum(irhe42nt)    = ratraw(irhe42nt) * sctt
       dratdumdt(irhe42nt) = dratrawdt(irhe42nt)*sctt + ratraw(irhe42nt)*scttdt
       dratdumdd(irhe42nt) = dratrawdd(irhe42nt)*sctt + ratraw(irhe42nt)*scttdd
       scfac(irhe42nt)     = sctt
       dscfacdt(irhe42nt)  = scttdt
       dscfacdd(irhe42nt)  = scttdd
      end if



! if we have he3
      if (ihe3 .ne. 0) then


! he3(he3,2p)he4
       i = i + 1
       ratdum(ir33)    = ratraw(ir33) * sc33
       dratdumdt(ir33) = dratrawdt(ir33) * sc33 + ratraw(ir33) * sc33dt
       dratdumdd(ir33) = dratrawdd(ir33) * sc33 + ratraw(ir33) * sc33dd
       scfac(ir33)     = sc33
       dscfacdt(ir33)  = sc33dt
       dscfacdd(ir33)  = sc33dd

       i = i + 1
       ratdum(ir33inv)    = ratraw(ir33inv) * sc33
       dratdumdt(ir33inv) = dratrawdt(ir33inv)*sc33 + ratraw(ir33inv)*sc33dt
       dratdumdd(ir33inv) = dratrawdd(ir33inv)*sc33 + ratraw(ir33inv)*sc33dd
       scfac(ir33inv)     = sc33
       dscfacdt(ir33inv)  = sc33dt
       dscfacdd(ir33inv)  = sc33dd

! he3(p,e-nu)he4
       i = i + 1
       ratdum(irhep)    = ratraw(irhep) * sche3p
       dratdumdt(irhep) = dratrawdt(irhep) * sche3p + ratraw(irhep) * sche3pdt
       dratdumdd(irhep) = dratrawdd(irhep) * sche3p + ratraw(irhep) * sche3pdd
       scfac(irhep)     = sche3p
       dscfacdt(irhep)  = sche3pdt
       dscfacdd(irhep)  = sche3pdd


! he3(n,g)he4
       i = i + 1
       ratdum(irhe3ng)    = ratraw(irhe3ng)
       dratdumdt(irhe3ng) = dratrawdt(irhe3ng)
       dratdumdd(irhe3ng) = dratrawdd(irhe3ng)
       scfac(irhe3ng)     = 1.0d0
       dscfacdt(irhe3ng)  = 0.0d0
       dscfacdd(irhe3ng)  = 0.0d0

       i = i + 1
       ratdum(irhe4gn)    = ratraw(irhe4gn)
       dratdumdt(irhe4gn) = dratrawdt(irhe4gn)
       dratdumdd(irhe4gn) = dratrawdd(irhe4gn)
       scfac(irhe4gn)     = 1.0d0
       dscfacdt(irhe4gn)  = 0.0d0
       dscfacdd(irhe4gn)  = 0.0d0
      end if


! if we have he3 and deuterium

      if (ihe3 .ne. 0 .and. ih2 .ne. 0) then

! he3(d,p)he4
       i = i + 1
       ratdum(irhe3dp)    = ratraw(irhe3dp) * sche3d
       dratdumdt(irhe3dp) = dratrawdt(irhe3dp) * sche3d + ratraw(irhe3dp) * sche3ddt
       dratdumdd(irhe3dp) = dratrawdd(irhe3dp) * sche3d + ratraw(irhe3dp) * sche3ddd
       scfac(irhe3dp)     = sche3d
       dscfacdt(irhe3dp)  = sche3ddt
       dscfacdd(irhe3dp)  = sche3ddd

       i = i + 1
       ratdum(irhe4pd)    = ratraw(irhe4pd) * sche3d
       dratdumdt(irhe4pd) = dratrawdt(irhe4pd) * sche3d + ratraw(irhe4pd) * sche3ddt
       dratdumdd(irhe4pd) = dratrawdd(irhe4pd) * sche3d + ratraw(irhe4pd) * sche3ddd
       scfac(irhe4pd)     = sche3d
       dscfacdt(irhe4pd)  = sche3ddt
       dscfacdd(irhe4pd)  = sche3ddd

! d(d,n)he3
       i = i + 1
       ratdum(irddn)    = ratraw(irddn)
       dratdumdt(irddn) = dratrawdt(irddn)
       dratdumdd(irddn) = dratrawdd(irddn)
       scfac(irddn)     = 1.0d0
       dscfacdt(irddn)  = 0.0d0
       dscfacdd(irddn)  = 0.0d0

       i = i + 1
       ratdum(irhe3nd)    = ratraw(irhe3nd)
       dratdumdt(irhe3nd) = dratrawdt(irhe3nd)
       dratdumdd(irhe3nd) = dratrawdd(irhe3nd)
       scfac(irhe3nd)     = 1.0d0
       dscfacdt(irhe3nd)  = 0.0d0
       dscfacdd(irhe3nd)  = 0.0d0
      end if


! if we have he3, deuterium, tritium

      if (ihe3 .ne. 0  .and. ih2 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,d)he4
       i = i + 1
       ratdum(irhe3td)    = ratraw(irhe3td) * sche3t
       dratdumdt(irhe3td) = dratrawdt(irhe3td) * sche3t + ratraw(irhe3td) * sche3tdt
       dratdumdd(irhe3td) = dratrawdd(irhe3td) * sche3t + ratraw(irhe3td) * sche3tdd
       scfac(irhe3td)     = sche3t
       dscfacdt(irhe3td)  = sche3tdt
       dscfacdd(irhe3td)  = sche3tdd

       i = i + 1
       ratdum(irhe4dt)    = ratraw(irhe4dt) * sche3t
       dratdumdt(irhe4dt) = dratrawdt(irhe4dt) * sche3t + ratraw(irhe4dt) * sche3tdt
       dratdumdd(irhe4dt) = dratrawdd(irhe4dt) * sche3t + ratraw(irhe4dt) * sche3tdd
       scfac(irhe4dt)     = sche3t
       dscfacdt(irhe4dt)  = sche3tdt
       dscfacdd(irhe4dt)  = sche3tdd
      end if


! if we have he3 and tritium

      if (ihe3 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,np)he4
       i = i + 1
       ratdum(irhe3tnp)    = ratraw(irhe3tnp) * sche3t
       dratdumdt(irhe3tnp) = dratrawdt(irhe3tnp) * sche3t + ratraw(irhe3tnp) * sche3tdt
       dratdumdd(irhe3tnp) = dratrawdd(irhe3tnp) * sche3t + ratraw(irhe3tnp) * sche3tdd
       scfac(irhe3tnp)     = sche3t
       dscfacdt(irhe3tnp)  = sche3tdt
       dscfacdd(irhe3tnp)  = sche3tdd

       i = i + 1
       ratdum(irhe4pnt)    = ratraw(irhe4pnt) * sche3t
       dratdumdt(irhe4pnt) = dratrawdt(irhe4pnt) * sche3t + ratraw(irhe4pnt) * sche3tdt
       dratdumdd(irhe4pnt) = dratrawdd(irhe4pnt) * sche3t + ratraw(irhe4pnt) * sche3tdd
       scfac(irhe4pnt)     = sche3t
       dscfacdt(irhe4pnt)  = sche3tdt
       dscfacdd(irhe4pnt)  = sche3tdd
      end if



! if we have li7

      if (ili7 .ne. 0) then

! li7(p,g)be8 and li7(p,a)he4
       i = i + 1
       ratdum(irli7pag)    = ratraw(irli7pag)
       dratdumdt(irli7pag) = dratrawdt(irli7pag)
       dratdumdd(irli7pag) = dratrawdd(irli7pag)
       scfac(irli7pag)     = 1.0d0
       dscfacdt(irli7pag)  = 0.0d0
       dscfacdd(irli7pag)  = 0.0d0

       i = i + 1
       ratdum(ir2he4ga)    = ratraw(ir2he4ga)
       dratdumdt(ir2he4ga) = dratrawdt(ir2he4ga)
       dratdumdd(ir2he4ga) = dratrawdd(ir2he4ga)
       scfac(ir2he4ga)     = 1.0d0
       dscfacdt(ir2he4ga)  = 0.0d0
       dscfacdd(ir2he4ga)  = 0.0d0
      end if

! if we have li7 and deuterium

      if (ili7 .ne. 0 .and. ih2 .ne. 0) then

! li7(d,n)2a
       i = i + 1
       ratdum(irli7dn)    = ratraw(irli7dn) * scli7d
       dratdumdt(irli7dn) = dratrawdt(irli7dn) * scli7d + ratraw(irli7dn) * scli7ddt
       dratdumdd(irli7dn) = dratrawdd(irli7dn) * scli7d + ratraw(irli7dn) * scli7ddd
       scfac(irli7dn)     = scli7d
       dscfacdt(irli7dn)  = scli7ddt
       dscfacdd(irli7dn)  = scli7ddd

       i = i + 1
       ratdum(ir2and)    = ratraw(ir2and) * scli7d
       dratdumdt(ir2and) = dratrawdt(ir2and) * scli7d + ratraw(ir2and) * scli7ddt
       dratdumdd(ir2and) = dratrawdd(ir2and) * scli7d + ratraw(ir2and) * scli7ddd
       scfac(ir2and)     = scli7d
       dscfacdt(ir2and)  = scli7ddt
       dscfacdd(ir2and)  = scli7ddd
      end if


! if we have li7 and tritium

      if (ili7 .ne. 0 .and. ih3 .ne. 0) then

! li7(t,2n)2a
       i = i + 1
       ratdum(irli7t2n)    = ratraw(irli7t2n) * scli7t
       dratdumdt(irli7t2n) = dratrawdt(irli7t2n) * scli7t + ratraw(irli7t2n) * scli7tdt
       dratdumdd(irli7t2n) = dratrawdd(irli7t2n) * scli7t + ratraw(irli7t2n) * scli7tdd
       scfac(irli7t2n)     = scli7t
       dscfacdt(irli7t2n)  = scli7tdt
       dscfacdd(irli7t2n)  = scli7tdd

       i = i + 1
       ratdum(ir2a2nt)    = ratraw(ir2a2nt) * scli7t
       dratdumdt(ir2a2nt) = dratrawdt(ir2a2nt) * scli7t + ratraw(ir2a2nt) * scli7tdt
       dratdumdd(ir2a2nt) = dratrawdd(ir2a2nt) * scli7t + ratraw(ir2a2nt) * scli7tdd
       scfac(ir2a2nt)     = scli7t
       dscfacdt(ir2a2nt)  = scli7tdt
       dscfacdd(ir2a2nt)  = scli7tdd
      end if


! if we have li7 and he3

      if (ili7 .ne. 0 .and. ihe3 .ne. 0) then

! li7(he3,np)2a
       i = i + 1
       ratdum(irli7he3np)    = ratraw(irli7he3np) * scli7he3
       dratdumdt(irli7he3np) = dratrawdt(irli7he3np) * scli7he3 + ratraw(irli7he3np) * scli7he3dt
       dratdumdd(irli7he3np) = dratrawdd(irli7he3np) * scli7he3 + ratraw(irli7he3np) * scli7he3dd
       scfac(irli7he3np)     = scli7he3
       dscfacdt(irli7he3np)  = scli7he3dt
       dscfacdd(irli7he3np)  = scli7he3dd

       i = i + 1
       ratdum(ir2apnhe3)    = ratraw(ir2apnhe3) * scli7he3
       dratdumdt(ir2apnhe3) = dratrawdt(ir2apnhe3) * scli7he3 + ratraw(ir2apnhe3) * scli7he3dt
       dratdumdd(ir2apnhe3) = dratrawdd(ir2apnhe3) * scli7he3 + ratraw(ir2apnhe3) * scli7he3dd
       scfac(ir2apnhe3)     = scli7he3
       dscfacdt(ir2apnhe3)  = scli7he3dt
       dscfacdd(ir2apnhe3)  = scli7he3dd
      end if




! if we have be7 and deutrium
! be7(d,p)2a
      if (ibe7 .ne. 0 .and. ih2 .ne. 0) then
       i = i + 1
       ratdum(irbe7dp)    = ratraw(irbe7dp) * scbe7d
       dratdumdt(irbe7dp) = dratrawdt(irbe7dp) * scbe7d + ratraw(irbe7dp) * scbe7ddt
       dratdumdd(irbe7dp) = dratrawdd(irbe7dp) * scbe7d + ratraw(irbe7dp) * scbe7ddd
       scfac(irbe7dp)     = scbe7d
       dscfacdt(irbe7dp)  = scbe7ddt
       dscfacdd(irbe7dp)  = scbe7ddd

       i = i + 1
       ratdum(ir2apd)    = ratraw(ir2apd) * scbe7d
       dratdumdt(ir2apd) = dratrawdt(ir2apd) * scbe7d + ratraw(ir2apd) * scbe7ddt
       dratdumdd(ir2apd) = dratrawdd(ir2apd) * scbe7d + ratraw(ir2apd) * scbe7ddd
       scfac(ir2apd)     = scbe7d
       dscfacdt(ir2apd)  = scbe7ddt
       dscfacdd(ir2apd)  = scbe7ddd
      end if


! if we have be7 and tritium
! be7(t,np)2a
      if (ibe7 .ne. 0 .and. ih3 .ne. 0) then
       i = i + 1
       ratdum(irbe7tnp)    = ratraw(irbe7tnp) * scbe7t
       dratdumdt(irbe7tnp) = dratrawdt(irbe7tnp) * scbe7t + ratraw(irbe7tnp) * scbe7tdt
       dratdumdd(irbe7tnp) = dratrawdd(irbe7tnp) * scbe7t + ratraw(irbe7tnp) * scbe7tdd
       scfac(irbe7tnp)     = scbe7t
       dscfacdt(irbe7tnp)  = scbe7tdt
       dscfacdd(irbe7tnp)  = scbe7tdd

       i = i + 1
       ratdum(ir2apnt)    = ratraw(ir2apnt) * scbe7t
       dratdumdt(ir2apnt) = dratrawdt(ir2apnt) * scbe7t + ratraw(ir2apnt) * scbe7tdt
       dratdumdd(ir2apnt) = dratrawdd(ir2apnt) * scbe7t + ratraw(ir2apnt) * scbe7tdd
       scfac(ir2apnt)     = scbe7t
       dscfacdt(ir2apnt)  = scbe7tdt
       dscfacdd(ir2apnt)  = scbe7tdd
      end if

! if we have be7 and he3
! be7(he3,2p)2a
      if (ibe7 .ne. 0 .and. ihe3 .ne. 0) then
       i = i + 1
       ratdum(irbe7he32p)    = ratraw(irbe7he32p) * scbe7he3
       dratdumdt(irbe7he32p) = dratrawdt(irbe7he32p) * scbe7he3 + ratraw(irbe7he32p) * scbe7he3dt
       dratdumdd(irbe7he32p) = dratrawdd(irbe7he32p) * scbe7he3 + ratraw(irbe7he32p) * scbe7he3dd
       scfac(irbe7he32p)     = scbe7he3
       dscfacdt(irbe7he32p)  = scbe7he3dt
       dscfacdd(irbe7he32p)  = scbe7he3dd

       i = i + 1
       ratdum(ir2a2phe3)    = ratraw(ir2a2phe3) * scbe7he3
       dratdumdt(ir2a2phe3) = dratrawdt(ir2a2phe3) * scbe7he3 + ratraw(ir2a2phe3) * scbe7he3dt
       dratdumdd(ir2a2phe3) = dratrawdd(ir2a2phe3) * scbe7he3 + ratraw(ir2a2phe3) * scbe7he3dd
       scfac(ir2a2phe3)     = scbe7he3
       dscfacdt(ir2a2phe3)  = scbe7he3dt
       dscfacdd(ir2a2phe3)  = scbe7he3dd
      end if



! if we have be9
      if (ibe9 .ne. 0) then

! a(an,g)be9
       i = i + 1
       ratdum(iraan)    = ratraw(iraan) * sc1a
       dratdumdt(iraan) = dratrawdt(iraan)*sc1a + ratraw(iraan)*sc1adt
       dratdumdd(iraan) = dratrawdd(iraan)*sc1a + ratraw(iraan)*sc1add
       scfac(iraan)     = sc1a
       dscfacdt(iraan)  = sc1adt
       dscfacdd(iraan)  = sc1add

       i = i + 1
       ratdum(irgaan)     = ratraw(irgaan) * sc1a
       dratdumdt(irgaan)  = dratrawdt(irgaan)*sc1a + ratraw(irgaan)*sc1adt
       dratdumdd(irgaan)  = dratrawdd(irgaan)*sc1a + ratraw(irgaan)*sc1add
       scfac(irgaan)      = sc1a
       dscfacdt(irgaan)   = sc1adt
       dscfacdd(irgaan)   = sc1add
      end if

! if we have be9 and deutrium
      if (ibe9 .ne. 0 .and. ih2 .ne. 0) then

! be9(p,d)be8 =>2a
       i = i + 1
       ratdum(irbe9pd)    = ratraw(irbe9pd)
       dratdumdt(irbe9pd) = dratrawdt(irbe9pd)
       dratdumdd(irbe9pd) = dratrawdd(irbe9pd)
       scfac(irbe9pd)     = 1.0d0
       dscfacdt(irbe9pd)  = 0.0d0
       dscfacdd(irbe9pd)  = 0.0d0

       i = i + 1
       ratdum(ir2adp)    = ratraw(ir2adp)
       dratdumdt(ir2adp) = dratrawdt(ir2adp)
       dratdumdd(ir2adp) = dratrawdd(ir2adp)
       scfac(ir2adp)     = 1.0d0
       dscfacdt(ir2adp)  = 0.0d0
       dscfacdd(ir2adp)  = 0.0d0
      end if


! if we have b8
      if (ib8 .ne. 0) then

! b8(p=>n)be8 =>2a reactions
       i = i + 1
       ratdum(irb8ep)    = ratraw(irb8ep)
       dratdumdt(irb8ep) = dratrawdt(irb8ep)
       dratdumdd(irb8ep) = dratrawdd(irb8ep)
       scfac(irb8ep)     = 1.0d0
       dscfacdt(irb8ep)  = 0.0d0
       dscfacdd(irb8ep)  = 0.0d0
      end if


! if we have b11
      if (ib11 .ne. 0) then

! b11(p,a)be8 => 2a reactions
       i = i + 1
       ratdum(irb11pa)    = ratraw(irb11pa)
       dratdumdt(irb11pa) = dratrawdt(irb11pa)
       dratdumdd(irb11pa) = dratrawdd(irb11pa)
       scfac(irb11pa)     = 1.0d0
       dscfacdt(irb11pa)  = 0.0d0
       dscfacdd(irb11pa)  = 0.0d0

       i = i + 1
       ratdum(ir3ap)    = ratraw(ir3ap)
       dratdumdt(ir3ap) = dratrawdt(ir3ap)
       dratdumdd(ir3ap) = dratrawdd(ir3ap)
       scfac(ir3ap)     = 1.0d0
       dscfacdt(ir3ap)  = 0.0d0
       dscfacdd(ir3ap)  = 0.0d0
      end if


! if we have c11
! c11(n,a)be8 => 2a
      if (ic11 .ne. 0) then
       i = i + 1
       ratdum(irc11na)     = ratraw(irc11na)
       dratdumdt(irc11na)  = dratrawdt(irc11na)
       dratdumdd(irc11na)  = dratrawdd(irc11na)
       scfac(irc11na)      = 1.0d0
       dscfacdt(irc11na)   = 0.0d0
       dscfacdd(irc11na)   = 0.0d0

       i = i + 1
       ratdum(ir3an)     = ratraw(ir3an)
       dratdumdt(ir3an)  = dratrawdt(ir3an)
       dratdumdd(ir3an)  = dratrawdd(ir3an)
       scfac(ir3an)      = 1.0d0
       dscfacdt(ir3an)   = 0.0d0
       dscfacdd(ir3an)   = 0.0d0
      end if



! bullet check the counting
      if (i .ne. nrat) then
       write(6,*)
       write(6,*) 'in screen_torch i=',i
       write(6,*) 'in screen_torch nrat=',nrat
       write(6,*) 'i is not equal to nrat'
       write(6,*) 'fatal counting error'
       write(6,*)
       stop 'fatal counting error in screen_torch'
      end if



! shut down the specials
       if (shut_down_specials .eq. 1) then
        do i=irpen,nrat
!        do i=irpen+24,nrat
!        do i=irpen+55,nrat
         ratraw(i)    = 0.0d0
         dratrawdt(i) = 0.0d0
         dratrawdd(i) = 0.0d0
         ratdum(i)    = 0.0d0
         dratdumdt(i) = 0.0d0
         dratdumdd(i) = 0.0d0
        end do
       end if



! debugs
!      do i=1,nrat
!       write(6,110) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
! 110   format(1x,i4,' ',a,' ',1p3e12.4)
!       if (ratraw(i) .lt. 0.0) stop 'negative rate'
!      enddo
!      read(5,*)


! write a specific set of rates into and out of a specified isotope

!108   format(1x,i4,' ',a,' ',1pe12.4)
!      i = 0
!      do j=ionbeg,ionend

! (n,g) and (g,n) reactions
!       k = nrr(1,j)
!       if (k .gt. 0) then
!        i            = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i            = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!       end if

! set up the (p,n) beta- beta+ decay components
!       k = nrr(2,j)
!       if (k .gt. 0) then
!        i            = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i            = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i            = i + 1
!        write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i            = i + 1
!        write(6,108) i, trim(ratnam(i)), ratraw(i)
!       end if

! (p,g) and (g,p) reactions
!       k = nrr(3,j)
!       if (k .gt. 0) then
!        i            = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i            = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!       end if

! (a,p) and (p,a) reactions
!       k = nrr(4,j)
!       if (k .gt. 0) then
!        i           = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i           = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!       end if

! (a,n) and (n,a) reactions
!       k = nrr(5,j)
!       if (k .gt. 0) then
!        i           = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i           = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!       end if

! (a,g) and (g,a) reactions
!       k = nrr(6,j)
!       if (k .gt. 0) then
!        i           = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!        i           = i + 1
!        if (ionam(j)(1:4) .eq. 'ca45') write(6,108) i, trim(ratnam(i)), ratraw(i)
!       end if
!      enddo
!      read(5,*)


      return
      end subroutine screen_torch

      subroutine screen6(jscreen, &
                         temp,den,zbar,abar,z2bar, &
                         z1,a1,z2,a2, &
                         scor,scordt,scordd)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'vector_eos.dek'

! this subroutine calculates screening factors and their derivatives
! for nuclear reaction rates in the weak, intermediate and strong regimes.
! based on dewitt, graboske, cooper apj, 181, 439, 1973 and 
! graboske, dewit, grossman and cooper apj 181 457 1973 for weak screening. 
!
! based on alastuey and jancovici apj 226 1034 1978 
! and jancovici jsp 17, 5. 357 1977 with plasma parameters 
! from itoh et al apj 234 1079 1979, for strong screening.

! vector version

! input:
! jscreen = length of vector
! temp    = temperature
! den     = density
! zbar    = mean charge per nucleus
! abar    = mean number of nucleons per nucleus
! z2bar   = mean square charge per nucleus
! z1 a1   = charge and number in the entrance channel
! z2 a2   = charge and number in the exit channel

! output:
! scor    = screening correction
! scordt  = derivative of screening correction with temperature
! scordd  = derivative of screening correction with density


! declare the pass
      integer          jscreen
      double precision temp,den,zbar,abar,z2bar, &
                       z1(jscreen),a1(jscreen),z2(jscreen),a2(jscreen), &
                       scor(jscreen),scordt(jscreen),scordd(jscreen)


! local variables
      integer          i, iscreen_type
      double precision bb,cc,dccdt,dccdd, &
                       pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd, &
                       ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd, &
                       vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni, &
                       qlam0z,qlam0zdt,qlam0zdd, &
                       h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd, &
                       h12x,dh12xdt,dh12xdd,alfa,beta, &
                       taufac,taufacdt,gamp,gampdt,gampdd, &
                       gamef,gamefdt,gamefdd, &
                       tau12,tau12dt,alph12,alph12dt,alph12dd, &
                       xlgfac,dxlgfacdt,dxlgfacdd, &
                       gamp14,gamp14dt,gamp14dd, &
                       xni,dxnidd,ytot 

     double precision  zfermim12,zfermi12,eta,detadt,detadd,theta_e,dtheta_e_dt,dtheta_e_dd, &
                       fdm12,fdetam12,fdthetam12,fdeta2m12,fdtheta2m12,fdetadthetam12, &
                       fdp12,fdetap12,fdthetap12,fdeta2p12,fdtheta2p12,fdetadthetap12


! screening variables
! zs13    = (z1+z2)**(1./3.)
! zhat    = combination of z1 and z2 raised to the 5/3 power
! zhat2   = combination of z1 and z2 raised to the 5/12 power
! lzav    = log of effective charge
! aznut   = combination of a1, z1, a2, z2 raised to 1/3 power


      integer          nscreen_max
      parameter        (nscreen_max = 2*abignet + 40)


      double precision zs13(nscreen_max),zhat(nscreen_max), &
                       zhat2(nscreen_max),lzav(nscreen_max), &
                       aznut(nscreen_max),zs13inv(nscreen_max), &
                       fac1(nscreen_max),fac2(nscreen_max), &
                       h12_vec(nscreen_max), &
                       dh12dt_vec(nscreen_max), &
                       dh12dd_vec(nscreen_max)


! parameter fact is the cube root of 2
      double precision  x13,x14,x53,x532,x512,fact,co2,gamefx,gamefs,dgamma,blend_frac
      parameter        (x13        = 1.0d0/3.0d0, &
                        x14        = 1.0d0/4.0d0, &
                        x53        = 5.0d0/3.0d0, &
                        x532       = 5.0d0/32.0d0, &
                        x512       = 5.0d0/12.0d0, &
                        fact       = 1.25992104989487d0, &
                        co2        = x13 * 4.248719d3, &
                        gamefx     = 0.3d0, &
                        gamefs     = 0.8d0, &
                        dgamma     = 1.0d0/(gamefs - gamefx), &
                        blend_frac = 0.05d0)


! compute and store the more expensive screening factors
      if (init_screen .eq. 1) then
       init_screen = 0
!       write(6,*) 'begin store of screening factors'

       if (jscreen .gt. nscreen_max) &
       stop 'jscreen > nscreen_max in screen6'

       do i=1,jscreen
        zs13(i)    = (z1(i) + z2(i))**x13
        zs13inv(i) = 1.0d0/zs13(i)
        zhat(i)    = (z1(i) + z2(i))**x53  - z1(i)**x53 - z2(i)**x53
        zhat2(i)   = (z1(i) + z2(i))**x512 - z1(i)**x512 -z2(i)**x512
        lzav(i)    = x53 * log(z1(i)*z2(i)/(z1(i) + z2(i)))
        aznut(i)   = (z1(i)**2*z2(i)**2*a1(i)*a2(i)/(a1(i)+a2(i)))**x13
        fac1(i)    = 0.896434d0 * zhat(i)
        fac2(i)    = 3.44740d0  * zhat2(i)
       enddo
!       write(6,*) 'done store of screening factors'
      endif


! calculate average plasma
      ytot     = 1.0d0/abar
      rr       = den * ytot
      tempi   = 1.0d0/temp
      dtempi  = -tempi*tempi
      deni    = 1.0d0/den


! etaele_row(1), detat_row(1),  detad_row(1) from torch_react

! call an eos to get the chemical potential
       temp_row(1) = temp
       den_row(1)  = den
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1
       call helmeos
      
! for the salpeter factor fprime/f 
       eta    = etaele_row(1)
       detadt = detat_row(1)
       detadd = detad_row(1)

       call dfermi(-0.5d0,eta,0.0d0,fdm12,fdetam12,fdthetam12,fdeta2m12,fdtheta2m12,fdetadthetam12)
       call dfermi(+0.5d0,eta,0.0d0,fdp12,fdetap12,fdthetap12,fdeta2p12,fdtheta2p12,fdetadthetap12)
       theta_e = fdm12/fdp12
       bb      = (fdetam12 - theta_e*fdetap12)/fdp12 
       dtheta_e_dt = bb * detadt
       dtheta_e_dd = bb * detadd

! no salpeter factor
!       theta_e     = 1.0d0
!       dtheta_e_dt = 0.0d0
!       dtheta_e_dd = 0.0d0

!       write(6,'(1x,a,1p3e12.4)') 'theta_e=',theta_e,dtheta_e_dt,dtheta_e_dd
    

! equation 4 of dewitt et al becomes ztilde = z2bar + zbar
      pp       = sqrt(rr*tempi*(z2bar + zbar*theta_e))
      qq       = 0.5d0/pp * (z2bar + zbar*theta_e)
      bb       = rr*tempi*zbar
      dppdt    = qq * (rr*dtempi  + bb*dtheta_e_dt)
      dppdd    = qq * (ytot*tempi + bb*dtheta_e_dd)


! equation 6 of dewitt et al
      qlam0z   = 1.88d8 * tempi * pp
      qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
      qlam0zdd = 1.88d8 * tempi * dppdd

      taufac   = co2 * tempi**x13
      taufacdt = -x13*taufac*tempi

      qq      = rr * zbar
      xni     = qq**x13
      dxnidd  = x13 * xni * deni

      gamp    = 2.27493d5 * tempi * xni
      gampdt  = 2.27493d5 * dtempi * xni
      gampdd  = 2.27493d5 * tempi * dxnidd



! calculate individual screening factors, start the pipeline
      do i=1,jscreen

       bb       = z1(i) * z2(i)
       qq       = fact * bb * zs13inv(i)
       gamef    = qq * gamp
       gamefdt  = qq * gampdt
       gamefdd  = qq * gampdd


! equation 6 of itoh et al
       tau12    = taufac * aznut(i)
       tau12dt  = taufacdt * aznut(i)

       qq       = 1.0d0/tau12
       alph12   = gamef * qq
       alph12dt = (gamefdt - alph12*tau12dt) * qq
       alph12dd = gamefdd * qq



! limit alph12 to 1.6 to prevent unphysical behavior.
! this should really be replaced by a pycnonuclear reaction rate formula
       if (alph12 .gt. 1.6) then
        alph12   = 1.6d0
        alph12dt = 0.0d0
        alph12dd = 0.0d0

        gamef    = 1.6d0 * tau12
        gamefdt  = 1.6d0 * tau12dt
        gamefdd  = 0.0d0

        qq       = zs13(i)/(fact * bb)
        gamp     = gamef * qq
        gampdt   = gamefdt * qq
        gampdd   = 0.0d0
       end if


! check gamef
!      write(6,'(1x,a,1p3e12.4)') 'gamef, gamefx, gamefs =',gamef, gamefx, gamefs



! weak screening regime
       h12w    = bb * qlam0z
       dh12wdt = bb * qlam0zdt
       dh12wdd = bb * qlam0zdd

       h12     = h12w
       dh12dt  = dh12wdt
       dh12dd  = dh12wdd



! intermediate and strong sceening regime
       if (gamef .gt. gamefx) then

        qq       = sqrt(gamp)
        gamp14   = sqrt(qq)
        rr       = 1.0d0/gamp
        qq       = 0.25d0*gamp14*rr
        gamp14dt = qq * gampdt
        gamp14dd = qq * gampdd


! equation 25 of alastuey
! coefficients from equation 16 of jancovici

        cc       =   gamp * fac1(i) - gamp14 * fac2(i) &
                   - 0.5551d0   * (log(gamp) + lzav(i)) &
                   - 2.996d0

        dccdt    =   gampdt * fac1(i)  - gamp14dt * fac2(i) &
                   - 0.5551d0*rr*gampdt

        dccdd    =   gampdd * fac1(i) - gamp14dd * fac2(i) &
                   - 0.5551d0*rr*gampdd

        qq     = alph12 * alph12
        a3     = qq * alph12
        da3    = 3.0d0 * qq


! equation 28 of alastuey
        qq     = 0.014d0 + 0.0128d0*alph12
        dqqdt  = 0.0128d0*alph12dt
        dqqdd  = 0.0128d0*alph12dd

        rr     = x532 - alph12*qq
        drrdt  = -(alph12dt*qq + alph12*dqqdt)
        drrdd  = -(alph12dd*qq + alph12*dqqdd)

        ss     = tau12*rr
        dssdt  = tau12dt*rr + tau12*drrdt
        dssdd  = tau12*drrdd


! equation 31 of alastuey
        tt     =  -0.0098d0 + 0.0048d0*alph12
        dttdt  = 0.0048d0*alph12dt
        dttdd  = 0.0048d0*alph12dd

        uu     =  0.0055d0 + alph12*tt
        duudt  = alph12dt*tt + alph12*dttdt
        duudd  = alph12dd*tt + alph12*dttdd

        vv   = gamef * alph12 * uu
        dvvdt= gamefdt*alph12*uu + gamef*(alph12dt*uu + alph12*duudt)
        dvvdd= gamefdd*alph12*uu + gamef*(alph12dd*uu + alph12*duudd)


! equation 32 of alastuey
        h12     = cc - a3 * (ss + vv)
        rr      = da3 * (ss + vv)
        dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
        dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)


! equation 34 of alastuey
        rr     =  1.0d0 - 0.0562d0*a3
        ss     =  -0.0562d0*da3
        drrdt  = ss*alph12dt
        drrdd  = ss*alph12dd


        if (rr .ge. 0.77d0) then
         xlgfac    = rr
         dxlgfacdt = drrdt
         dxlgfacdd = drrdd
        else
         xlgfac    = 0.77d0
         dxlgfacdt = 0.0d0
         dxlgfacdd = 0.0d0
        end if


        h12    = log(xlgfac) + h12
        rr     = 1.0d0/xlgfac
        dh12dt = rr*dxlgfacdt + dh12dt
        dh12dd = rr*dxlgfacdd + dh12dd


        if (gamef .le. gamefs) then
         iscreen_type = intermediate_screening

         rr     =  dgamma*(gamefs - gamef)
         drrdt  = -dgamma*gamefdt
         drrdd  = -dgamma*gamefdd

         ss     = dgamma*(gamef-gamefx)
         dssdt  = dgamma*gamefdt
         dssdd  = dgamma*gamefdd


! store current values for possible blending with weak regime
         h12x    = h12
         dh12xdt = dh12dt
         dh12xdd = dh12dd

         vv     = h12
         h12    = h12w*rr + vv*ss
         dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
         dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd


! blend the transition region - from bill paxton
!        if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
!          alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
!          alfa   = 0.5d0 * (1d0 - cos(pi*alfa))
!          beta   = 1.0d0 - alfa
!          h12    = alfa * h12 + beta * h12x
!          dh12dt = alfa * dh12dt + beta * dh12xdt
!          dh12dd = alfa * dh12dd + beta * dh12xdd
!         end if

        else
         iscreen_type = strong_screening
        end if


! end of intermediate and strong screening if
       else
       iscreen_type = weak_screening
!        write(6,'(1x,a,1p2e12.4)') 'weak screening gamef gamefx =',gamef,gamefx
       end if


! limit h12 to avoid the pyconuclear regime
       h12           = max(min(h12,30.0d0),0.0d0)
       if (h12 .eq. 30.0d0) then
        dh12dt = 0.0d0
        dh12dd = 0.0d0
       end if

! store what we got
       scor(i)          = exp(h12)
       scordt(i)        = scor(i) * dh12dt
       scordd(i)        = scor(i) * dh12dd
       screen_regime(i) = iscreen_type


!       h12_vec(i)    = h12
!       dh12dt_vec(i) = dh12dt
!       dh12dt_vec(i) = dh12dd

!       scor(i)   = exp(h12_vec(i))
!       scordt(i) = scor(i) * dh12dt_vec(i)
!       scordd(i) = scor(i) * dh12dd_vec(i)


! end of individual screening pipeline
      end do



! crank the exponential
!       scor(1:jscreen)   = exp(h12_vec(1:jscreen))
!       scordt(1:jscreen) = scor(1:jscreen) * dh12dt_vec(1:jscreen)
!       scordd(1:jscreen) = scor(1:jscreen) * dh12dd_vec(1:jscreen)

      return
    end subroutine screen6
       
      subroutine rhstorch(y,rte,rate,dydt)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'


! evaluates the right hand side of the aprox13 odes
! evaluating the running sum in real*16 precision significantly
! improves the efficiency of the integration and the mass conservation.

! declare the pass
      real*8 :: y(*),rte(14,abignet),rate(*),dydt(*)


! local variables
      integer           :: i,j,k,in,ip,ia
      real*8            :: aaa,aap,aan
      real*8, parameter :: sixth = 1.0d0/6.0d0


! double precision dydt sums
!      real*8           :: da0,da1,da2,da3,da4,db1,db2,db3,db4,db5
!      real*8           :: dray(abignet)


! quad precision dydt sums
      real*8           :: a0,a1,a2,a3,a4,b1,b2,b3,b4,b5
      real*8           :: qray(abignet)


! for storing the terms of the sum
!      integer, parameter :: cdim2 = 4 * abignet
!      integer :: icasc(abignet),indx(abignet)
!      real*8  :: casc(abignet,cdim2),cdum(abignet)
!      real*8  :: ksum,kdum(abignet)
!      real*16 :: mpdum(abignet)
!      real*16 :: value16
!      real*16 :: qasc(abignet,cdim2),qdum(abignet)



! for the statistics
!      integer           :: icount=0
!      real*16           :: bigdiff
!      integer           :: ibigdiff(1)


! declare for brent mp
!      integer       :: ifirst = 1
!      integer       :: mm
!      character*110 :: string
!      integer       :: iterm(64), XS(64), C(110)


! communicate with brent mp
!      INTEGER :: B, T, M, LUN, MXR, SPTR, MXSPTR, DUMMY(16), R(500)
!      COMMON R
!      COMMON /MPCOM/ B, T, M, LUN, MXR, SPTR, MXSPTR, DUMMY


! initialize brent mp
!      if (ifirst .eq. 1) then
!       ifirst = 0
!       call mpset2 (6, 110, 64, 1, 500)
!       call mpparb(4,'EXWID')
!     end if




! zero the odes
       dydt(1:neqs) = 0.0d0
!       dray(1:neqs) = 0.0d0
       qray(1:neqs) = 0.0q0
!       icasc(1:ionmax) = 0
!       casc(1:ionmax,1:cdim2) = 0.0d0
!       qasc(1:ionmax,1:cdim2) = 0.0q0




! set the pointers for n,p and a reactions
      in  = ineut
      aan = y(ineut)
      ip  = iprot
      aap = y(iprot)
      ia  = ihe4
      aaa = y(ihe4)


! for every isotope in the network
      do  j=ionbeg,ionend

! set up the y(j)(n,g)y(k) and y(k)(g,n)y(j) components
       k = nrr(1,j)
       if (k .gt. 0) then

!        da1       = -aan*rte(1,j)*y(j)
!        da2       = rte(2,j)*y(k)
!        db1       = da1 + da2
!        dray(j)  = dray(j)  + db1
!        dray(k)  = dray(k)  - db1
!        dray(in) = dray(in) + db1


        a1       = -aan*rte(1,j)*y(j)
        a2       = rte(2,j)*y(k)
        b1       = a1 + a2

        qray(j)  = qray(j)  + b1
        qray(k)  = qray(k)  - b1
        qray(in) = qray(in) + b1


!        icasc(J) = icasc(j) + 1
!        casc(j,icasc(j)) = db1
!        qasc(j,icasc(j)) = b1
!        icasc(k) = icasc(k) + 1
!        casc(k,icasc(k)) = -db1
!        qasc(k,icasc(k)) = -b1
!        icasc(in) = icasc(in) + 1
!        casc(in,icasc(in)) = db1
!        qasc(in,icasc(in)) = b1
       end if


! set up the (p,n) (n,p) and beta- beta+ decay components
       k = nrr(2,j)
       if (k .gt. 0) then

!        da1       = -aap*rte(3,j)*y(j)
!        da2       = aan*rte(4,j)*y(k)
!        da3       = -y(j)*rte(5,j) 
!        da4       = y(k)*rte(6,j)
!        db1       = da1 + da2
!        db2       = da3 + da4 
!        db3       = db1 + db2
!        dray(j)  = dray(j)  + db3
!        dray(k)  = dray(k)  - db3
!        dray(in) = dray(in) - db1
!        dray(ip) = dray(ip) + db1


        a1       = -aap*rte(3,j)*y(j)
        a2       = aan*rte(4,j)*y(k)
        a3       = -y(j)*rte(5,j) 
        a4       = y(k)*rte(6,j)
        b1       = a1 + a2
        b2       = a3 + a4 
        b3       = b1 + b2

        qray(j)  = qray(j)  + b3
        qray(k)  = qray(k)  - b3
        qray(in) = qray(in) - b1
        qray(ip) = qray(ip) + b1

!        icasc(J) = icasc(j) + 1
!        casc(j,icasc(j)) = db3
!        qasc(j,icasc(j)) = b3
!        icasc(k) = icasc(k) + 1
!        casc(k,icasc(k)) = -db3
!        qasc(k,icasc(k)) = -b3
!        icasc(in) = icasc(in) + 1
!        casc(in,icasc(in)) = -db1
!        qasc(in,icasc(in)) = -b1
!        icasc(ip) = icasc(ip) + 1
!        casc(ip,icasc(ip)) = db1
!        qasc(ip,icasc(ip)) = b1
       end if


! set up the (p,g) and (g,p) components
       k = nrr(3,j)
       if (k .gt. 0) then

!        da1       = -aap*rte(7,j)*y(j)
!        da2       = rte(8,j)*y(k)
!        db1       = da1 + da2 
!        dray(j)  = dray(j)  + db1
!        dray(k)  = dray(k)  - db1
!        dray(ip) = dray(ip) + db1


        a1       = -aap*rte(7,j)*y(j)
        a2       = rte(8,j)*y(k)
        b1       = a1 + a2 

        qray(j)  = qray(j)  + b1
        qray(k)  = qray(k)  - b1
        qray(ip) = qray(ip) + b1


!        icasc(J) = icasc(j) + 1
!        casc(j,icasc(j)) = db1
!        qasc(j,icasc(j)) = b1
!        icasc(k) = icasc(k) + 1
!        casc(k,icasc(k)) = -db1
!        qasc(k,icasc(k)) = -b1
!        icasc(ip) = icasc(ip) + 1
!        casc(ip,icasc(ip)) = db1
!        qasc(ip,icasc(ip)) = b1
       end if


! set up the (a,p) and (p,a) reactions
       k = nrr(4,j)
       if (k .gt. 0) then

!        da1       = -aaa*rte(9,j)*y(j) 
!        da2       = aap*rte(10,j)*y(k)
!        db1       = da1 + da2 
!        dray(j)  = dray(j)  + db1
!        dray(k)  = dray(k)  - db1
!        dray(ip) = dray(ip) - db1
!        dray(ia) = dray(ia) + db1


        a1       = -aaa*rte(9,j)*y(j) 
        a2       = aap*rte(10,j)*y(k)
        b1       = a1 + a2 

        qray(j)  = qray(j)  + b1
        qray(k)  = qray(k)  - b1
        qray(ip) = qray(ip) - b1
        qray(ia) = qray(ia) + b1

!        icasc(J) = icasc(j) + 1
!        casc(j,icasc(j)) = db1
!        qasc(j,icasc(j)) = b1
!        icasc(k) = icasc(k) + 1
!        casc(k,icasc(k)) = -db1
!        qasc(k,icasc(k)) = -b1
!        icasc(ip) = icasc(ip) + 1
!        casc(ip,icasc(ip)) = -db1
!        qasc(ip,icasc(ip)) = -b1
!        icasc(ia) = icasc(ia) + 1
!        casc(ia,icasc(ia)) = db1
!        qasc(ia,icasc(ia)) = b1
       end if

! set up the (a,n) and (n,a) components
       k = nrr(5,j)
       if (k .gt. 0) then

!        da1       = -aaa*rte(11,j)*y(j)
!        da2       = aan*rte(12,j)*y(k)
!        db1       = da1 + da2 
!        dray(j)  = dray(j)  + db1
!        dray(k)  = dray(k)  - db1
!        dray(in) = dray(in) - db1
!        dray(ia) = dray(ia) + db1


        a1       = -aaa*rte(11,j)*y(j)
        a2       = aan*rte(12,j)*y(k)
        b1       = a1 + a2 

        qray(j)  = qray(j)  + b1
        qray(k)  = qray(k)  - b1
        qray(in) = qray(in) - b1
        qray(ia) = qray(ia) + b1


 !       icasc(J) = icasc(j) + 1
 !       casc(j,icasc(j)) = db1
 !       qasc(j,icasc(j)) = b1
 !       icasc(k) = icasc(k) + 1
 !       casc(k,icasc(k)) = -db1
 !       qasc(k,icasc(k)) = -b1
 !       icasc(in) = icasc(in) + 1
 !       casc(in,icasc(in)) = -db1
 !       qasc(in,icasc(in)) = -b1
 !       icasc(ia) = icasc(ia) + 1
 !       casc(ia,icasc(ia)) = db1
 !       qasc(ia,icasc(ia)) = b1

       end if

! and the (a,g) and (g,a) components
       k = nrr(6,j)
       if (k .gt. 0) then

!        da1       = -aaa*rte(13,j)*y(j)
!        da2       = rte(14,j)*y(k)
!        db1       = da1 + da2 
!        dray(j)  = dray(j)  + db1
!        dray(k)  = dray(k)  - db1
!        dray(ia) = dray(ia) + db1


        a1       = -aaa*rte(13,j)*y(j)
        a2       = rte(14,j)*y(k)
        b1       = a1 + a2 

        qray(j)  = qray(j)  + b1
        qray(k)  = qray(k)  - b1
        qray(ia) = qray(ia) + b1


 !       icasc(J) = icasc(j) + 1
 !       casc(j,icasc(j)) = db1
 !       qasc(j,icasc(j)) = b1
 !       icasc(k) = icasc(k) + 1
 !       casc(k,icasc(k)) = -db1
 !       qasc(k,icasc(k)) = -b1
 !       icasc(ia) = icasc(ia) + 1
 !       casc(ia,icasc(ia)) = db1
 !       qasc(ia,icasc(ia)) = b1
       end if
      enddo



! now the special matrix elements
! for p(e-,nu)n and n(e+,nub)p reactions

!      da1 = aan*rate(irnep)
!      da2 = aap*rate(irpen)
!      da3 = da2 - da1
!      da4 = da1 - da2
!      dray(in)    = dray(in) + da3
!      dray(ip)    = dray(ip) + da4


      a1 = aan*rate(irnep)
      a2 = aap*rate(irpen)
      a3 = a2 - a1
      a4 = a1 - a2

      qray(in)    = qray(in) + a3
      qray(ip)    = qray(ip) + a4


!      icasc(in) = icasc(in) + 1
!      casc(in,icasc(in)) = da3
!      qasc(in,icasc(in)) = a3
!      icasc(ip) = icasc(ip) + 1
!      casc(ip,icasc(ip)) = da4
!      qasc(ip,icasc(ip)) = a4


! triple alpha reactions; must have c12 in the network
      if (ic12 .ne. 0) then

!       da1 = sixth*rate(ir3a)*aaa*aaa*aaa
!       da2 = rate(irg3a)*y(ic12)
!       da3 = da1 - da2
!       dray(ia)   = dray(ia) - 3.0d0*da3
!       dray(ic12) = dray(ic12) + da3


       a1 = sixth*rate(ir3a)*aaa*aaa*aaa
       a2 = rate(irg3a)*y(ic12)
       a3 = a1 - a2

       qray(ia)   = qray(ia) - 3.0d0*a3
       qray(ic12) = qray(ic12) + a3


!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -3.0d0*da3
!       qasc(ia,icasc(ia)) = -3.0d0*a3
!       icasc(ic12) = icasc(ic12) + 1
!       casc(ic12,icasc(ic12)) = da3
!       qasc(ic12,icasc(ic12)) = a3


! c12+c12 reactions; must have ne20, na23, mg23 in the network
       if (ine20 .ne. 0  .and. ina23 .ne. 0  .and.  img23 .ne. 0) then


!        da1 = 0.5d0*y(ic12)*y(ic12)
!        da2 = aaa*y(ine20)*rate(irne20ac)
!        da3 = aap*y(ina23)*rate(irna23pc)
!        da4 = aan*y(img23)*rate(irmg23nc)

!        db1          = 2.0d0*( (da2 + da3 + da4) - da1*(rate(ir1212n) + rate(ir1212p) + rate(ir1212a)))
!        dray(ic12)  = dray(ic12)  + db1

!        db2 = da1*rate(ir1212a) - da2
!        dray(ine20) = dray(ine20) + db2
!        dray(ia)    = dray(ia)    + db2

!        db3 = da1*rate(ir1212p) - da3 
!        dray(ina23) = dray(ina23) + db3
!        dray(ip)    = dray(ip)    + db3

!        db4 = da1*rate(ir1212n) - da4
!        dray(img23) = dray(img23) + db4
!        dray(in)    = dray(in)    + db4


        a1 = 0.5d0*y(ic12)*y(ic12)
        a2 = aaa*y(ine20)*rate(irne20ac)
        a3 = aap*y(ina23)*rate(irna23pc)
        a4 = aan*y(img23)*rate(irmg23nc)

        b1          = 2.0d0*( (a2 + a3 + a4) - a1*(rate(ir1212n) + rate(ir1212p) + rate(ir1212a)))
        qray(ic12)  = qray(ic12)  + b1

        b2 = a1*rate(ir1212a) - a2
        qray(ine20) = qray(ine20) + b2
        qray(ia)    = qray(ia)    + b2

        b3 = a1*rate(ir1212p) - a3 
        qray(ina23) = qray(ina23) + b3
        qray(ip)    = qray(ip)    + b3

        b4 = a1*rate(ir1212n) - a4
        qray(img23) = qray(img23) + b4
        qray(in)    = qray(in)    + b4



!       icasc(ic12) = icasc(ic12) + 1
!       casc(ic12,icasc(ic12)) = db1
!       qasc(ic12,icasc(ic12)) = b1

!       icasc(ine20) = icasc(ine20) + 1
!       casc(ine20,icasc(ine20)) = db2
!       qasc(ine20,icasc(ine20)) = b2
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db2
!       qasc(ia,icasc(ia)) = b2

!       icasc(ina23) = icasc(ina23) + 1
!       casc(ina23,icasc(ina23)) = db3
!       qasc(ina23,icasc(ina23)) = b3
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db3
!       qasc(ip,icasc(ip)) = b3

!       icasc(img23) = icasc(img23) + 1
!       casc(img23,icasc(img23)) = db4
!       qasc(img23,icasc(img23)) = b4
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db4
!       qasc(in,icasc(in)) = b4

       end if
      end if



! o16+o16 reactions; must have si28, p30, p31 and s31 in the network
      if (io16 .ne. 0 .and. isi28 .ne. 0 .and. ip30 .ne. 0 .and. &
          ip31 .ne. 0 .and. is31 .ne. 0) then

!       da0 = 0.5d0*y(io16)*y(io16)
!       da1 = aaa*y(isi28)*rate(irsi28ao)
!       da2 = aap*y(ip31)*rate(irp31po)
!       da3 = aan*y(is31)*rate(irs31no)

!       db1 = 2.0d0*(da1 + da2 + da3 - da0*(rate(ir1616n) + rate(ir1616p) + rate(ir1616a) + rate(ir1616d)))
!       dray(io16)  = dray(io16)  + db1 

!       db2 = da0*rate(ir1616a) - da1
!       dray(isi28) = dray(isi28) + db2
!       dray(ia)    = dray(ia)    + db2

!       db3 = da0*rate(ir1616p) - da2
!       dray(ip31)  = dray(ip31)  + db3
!       dray(ip)    = dray(ip)    + db3

!       db4 = da0*rate(ir1616n) - da3
!       dray(is31)  = dray(is31)  + db4
!       dray(in)    = dray(in)    + db4

!       db5 = da0*rate(ir1616d)
!       dray(ip30)  = dray(ip30)  + db5
!       dray(in)    = dray(in)    + db5
!       dray(ip)    = dray(ip)    + db5



       a0 = 0.5d0*y(io16)*y(io16)
       a1 = aaa*y(isi28)*rate(irsi28ao)
       a2 = aap*y(ip31)*rate(irp31po)
       a3 = aan*y(is31)*rate(irs31no)

       b1 = 2.0d0*(a1 + a2 + a3 - a0*(rate(ir1616n) + rate(ir1616p) + rate(ir1616a) + rate(ir1616d)))
       qray(io16)  = qray(io16)  + b1 

       b2 = a0*rate(ir1616a) - a1
       qray(isi28) = qray(isi28) + b2
       qray(ia)    = qray(ia)    + b2

       b3 = a0*rate(ir1616p) - a2
       qray(ip31)  = qray(ip31)  + b3
       qray(ip)    = qray(ip)    + b3

       b4 = a0*rate(ir1616n) - a3
       qray(is31)  = qray(is31)  + b4
       qray(in)    = qray(in)    + b4

       b5 = a0*rate(ir1616d)
       qray(ip30)  = qray(ip30)  + b5
       qray(in)    = qray(in)    + b5
       qray(ip)    = qray(ip)    + b5



!       icasc(io16) = icasc(io16) + 1
!       casc(io16,icasc(io16)) = db1
!       qasc(io16,icasc(io16)) = b1

!       icasc(isi28) = icasc(isi28) + 1
!       casc(isi28,icasc(isi28)) = db2
!       qasc(isi28,icasc(isi28)) = b2
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db2
!       qasc(ia,icasc(ia)) = b2

!       icasc(ip31) = icasc(ip31) + 1
!       casc(ip31,icasc(ip31)) = db3
!       qasc(ip31,icasc(ip31)) = b3
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db3
!       qasc(ip,icasc(ip)) = b3

!       icasc(is31) = icasc(is31) + 1
!       casc(is31,icasc(is31)) = db4
!       qasc(is31,icasc(is31)) = b4
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db4
!       qasc(in,icasc(in)) = b4

!       icasc(ip30) = icasc(ip30) + 1
!       casc(ip30,icasc(ip30)) = db5
!       qasc(ip30,icasc(ip30)) = b5
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db5
!       qasc(in,icasc(in)) = b5
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db5
!       qasc(ip,icasc(ip)) = b5

      end if


! c12+o16 reactions; must have mg24, al27, si27 in the network
      if (ic12 .ne. 0 .and. io16 .ne. 0 .and. &
          img24 .ne. 0 .and. ial27 .ne. 0 .and. isi27 .ne. 0) then

!       da0 = y(ic12)*y(io16)
!       da1 = aap*y(ial27)*rate(iral27po)
!       da2 = aan*y(isi27)*rate(irsi27no)
!       da3 = aaa*y(img24)*rate(irmg24ao)

!       db1 = da1 + da2 + da3 - da0*(rate(ir1216n) + rate(ir1216p) + rate(ir1216a))
!       dray(ic12)  = dray(ic12)  + db1
!       dray(io16)  = dray(io16)  + db1

!       db2 = da0*rate(ir1216n) - da2
!       dray(isi27) = dray(isi27) + db2
!       dray(in)    = dray(in)    + db2

!       db3 = da0*rate(ir1216p) - da1
!       dray(ial27) = dray(ial27) + db3
!       dray(ip)    = dray(ip)    + db3

!       db4 = da0*rate(ir1216a) - da3 
!       dray(img24) = dray(img24) + db4
!       dray(ia)    = dray(ia)    + db4



       a0 = y(ic12)*y(io16)
       a1 = aap*y(ial27)*rate(iral27po)
       a2 = aan*y(isi27)*rate(irsi27no)
       a3 = aaa*y(img24)*rate(irmg24ao)

       b1 = a1 + a2 + a3 - a0*(rate(ir1216n) + rate(ir1216p) + rate(ir1216a))
       qray(ic12)  = qray(ic12)  + b1
       qray(io16)  = qray(io16)  + b1

       b2 = a0*rate(ir1216n) - a2
       qray(isi27) = qray(isi27) + b2
       qray(in)    = qray(in)    + b2

       b3 = a0*rate(ir1216p) - a1
       qray(ial27) = qray(ial27) + b3
       qray(ip)    = qray(ip)    + b3

       b4 = a0*rate(ir1216a) - a3 
       qray(img24) = qray(img24) + b4
       qray(ia)    = qray(ia)    + b4


!       icasc(ic12) = icasc(ic12) + 1
!       casc(ic12,icasc(ic12)) = db1
!       qasc(ic12,icasc(ic12)) = b1
!       icasc(io16) = icasc(io16) + 1
!       casc(io16,icasc(io16)) = db1
!       qasc(io16,icasc(io16)) = b1

!       icasc(isi27) = icasc(isi27) + 1
!       casc(isi27,icasc(isi27)) = db2
!       qasc(isi27,icasc(isi27)) = b2
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db2
!       qasc(in,icasc(in)) = b2

!       icasc(ial27) = icasc(ial27) + 1
!       casc(ial27,icasc(ial27)) = db3
!       qasc(ial27,icasc(ial27)) = b3
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db3
!       qasc(ip,icasc(ip)) = b3

!       icasc(img24) = icasc(img24) + 1
!       casc(img24,icasc(img24)) = db4
!       qasc(img24,icasc(img24)) = b4
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db4
!       qasc(ia,icasc(ia)) = b4
      end if




! if we have deuterium
      if (ih2 .ne. 0) then

! pp and pep
!       da1          = aap*aap*rate(irpp)
!       da2          = aap*aap*rate(irpep)
!       db1          = da1 + da2 
!       dray(ip)    = dray(ip)  - db1
!       dray(ih2)   = dray(ih2) + 0.5d0 * db1


       a1          = aap*aap*rate(irpp)
       a2          = aap*aap*rate(irpep)
       b1          = a1 + a2 

       qray(ip)    = qray(ip)  - b1
       qray(ih2)   = qray(ih2) + 0.5d0 * b1

!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = -db1
!       qasc(ip,icasc(ip)) = -b1
!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = 0.5d0*db1
!       qasc(ih2,icasc(ih2)) = 0.5d0*b1



! p(n,g)d
!       da1          = -aan*rate(irpng)*aap 
!       da2          = rate(irdgn)*y(ih2)
!       db1          = da1 + da2 
!       dray(ip)    = dray(ip)  + db1
!       dray(in)    = dray(in)  + db1
!       dray(ih2)   = dray(ih2) - db1


       a1          = -aan*rate(irpng)*aap 
       a2          = rate(irdgn)*y(ih2)
       b1          = a1 + a2 

       qray(ip)    = qray(ip)  + b1
       qray(in)    = qray(in)  + b1
       qray(ih2)   = qray(ih2) - b1


!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -db1
!       qasc(ih2,icasc(ih2)) = -b1


! d(p,n)2p
!       da1          = aap*y(ih2)*rate(irdpn) 
!       da2          = -0.5d0*aan*aap*aap*rate(ir2pnp)
!       db1          = da1 + da2
!       dray(in)    = dray(in)  + db1
!       dray(ip)    = dray(ip)  + db1
!       dray(ih2)   = dray(ih2) - db1


       a1          = aap*y(ih2)*rate(irdpn) 
       a2          = -0.5d0*aan*aap*aap*rate(ir2pnp)
       b1          = a1 + a2

       qray(in)    = qray(in)  + b1
       qray(ip)    = qray(ip)  + b1
       qray(ih2)   = qray(ih2) - b1


!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -db1
!       qasc(ih2,icasc(ih2)) = -b1



! d(d,g)he4
!       da1          = -y(ih2)*y(ih2)*rate(irddg) 
!       da2          = 2.0d0*aaa*rate(irhe4gd)
!       db1          = da1 + da2 
!       dray(ih2)   = dray(ih2) + db1
!       dray(ia)    = dray(ia)  - 0.5d0 * db1


       a1          = -y(ih2)*y(ih2)*rate(irddg) 
       a2          = 2.0d0*aaa*rate(irhe4gd)
       b1          = a1 + a2 

       qray(ih2)   = qray(ih2) + b1
       qray(ia)    = qray(ia)  - 0.5d0 * b1

!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = db1
!       qasc(ih2,icasc(ih2)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -0.5d0*db1
!       qasc(ia,icasc(ia)) = -0.5d0*b1

      end if


! if we have deuterium and tritium
      if (ih2 .ne. 0 .and. ih3 .ne. 0) then

! d(t,n)he4, dt reaction
!       da1         = y(ih2)*y(ih3)*rate(irtdn) 
!       da2         = -aan*aaa*rate(irhe4nd)
!       db1         = da1 + da2 
!       dray(ih2)  = dray(ih2) - db1
!       dray(ih3)  = dray(ih3) - db1
!       dray(ia)   = dray(ia)  + db1
!       dray(in)   = dray(in)  + db1


       a1         = y(ih2)*y(ih3)*rate(irtdn) 
       a2         = -aan*aaa*rate(irhe4nd)
       b1         = a1 + a2 

       qray(ih2)  = qray(ih2) - b1
       qray(ih3)  = qray(ih3) - b1
       qray(ia)   = qray(ia)  + b1
       qray(in)   = qray(in)  + b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -db1
!       qasc(ih2,icasc(ih2)) = -b1
!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = -db1
!       qasc(ih3,icasc(ih3)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db1
!       qasc(ia,icasc(ia)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1


! d(d,p)t
!       da1           = 0.5d0*y(ih2)*y(ih2)*rate(irddp) 
!       da2           = -aap*y(ih3)*rate(irtpd)
!       db1           = da1 + da2
!       dray(ih2)   = dray(ih2) - 2.0d0 * db1
!       dray(ih3)   = dray(ih3) + db1
!       dray(ip)    = dray(ip)  + db1


       a1           = 0.5d0*y(ih2)*y(ih2)*rate(irddp) 
       a2           = -aap*y(ih3)*rate(irtpd)
       b1           = a1 + a2

       qray(ih2)   = qray(ih2) - 2.0d0 * b1
       qray(ih3)   = qray(ih3) + b1
       qray(ip)    = qray(ip)  + b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -2.0d0*db1
!       qasc(ih2,icasc(ih2)) = -2.0d0*b1
!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = db1
!       qasc(ih3,icasc(ih3)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
      end if


! if we have tritium

      if (ih3 .ne. 0) then

! t(p,g)he4
!       da1         = -aap*rate(irh3pg)*y(ih3) 
!       da2         = rate(irhe4gp)*aaa
!       db1         = da1 + da2
!       dray(ih3)  = dray(ih3)  + db1
!       dray(ia)   = dray(ia)   - db1
!       dray(ip)   = dray(ip)   + db1



       a1         = -aap*rate(irh3pg)*y(ih3) 
       a2         = rate(irhe4gp)*aaa
       b1         = a1 + a2

       qray(ih3)  = qray(ih3)  + b1
       qray(ia)   = qray(ia)   - b1
       qray(ip)   = qray(ip)   + b1


!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = db1
!       qasc(ih3,icasc(ih3)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -db1
!       qasc(ia,icasc(ia)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1


! t(t,2n)he4
!       da1          = y(ih3)*y(ih3)*rate(irtt2n) 
!       da2          = -aaa*aan*aan*rate(irhe42nt)
!       db1          = da1 + da2 
!       dray(ih3)  = dray(ih3)  - db1
!       dray(in)   = dray(in)   + db1
!       dray(ia)   = dray(ia)   + 0.5d0 * db1


       a1          = y(ih3)*y(ih3)*rate(irtt2n) 
       a2          = -aaa*aan*aan*rate(irhe42nt)
       b1          = a1 + a2 

       qray(ih3)  = qray(ih3)  - b1
       qray(in)   = qray(in)   + b1
       qray(ia)   = qray(ia)   + 0.5d0 * b1


!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = -db1
!       qasc(ih3,icasc(ih3)) = -b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 0.5d0*db1
!       qasc(ia,icasc(ia)) = 0.5d0*b1
      end if


! if we have he3

      if (ihe3 .ne. 0) then

! he3(he3,2p)he4
!       da1           = y(ihe3)*y(ihe3)*rate(ir33) 
!       da2           = -aaa*aap*aap*rate(ir33inv)
!       db1           = da1 + da2 
!       dray(ihe3)  = dray(ihe3) - db1
!       dray(ip)    = dray(ip)   + db1
!       dray(ia)    = dray(ia)   + 0.5d0 * db1


       a1           = y(ihe3)*y(ihe3)*rate(ir33) 
       a2           = -aaa*aap*aap*rate(ir33inv)
       b1           = a1 + a2 

       qray(ihe3)  = qray(ihe3) - b1
       qray(ip)    = qray(ip)   + b1
       qray(ia)    = qray(ia)   + 0.5d0 * b1


!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 0.5d0*db1
!       qasc(ia,icasc(ia)) = 0.5d0*b1


! he3(p,e+nu)he4
!       db1 = y(ihe3)*aap*rate(irhep)
!       dray(ihe3)  = dray(ihe3) - db1
!       dray(ip)    = dray(ip)   - db1
!       dray(ia)    = dray(ia)   + db1


       b1 = y(ihe3)*aap*rate(irhep)

       qray(ihe3)  = qray(ihe3) - b1
       qray(ip)    = qray(ip)   - b1
       qray(ia)    = qray(ia)   + b1

!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = -db1
!       qasc(ip,icasc(ip)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db1
!       qasc(ia,icasc(ia)) = b1


! he3(n,g)he4
!       da1          = -aan*rate(irhe3ng)*y(ihe3) 
!       da2          = rate(irhe4gn)*aaa
!       db1          = da1 + da2 
!       dray(ihe3)  = dray(ihe3) + db1
!       dray(in)    = dray(in)   + db1
!       dray(ia)    = dray(ia)   - db1


       a1          = -aan*rate(irhe3ng)*y(ihe3) 
       a2          = rate(irhe4gn)*aaa
       b1          = a1 + a2 

       qray(ihe3)  = qray(ihe3) + b1
       qray(in)    = qray(in)   + b1
       qray(ia)    = qray(ia)   - b1


!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = db1
!       qasc(ihe3,icasc(ihe3)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -db1
!       qasc(ia,icasc(ia)) = -b1
      end if


! if we have he3 and deuterium

      if (ihe3 .ne. 0 .and. ih2 .ne. 0) then

! he3(d,p)he4
!       da1          = y(ihe3)*y(ih2)*rate(irhe3dp) 
!       da2          = -aap*aaa*rate(irhe4pd)
!       db1           = da1 + da2 
!       dray(ih2)   = dray(ih2)  - db1
!       dray(ihe3)  = dray(ihe3) - db1
!       dray(ia)    = dray(ia)   + db1
!       dray(ip)    = dray(ip)   + db1


       a1          = y(ihe3)*y(ih2)*rate(irhe3dp) 
       a2          = -aap*aaa*rate(irhe4pd)
       b1           = a1 + a2 

       qray(ih2)   = qray(ih2)  - b1
       qray(ihe3)  = qray(ihe3) - b1
       qray(ia)    = qray(ia)   + b1
       qray(ip)    = qray(ip)   + b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -db1
!       qasc(ih2,icasc(ih2)) = -b1
!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db1
!       qasc(ia,icasc(ia)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1


! d(d,n)he3
!       da1          = 0.5d0*y(ih2)*y(ih2)*rate(irddn) 
!       da2          = -aan*y(ihe3)*rate(irhe3nd)
!       db1          = da1 + da2 
!       dray(ih2)   = dray(ih2)  - 2.0d0 * db1
!       dray(ihe3)  = dray(ihe3) + db1
!       dray(in)    = dray(in)   + db1


       a1          = 0.5d0*y(ih2)*y(ih2)*rate(irddn) 
       a2          = -aan*y(ihe3)*rate(irhe3nd)
       b1          = a1 + a2 

       qray(ih2)   = qray(ih2)  - 2.0d0 * b1
       qray(ihe3)  = qray(ihe3) + b1
       qray(in)    = qray(in)   + b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -2.0d0*db1
!       qasc(ih2,icasc(ih2)) = -2.0d0*b1
!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = db1
!       qasc(ihe3,icasc(ihe3)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
      end if



! if we have he3, deuterium, tritium

      if (ihe3 .ne. 0  .and. ih2 .ne. 0  .and. ih3 .ne. 0) then

! he3(t,d)he4
!       da1          = y(ihe3)*y(ih3)*rate(irhe3td) 
!       da2          = -aaa*y(ih2)*rate(irhe4dt)
!       db1          = da1 + da2 
!       dray(ih2)  = dray(ih2)   + db1
!       dray(ih3)  = dray(ih3)   - db1
!       dray(ihe3) = dray(ihe3)  - db1
!       dray(ia)   = dray(ia)    + db1


       a1          = y(ihe3)*y(ih3)*rate(irhe3td) 
       a2          = -aaa*y(ih2)*rate(irhe4dt)
       b1          = a1 + a2 

       qray(ih2)  = qray(ih2)   + b1
       qray(ih3)  = qray(ih3)   - b1
       qray(ihe3) = qray(ihe3)  - b1
       qray(ia)   = qray(ia)    + b1

!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = db1
!       qasc(ih2,icasc(ih2)) = b1
!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = -db1
!       qasc(ih3,icasc(ih3)) = -b1
!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db1
!       qasc(ia,icasc(ia)) = b1
      end if


! if we have he3 and tritium
      if (ihe3 .ne. 0 .and. ih3 .ne. 0) then

! he3(t,np)he4
!       da1         = y(ih3)*y(ihe3)*rate(irhe3tnp) 
!       da2         = aan*aap*aaa*rate(irhe4pnt) 
!       db1         = da1 - da2
!       dray(ihe3) = dray(ihe3) - db1
!       dray(ih3)  = dray(ih3)  - db1
!       dray(ia)   = dray(ia)   + db1
!       dray(ip)   = dray(ip)   + db1
!       dray(in)   = dray(in)   + db1


       a1         = y(ih3)*y(ihe3)*rate(irhe3tnp) 
       a2         = aan*aap*aaa*rate(irhe4pnt) 
       b1         = a1 - a2

       qray(ihe3) = qray(ihe3) - b1
       qray(ih3)  = qray(ih3)  - b1
       qray(ia)   = qray(ia)   + b1
       qray(ip)   = qray(ip)   + b1
       qray(in)   = qray(in)   + b1


!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = -db1
!       qasc(ih3,icasc(ih3)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = db1
!       qasc(ia,icasc(ia)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
      end if




! if we have li7

      if (ili7 .ne. 0) then

! li7(p,g)be8=>2a + li7(p,a)a reactions
!       da1          = -y(ili7)*aap*rate(irli7pag) 
!       da2          = 0.5d0*aaa*aaa*rate(ir2he4ga)
!       db1          = da1 + da2 
!       dray(ili7)  = dray(ili7) + db1
!       dray(ip)    = dray(ip)   + db1
!       dray(ia)    = dray(ia)   - 2.0d0*db1


       a1          = -y(ili7)*aap*rate(irli7pag) 
       a2          = 0.5d0*aaa*aaa*rate(ir2he4ga)
       b1          = a1 + a2 

       qray(ili7)  = qray(ili7) + b1
       qray(ip)    = qray(ip)   + b1
       qray(ia)    = qray(ia)   - 2.0d0*b1

!       icasc(ili7) = icasc(ili7) + 1
!       casc(ili7,icasc(ili7)) = db1
!       qasc(ili7,icasc(ili7)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -2.0d0*db1
!       qasc(ia,icasc(ia)) = -2.0d0*b1
      end if


! if we have li7 and deuterium

      if (ili7 .ne. 0 .and. ih2 .ne. 0) then

! li7(d,n)2a
!       da1           = y(ili7)*y(ih2)*rate(irli7dn)
!       da2           = 0.5d0 * aan*aaa*aaa*rate(ir2and)
!       db1           = da1 - da2
!       dray(ih2)   = dray(ih2)  - db1
!       dray(ili7)  = dray(ili7) - db1
!       dray(in)    = dray(in)   + db1
!       dray(ia)    = dray(ia)   + 2.0d0*db1


       a1           = y(ili7)*y(ih2)*rate(irli7dn)
       a2           = 0.5d0 * aan*aaa*aaa*rate(ir2and)
       b1           = a1 - a2

       qray(ih2)   = qray(ih2)  - b1
       qray(ili7)  = qray(ili7) - b1
       qray(in)    = qray(in)   + b1
       qray(ia)    = qray(ia)   + 2.0d0*b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -db1
!       qasc(ih2,icasc(ih2)) = -b1
!       icasc(ili7) = icasc(ili7) + 1
!       casc(ili7,icasc(ili7)) = -db1
!       qasc(ili7,icasc(ili7)) = -b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if


! if we have li7 and tritium

      if (ili7 .ne. 0 .and. ih3 .ne. 0) then

! li7(t,2n)2a
!       da1 = y(ili7)*y(ih3)*rate(irli7t2n)
!       da2 = 0.25d0*aan*aan*aaa*aaa*rate(ir2a2nt)
!       db1 = da1 - da2
!       dray(ih3)   = dray(ih3)  - db1
!       dray(ili7)  = dray(ili7) - db1
!       dray(in)    = dray(in)   + 2.0d0*db1
!       dray(ia)    = dray(ia)   + 2.0d0*db1



       a1 = y(ili7)*y(ih3)*rate(irli7t2n)
       a2 = 0.25d0*aan*aan*aaa*aaa*rate(ir2a2nt)
       b1 = a1 - a2

       qray(ih3)   = qray(ih3)  - b1
       qray(ili7)  = qray(ili7) - b1
       qray(in)    = qray(in)   + 2.0d0*b1
       qray(ia)    = qray(ia)   + 2.0d0*b1

!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = -db1
!       qasc(ih3,icasc(ih3)) = -b1
!       icasc(ili7) = icasc(ili7) + 1
!       casc(ili7,icasc(ili7)) = -db1
!       qasc(ili7,icasc(ili7)) = -b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = 2.0d0*db1
!       qasc(in,icasc(in)) = 2.0d0*b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if



! if we have li7 and he3 

      if (ili7 .ne. 0 .and. ihe3 .ne. 0) then

! li7(he3,np)2a
!       da1 = y(ili7)*y(ihe3)*rate(irli7he3np)
!       da2 = 0.5d0*aan*aap*aaa*aaa*rate(ir2apnhe3)
!       db1 = da1 - da2
!       dray(ihe3)  = dray(ihe3) - db1
!       dray(ili7)  = dray(ili7) - db1
!       dray(in)    = dray(in)   + db1
!       dray(ip)    = dray(ip)   + db1
!       dray(ia)    = dray(ia)   + 2.0d0*db1


       a1 = y(ili7)*y(ihe3)*rate(irli7he3np)
       a2 = 0.5d0*aan*aap*aaa*aaa*rate(ir2apnhe3)
       b1 = a1 - a2

       qray(ihe3)  = qray(ihe3) - b1
       qray(ili7)  = qray(ili7) - b1
       qray(in)    = qray(in)   + b1
       qray(ip)    = qray(ip)   + b1
       qray(ia)    = qray(ia)   + 2.0d0*b1


!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ili7) = icasc(ili7) + 1
!       casc(ili7,icasc(ili7)) = -db1
!       qasc(ili7,icasc(ili7)) = -b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if



! if we have be7 and deuterium

      if (ibe7 .ne. 0 .and. ih2 .ne. 0) then

! be7(d,p)2a
!       da1 = y(ibe7)*y(ih2)*rate(irbe7dp)
!       da2 = 0.5d0*aap*aaa*aaa*rate(ir2apd)
!       db1 = da1 - da2
!       dray(ih2)   = dray(ih2)  - db1
!       dray(ibe7)  = dray(ibe7) - db1
!       dray(ip)    = dray(ip)   + db1
!       dray(ia)    = dray(ia)   + 2.0d0*db1


       a1 = y(ibe7)*y(ih2)*rate(irbe7dp)
       a2 = 0.5d0*aap*aaa*aaa*rate(ir2apd)
       b1 = a1 - a2

       qray(ih2)   = qray(ih2)  - b1
       qray(ibe7)  = qray(ibe7) - b1
       qray(ip)    = qray(ip)   + b1
       qray(ia)    = qray(ia)   + 2.0d0*b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = -db1
!       qasc(ih2,icasc(ih2)) = -b1
!       icasc(ibe7) = icasc(ibe7) + 1
!       casc(ibe7,icasc(ibe7)) = -db1
!       qasc(ibe7,icasc(ibe7)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if


! if we have be7 and tritium 

      if (ibe7 .ne. 0 .and. ih3 .ne. 0) then

! be7(t,np)2a
!       da1 = y(ibe7)*y(ih3)*rate(irbe7tnp)
!       da2 = 0.5d0*aan*aap*aaa*aaa*rate(ir2apnt)
!       db1 = da1 - da2 
!       dray(ih3)   = dray(ih3)  - db1
!       dray(ibe7)  = dray(ibe7) - db1
!       dray(ip)    = dray(ip)   + db1
!       dray(in)    = dray(in)   + db1
!       dray(ia)    = dray(ia)   + 2.0d0*db1


       a1 = y(ibe7)*y(ih3)*rate(irbe7tnp)
       a2 = 0.5d0*aan*aap*aaa*aaa*rate(ir2apnt)
       b1 = a1 - a2 

       qray(ih3)   = qray(ih3)  - b1
       qray(ibe7)  = qray(ibe7) - b1
       qray(ip)    = qray(ip)   + b1
       qray(in)    = qray(in)   + b1
       qray(ia)    = qray(ia)   + 2.0d0*b1


!       icasc(ih3) = icasc(ih3) + 1
!       casc(ih3,icasc(ih3)) = -db1
!       qasc(ih3,icasc(ih3)) = -b1
!       icasc(ibe7) = icasc(ibe7) + 1
!       casc(ibe7,icasc(ibe7)) = -db1
!       qasc(ibe7,icasc(ibe7)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = db1
!       qasc(in,icasc(in)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if


! if we have be7 and he3 

      if (ibe7 .ne. 0 .and. ihe3 .ne. 0) then

! be7(he3,2p)2a
!       da1          = y(ibe7)*y(ihe3)*rate(irbe7he32p)
!       da2          = 0.25d0*aap*aap*aaa*aaa*rate(ir2a2phe3)
!       db1          = da1 - da2
!       dray(ihe3)  = dray(ihe3) - db1
!       dray(ibe7)  = dray(ibe7) - db1
!       dray(ip)    = dray(ip)   + 2.0d0*db1
!       dray(ia)    = dray(ia)   + 2.0d0*db1


       a1          = y(ibe7)*y(ihe3)*rate(irbe7he32p)
       a2          = 0.25d0*aap*aap*aaa*aaa*rate(ir2a2phe3)
       b1          = a1 - a2

       qray(ihe3)  = qray(ihe3) - b1
       qray(ibe7)  = qray(ibe7) - b1
       qray(ip)    = qray(ip)   + 2.0d0*b1
       qray(ia)    = qray(ia)   + 2.0d0*b1


!       icasc(ihe3) = icasc(ihe3) + 1
!       casc(ihe3,icasc(ihe3)) = -db1
!       qasc(ihe3,icasc(ihe3)) = -b1
!       icasc(ibe7) = icasc(ibe7) + 1
!       casc(ibe7,icasc(ibe7)) = -db1
!       qasc(ibe7,icasc(ibe7)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = 2.0d0*db1
!       qasc(ip,icasc(ip)) = 2.0d0*b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if


! if we have be9
      if (ibe9 .ne. 0) then

! a(an,g)be9
!       da1         = 0.5d0*aaa*aaa*aan*rate(iraan) 
!       da2         = y(ibe9)*rate(irgaan)
!       db1         = da1 - da2 
!       dray(ibe9) = dray(ibe9) + db1
!       dray(in)   = dray(in)   - db1
!       dray(ia)   = dray(ia) - 2.0d0 * db1


       a1         = 0.5d0*aaa*aaa*aan*rate(iraan) 
       a2         = y(ibe9)*rate(irgaan)
       b1         = a1 - a2 

       qray(ibe9) = qray(ibe9) + b1
       qray(in)   = qray(in)   - b1
       qray(ia)   = qray(ia) - 2.0d0 * b1


!       icasc(ibe9) = icasc(ibe9) + 1
!       casc(ibe9,icasc(ibe9)) = db1
!       qasc(ibe9,icasc(ibe9)) = b1
!       icasc(in) = icasc(in) + 1
!       casc(in,icasc(in)) = -db1
!       qasc(in,icasc(in)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -2.0d0*db1
!       qasc(ia,icasc(ia)) = -2.0d0*b1
      end if


! if we have be9 and deuterium

      if (ibe9 .ne. 0  .and. ih2 .ne. 0) then

! be9(p,d)be8 =>2a
!       da1           = y(ibe9)*aap*rate(irbe9pd)
!       da2           = 0.5d0*y(ih2)*aaa*aaa*rate(ir2adp)
!       db1           = da1 - da2
!       dray(ibe9)  = dray(ibe9)  - db1
!       dray(ip)    = dray(ip)    - db1
!       dray(ih2)   = dray(ih2)   + db1
!       dray(ia)    = dray(ia)    + 2.0d0*db1


       a1           = y(ibe9)*aap*rate(irbe9pd)
       a2           = 0.5d0*y(ih2)*aaa*aaa*rate(ir2adp)
       b1           = a1 - a2

       qray(ibe9)  = qray(ibe9)  - b1
       qray(ip)    = qray(ip)    - b1
       qray(ih2)   = qray(ih2)   + b1
       qray(ia)    = qray(ia)    + 2.0d0*b1


!       icasc(ih2) = icasc(ih2) + 1
!       casc(ih2,icasc(ih2)) = db1
!       qasc(ih2,icasc(ih2)) = b1
!       icasc(ibe9) = icasc(ibe9) + 1
!       casc(ibe9,icasc(ibe9)) = -db1
!       qasc(ibe9,icasc(ibe9)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = -db1
!       qasc(ip,icasc(ip)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if


! if we have b8

      if (ib8 .ne. 0) then

! b8(p=>n)be8 =>2a reactions
!       db1 = y(ib8)*rate(irb8ep)
!       dray(ib8)   = dray(ib8) - db1
!       dray(ia)    = dray(ia) + 2.0d0*db1


       b1 = y(ib8)*rate(irb8ep)

       qray(ib8)   = qray(ib8) - b1
       qray(ia)    = qray(ia) + 2.0d0*b1


!       icasc(ib8) = icasc(ib8) + 1
!       casc(ib8,icasc(ib8)) = -db1
!       qasc(ib8,icasc(ib8)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 2.0d0*db1
!       qasc(ia,icasc(ia)) = 2.0d0*b1
      end if


! if we have b11
      if (ib11 .ne. 0) then

! b11(p,a)be8 => 2a reactions
!       da1          = -y(ib11)*aap*rate(irb11pa) 
!       da2          = sixth*aaa*aaa*aaa*rate(ir3ap)
!       db1          = da1 + da2 
!       dray(ib11)  = dray(ib11) + db1
!       dray(ip)    = dray(ip)   + db1
!       dray(ia)    = dray(ia)   - 3.0d0*db1

       a1          = -y(ib11)*aap*rate(irb11pa) 
       a2          = sixth*aaa*aaa*aaa*rate(ir3ap)
       b1          = a1 + a2 

       qray(ib11)  = qray(ib11) + b1
       qray(ip)    = qray(ip)   + b1
       qray(ia)    = qray(ia)   - 3.0d0*b1


!       icasc(ib11) = icasc(ib11) + 1
!       casc(ib11,icasc(ib11)) = db1
!       qasc(ib11,icasc(ib11)) = b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = db1
!       qasc(ip,icasc(ip)) = b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = -3.0d0*db1
!       qasc(ia,icasc(ia)) = -3.0d0*b1
      end if


! if we have c11
      if (ic11 .ne. 0) then

! c11(n,a)be8 => 2a
!       db1          = aan*y(ic11)*rate(irc11na)
!       dray(ic11)  = dray(ic11) - db1
!       dray(in)    = dray(in)   - db1
!       dray(ia)    = dray(ia)   + 3.0d0*db1


       b1          = aan*y(ic11)*rate(irc11na)

       qray(ic11)  = qray(ic11) - b1
       qray(in)    = qray(in)   - b1
       qray(ia)    = qray(ia)   + 3.0d0*b1


!       icasc(ic11) = icasc(ic11) + 1
!       casc(ic11,icasc(ic11)) = -db1
!       qasc(ic11,icasc(ic11)) = -b1
!       icasc(ip) = icasc(ip) + 1
!       casc(ip,icasc(ip)) = -db1
!       qasc(ip,icasc(ip)) = -b1
!       icasc(ia) = icasc(ia) + 1
!       casc(ia,icasc(ia)) = 3.0d0*db1
!       qasc(ia,icasc(ia)) = 3.0d0*b1
      end if


! cascade sums
! double, sorted
!       do i=1,ionmax
!        if (icasc(i) .gt. cdim2) stop 'cdim2 too small'
!        do j=1,icasc(i)
!         cdum(i) = abs(casc(i,j))
!        end do
!        call indexx(icasc(i),cdum,indx)
!        cdum(i) = 0.0d0
!        do j=1,icasc(i)
!         cdum(i) = cdum(i) + casc(i,indx(j))
!!         cdum(i) = cdum(i) + casc(i,j)
!        end do
!       end do


! double, kaham sorted
!       do i=1,ionmax
!        if (icasc(i) .gt. cdim2) stop 'cdim2 too small'
!        do j=1,icasc(i)
!         kdum(i) = abs(casc(i,j))
!        end do
!        call indexx(icasc(i),kdum,indx)
!        kdum(i) = 0.0d0
!        do j=1,icasc(i)
!           zwork1(j) = casc(i,indx(j))
!           zwork1(j) = casc(i,j)
!        end do
!        kdum(i) = ksum(icasc(i),zwork1)
!       end do


! quad, sorted
!       do i=1,ionmax
!        if (icasc(i) .gt. cdim2) stop 'cdim2 too small'
!        do j=1,icasc(i)
!         qdum(i) = abs(qasc(i,j))
!        end do
!        call indexx16(icasc(i),qdum,indx)
!        qdum(i) = 0.0q0
!        do j=1,icasc(i)
!         qdum(i) = qdum(i) + qasc(i,indx(j))
!!         qdum(i) = qdum(i) + qasc(i,j)
!        end do
!       end do


! mp, sorted
!       do i=1,ionmax
!!       do i=ionmax-2,ionmax
!        if (icasc(i) .gt. cdim2) stop 'cdim2 too small'
!        do j=1,icasc(i)
!         qdum(i) = abs(qasc(i,j))
!        end do
!        call indexx16(icasc(i),qdum,indx)
!        xs(1:64) = 0
!        iterm(1:64) = 0
!        do j=1,icasc(i)
!         write(string,01) qasc(i,indx(j)),'$'
!01       format(1pe38.30,a1)
!         call mpcam(string,iterm)
!         call mpadd(xs, iterm, xs)
!        end do
!
!        call mpcmd(xs,mpdum(i))
!
!        if (izwork1(1) .eq. -999) then
!         if (i .eq. in-1) then
!          write(6,332) 
!332       format(4x,'12345678901234567890123456789012345678901234567890123456789012345678901234567890')
!         end if
!
!         if (i .eq. in  .or. i .eq. ip .or. i .eq. ia) then
!!          call mpfout(xs,64)
!          call mpout (xs, c, 110, 90)
!!          write (lun, 333) c
!!   333    format (1x, 110a1)
!
!          write(string,334) c
!334       format(110a1)
!          call sqeeze(string)
!
!          mpdum(i) = value16(string)         
!          write(6,'(1pe40.32)') mpdum(i)
!
!         end if
!        end if
!       end do



! now set the double precision return argument dydt

!      dydt(1:ionmax) = dray(1:ionmax)
!      dydt(1:ionmax) = cdum(1:ionmax)
!      dydt(1:ionmax) = kdum(1:ionmax)
!      dydt(1:ionmax) = qdum(1:ionmax)
!      dydt(1:ionmax) = mpdum(1:ionmax)
      dydt(1:ionmax) = qray(1:ionmax)


! n,p,a dominate sums, but here is the checker
!      bigdiff = maxval(abs(qray(1:ionmax)-dray(1:ionmax)))
!      ibigdiff = maxloc(abs(qray(1:ionmax)-dray(1:ionmax)))
!      write(6,466) 'max difference: ',ionam(ibigdiff),bigdiff,qray(ibigdiff),dray(ibigdiff)
!466   format(1x,a,a,1p4e38.30)


!      if (izwork1(1) .eq. -999) then
!      icount = icount + 1

!!       write(6,466) icount
!!466    format(1x,i6,'    12345678901234567890123456789012',8x,'12345678901234567890123456789012',8x,'12345678901234567890123456789012')

!       write(6,467) icount,qray(in),qray(ip),qray(ia)
!467    format(1x,i6,1p3e40.32)

!!       write(6,467) icount,dray(in),dray(ip),dray(ia)
!!       write(6,467) icount,cdum(in),cdum(ip),cdum(ia)
!!       write(6,467) icount,kdum(in),kdum(ip),kdum(ia)
!        write(6,467) icount,mpdum(in),mpdum(ip),mpdum(ia)
!!468    format(1x,i6,1pe24.16,16x,1pe24.16,16x,1pe24.16)


!!       write(6,467) icount, &
!!                    abs(qray(in)-dray(in)), &
!!                    abs(qray(ip)-dray(ip)), & 
!!                    abs(qray(ia)-dray(ia))

!!       write(6,467) icount, &
!!                    abs((qray(in)-dray(in))/qray(in)), &
!!                    abs((qray(ip)-dray(ip))/qray(ip)), &
!!                    abs((qray(ia)-dray(ia))/qray(ia))

!       write(6,469) icount, &
!                    int(abs(log10(abs((qray(in)-mpdum(in))/qray(in))))), &
!                    int(abs(log10(abs((qray(ip)-mpdum(ip))/qray(ip))))), &
!                    int(abs(log10(abs((qray(ia)-mpdum(ia))/qray(ia)))))
!469    format(1x,i6,i4,i4,i4)


!       write(16,469) icount, &
!                    int(abs(log10(abs((qray(in)-mpdum(in))/qray(in))))), &
!                    int(abs(log10(abs((qray(ip)-mpdum(ip))/qray(ip))))), &
!                    int(abs(log10(abs((qray(ia)-mpdum(ia))/qray(ia)))))

!        write(6,*) 
!       end if


      return
      end subroutine rhstorch

      subroutine read_partition_functions
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
!
! this routine reads and initializes the rauscher & thielemann partition function data
!
! declare

      character*132    :: string
      integer          :: i,j,izpf,iapf,inpf
      double precision :: xx


! initialize
      pf(1:pf_tmax,1:pf_zmax,1:pf_nmax) = -101.0d0 
      pf_gspin(1:pf_zmax,1:pf_nmax)     = -101.0d0

! open the partition function data file and read the header
      open(unit=22, file='partition_function_frdm.dat', status='old')
      do i=1,4
       read(22,11) string
11     format(a)
      enddo

! read the file
      do i=1,40000
       read(22,11,end=30) string
       if (string(1:12) .eq. 'END OF TABLE') goto 30
       read(22,*) izpf,iapf,xx
       inpf = iapf - izpf
       if (izpf .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (inpf .lt. 1) stop 'inpf < 1 during read of partition functions'
       if (izpf .gt. pf_zmax) stop 'izpf > pf_izmax during read of partition functions'
       if (inpf .gt. pf_nmax) stop 'inpf > pf_inmax during read of partition functions'
       pf_gspin(izpf,inpf) = xx
       read(22,*) (pf(j,izpf,inpf), j=1,24)
      end do

! close up
30    close(unit=22)


! set the temperatures on which the partition functions are calculated
      pf_t9 = (/0.01d0, 0.15d0,  0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, &
                 0.8d0,  0.9d0,  1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, & 
                 4.0d0,  4.5d0,  5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0, 10.0d0/)

      return
    end subroutine read_partition_functions

      subroutine get_partition_functions(t9)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'vector_eos.dek'
      include 'network.dek'


! this routine generates rauscher & thielemann partition function values

! declare the pass
     double precision :: t9


! declare locals
! for the partition function interpolation
      integer            :: i,jat
      integer, parameter :: norder=2
      double precision   :: x,x1,x2,a,b,e,alfa,beta,stat_weight


! waprt(i) is the temperature dependent partition functions
! zwork1(i) is the ratio of wpart(i) to the ground state partition function
! zwork2(i) is the temperature derivative of zwork1
! initialize


      zwork1(ionbeg:ionend) = 1.0d0
      zwork2(ionbeg:ionend) = 0.0d0
      wpart(ionbeg:ionend)  = 1.0d0

      if (ineut .ne. 0) then
       zwork1(ineut) = 1.0d0
       zwork2(ineut) = 0.0d0
       wpart(ineut)  = 2.0d0
      end if

      if (iprot .ne. 0) then
       zwork1(iprot) = 1.0d0
       zwork2(iprot) = 0.0d0
       wpart(iprot)  = 2.0d0
      end if

      if (ih1 .ne. 0) then
       zwork1(ih1) = 1.0d0
       zwork2(ih1) = 0.0d0
       wpart(ih1)  = 2.0d0
      end if

      if (ihe4 .ne. 0) then
       zwork1(ihe4) = 1.0d0
       zwork2(ihe4) = 0.0d0
       wpart(ihe4)  = 1.0d0
      end if



! for every isotope not named n, p or alpha

      do i=ionbeg,ionend

! bullet check
       if (int(zion(i)) .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (int(nion(i)) .lt. 1) stop 'izpf < 1 during read of partition functions'
       if (int(zion(i)) .gt. pf_zmax) stop 'izpf > pf_izmax during read of partition functions'
       if (int(nion(i)) .gt. pf_nmax) stop 'inpf > pf_inmax during read of partition functions'


! if isotope is in the network 
       if (pf_gspin(int(zion(i)),int(nion(i))) .lt. 0.0d0) cycle

!       write(6,*) int(zion(i)),int(nion(i)),pf_gspin(int(zion(i)),int(nion(i)))


        stat_weight = 2.0d0 * pf_gspin(int(zion(i)),int(nion(i))) + 1.0d0

! clamp at both ends of the temperature grid
       if (t9 .le. pf_t9(1)) then
        zwork1(i) = pf(1,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

       else if (t9 .ge. pf_t9(pf_tmax)) then
        zwork1(i) = pf(pf_tmax,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

! linear larange interpolation

       else
        call locate(pf_t9,pf_tmax,t9,jat)
        jat = max(1,min(jat - norder/2 + 1,pf_tmax - norder + 1))

!        write(6,*) jat,t9
!        write(6,*) pf_t9(jat),pf_t9(jat+1)
!        write(6,*) pf(jat,int(zion(i)),int(nion(i))),pf(jat+1,int(zion(i)),int(nion(i)))

        x  = t9
        x1 = pf_t9(jat)
        x2 = pf_t9(jat+1)
        a  = x - x1
        b  = x - x2
        e  = x1 - x2
        alfa =  b/e
        beta = -a/e
        zwork1(i) = alfa * pf(jat,int(zion(i)),int(nion(i))) + beta * pf(jat+1,int(zion(i)),int(nion(i)))
        zwork2(i) = 0.0d0
        wpart(i)  = zwork1(i) * stat_weight

! no partition functions
!       zwork1(i) = 1.0d0
!       zwork2(i) = 0.0d0
!       wpart(i)  = 1.0d0



       end if

!       write(6,156) i, ionam(i), zwork1(i), wpart(i)
!156    format(1x,i4,' ',a,' ',1p2e12.4)
!       read(5,*)

      enddo
      return
    end subroutine get_partition_functions

      subroutine bdat_partition_functions
      include 'implno.dek'
      include 'burn_common.dek'
      include 'vector_eos.dek'
      include 'network.dek'
      include 'tfactors.dek'


! this routine generates partition function values from bdat data

! declare locals
      integer            :: i,jd1,jd2,jd3,jd4,jd5
      double precision   :: g0,dg0,aa,daa


! wpart(i) is the temperature dependent partition functions
! zwork1(i) is the ratio of wpart(i) to the ground state partition function
! zwork2(i) is the temperature derivative of zwork1


! initialize
      zwork1(ionbeg:ionend) = 1.0d0
      zwork2(ionbeg:ionend) = 0.0d0
      wpart(ionbeg:ionend)  = 1.0d0


! n,p,a 
      if (ineut .ne. 0) then
       zwork1(ineut) = 1.0d0
       zwork2(ineut) = 0.0d0
       wpart(ineut)  = 2.0d0
      end if

      if (iprot .ne. 0) then
       zwork1(iprot) = 1.0d0
       zwork2(iprot) = 0.0d0
       wpart(iprot)  = 2.0d0
      end if

      if (ih1 .ne. 0) then
       zwork1(ih1) = 1.0d0
       zwork2(ih1) = 0.0d0
       wpart(ih1)  = 2.0d0
      end if

      if (ihe4 .ne. 0) then
       zwork1(ihe4) = 1.0d0
       zwork2(ihe4) = 0.0d0
       wpart(ihe4)  = 1.0d0
      end if



! for every isotope not named n, p or alpha

      do i=ionbeg,ionend
       jd1  = 5*(i-1) + 1
       jd2  = jd1 + 1
       jd3  = jd2 + 1
       jd4  = jd3 + 1
       jd5  = jd4 + 1
       g0   = 1.0d0
       dg0  = 0.0d0

       if (as(jd2) .ne. 0.0) then
        aa        = as(jd2)*t9i + as(jd3) + as(jd4)*t9 + as(jd5)*t92
        daa       = -as(jd2)*t9i2 + as(jd4) + 2.0d0*as(jd5)*t9
        zwork1(i) = exp(aa)
        zwork2(i) = zwork1(i)*daa
       end if

       zwork1(i) = g0  + zwork1(i)
       zwork2(i) = dg0 + zwork2(i)
       wpart(i)  = as(jd1) * zwork1(i)

!       write(6,156) i, ionam(i), zwork1(i), wpart(i)
!156    format(1x,i4,' ',a,' ',1p2e12.4)
!       read(5,*)

      enddo
      return
    end subroutine bdat_partition_functions
      


      subroutine rate_aan(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd

! he4(an,g)be9
      aa  = 1.0d0 + 0.344*t9
      bb  = t92 * aa
      dbb = 2.0d0 * t9 * aa + t92*0.344

      cc  = 1.0d0/bb
      dcc = -cc*cc*dbb

      dd  = 2.59e-6 * exp(-1.062*t9i)
      ddd = dd*1.062*t9i2

      term    = cc * dd
      dtermdt = dcc*dd + cc*ddd

! rates

      fr    = den * den * term
      dfrdt = den * den * dtermdt * 1.0d-9
      dfrdd = 2.0d0 * den * term

      rev      = 5.84e19 * t93 * exp(-18.260*t9i)
      drevdt   = rev*(3.0d0*t9i + 18.260*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9pd(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/0.2704d0)


! be9(p,d)be8 =>2a
      aa  = 2.11e+11 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0  + 0.04*t913 + 1.09*t923 + 0.307*t9 &
            + 3.21*t943 + 2.30*t953
      dbb = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307 &
            + fourth*3.21*t913 + fiveth*2.30*t923

      cc  = 5.79e+08 * t9i * exp(-3.046*t9i)
      dcc = cc*(-t9i + 3.046*t9i2)

      dd  = 8.50e+08 * t9i34 * exp(-5.800*t9i)
      ddd = dd*(-0.75d0*t9i + 5.800*t9i2)

      term    = aa*bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.07e-11 * t9i32 *exp(-7.555*t9i)
      drevdt   = rev*(-1.5d0*t9i + 7.555*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end






      subroutine rate_be9pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


! be9(p,n)b9
      aa  = 5.58e7*(1.0d0 + 0.042*t912 + 0.985*t9)
      daa = 5.58e7*(0.5d0*0.042*t9i12 + 0.985)

      zz  = exp(-21.473*t9i)
      dzz = zz*21.473*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 1.02e+09 * t9i32 * exp(-26.725*t9i)
      dcc = cc*(-1.5d0*t9i + 26.725*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      bb  = 0.998 * aa
      dbb = 0.998 * daa

      cc  = 0.998 * 1.02e+09 * t9i32 * exp(-5.252*t9i)
      dcc = cc*(-1.5d0*t9i + 5.252*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b10an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,q1
      parameter        (q1 = 1.0d0/91.948921d0)


! b10(a,n)n13
      term    = 1.20e+13 * t9i23 * exp(-27.989*t9i13 - t92*q1)
      dtermdt = -twoth*term*t9i &
                + term*(oneth*27.989*t9i43 - 2.0d0*t9*q1)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.34 * exp(-12.287*t9i)
      drevdt   = rev*12.287*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b11an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.0196d0)


! b11(a,n)n14
      aa  = 6.97e+12 * t9i23 * exp(-28.234*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.015*t913 + 8.115*t923 + 0.838*t9 &
            + 39.804*t943 + 10.456*t953
      dbb = oneth*0.015*t9i23 + twoth*8.115*t9i13 + 0.838 &
            + fourth*39.804*t913 + fiveth*10.456*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.79 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 1.71e+03 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 4.49e+06 * t935 * exp(-8.596*t9i)
      dff  = ff*(0.6d0*t9i + 8.596*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.67 * exp(-1.835*t9i)
      drevdt   = rev*1.835*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_n14pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! n14(p,a)b11
       aa     = 1.0d0 + 0.0478*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0478

       zz     = 1.0d0/bb
       cc     = aa + 7.56e-03*t953*zz
       dcc    = 0.0478 + (fiveth*7.56e-3*t923 - 7.56e-3*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix*t9a56*zz

       dd      = 2.63e+16 * t9a56 * t9i32 * exp(-31.883/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 31.883/t9a13**2 * dt9a13)

       term    = dd * exp(-33.915*t9i)
       dtermdt = term*(ddd/dd + 33.915*t9i2)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.272 * dd
      drevdt   = 0.272 * ddd

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_c11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/2.647129d0)


! c11(p,g)n12
      aa  = 4.24e+04 * t9i23 * exp(-13.658*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*13.658*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0  + 0.031*t913 + 3.11*t923 + 0.665*t9 &
            + 4.61*t943 + 2.50*t953
      dbb = oneth*0.031*t9i23 + twoth*3.11*t9i13 + 0.665 &
            + fourth*4.61*t913 + fiveth*2.50*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 8.84e+03 * t9i32 * exp(-7.021*t9i)
      ddd  = dd*(-1.5d0*t9i + 7.021*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.33e+10 * t932 * exp(-6.975*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.975*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! c12(a,n)o15
      aa  = 2.48e7 * (1.0d0 + 0.188*t912 + 0.015*t9)
      daa = 2.48e7 * (0.5d0*0.188*t9i12 + 0.015)

      bb  = exp(-98.661*t9i)
      dbb = bb*98.661*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 1.41 * aa
      drrdt = den * 1.41 * daa  * 1.0d-9
      drrdd = 1.41 * aa

      return
      end





      subroutine rate_c13an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.648656d0)


! c13(a,n)o16
      aa  = 6.77e+15 * t9i23 * exp(-32.329*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*32.329*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.013*t913 + 2.04*t923 + 0.184*t9
      dbb = oneth*0.013*t9i23 + twoth*2.04*t9i13 + 0.184

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.82e+05 * t9i32 * exp(-9.373*t9i)
      ddd  = dd*(-1.5d0*t9i + 9.373*t9i2)

      ee   = 1.41e+06 * t9i32 * exp(-11.873*t9i)
      dee  = ee*(-1.5d0*t9i + 11.873*t9i2)

      ff   = 2.0e+09 * t9i32 * exp(-20.409*t9i)
      dff  = ff*(-1.5d0*t9i + 20.409*t9i2)

      gg   = 2.92e+09 * t9i32 * exp(-29.283*t9i)
      dgg  = gg*(-1.5d0*t9i + 29.283*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.79e+00 * exp(-25.711*t9i)
      drevdt = rev*25.711*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_c13pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! c13(p,n)n13
      aa  = 1.88e+08*(1.0d0 - 0.167*t912 + 0.037*t9)
      daa = 1.88e+08*(0.037 - 0.5d0*0.167*t9i12)

      bb  = exp(-34.846*t9i)
      dbb = bb*34.846*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 0.998 * aa
      drrdt = den * 0.998 * daa * 1.0d-9
      drrdd = 0.998 * aa

      return
      end






      subroutine rate_c14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/7.086244d0)


! c14(a,g)o18
      aa  = 1.528e+09 * t9i23 * exp(-32.513*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*32.513*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.0128*t913 - 0.869*t923 - 0.0779*t9 &
            + 0.321*t943 + 0.0732*t953
      dbb = oneth*0.0128*t9i23 - twoth*0.869*t9i13 - 0.0779 &
            + fourth*0.321*t913 + fiveth*0.0732*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.375e+08 * t9i2 * exp(-32.513*t9i13)
      ddd  = dd*(-2.0d0*t9i + oneth*32.513*t9i43)

      ee   = 9.29e-08 * t9i32 * exp(-2.048*t9i)
      dee  = ee*(-1.5d0*t9i + 2.048*t9i2)

      ff   = 2.77e+03 * t9i45 * exp(-9.876*t9i)
      dff  = ff*(-0.8d0*t9i + 9.876*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-72.262*t9i)
      drevdt   = rev*(1.5d0*t9i + 72.262*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


! c14(p,n)n14
      aa  = 7.19e+05*(1.0d0 + 0.361*t912 + 0.502*t9)
      daa = 7.19e+05*(0.5d0*0.361*t9i12 + 0.502)

      zz  = exp(-7.263*t9i)
      dzz = zz*7.263*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 3.34e+08 * t9i12 * exp(-12.246*t9i)
      dcc = cc*(-0.5d0*t9i + 12.246*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      cc  = 3.34e+08 * t9i12 * exp(-4.983*t9i)
      dcc = cc*(-0.5d0*t9i + 4.983*t9i2)

      rr    = den * 0.333 * (aa + cc)
      drrdt = den * 0.333 * (daa + dcc) * 1.0d-9
      drrdd = 0.333 * (aa + cc)

      return
      end






      subroutine rate_c14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/32.729841d0)


! c14(p,g)n14
      aa  = 6.80e+06 * t9i23 * exp(-13.741*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*13.741*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.03*t913 + 0.503*t923 + 0.107*t9 &
            + 0.213*t943 + 0.115*t953
      dbb = oneth*0.03*t9i23 + twoth*0.503*t9i13 + 0.107 &
            + fourth*0.213*t913 + fiveth*0.115*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.36e+03 * t9i32 * exp(-3.811*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.811*t9i2)

      ee   = 9.82e+04 * t9i13 * exp(-4.739*t9i)
      dee  = ee*(-oneth*t9i + 4.739*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.00e+09 * t932 * exp(-118.452*t9i)
      drevdt   = rev*(1.5d0*t9i + 118.452*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! o16(p,a)n13
       aa     = 1.0d0 + 0.0776*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0776

       zz     = 1.0d0/bb
       cc     = aa + 0.0264*t953*zz
       dcc    = 0.0776 + (fiveth*0.0264*t923 - 0.0264*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd      = 1.88e+18 * t9a56 * t9i32 * exp(-35.829/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 35.829/t9a13**2 * dt9a13)

       term    = dd * exp(-60.561*t9i)
       dtermdt = term*(ddd/dd + 60.561*t9i2)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.172 * dd
      drevdt   = 0.172 * ddd

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_n14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! n14(p,n)o14
      aa  = 6.74e+07 * (1.0d0 + 0.658*t912 + 0.379*t9)
      daa = 6.74e+07 * (0.5d0*0.658*t9i12 + 0.379)

      bb  = exp(-68.762*t9i)
      dbb = bb*68.762*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 2.99 * aa
      drrdt = den * 2.99 * daa * 1.0d-9
      drrdd = 2.99 * aa

      return
      end





      subroutine rate_n14an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz,q1
      parameter        (q1 = 1.0d0/7.828804d0)


! n14(a,n)f17
      aa  = 5.24e9*(1.0d0 - 1.15*t912 + 0.365*t9)
      daa = 5.24e9*(0.365 - 0.5d0*1.15*t9i12)

      zz  = exp(-t92*q1)
      dzz = -zz*2.0d0*t9*q1

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc   = 3.28e10 * t9i32 * exp(-1.5766e1*t9i)
      dcc  = cc*(-1.5d0*t9i + 1.5766e1*t9i2)

      term     = (bb + cc) * exp(-54.942*t9i)
      dtermdt  = term*((dbb+dcc)/(bb+cc) + 54.942*t9i2)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 1.48 * (bb + cc)
      dtermdt  = 1.48 * (dbb + dcc)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,aa,daa,bb,dbb


! n15(p,n)o15
      t9a = min(t9,10.0d0)
      aa  = 3.51e+08 * (1.0d0 + 0.452*t912 - 0.191*t9a)
      if (t9a .eq. 10.0) then
       daa = 3.51e+08 * 0.5d0*0.452*t9i12
      else
       daa = 3.51e+08 * (0.5d0*0.452*t9i12 - 0.191)
      end if

      bb  = exp(-41.032*t9i)
      dbb = bb*41.032*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! n15(a,n)f18
      aa  = 3.14e8 * (1.0d0 - 0.641*t912 + 0.108*t9)
      daa = 3.14e8 * (0.108 - 0.5d0*0.641*t9i12)

      bb  = exp(-74.479*t9i)
      dbb = bb*74.479*t9i2

      term     = aa * bb
      dtermdt  = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 2.0d0 * aa
      dtermdt  = 2.0d0 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_n15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.379456d0)


! n15(a,g)f19
      aa  = 2.54e+10 * t9i23 * exp(-36.211*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*36.211*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.012*t913 + 1.69*t923 + 0.136*t9 &
            + 1.91*t943 + 0.391*t953
      dbb = oneth*0.012*t9i23 + twoth*1.69*t9i13 + 0.136 &
            + fourth*1.91*t913 + fiveth*0.391*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.83e-03 * t9i32 * exp(-4.232*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.232*t9i2)

      ee   = 1.52e+03 * t9 * exp(-9.747*t9i)
      dee  = ee*(t9i + 9.747*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.54e+10 * t932 * exp(-46.578*t9i)
      drevdt = rev*(1.5d0*t9i + 46.578*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.514089d0)


! o14(a,g)ne18
      aa  = 9.47e+08 * t9i23 * exp(-39.388*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*39.388*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 + 1.974*t923 + 0.146*t9 &
            + 3.036*t943 + 0.572*t953
      dbb = oneth*0.011*t9i23 + twoth*1.974*t9i13 + 0.146 &
            + fourth*3.036*t913 + fiveth*0.572*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.16e-01 * t9i32 * exp(-11.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 11.733*t9i2)

      ee   = 3.39e+01 * t9i32 * exp(-22.609*t9i)
      dee  = ee*(-1.5d0*t9i + 22.609*t9i2)

      ff   = 9.10e-03 * t95 * exp(-12.159*t9i)
      dff  = ff*(5.0d0*t9i + 12.159*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.42e+10 * t932 * exp(-59.328*t9i)
      drevdt = rev*(1.5d0*t9i + 59.328*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x,zz


! o17(a,g)ne21
       aa    = 1.0d0 + 0.1646*t9
       zz    = 1.0d0/aa
       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.1646)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.786/t9a
       daa    = -aa*zz
       bb     = aa**3.51
       dbb    = 3.51*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/1.084
       bb     = aa**1.69
       dbb    = 1.69*bb/aa * dt9a/1.084
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = oneth*exp(-10.106*t9i)
       daa    = aa*10.106*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 1.73e17 * fpt9a*zz
       daa    = (1.73e17*dfpt9a - aa*dgt9x)*zz

       bb     = 3.50e15 * ft9a*zz
       dbb    = (3.50e15*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 &
                 - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.63e+10 * t932 * exp(-85.305*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.305*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9x,dgt9x,zz


! o17(a,n)ne20
       aa     = 1.0d0 + 0.0268*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0268

       zz     = 1.0d0/bb
       cc     = aa + 0.0232*t953*zz
       dcc    = 0.0268 + (fiveth*0.0232*t923 - 0.0232*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd     = oneth*exp(-10.106*t9i)
       gt9x   = 1.0d0 + dd
       dgt9x  = dd*10.106*t9i2

       term      = 1.03e+18/gt9x * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt   = term*(-dgt9x/gt9x + dt9a56/t9a56 &
                   - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.86e+01 * exp(-6.852*t9i)
      drevdt   = rev*6.852*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_o18ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.117649d0)


! o18(a,g)ne22
! giessen et al 1994 nuc phys a 567, 146 for t9 less than 0.3
! cf88 otherwise

      if (t9.lt.0.3) then
       aa   = 1.066d-41 * t9i32 * exp(-5.507d-01*t9i)
       daa  = aa*(-1.5d0*t9i + 5.507d-1*t9i2)

       bb   = 1.852d-13 * t9i32 * exp(-2.070*t9i)
       dbb  = bb*(-1.5d0*t9i + 2.070*t9i2)

       cc   = 1.431d-02 * t9i32 * exp(-4.462*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.462*t9i2)

       dd   = 2.055d-04 * t9i32 * exp(-5.374*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.374*t9i2)

       ee   = 5.332d+00 * t9i32 * exp(-6.285*t9i)
       dee  = ee*(-1.5d0*t9i + 6.285*t9i2)

       ff   = 1.457d+00 * t9i32 * exp(-7.121*t9i)
       dff  = ff*(-1.5d0*t9i + 7.121*t9i2)

       gg   = 3.121d-02 * t9i32 * exp(-7.292*t9i)
       dgg  = gg*(-1.5d0*t9i + 7.292*t9i2)

       hh   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dhh  = hh*(t9i + 16.987*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg + hh
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


      else
       aa  = 1.82d+12 * t9i23 * exp(-40.057*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*40.057*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.01*t913 + 0.988*t923 + 0.072*t9 &
             + 3.17*t943 + 0.586*t953
       dbb = oneth*0.01*t9i23 + twoth*0.988*t9i13 + 0.072 &
            + fourth*3.17*t913 + fiveth*0.586*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 7.54 * t9i32 * exp(-6.228*t9i)
       ddd  = dd*(-1.5d0*t9i + 6.228*t9i2)

       ee   = 34.8 * t9i32 * exp(-7.301*t9i)
       dee  = ee*(-1.5d0*t9i + 7.301*t9i2)

       ff   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dff  = ff*(t9i + 16.987*t9i2)

       gg   = theta * 1.0d-11 * t9i32 * exp(-1.994*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.994*t9i2)

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg
      end if

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.85d+10 * t932 * exp(-112.208*t9i)
      drevdt   = rev*(1.5d0*t9i + 112.208*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_o18an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,bb,dbb,cc,dcc,dd, &
                       ee,dee,ff,dff,gg,dgg,hh,dhh,ft9a,dft9a,gt9,dgt9, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9i,zz


! o18(a,n)ne21
      aa     = 1.0d0 + 0.0483*t9
      bb     = aa**twoth
      dbb    = twoth*bb/aa*0.0483

      zz     = 1.0d0/bb
      cc     = aa + 0.00569*t953*zz
      dcc    = 0.0483 + (fiveth*0.00569*t923 - 0.00569*t953*zz*dbb)*zz

      zz     = 1.0d0/cc
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*dcc)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix * t9a56*zz

      dd     = 5.0d0 * exp(-23.002*t9i)
      gt9    = 1.0d0 + dd
      gt9i   = 1.0d0/gt9
      dgt9   = dd*23.002*t9i2

      ee     = 0.431/t9a
      dee    = -ee*zz
      ff     = ee**3.89
      dff    = 3.89*ff/ee*dee
      ft9a   = exp(-ff)
      dft9a  = -ft9a*dff

      gg     = 7.22e+17 * ft9a*gt9i * t9a56 * t9i32 * exp(-40.056/t9a13)
      dgg    = gg*(-dff - gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i &
               + 40.056/t9a13**2 *dt9a13)

      hh     = 150.31 / gt9 * exp(-8.045*t9i)
      dhh    = hh*(-gt9i*dgt9 + 8.045*t9i2)

      term     = gg + hh
      dtermdt  = dgg + dhh

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

! must protect the 8.045*t9i from overflow, so write it this way

      gg     = 7.22e+17*gt9i * t9a56 * t9i32 &
               * exp(-ff - 40.056/t9a13 + 8.045*t9i)
      dgg    = gg*(-gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i &
                - dff + 40.056/t9a13**2*dt9a13 - 8.045*t9i2)

      hh     = 150.31 * gt9i
      dhh    = -hh*gt9i*dgt9

      term    = 0.784 * (gg + hh)
      dtermdt = 0.784 * (dgg + dhh)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_ne20pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       ee,dee,ff,dff,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       zz,t9b


! ne20(p,a)f17
       aa     = 1.0d0 + 0.0612*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0612

       zz     = 1.0d0/bb
       cc     = aa + 0.013*t953*zz
       dcc    = 0.0612 + (fiveth*0.013*t923 - 0.013*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       t9b     = min(t9,10.0d0)
       dd      = 5.31 + 0.544*t9b - 0.0523*t9b*t9b
       ddd     = 0.544 - 2.0d0*0.0523*t9b
       if (t9b .eq. 10.0) ddd     = 0.0d0

       ee      = 3.25e19 * dd * t9a56 * t9i32 * exp(-43.176/t9a13)
       dee     = ee*(ddd/dd + dt9a56/t9a56 - 1.5d0*t9i &
                 +   43.176/t9a13**2 * dt9a13)

       ff      = exp(-47.969*t9i)
       dff     = ff*47.969*t9i2

       term    = ee * ff
       dtermdt = dee*ff + ee*dff


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0537 * ee
      drevdt   = 0.0537 * dee

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_f18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! f18(p,g)ne19
! wiescher and kettner, apj 263, 891 1982
      aa  = 1.658e+7 * t9i23 * exp(-18.06*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.06*t9i43)

      bb  = 4.604 + 0.106*t913 + 0.053*t923 + 0.009*t9 &
            - 0.036*t943 - 0.015*t953
      dbb = oneth*0.106*t9i23 + twoth*0.053*t9i13 + 0.009 &
            - fourth*0.036*t913 - fiveth*0.015*t923

! for temps greater than about t9 = 20, bb goes negative
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.55e-14 * t9i32* exp(-0.302*t9i)
      ddd  = dd*(-1.5d0*t9i + 0.302*t9i2)

      ee   = 327.0 * t9i32 * exp(-3.84*t9i)
      dee  = ee*(-1.5d0*t9i + 3.84*t9i2)

      ff   = 1.32e+04 * t9i32 * exp(-5.22*t9i)
      dff  = ff*(-1.5d0*t9i + 5.22*t9i2)

      gg   = 93.0 * t9i32 * exp(-4.29*t9i)
      dgg  = gg*(-1.5d0*t9i + 4.29*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.73e+10 * t932 * exp(-74.396*t9i)
      drevdt   = rev*(1.5d0*t9i + 74.396*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_f19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,zz,q1
      parameter        (q1 = 1.0d0/0.173056d0)


! f19(p,g)ne20
      aa  = 6.04e+07 * t9i23 * exp(-18.113*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*18.113*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.023*t913 + 2.06*t923 + 0.332*t9 &
            + 3.16*t943 + 1.30*t953
      dbb = oneth*0.023*t9i23 + twoth*2.06*t9i13 + 0.332 &
            + fourth*3.16*t913 + fiveth*1.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.32e+02 * t9i32 * exp(-3.752*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.752*t9i2)

      ee   = 7.56e+04 * t9i27 * exp(-5.722*t9i)
      dee  = ee*(-twosev*t9i + 5.722*t9i2)

      ff   = 7.0*exp(-16.44*t9i)
      dff  = ff*16.44*t9i2

      gg   = 4.0 * exp(-2.09*t9i)
      dgg  = gg*2.09*t9i2

      hh   = 1.0d0 + ff + gg
      dhh  = dff + dgg

      zz      = 1.0d0/hh
      term    = (cc + dd + ee)*zz
      dtermdt = (dcc + ddd + dee - term*dhh)*zz

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.7e+10 * t932 * exp(-149.093*t9i)
      drevdt   = rev*(1.5d0*t9i + 149.093*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! f19(p,n)ne19
      aa  = 1.27e+08 * (1.0d0 - 0.147*t912 + 0.069*t9)
      daa = 1.27e+08 * (0.069 - 0.5d0*0.147*t9i12)

      bb  = exp(-46.659*t9i)
      dbb = bb*46.659*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_f19ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,q1
      parameter        (q1 = 1.0d0/0.405769)


! f19(a,p)ne22
      aa  = 4.50e+18 * t9i23 * exp(-43.467*t9i13  - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*43.467*t9i43 - 2.0d0*t9*q1)

      bb   = 7.98e+04 * t932 * exp(-12.760*t9i)
      dbb  = 1.5d0*bb*t9i + bb*12.760*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.36 * exp(-19.439*t9i)
      drevdt   = rev*19.439*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_na22na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa, &
                       t9b,t9b2,t9b3


! na22(n,a)f19
      t9b  = min(t9,10.0d0)
      t9b2 = t9b*t9b
      t9b3 = t9b2*t9b
      aa  = 1.0d0 + 0.8955*t9b - 0.05645*t9b2 + 7.302e-04*t9b3
      daa = 0.8955 - 2.0d0*0.05645*t9b + 3.0d0*7.302e-4*t9b2
      if (t9b .eq. 10.0) daa = 0.0d0

      term     = 1.21e6 * exp(aa)
      dtermdt  = term*daa

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10 * exp(-22.620*t9i)
      drevdt   = rev*22.620*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne20pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ff,gg,dgg,zz


! ne20(p,g)na21
      aa  = 9.55e+06 * exp(-19.447*t9i13)
      daa = aa*oneth*19.447*t9i43

      bb  = 1.0d0 + 0.0127*t9i23
      dbb = -twoth*0.0127*t9i53

      cc  = t92 * bb * bb
      dcc = 2.0d0*cc*t9i + 2.0d0*t92*bb*dbb

      zz  = 1.0d0/cc
      dd  = aa*zz
      ddd = (daa - dd*dcc)*zz

      aa  = 2.05e+08 * t9i23 * exp(-19.447*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.447*t9i43)

      bb  = sqrt (t9/0.21)
      dbb = 0.5d0/(bb * 0.21)

      cc  = 2.67 * exp(-bb)
      dcc = -cc*dbb

      ff  = 1.0d0 + cc

      gg  = aa*ff
      dgg = daa*ff + aa*dcc


      aa  = 18.0 * t9i32 * exp(-4.242*t9i)
      daa = aa*(-1.5d0*t9i + 4.242*t9i2)

      bb  = 10.2 * t9i32 * exp(-4.607*t9i)
      dbb = bb*(-1.5d0*t9i + 4.607*t9i2)

      cc  = 3.6e+04 * t9i14 * exp(-11.249*t9i)
      dcc = cc*(-0.25d0*t9i + 11.249*t9i2)

      term    = dd + gg + aa + bb + cc
      dtermdt = ddd + dgg + daa + dbb + dcc


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.63e+09 * t932 * exp(-28.216*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.216*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta,q1,q2
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.0169d0, &
                        q2    = 1.0d0/0.017161d0)


! na23(p,a)ne20
! el eid & champagne 1995
      if (t9 .le. 2.0) then
       aa  = 1.26d+10 * t9i23 * exp(-20.758*t9i13 - t92*q1)
       daa = -twoth*aa*t9i + aa*(oneth*20.758*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0  + 0.02*t913 - 13.8*t923 - 1.93*t9 &
             + 234.0*t943 + 83.6*t953
       dbb = oneth*0.02*t9i23 - twoth*13.8*t9i13 - 1.93 &
             + fourth*234.0*t913 + fiveth*83.6*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.38 * t9i32 * exp(-1.979*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.979*t9i2

       ee   = 6.50d+06 * (t9**(-1.366d0)) * exp(-6.490*t9i)
       dee  = -1.366d0*ee*t9i + ee*6.490*t9i2

       ff   = 1.19d+08 * (t9**(-1.055d0)) * exp(-11.411*t9i)
       dff  = -1.055d0*ff*t9i + ff*11.411*t9i2

       gg   = theta * 9.91d-14 * t9i32 * exp(-0.418*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.418*t9i2

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg



! cf88 + one term from gorres, wiesher & rolfs 1989, apj 343, 365
      else
       aa  = 8.56d+09 * t9i23 * exp(-20.766*t9i13 - t92*q2)
       daa = -twoth*aa*t9i + aa*(oneth*20.766*t9i43 - 2.0d0*t9*q2)

       bb  = 1.0d0  + 0.02*t913 + 8.21*t923 + 1.15*t9 &
             + 44.36*t943 + 15.84*t953
       dbb = oneth*0.02*t9i23 + twoth*8.21*t9i13 + 1.15 &
             + fourth*44.36*t913 + fiveth*15.84*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.02 * t9i32 * exp(-1.99*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.99*t9i2

       ee   = 1.18d+04 * t9i54 * exp(-3.148*t9i)
       dee  = -1.25d0*ee*t9i + ee*3.148*t9i2

       ff   = 8.59d+05 * t943 * exp(-4.375*t9i)
       dff  = fourth*ff*t9i + ff*4.375*t9i2

       gg   = theta * 3.06d-12 * t9i32 * exp(-0.447*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.447*t9i2

       hh   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dhh  = -1.5d0*hh*t9i + hh*1.601*t9i2

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.25 * exp(-27.606*t9i)
      drevdt   = rev*27.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne20ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! ne20(n,g)ne21
! wm88 Apj 239, 943; fit over range of experimental data, constant otherwise

      if (t9 .lt. 5.8025d-2) then
       term    = 5.449d+03
       dtermdt = 0.0d0

      else if (t9 .gt. 1.1605) then
       term    = 6.977d+04
       dtermdt = 0.0d0

      else if (t9 .ge. 5.8025d-2 .and. t9 .le. 2.9012d-1) then
       term    = 4.7219d+3 + 2.5248d+4*t9 - 2.7448d+5*t92 &
                 + 9.2848d+5*t93
       dtermdt = 2.5248d+4 - 2.0d0*2.7448d+5*t9 &
                 + 3.0d0*9.2848d+5*t92

      else

       aa  = 1.802d+04 * (t9/0.348)**4.43
       daa = 4.43 * aa * t9i

       bb  = -5.931 * (t9-0.348) + 1.268 * (t9-0.348)**2
       dbb = -5.931 + 2.0d0*1.268*(t9 - 0.348)

       cc  = exp(bb)
       dcc = cc*dbb

       term    = aa * cc
       dtermdt = daa*cc + aa*dcc

      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      =  4.650d+09 * t932 * exp(-78.46*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.46*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,xx,dxx,zz, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.003364d0)


! ne21(p,g)na22

! el eid & champagne 1995

      if (t9.le.2.0) then
       aa  = 3.4d+08 * t9i23 * exp(-19.41*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.41*t9i43)

       bb  = (16.7*t9 - 1.0)**2
       dbb = 2.0d0*(16.7*t9 - 1.0)*16.7

       cc  = 0.56 * exp(-bb)
       dcc = -cc*dbb

       dd  = 1.0d0 + cc
       ddd = dcc

       ee   = aa * dd
       dee  = daa*dd + aa*ddd

       ff   = 6.12 * t9i32 * exp(-1.403*t9i)
       dff  = ff*(-1.5d0*t9i + 1.403*t9i2)

       gg   = 1.35d+04 * t9i32 * exp(-3.008*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.008*t9i2)

       aa   = t9**0.67
       daa  = 0.67*aa*t9i
       zz   = 1.0d0/aa

       hh   = 3.12d+06 * t9**(-0.72) * exp(-8.268*zz)
       dhh  = hh*(-0.72d0*t9i + 8.268*zz*zz*daa)

       xx   = theta * 1.1d-03 * t9i32 * exp(-1.114*t9i)
       dxx  = xx*(-1.5d0*t9i + 1.114*t9i2)

       term    = ee + ff + gg + hh + xx
       dtermdt = dee + dff + dgg + dhh + dxx


! cf88
      else

       aa  = theta * 2.95d+08 * t9i23 * exp(-19.462*t9i13 -t92*q1)
       daa = aa*(-twoth*t9i + oneth*19.462*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.021*t913 + 13.29*t923 + 1.99*t9 &
             + 124.1*t943 + 47.29*t953
       dbb = oneth*0.021*t9i23 + twoth*13.29*t9i13 + 1.99 &
             + fourth*124.1*t913 + fiveth*47.29*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = theta * 7.80d-01 * t9i32 * exp(-1.085*t9i)
       ddd  = dd*(-1.5d0*t9i + 1.085*t9i2)

       ee   = 4.37d+08 * t9i23 * exp(-19.462*t9i13)
       dee  = ee*(-twoth*t9i + oneth*19.462*t9i43)

       ff   = 5.85 * t9i32 * exp(-1.399*t9i)
       dff  = ff*(-1.5d0*t9i + 1.399*t9i2)

       gg   = 1.29d+04 * t9i32 * exp(-3.009*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.009*t9i2)

       hh   = 3.15d+05 * t9i35 * exp(-5.763*t9i)
       dhh  = hh*(-0.6d0*t9i + 5.763*t9i2)

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.06d+10 * t932 * exp(-78.194*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.194*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ne21ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,dgg,hh,dhh,zz, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! ne21(a,g)mg25
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 8.72e-03*t9 - 6.87e-04*t92 + 2.15e-05*t93
       daa    = 8.72e-3 - 2.0d0*6.87e-4*t9 + 3.0d0*2.15e-5*t92

       bb     = 1.52e-04 * exp(-46.90*t9i13*aa)
       dbb    = bb*46.90*(oneth*t9i43*aa - t9i13*daa)

       cc     = 1.5*exp(-4.068*t9i)
       dcc    =  cc*4.068*t9i2

       gg     = 2.0 * exp(-20.258*t9i)
       dgg    = gg*20.258*t9i2

       hh     = 1.0d0 + cc + gg
       dhh    = dcc + dgg

       zz     = 1.0d0/hh
       dd     = bb*zz
       ddd    = (dbb - dd*dhh)*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = aa + bb
       dcc    = daa + dbb

       term    = dd * cc
       dtermdt = ddd*cc + dd*dcc


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 4.06e+10 * t932 * exp(-114.676*t9i)
      drevdt = rev*(1.5d0*t9i + 114.676*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_ne21an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,zz, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! ne21(a,n)mg24
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = 2.0d0*exp(-20.258*t9i)
       dcc    = cc*20.258*t9i2

       dd     = 1.5*exp(-4.068*t9i)
       ddd    = dd*4.068*t9i2

       ee     = 1.0d0 + cc + dd
       dee    = dcc + ddd

       zz      = 1.0d0/ee
       term    = (aa + bb)*zz
       dtermdt = (daa + dbb - term*dee)*zz


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 12.9 * exp(-29.606*t9i)
      drevdt   = rev*29.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


! ne22(p,g)na23

! el eid & champagne 1995

      if (t9.le.2.0) then
       aa  = 1.05d+09 * t9i23 * exp(-19.431*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

       bb  = 1.24d-09 * t9i32 * exp(-0.414*t9i)
       dbb = bb*(-1.5d0*t9i + 0.414*t9i2)

       cc  = 2.90d-02 * t9i32 * exp(-1.752*t9i)
       dcc = cc*(-1.5d0*t9i + 1.752*t9i2)

       dd  = 9.30d+04 * t9**(-1.174) * exp(-5.100*t9i)
       ddd = dd*(-1.174*t9i + 5.100*t9i2)

       ee   = 5.71d+05 * t9**(0.249) * exp(-7.117*t9i)
       dee  = ee*(0.249*t9i + 7.117*t9i2)

       ff   = theta * 3.25d-04 * t9i32 * exp(-0.789*t9i)
       dff  = ff*(-1.5d0*t9i + 0.789*t9i2)

       gg   = theta * 0.10 * t9i32 * exp(-1.161*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.161*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! cf88
      else

       aa  = 1.15d+09 * t9i23 * exp(-19.475*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.475*t9i43)

       bb  = 9.77d-12 * t9i32 * exp(-0.348*t9i)
       dbb = bb*(-1.5d0*t9i + 0.348*t9i2)

       cc   = 8.96d+03 * t9i32 * exp(-4.84*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.84*t9i2)

       dd   = 6.52d+04 * t9i32 * exp(-5.319*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.319*t9i2)

       ee   = 7.97d+05 * t9i12 * exp(-7.418*t9i)
       dee  = ee*(-0.5d0*t9i + 7.418*t9i2)

       ff   = theta * 1.63d-01 * t9i32 * exp(-1.775*t9i)
       dff  = ff*(-1.5d0*t9i + 1.775*t9i2)

       term    = aa + bb + cc + dd + ee + ff
       dtermdt = daa + dbb + dcc + ddd + dee + dff

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.67d+09 * t932 * exp(-102.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne22ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,res1,dres1, &
                       ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       rdmass,res2,zz
      parameter        (rdmass = 22.0d0*4.0d0/26.0d0, &
                        res2   = -11.604d0 * 22.0d0/26.0d0)


! ne22(a,g)mg26
! kappeler 1994 apj 437, 396

      if (t9 .lt. 1.25) then

       res1 = 1.54d-01*(t9*rdmass)**(-1.5)
       dres1 = -1.5d0 * res1 * t9i

       aa    = 1.7d-36 * res1 * exp(res2*t9i*0.097)
       daa   = aa/res1*dres1 - aa*res2*0.097*t9i2

       bb    = 1.5d-7 * res1 * exp(res2*t9i*0.400)
       dbb   = bb/res1*dres1 - bb * res2 * 0.400 * t9i2

       cc    = 0.5 * res1 * 3.7d-2 * exp(res2*t9i*0.633)
       dcc   = cc/res1*dres1 - cc*res2*0.633*t9i2

       dd    = res1 * 3.6d+1 * exp(res2*t9i*0.828)
       ddd   = dd/res1*dres1 - dd*res2*0.828*t9i2

       term    = aa + bb + cc + dd
       dtermdt = daa + dbb + dcc + ddd


! cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/0.249
       bb     = aa**2.31
       dbb    = 2.31*bb/aa * dt9a/0.249
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 4.16e19 * fpt9a*zz
       daa    = (4.16e19*dfpt9a - aa*dgt9x)*zz

       bb     = 2.08e16 * ft9a*zz
       dbb    = (2.08e16*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-47.004/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 &
                       - 1.5d0*t9i + 47.004/t9a13**2 * dt9a13)
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.15d+10 * t932 * exp(-123.151*t9i)
      drevdt = rev*(1.5d0*t9i + 123.151*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na22np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! na22(n,p)ne22
      aa  = 1.0d0 - 3.037e-02*t9 + 8.380e-03*t92 - 7.101e-04*t93
      daa =  -3.037e-02 + 2.0d0*8.380e-03*t9 - 3.0d0*7.101e-04*t92

      term    = 1.24e+08 * exp(aa)
      dtermdt = term*daa

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev     = 7.01*exp(-42.059*t9i)
      drevdt  = rev*42.059*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne22an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,dgg,ft9a,dft9a,gt9x,dgt9x, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,res1,res2, &
                       zz
      parameter        (res1 = 2.4731857150793075d-2, &
                        res2 = -9.8187694549560547d0)

! note: res1=1.54d-1*(88./26.)**(-1.5)   res2=-11.604*(22./26.)

! ne22(a,n)mg25
! kappeler 1994 apj 437, 396 ; wiescher suggest only 828 kev, ignore 633 kev

      if (t9 .lt. 0.6) then
       term = res1*1.64d+02 * t9i32 * exp(t9i*0.828*res2)
       dtermdt = -1.5d0*term*t9i - term*res2*0.828*t9i2

! cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       gg     = bb
       dgg    = dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = ft9a*zz
       daa    = (dft9a - aa*dgt9x)*zz

       bb     = 4.16e+19 * t9a56 * t9i32 * exp(-47.004/t9a13)
       dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 47.004/t9a13**2 * dt9a13)

       cc     = aa*bb
       dcc    = daa*bb + aa*dbb

       dd      = 1.44e-04*zz * exp(-5.577*t9i)
       ddd     = -dd*zz*dgt9x + dd*5.577*t9i2

       term    = cc + dd
       dtermdt = dcc + ddd
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 7.833d-5
      drevdt = 0.0d0
      if (t9 .gt. 0.008) then
       rev    = 0.544 * exp(5.577*t9i)
       drevdt = -rev*5.577*t9i2
      end if

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_na21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/0.133956d0)


! na21(p,g)mg22
      aa  = 1.41e+05 * t9i23 * exp(-20.739*t9i13 -  t92*q1)
      daa = aa*(-twoth*t9i + oneth*20.739*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.020*t913 + 4.741*t923 + 0.667*t9 &
            + 16.380*t943 + 5.858*t953
      dbb = oneth*0.020*t9i23 + twoth*4.741*t9i13 + 0.667 &
            + fourth*16.380*t913 + fiveth*5.858*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.72e+02 * t9i34 * exp(-2.436*t9i)
      ddd  = dd*(-0.75d0*t9i + 2.436*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.44e+10 * t932 * exp(-63.790*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.790*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! mg24(p,a)na21
       aa     = 1.0d0 + 0.127*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.127)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       gg      = min(t9,12.0d0)
       aa      = 4.43 + 3.31*gg - 0.229*gg*gg
       daa     = 3.31 - 2.0d0*0.229*gg
       if (gg .eq. 12.0) daa = 0.0d0

       bb      = 1.81e21 * t9a56 * t9i32 * exp(-49.967/t9a13)
       dbb     = bb*(dt9a56/t9a56 - 1.5d0*t9i &
                     + 49.967/t9a13**2 * dt9a13)

       cc      = aa*bb
       dcc     = daa*bb + aa*dbb

       dd      = exp(-79.843*t9i)
       ddd     = dd*79.843*t9i2

       term    = cc * dd
       dtermdt = dcc*dd + cc*ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0771 * cc
      drevdt   = 0.0771 * dcc

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end






      subroutine rate_na22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! na22(p,g)mg23
      aa  = 9.63e-05 * t932 * exp(-0.517*t9i)
      daa = aa*(1.5d0*t9i + 0.517*t9i2)

      bb  = 2.51e+04 * t9 * exp(-2.013*t9i)
      dbb = bb*(t9i + 2.013*t9i2)

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.27e+10 * t932 * exp(-87.933*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_na23pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,hhi,xx,dxx, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.088209d0)


! na23(p,g)mg24

! el eid & champagne 1995
      if (t9 .le. 2.0) then
       aa  = 2.47d+09 * t9i23 * exp(-20.758*t9i13)
       daa = aa*(-twoth*t9i + oneth*20.758*t9i43)

       bb  = 9.19d+01 * t9i32 * exp(-2.789*t9i)
       dbb = bb*(-1.5d0*t9i + 2.789*t9i2)

       cc  = 1.72d+04 * t9i32 * exp(-3.433*t9i)
       dcc = cc*(-1.5d0*t9i + 3.433*t9i2)

       dd  = 3.44d+04 * t9**0.323 * exp(-5.219*t9i)
       ddd = dd*(0.323*t9i + 5.219*t9i2)

       ee   = theta * 2.34d-04 * t9i32 * exp(-1.590*t9i)
       dee  = ee*(-1.5d0*t9i + 1.590*t9i2)

       term    = aa + bb + cc + dd + ee
       dtermdt = daa + dbb + dcc + ddd + dee


! cf88 + gorres, wiesher & rolfs 1989, apj 343, 365
      else

       aa  = 2.93d+08 * t9i23 * exp(-20.766*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*20.766*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.02*t913 + 1.61*t923 + 0.226*t9 &
             + 4.94*t943 + 1.76*t953
       dbb = oneth*0.02*t9i23 + twoth*1.61*t9i13 + 0.226 &
            + fourth*4.94*t913 + fiveth*1.76*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 9.34d+01 * t9i32 * exp(-2.789*t9i)
       dcc  = cc*(-1.5d0*t9i + 2.789*t9i2)

       dd   = 1.89d+04 * t9i32 * exp(-3.434*t9i)
       ddd  = dd*(-1.5d0*t9i + 3.434*t9i2)

       ee   = 5.1d+04 * t915 * exp(-5.51*t9i)
       dee  = ee*(0.2d0*t9i + 5.51*t9i2)

       ff   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dff  = ff*(-1.5d0*t9i + 1.601*t9i2)

       gg   = 1.5 * exp(-5.105*t9i)
       dgg  = gg*5.105*t9i2

       hh   = 1.0d0 + gg
       hhi  = 1.0d0/hh

       term    = (xx + cc + dd + ee + ff) * hhi
       dtermdt = (dxx + dcc + ddd + dee + dff - term*dgg)*hhi

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.49d+10 * t932 * exp(-135.665*t9i)
      drevdt   = rev*(1.5d0*t9i + 135.665*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,bb,dbb,cc,dcc, &
                       t9a,dt9a,t9a32,dt9a32,zz


! na23(p,n)mg24
      aa  = 1.0d0 + 0.141*t9
      zz  = 1.0d0/aa

      t9a = t9*zz
      dt9a = (1.0d0 - t9a*0.141)*zz

      aa    = sqrt(t9a)
      t9a32 = t9a * aa
      dt9a32 = 1.5d0 * aa * dt9a

      bb   = 9.29d8 * (1.0d0 - 0.881d0 * t9a32 * t9i32)
      dbb  = -9.29d8 * 0.881d0 * t9i32*(dt9a32 - 1.5d0*t9a32*t9i)

      cc   = exp(-56.173*t9i)
      dcc  = cc*56.173*t9i2

      term    = bb * cc
      dtermdt = dbb*cc + bb*dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * bb
      dtermdt = 0.998 * dbb

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_mg24pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,ggi


! mg24(p,g)al25
      aa  = 5.60e+08 * t9i23 * exp(-22.019*t9i13)
      daa = aa*(-twoth*t9i + oneth*22.019*t9i43)

      bb  = 1.0d0 + 0.019*t913 - 0.173*t923 - 0.023*t9
      dbb = oneth*0.019*t9i23 - twoth*0.173*t9i13 - 0.023

! stop negative rates above t9 = 10
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd   = 1.48e+03 * t9i32 * exp(-2.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.484*t9i2)

      ee   = 4.00e+03 * exp(-4.180*t9i)
      dee  = ee*4.180*t9i2

      ff   = 5.0 * exp(-15.882*t9i)
      dff  = ff*15.882*t9i2

      gg   = 1.0d0 + ff
      ggi  = 1.0d0/gg

      term    = (cc + dd + ee) * ggi
      dtermdt = (dcc + ddd + dee - term*dff)*ggi


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.13e+09 * t932 * exp(-26.358*t9i)
      drevdt   = rev*(1.5d0*t9i + 26.358*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


! al27(p,a)mg24
! champagne 1996
      aa  = 4.71d+05 * t9i23 * exp(-23.25*t9i13 - 3.57*t92)
      daa = -twoth*aa*t9i + aa*(oneth*23.25*t9i43 - 2.0d0*3.57*t9)

      bb  = 1.0d0 + 0.018*t913 - 7.29*t923 - 0.914*t9 &
            + 77.2*t943 + 24.6*t953
      dbb = oneth*0.018*t9i23 - twoth*7.29*t9i13 - 0.914 &
            + fourth*77.2*t913 + fiveth*24.6*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 2.23d+04 * (t9**3.989) * exp(-2.148 * t9**(-1.293))
      ddd = 3.989*dd*t9i + 1.293*dd*2.148*t9**(-2.293)

      ee   = 0.17 * 1.29d-09 * t9i32 * exp(-0.836*t9i)
      dee  = -1.5d0*ee*t9i + ee*0.836*t9i2

      ff   = theta * 2.73d-03 * t9i32 * exp(-2.269*t9i)
      dff  = -1.5d0*ff*t9i + ff*2.269*t9i2

      gg   = theta * 2.60d-02 * t9i32 * exp(-2.492*t9i)
      dgg  = -1.5d0*gg*t9i + gg*2.492*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.81*exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mg25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,q1,q2
      parameter        (q1 = 1.0d0/0.0036d0, &
                        q2 = 1.0d0/169.0d0)


! mg25(p,g)al26
      aa  = 3.57e+09 * t9i23 * exp(-22.031*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*22.031*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.019*t913 + 7.669*t923 + 1.015*t9 &
            + 167.4*t943 + 56.35*t953
      dbb = oneth*0.019*t9i23 + twoth*7.669*t9i13 + 1.015 &
            + fourth*167.4*t913 + fiveth*56.35*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.07e-13 * t9i32 * exp(-0.435*t9i)
      ddd = dd*(-1.5d0*t9i + 0.435*t9i2)

      ee   = 1.94e-07 * t9i32 * exp(-0.673*t9i)
      dee  = ee*(-1.5d0*t9i + 0.673*t9i2)

      ff   = 3.15e-05 * t9**(-3.40)* exp(-1.342*t9i - t92*q2)
      dff  = ff*(-3.40d0*t9i + 1.342*t9i2 - 2.0d0*t9*q2)

      gg   = 1.77e+04 * t958 * exp(-3.049*t9i - t92*q2)
      dgg  = gg*(0.625*t9i + 3.049*t9i2 - 2.0d0*t9*q2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.03e+10 * t932 * exp(-73.183*t9i)
      drevdt   = rev*(1.5d0*t9i + 73.183*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_mg25ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc


! mg25(a,p)si28
      aa  = -23.271*t9i13 + 6.46*t9 - 2.39*t92 + 0.506*t93 &
            - 6.04e-2*t94 + 3.75e-3*t95 - 9.38e-5*t96

      daa = oneth*23.271*t9i43 + 6.46 - 2.0d0*2.39*t9 + 3.0d0*0.506*t92 &
            - 4.0d0*6.04e-2*t93 + 5.0d0*3.75e-3*t94 - 6.0d0*9.38e-5*t95

      bb  = 3.23e8 * t9i23 * exp(aa)
      dbb  = -twoth*bb*t9i + bb*daa

! dbb/bb
      cc   = -twoth*t9i + daa

      term    = bb * exp(-13.995*t9i)
      dtermdt = term*cc + term*13.995*t9i2

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.86 * bb
      drevdt   = 2.86 * dbb

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_mg25ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,zz


! mg25(a,g)si29
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                   + 53.41/t9a13**2 * dt9a13)

      dd     = 0.0156*t9 - 1.79e-03*t92 + 9.08e-05*t93
      ddd    = 0.0156 - 2.0d0*1.79e-03*t9 + 3.0d0*9.08e-05*t92

      ee     = 5.87e-04*exp(-53.42*t9i13*dd)
      dee    = ee*53.42*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.90e+11 * t932 * exp(-129.128*t9i)
      drevdt = rev*(1.5d0*t9i + 129.128*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg25an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz, &
                       gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! mg25(a,n)si28
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term    = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dtermdt = term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                + 53.41/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 20.0*exp(-30.792*t9i)
      drevdt   = rev*30.792*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_mg26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta
      parameter        (theta = 0.1d0)


! mg26(p,g)al27
! champagne 1996
      aa  = 8.54d-12 * t9i32 * exp(-0.605*t9i)
      daa = aa*(-1.5d0*t9i + 0.605*t9i2)

      bb  = 2.75d-06 * t9i32 * exp(-1.219*t9i)
      dbb = bb*(-1.5d0*t9i + 1.219*t9i2)

      cc  = 1.30d-02 * t9i32 * exp(-1.728*t9i)
      dcc = cc*(-1.5d0*t9i + 1.728*t9i2)

      dd  = 8.06d+00 * t9i32 * exp(-2.537*t9i)
      ddd = dd*(-1.5d0*t9i + 2.537*t9i2)

      ee   = 1.45d+03 * t9i32 * exp(-3.266*t9i)
      dee  = ee*(-1.5d0*t9i + 3.266*t9i2)

      ff   = 4.03d+04 * t9i32 * exp(-3.784*t9i)
      dff  = ff*(-1.5d0*t9i + 3.784*t9i2)

      gg   = 8.82d+04 * t9**(-0.21) * exp(-4.194*t9i)
      dgg  = gg*(-0.21*t9i + 4.194*t9i2)

      hh   = theta * 1.93d-05 * t9i32 * exp(-1.044*t9i)
      dhh  = hh*(-1.5d0*t9i + 1.044*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg + hh
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.14d+09 * t932 * exp(-95.99*t9i)
      drevdt   = rev*(1.5d0*t9i + 95.99*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,zz


! mg26(a,g)si30
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                   + 53.505/t9a13**2 * dt9a13)

      dd     = 0.0751*t9 - 0.0105*t92 + 5.57e-04*t93
      ddd    = 0.0751 - 2.0d0*0.0105*t9 + 3.0d0*5.57e-04*t92

      ee     = 4.55e-2 * exp(-53.51*t9i13*dd)
      dee    = ee*53.51*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.38e+10 * t932 * exp(-123.52*t9i)
      drevdt = rev*(1.5d0*t9i + 123.52*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz, &
                       gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! mg26(a,n)si29
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term   = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dtermdt= term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
               + 53.505/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.68*exp(-0.401*t9i)
      drevdt   = rev*0.401*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_al25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee

! al25(p,g)si26
! coc et al 1995 a&a 299, 479 , case b

      aa  = 8.98d+1 * t9i32 * exp(-4.874*t9i)
      daa = aa*(-1.5d0*t9i + 4.874*t9i2)

      bb  = 1.568d+3 * t9i32 * exp(-9.632*t9i)
      dbb = bb*(-1.5d0*t9i + 9.632*t9i2)

      cc  = 2.42d+8 * t9i23 * exp(-23.18*t9i13)
      dcc = cc*(-twoth*t9i + oneth*23.18*t9i43)

      dd  = 4.10d-02 * t9i32 * exp(-1.741*t9i)
      ddd = dd*(-1.5d0*t9i + 1.741*t9i2)

      ee  = 2.193d+3 * t9i32 * exp(-4.642*t9i)
      dee = ee*(-1.5d0*t9i + 4.642*t9i2)

      term    = aa + bb + cc + dd + ee
      dtermdt = daa + dbb + dcc + ddd + dee


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.117d+11 * t932 * exp(-64.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 64.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,theta
      parameter        (theta = 0.1d0)


! al26(p,g)si27
! coc et al 1995 a&a 299, 479

      aa  = 1.53d+9 * t9**(-1.75) * exp(-23.19*t9i13)
      daa = aa*(-1.75*t9i + oneth*23.19*t9i43)

      bb  = theta*8.7d-7 * t9i32 * exp(-0.7845*t9i)
      dbb = bb*(-1.5d0*t9i + 0.7845*t9i2)

      cc  = theta*1.00d-3 * t9i32 * exp(-1.075*t9i)
      dcc = cc*(-1.5d0*t9i + 1.075*t9i2)

      dd  = 9.00d+00 * t9i32 * exp(-2.186*t9i)
      ddd = dd*(-1.5d0*t9i + 2.186*t9i2)

      ee  = 5.05d+02 * t9i32 * exp(-3.209*t9i)
      dee = ee*(-1.5d0*t9i + 3.209*t9i2)

      ff  = 9.45d+03 * t9i * exp(-4.008*t9i)
      dff = ff*(-t9i + 4.008*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.46d+10 * t932 * exp(-86.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 86.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! al27(a,n)p30

      aa    = 8.2e+04*exp(-30.588*t9i)
      daa   = aa*30.588*t9i2

      bb    = 5.21e+05 * t974 * exp(-33.554*t9i)
      dbb   = 1.75d0*bb*t9i + bb*33.554*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      aa  = 5.21e+05 * t974 * exp(-2.966*t9i)
      daa = aa*(1.75d0*t9i + 2.966*t9i2)

      rev      = 6.75d0 * (8.20e4 + aa)
      drevdt   = 6.75d0 * daa

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_si27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd


! si27(p,g)p28

      aa  = 1.64e+09 * t9i23 * exp(-24.439*t9i13)
      daa = aa*(-twoth*t9i + oneth*24.439*t9i43)

      bb  = 2.00e-08 * t9i32 * exp(-0.928*t9i)
      dbb = bb*(-1.5d0*t9i + 0.928*t9i2)

      cc  = 1.95e-02 * t9i32 * exp(-1.857*t9i)
      dcc = cc*(-1.5d0*t9i + 1.857*t9i2)

      dd  = 3.70e+02 * t9i47 * exp(-3.817*t9i)
      ddd = dd*(-foursev*t9i + 3.817*t9i2)

      term    = aa + bb + cc + dd
      dtermdt = daa + dbb + dcc + ddd


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.62e+10 * t932 * exp(-23.960*t9i)
      drevdt   = rev*(1.5d0*t9i + 23.960*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/8.4681d0)



! si28(p,g)p29

! champagne et al 96

      if (t9.le.5.0) then

       aa  = 8.44d+08 * t9i23 * exp(-24.389*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*24.389*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.17*t913 + 0.113*t923 + 0.0135*t9 &
             + 0.194*t943 + 0.0591*t953
       dbb = oneth*0.17*t9i23 + twoth*0.113*t9i13 + 0.0135 &
            + fourth*0.194*t913 + fiveth*0.0591*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 2.92d+02 * t9i32 * exp(-4.157*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.157*t9i2)

       dd   = 4.30d+05 * t9i32 * exp(-18.51*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.51*t9i2)

       ee   = 6.05d+03 * t9i32 * exp(-18.17*t9i)
       dee  = ee*(-1.5d0*t9i + 18.17*t9i2)

       term    = xx + cc + dd + ee
       dtermdt = dxx + dcc + ddd + dee


! cf88
      else

       aa  = 1.64d+08 * t9i23 * exp(-24.449*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*24.449*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.017*t913 - 4.11*t923 - 0.491*t9 &
             + 5.22*t943 + 1.58*t953
       dbb = oneth*0.017*t9i23 - twoth*4.11*t9i13 - 0.491 &
            + fourth*5.22*t913 + fiveth*1.58*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 3.52d+02 * t9i32 * exp(-4.152*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.152*t9i2)

       dd   = 6.3d+05 * t9i32 * exp(-18.505*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.505*t9i2)

       ee   = 1.69d+03 * exp(-14.518*t9i)
       dee  = ee*14.518*t9i2

       term    = xx + cc + dd + ee
       dtermdt = dxx + dcc + ddd + dee

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.46d+09 * t932 * exp(-31.879*t9i)
      drevdt   = rev*(1.5d0*t9i + 31.879*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si29pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/0.065536d0)



! si29(p,g)p30

      aa  = 3.26e+09 * t9i23 * exp(-24.459*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*24.459*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.017*t913 + 4.27*t923 + 0.509*t9 &
            + 15.40*t943 + 4.67*t953
      dbb = oneth*0.017*t9i23 + twoth*4.27*t9i13 + 0.509 &
           + fourth*15.40*t913 + fiveth*4.67*t923

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 2.98e+03 * t9i32 * exp(-3.667*t9i)
      dcc  = cc*(-1.5d0*t9i + 3.667*t9i2)

      dd   = 3.94e+04 * t9i32 * exp(-4.665*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.665*t9i2)

      ee   = 2.08e+04 * t912 * exp(-8.657*t9i)
      dee  = ee*(0.5d0*t9i + 8.657*t9i2)

      term    = xx + cc + dd + ee
      dtermdt = dxx + dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.26e+10 * t932 * exp(-65.002*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.002*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si30pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/0.4489d0)



! si30(p,g)p31

      aa  = 4.25e8 * t9i23 * exp(-24.468*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*24.468*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.017*t913 + 0.150*t923 + 0.018*t9 &
            + 5.53*t943 + 1.68*t953
      dbb = oneth*0.017*t9i23 + twoth*0.150*t9i13 + 0.018 &
           + fourth*5.53*t913 + fiveth*1.68*t923

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 1.86e4 * t9i32 * exp(-5.601*t9i)
      dcc  = cc*(-1.5d0*t9i + 5.601*t9i2)

      dd   = 3.15e5 * t9i32 * exp(-6.961*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.961*t9i2)

      ee   = 2.75e5 * t9i12 * exp(-10.062*t9i)
      dee  = ee*(-0.5d0*t9i + 10.062*t9i2)

      term    = xx + cc + dd + ee
      dtermdt = dxx + dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.50e9 * t932 * exp(-84.673*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.673*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end








      subroutine rate_weaknp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,cc,dcc, &
                       zm1,zm2,zm3,zm4,zm5, &
                       c1,c2
      parameter        (c1 = 1.0d0/5.93d0, &
                        c2 = 0.98d0/885.7d0)


! free decay of neutrons, n(e-nu)p and p(e-,nu)n
! fit formula from schramm and wagoner annual review 1977
! currently accepted best value for the neutron lifetime,
! 886.7 (+/- 1.9) seconds. P.R. Huffman et al., Nature, 6 January 2000.
! world average 885.7 +/- 0.8.

      zm1   = t9 * c1
      zm2   = zm1*zm1
      zm3   = zm1*zm2
      zm4   = zm1*zm3
      zm5   = zm1*zm4

      aa   = 27.512*zm5 + 36.492*zm4 + 11.108*zm3 &
             - 6.382*zm2 + 0.565*zm1 + 1.0d0
      daa  = (5.0d0*27.512*zm4 + 4.0d0*36.492*zm3 + 3.0d0*11.108*zm2 &
             - 2.0d0*6.382*zm1 + 0.565)*c1


! n=>p
      fr    = c2 * aa
      dfrdt = c2 * daa * 1.0d-9
      dfrdd = 0.0d0


      aa  = 27.617*zm5 + 34.181*zm4 + 18.059*zm3 &
            - 16.229*zm2 + 5.252*zm1
      daa = (5.0d0*27.617*zm4 + 4.0d0*34.181*zm3 + 3.0d0*18.059*zm2 &
            - 2.0d0*16.229*zm1 + 5.252)*c1

      bb = exp(-2.531d0/zm1)
      dbb = bb*2.531d0/zm2*c1

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

! p=>n
      rr    = c2 * cc
      drrdt = c2 * dcc * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_dpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! d(p,n)2p
       aa   = 3.35e7 * exp(-3.720*t9i13)
       daa  = aa*oneth*3.720d0*t9i43

       bb   = 1.0d0 + 0.784*t913 + 0.346*t923 + 0.690*t9
       dbb  = oneth*0.784*t9i23 + twoth*0.346*t9i13 + 0.690

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

! rate
      cc = exp(-25.815*t9i)
      dcc = cc*25.815*t9i2

      fr    = den * cc * term
      dfrdt = den * (dcc*term + cc*dtermdt) * 1.0d-9
      dfrdd = cc * term

      rev      =  4.24e-10 * t9i32
      drevdt   = -1.5d0*rev*t9i

      rr    = den**2 * rev * term
      drrdt = den**2 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0*den * rev * term

      return
      end




      subroutine rate_dng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,c1
      parameter        (c1 = 66.2d0*18.9d0)


! d(n,g)t
      term    = 66.2 * (1.0d0  + 18.9*t9)
      dtermdt = c1

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      =  1.63e+10 * t9i32 * exp(-72.62*t9i)
      drevdt   = rev*(-1.5d0*t9i + 72.62*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ddp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



! d(d,p)t
      aa  = 4.13e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 4.39e-2*t923 + 3.01e-2*t9 &
            + 0.543*t943 + 0.946*t953
      dbb = oneth*0.098*t9i23 + twoth*4.39e-2*t9i13 + 3.01e-2 &
            + fourth*0.543*t913 + fiveth*0.946*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.73 * exp(-46.798*t9i)
      drevdt   = rev*46.798*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ddn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



! d(d,n)he3
      aa  = 3.88e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 0.418*t923 + 0.287*t9 &
            + 0.638*t943 + 1.112*t953
      dbb = oneth*0.098*t9i23 + twoth*0.418*t9i13 + 0.287 &
            + fourth*0.638*t913 + fiveth*1.112*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.730 * exp(-37.935*t9i)
      drevdt   = rev*37.935*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_tpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! t(p,n)he3
      term    = 7.07e8 * (1.0d0 - 0.15*t912 + 0.098*t9)
      dtermdt = 7.07e8 * (-0.5d0*0.15*t9i12 + 0.098)

      aa  = exp(-8.863*t9i)
      daa = aa*8.863*t9i2

! rate
      fr    = den * aa * term
      dfrdt = den * (daa*term + aa*dtermdt) * 1.0d-9
      dfrdd = aa * term

      rev      = 0.998
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_ddg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! d(d,g)he4
      aa  = 4.84e+01 * t9i23 * exp(-4.258*t9i13)
      daa = aa*(-twoth*t9i + oneth*4.258*t9i43)

      bb  = 1.0d0 + 0.098*t913 - 0.203*t923 - 0.139*t9 &
            + 0.106*t943 + 0.185*t953
      dbb = oneth*0.098*t9i23 - twoth*0.203*t9i13 - 0.139 &
            + fourth*0.106*t913 + fiveth*0.185*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.53e+10 * t932 * exp(-276.729*t9i)
      drevdt   = rev*(1.5d0*t9i + 276.729*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_tpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! t(p,g)he4
      aa  = 2.20e+04 * t9i23 * exp(-3.869*t9i13)
      daa = aa*(-twoth*t9i + oneth*3.869*t9i43)

      bb  = 1. + 0.108*t913 + 1.68*t923 + 1.26*t9 &
            + 0.551*t943 + 1.06*t953
      dbb = oneth*0.108*t9i23 + twoth*1.68*t9i13 + 1.26 &
            + fourth*0.551*t913 + fiveth*1.06*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61e+10 * t932 * exp(-229.932*t9i)
      drevdt   = rev*(1.5d0*t9i + 229.932*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_tdn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,q1
      parameter        (q1 = 1.0d0/0.0144d0)


! t(d,n)he4 ; the "dt" reaction
      aa  = 8.09e+10 * t9i23 * exp(-4.524*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*4.524*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.092*t913 + 1.80*t923 + 1.16*t9 &
            + 10.52*t943 + 17.24*t953
      dbb = oneth*0.092*t9i23 + twoth*1.80*t9i13 + 1.16 &
            + fourth*10.52*t913 + fiveth*17.24*t923

      cc  = 8.73e+08 * t9i23 * exp(-0.523*t9i)
      dcc = -twoth*cc*t9i + cc*0.523*t9i2

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54*exp(-204.117*t9i)
      drevdt   = rev*204.117*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_tt2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! t(t,2n)he4
      aa  = 1.67e+09 * t9i23 * exp(-4.872*t9i13)
      daa = aa*(-twoth*t9i + oneth*4.872*t9i43)

      bb  = 1.0d0 + 0.086*t913 - 0.455*t923 - 0.272*t9 &
            + 0.148*t943 + 0.225*t953
      dbb = oneth*0.086*t9i23 - twoth*0.455*t9i13 - 0.272 &
            + fourth*0.148*t913 + fiveth*0.225*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.38e-10 * t9i32 * exp(-131.504*t9i)
      drevdt   = rev*(-1.5d0*t9i + 131.504*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end





      subroutine rate_he3dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,q1
      parameter        (q1 = 1.0d0/0.099225d0)


! he3(d,p)he4
      aa  = 5.86e+10 * t9i23 * exp(-7.181*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*7.181*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.058*t913 + 0.142*t923 + 0.0578*t9 &
            + 2.25*t943 + 2.32*t953
      dbb = oneth*0.058*t9i23 + twoth*0.142*t9i13 + 0.0578 &
            + fourth*2.25*t913 + fiveth*2.32*t923

      cc  = 4.36e+08 * t9i12 * exp(-1.72*t9i)
      dcc = -0.5d0*cc*t9i + cc*1.72*t9i2

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.55*exp(-212.980*t9i)
      drevdt   = rev*212.980*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3td(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(t,d)he4
      aa       = 1.0d0 + 0.128*t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a*0.128)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 5.46e+09 * t9a56 * t9i32 * exp(-7.733/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 7.733/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.60*exp(-166.182*t9i)
      drevdt   = rev*166.182*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(t,np)he4
      aa       = 1.0d0 + 0.115d0 * t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a * 0.115d0)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 7.71d9 * t9a56 * t9i32 * exp(-7.733d0/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 7.733d0/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39d-10 * t9i32 * exp(-140.367*t9i)
      drevdt   = rev*(-1.5d0 * t9i + 140.367d0 * t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_he4npg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he4(np,g)li6
      aa  = 4.62e-6 * t9i2 * exp(-19.353*t9i)
      daa = aa*(-2.0d0*t9i + 19.353*t9i2)

      bb  = 1.0d0 + 0.075*t9
      dbb = 0.075

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.22e19 * t93 * exp(-42.933*t9i)
      drevdt   = rev*(3.0d0*t9i + 42.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4dg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! he4(d,g)li6
      aa  = 3.01e1 * t9i23 * exp(-7.423*t9i13)
      daa = aa*(-twoth*t9i + oneth*7.423*t9i43)

      bb  = 1.0d0 + 0.056*t913 - 4.85*t923 + 8.85*t9 &
            - 0.585*t943 - 0.584*t953
      dbb = oneth*0.056*t9i23 - twoth*4.85*t9i13 + 8.850 &
            - fourth*0.585*t913 - fiveth*0.584*t923

! rate goes negative for t9 greater than about 15, so try this
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc =  8.55e1 * t9i32 * exp(-8.228*t9i)
      dcc = cc*(-1.5d0*t9i + 8.228*t9i2)

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.53e10 * t932 * exp(-17.1180*t9i)
      drevdt   = rev*(1.5d0*t9i + 17.1180*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he4tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,t9a,dt9a,t9a32,dt9a32,zz


! he4(t,n)li6
      aa     = 1.0d0 + 49.180*t9
      zz     = 1.0d0/aa

      t9a   = t9*zz
      dt9a   = (1.0d0 - t9a*49.180)*zz

      t9a32  = t9a * sqrt(t9a)
      dt9a32 = 1.5d0*t9a32/t9a * dt9a

      aa     = 1.80e8 * exp(-55.4940*t9i)
      daa    = aa*55.4940*t9i2

      bb     = 1.0d0 - 0.2610 * t9a32 * t9i32
      dbb    = -0.2610*(-1.5d0*t9a32*t9i52 + dt9a32*t9i32)

      cc     = aa*bb
      dcc    = daa*bb + aa*dbb

      dd     = 2.72e9 * t9i32 * exp(-57.8840*t9i)
      ddd    = dd*(-1.5d0*t9i + 57.8840*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.72e9 * t9i32 * exp(-2.39*t9i)
      drevdt   = rev*(-1.5d0*t9i + 2.39*t9i2)

      term = 0.935*(1.80e8*bb + rev)
      dtermdt = 0.935*(1.80e8*dbb + drevdt)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li6phe3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/30.25d0)


! li6(p,he3)he4
      aa  = 3.73e10 * t9i23 * exp(-8.413*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*8.413*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.050*t913 - 0.061*t923 - 0.0210*t9 &
            + 0.0060*t943 + 0.0050*t953
      dbb = oneth*0.050*t9i23 - twoth*0.061*t9i13 - 0.0210 &
            + fourth*0.0060*t913 + fiveth*0.0050*t923

      cc  = 1.33e10 * t9i32 * exp(-17.7630*t9i)
      dcc = -1.5d0*cc*t9i + cc*17.7630*t9i2

      dd  =  1.29e9 * t9i * exp(-21.82*t9i)
      ddd = -dd*t9i + dd*21.82*t9i2

      term    = aa * bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.070d0 * exp(-46.6310*t9i)
      drevdt   = rev*46.6310*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li6ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li6(n,g)li7
! malaney-fowler 1989

      term    = 5.10e3
      dtermdt = 0.0d0

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e10 * t932 * exp(-84.17*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.17*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4tg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he4(t,g)li7
      aa  = 8.67e5 * t9i23 * exp(-8.08*t9i13)
      daa = aa*(-twoth*t9i + oneth*8.08*t9i43)

      bb  = 1.0d0 + 0.052*t913 - 0.448*t923 - 0.165*t9 &
            + 0.144*t943 + 0.134*t953
      dbb = oneth*0.052*t9i23 - twoth*0.448*t9i13 - 0.165 &
            + fourth*0.144*t913 + fiveth*0.134*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e10 * t932 * exp(-28.64*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.64*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_li7dn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(d,n)2a
      term    = 2.92e11 * t9i23 * exp(-10.259*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*10.259*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.95e-10 * t9i32 * exp(-175.476*t9i)
      drevdt   = rev*(-1.5d0*t9i + 175.476*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_li7tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(t,n)be9
! malaney and fowler (apjl, 345, l5, 1989)
      term    = 1.46d+11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = -twoth*term*t9i + oneth*term*11.333*t9i43

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li7t2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(t,2n)2a
      term    = 8.81e11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*11.333*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.22e-19 * t9i3 * exp(-102.864*t9i)
      drevdt   = rev*(-3.0d0*t9i + 102.864*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_li7he3np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(he3,np)2a
      term    = 1.11e13 * t9i23 * exp(-17.989*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*17.989*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.09e-20 * t9i3 * exp(-111.727*t9i)
      drevdt   = rev*(-3.0d0*t9i + 111.727*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_li6pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! li6(p,g)be7
      if (t9 .gt. 10.0) then
       t9a  = 1.0d0
       dt9a = 0.0d0
      else
       aa   = 1.0d0 - 0.0969*t9

       bb   = aa**(-twoth)
       dbb  = twoth*bb/aa*0.0969

       cc   = aa + 0.0284*t953*bb
       dcc  = -0.0969 + 0.0284*(fiveth*t923*bb + t953*dbb)

       zz   = 1.0d0/cc
       t9a  = t9*zz
       dt9a = (1.0d0 - t9a*dcc)*zz
      end if

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      term    = 6.69e+05 * t9a56 * t9i32 * exp(-8.413/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 8.413/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-65.054*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.054*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0


      return
      end




      subroutine rate_li7pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! li7(p,n)be7
      aa  = 5.15e+09 * exp(-1.167*t913 - 19.081*t9i)
      daa = aa*(-oneth*1.167*t9i23 + 19.081*t9i2)

      bb  = 7.84e+09 * t9i32 * exp(-22.832*t9i)
      dbb = -1.5d0*bb*t9i + bb*22.832*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      aa  = 5.15e+09 * exp(-1.167*t913)
      daa = -aa*oneth*1.167*t9i23

      bb  = 0.998 * 7.84e+09 * t9i32 * exp(-3.751*t9i)
      dbb = -1.5d0*bb*t9i + bb*3.751*t9i2

      term    = aa + bb
      dtermdt = daa + dbb
      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li7ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! li7(n,g)li8
! apj 372, 1
      aa      = 4.26d+03 * t9i32 * exp(-2.576*t9i)
      daa     = aa*(-1.5d0*t9i + 2.576*t9i2)

      term    = 3.144d+03 + aa
      dtermdt = daa

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.2923d+10 * t932 * exp(-2.359d+01*t9i)
      drevdt = rev*(1.5d0*t9i + 2.359d+01*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be7dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(d,p)2a
      term    = 1.07e12 * t9i23 * exp(-12.428*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*12.428*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.97e-10 * t9i32 * exp(-194.557*t9i)
      drevdt   = rev*(-1.5d0*t9i + 194.557*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_be7tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(t,np)2a
      term    = 2.91e12 * t9i23 * exp(-13.729*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*13.729*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.09e-20 * t9i3 * exp(-121.944*t9i)
      drevdt   = rev*(-3.0d0*t9i + 121.944*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end



      subroutine rate_be7he32p(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(he3,2p)2a
      term    = 6.11e13 * t9i23 * exp(-21.793*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*21.793*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.22e-19 * t9i3 * exp(-130.807*t9i)
      drevdt   = rev*(-3.0d0*t9i + 130.807*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_be9pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.2704d0)


! be9(p,a)li6
      aa  = 2.11e11 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.04*t913 + 1.09*t923 + 0.307*t9 &
            + 3.21*t943 + 2.30*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307 &
             + fourth*3.21*t913 + fiveth*2.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.51e8 * t9i * exp(-3.046*t9i)
      ddd  = -dd*t9i + dd*3.046*t9i2

      ee   = 6.70e8 * t9i34 * exp(-5.160*t9i)
      dee  = -0.75d0*ee*t9i + ee*5.160*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.18e-1 * exp(-24.674*t9i)
      drevdt   = rev*24.674*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li6ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/1.758276d0)


! li6(a,g)b10
      aa  = 4.06e6 * t9i23 * exp(-18.790*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*18.790*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.022*t913 + 1.54*t923 + 0.239*t9 &
             +  2.20*t943 + 0.869*t953
      dbb  = oneth*0.022*t9i23 + twoth*1.54*t9i13 + 0.239 &
             + fourth*2.20*t913 + fiveth*0.869*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1910.0 * t9i32 * exp(-3.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.484*t9i2)

      ee   = 1.01e4 * t9i * exp(-7.269*t9i)
      dee  = ee*(-t9i + 7.269*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.58e10 * t932 * exp(-51.753*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.753*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_li7an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt

! li7(a,n)b10
      term    = 3.84e8 * exp(-32.382*t9i)
      dtermdt = term*32.382*t9i2

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.84e8*1.32
      drevdt   = 0.0d0

      rr    = den * rev
      drrdt = 0.0d0
      drrdd = rev

      return
      end





      subroutine rate_be9pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.715716d0)


! be9(p,g)b10
      aa  = 1.33e7 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.04*t913 + 1.52*t923 + 0.428*t9 &
             + 2.15*t943 + 1.54*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.52*t9i13 + 0.428 &
             + fourth*2.15*t913 + fiveth*1.54*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.64e4 * t9i32 * exp(-3.445*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.445*t9i2)

      ee   = 2.72e6 * t9i32 * exp(-10.62*t9i)
      dee  = ee*(-1.5d0*t9i + 10.62*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.73e9 * t932 * exp(-76.427*t9i)
      drevdt   = rev*(1.5d0*t9i + 76.427*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b10pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/19.377604d0)

! b10(p,a)li7
      aa  = 1.26e11 * t9i23 * exp(-12.062*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*12.062*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 - 0.498*t923 - 0.121*t9 &
             + 0.3*t943 + 0.184*t953
      dbb  = oneth*0.035*t9i23 - twoth*0.498*t9i13 - 0.121 &
             + fourth*0.3*t913 + fiveth*0.184*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.59e9 * t9i * exp(-12.260*t9i)
      ddd  = -dd*t9i + dd*12.260*t9i2

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.54e-01 * exp(-13.301*t9i)
      drevdt   = rev*13.301*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/17.598025d0)


! li7(a,g)b11
      aa  = 3.55e7 * t9i23 * exp(-19.161*t9i13 -t92*q1)
      daa = aa*(-twoth*t9i + oneth*19.161*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.022*t913 + 0.775*t923 + 0.118*t9 &
             + 0.884*t943 + 0.342*t953
      dbb  = oneth*0.022*t9i23 + twoth*0.775*t9i13 + 0.118 &
             + fourth*0.884*t913 + fiveth*0.342*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.33e2 * t9i32 * exp(-2.977*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.977*t9i2)

      ee   = 4.10e4 * t9i * exp(-6.227*t9i)
      dee  = ee*(-t9i + 6.227*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e10 * t932 * exp(-100.538*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.538*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b11pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/2.702736d0)

! b11(p,a)be8=>2a
      aa  = 2.20e12 * t9i23 * exp(-12.095*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0  + 0.034*t913 + 0.14*t923 + 0.034*t9 &
             + 0.19*t943 + 0.116*t953
      dbb  = oneth*0.034*t9i23 + twoth*0.14*t9i13 + 0.034 &
             + fourth*0.19*t913 + fiveth*0.116*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb


      dd   = 4.03e6 * t9i32 * exp(-1.734*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.734*t9i2)

      ee   = 6.73e9 * t9i32 * exp(-6.262*t9i)
      dee  = ee*(-1.5d0*t9i + 6.262*t9i2)

      ff   = 3.88e9*t9i * exp(-14.154*t9i)
      dff  = ff*(-t9i + 14.154*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd +dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.5e-10* t9i32 *exp(-100.753*t9i)
      drevdt   = rev*(-1.5d0*t9i + 100.753*t9i2)

      rr    = den * den * rev * term 
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term 

      return
      end




      subroutine rate_be7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/22.743361d0)

! be7(a,g)c11
      aa  = 8.45e+07 * t9i23 * exp(-23.212*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*23.212*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.018*t913 + 0.488*t923 + 0.061*t9 &
             + 0.296*t943 + 0.095*t953
      dbb  = oneth*0.018*t9i23 + twoth*0.488*t9i13 + 0.061 &
             + fourth*0.296*t913 + fiveth*0.095*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.25e+04 * t9i32 * exp(-6.510*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.510*t9i2)

      ee   = 1.29e+05 * t9i54 * exp(-10.039*t9i)
      dee  = ee*(-1.25d0*t9i + 10.039*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e+10 * t932 * exp(-87.539*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.539*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_b11pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! b11(p,n)c11
      aa  = 1.69e8*(1.0d0 - 0.048*t912 + 0.010*t9)
      daa = 1.69e8*(-0.5d0*0.048*t9i12 + 0.010)

      bb  = exp(-32.080*t9i)
      dbb = bb*32.080*t9i2

      term    = aa*bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 0.998*aa
      dtermdt  = 0.998*daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b8ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh, &
                       pp,dpp,qq,dqq,drr


! b8(a,p)c11
      aa  = 1.67e-09 * t9i32 * exp(-1.079*t9i)
      daa = -1.5d0*aa*t9i + aa*1.079*t9i2

      bb  = 9.55e-08 * t9i32 * exp(-1.311*t9i)
      dbb = -1.5d0*bb*t9i + bb*1.311*t9i2

      cc  = 1.98e-01 * t9i32 * exp(-2.704*t9i)
      dcc = -1.5d0*cc*t9i + cc*2.704*t9i2

      dd  = 1.34e+00 * t9i32 * exp(-4.282*t9i)
      ddd = -1.5d0*dd*t9i + dd*4.282*t9i2

      ee  = 3.22e+04 * t9i32 * exp(-6.650*t9i)
      dee = -1.5d0*ee*t9i + ee*6.650*t9i2

      ff  = 2.33e+05 * t9i32 * exp(-8.123*t9i)
      dff = -1.5d0*ff*t9i + ff*8.123*t9i2

      gg  = 2.55e+06 * t9i32 * exp(-11.99*t9i)
      dgg = -1.5d0*gg*t9i + gg*11.99*t9i2

      hh  = 9.90e+06 * t9i32 * exp(-13.50*t9i)
      dhh = -1.5d0*hh*t9i + hh*13.50*t9i2

      pp  = 1.41e+06 * t9i32 * exp(-16.51*t9i)
      dpp = -1.5d0*pp*t9i + pp*16.51*t9i2

      qq  = 1.99e+07 * t9i32 * exp(-18.31*t9i)
      dqq = -1.5d0*qq*t9i + qq*18.31*t9i2

      rr  = 6.01e+07 * t9i32 * exp(-20.63*t9i)
      drr = -1.5d0*rr*t9i + rr*20.63*t9i2

      term    = aa + bb + cc + dd + ee + ff + gg + hh + pp + qq + rr
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh &
                + dpp + dqq + drr

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.101 * exp(-85.95*t9i)
      drevdt   = rev*85.95*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b10pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/19.377604d0)


! b10(p,g)c11
      aa  = 4.61e+05 * t9i23 * exp(-12.062*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.062*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 + 0.426*t923 + 0.103*t9 &
             + 0.281*t943 + 0.173*t953
      dbb  = oneth*0.035*t9i23 + twoth*0.426*t9i13 + 0.103 &
             + fourth*0.281*t913 + fiveth*0.173*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.93e+05 * t9i32 * exp(-12.041*t9i)
      ddd  = dd*(-1.5d0*t9i + 12.041*t9i2)

      ee   = 1.14e+04 * t9i32 * exp(-16.164*t9i)
      dee  = ee*(-1.5d0*t9i + 16.164*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+10 * t932 * exp(-100.840*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.840*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_c11na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd


! c11(n,a)be8=>2a   get a real rate in here

!      fr    = den * 7.0e4
!      dfrdt = 0.0d0
!      dfrdd = 7.0e4

      fr    = 0.0d0
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh


! be9(a,n)c12
! Wrean 94 Phys Rev C (1994) vol 49, #2, 1205
      aa  = 6.476d+13 * t9i23 * exp(-23.8702*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*23.8702*t9i43

      bb  = (1.0d0 - 0.3270*t913)
      dbb = -oneth*0.3270*t9i23

! rate goes negative for t9 greater than about 15, so try this
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 6.044d-3*t9i32*exp(-1.401*t9i)
      ddd = -1.5d0*dd*t9i + dd*1.401*t9i2

      ee  = 7.268*t9i32*exp(-2.063*t9i)
      dee = -1.5d0*ee*t9i + ee*2.063*t9i2

      ff  = 3.256d+4*t9i32*exp(-3.873*t9i)
      dff = -1.5d0*ff*t9i + ff*3.873*t9i2

      gg  = 1.946d+5*t9i32*exp(-4.966*t9i)
      dgg = -1.5d0*gg*t9i + gg*4.966*t9i2

      hh  = 1.838e9*t9i32*exp(-15.39*t9i)
      dhh = -1.5d0*hh*t9i + hh*15.39*t9i2

      term    = cc + dd + ee + ff + gg + hh
      dtermdt = dcc + ddd + dee + dff + dgg + dhh


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 10.3d0 * exp(-66.160*t9i)
      drevdt   = rev*66.160d0*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term


      return
      end




      subroutine rate_b11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.57121d0)


! b11(p,g)c12
      aa  = 4.62e+07 * t9i23 * exp(-12.095*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 + 3.0*t923 + 0.723*t9 &
             + 9.91*t943 + 6.07*t953
      dbb  = oneth*0.035*t9i23 + twoth*3.0*t9i13 + 0.723 &
             + fourth*9.91*t913 + fiveth*6.07*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 7.89e+03 * t9i32 * exp(-1.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.733*t9i2)

      ee   = 9.68e+04 * t9i15 * exp(-5.617*t9i)
      dee  = ee*(-0.2d0*t9i + 5.617*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.01e+10 * t932 * exp(-185.173*t9i)
      drevdt   = rev*(1.5d0*t9i + 185.173*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_b11ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.120409d0)


! b11(a,p)c14
      aa  = 5.37e+11 * t9i23 * exp(-28.234*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.015*t913 + 5.575*t923 + 0.576*t9 &
             + 15.888*t943 + 4.174*t953
      dbb  = oneth*0.015*t9i23 + twoth*5.575*t9i13 + 0.576 &
             + fourth*15.888*t913 + fiveth*4.174*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.44e-03 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 3.36e+02 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 5.32e+06 * t9i38 * exp(-11.617*t9i)
      dff  = ff*(-0.375*t9i + 11.617*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10e+01 * exp(-9.098*t9i)
      drevdt   = rev*9.098*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_pp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! p(p,e+nu)d
      if (t9 .le. 3.0) then
       aa   = 4.01d-15 * t9i23 * exp(-3.380d0*t9i13)
       daa  = aa*(-twoth*t9i + oneth*3.380d0*t9i43)

       bb   = 1.0d0 + 0.123d0*t913 + 1.09d0*t923 + 0.938d0*t9
       dbb  = oneth*0.123d0*t9i23 + twoth*1.09d0*t9i13 + 0.938d0

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 1.1581136d-15
       dtermdt = 0.0d0
      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_pep(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! p(e-p,nu)d
      if (t9 .le. 3.0) then

       aa   = 1.36e-20 * t9i76 * exp(-3.380*t9i13)
       daa  = aa*(-sevsix*t9i + oneth*3.380d0*t9i43)

       bb   = (1.0d0 - 0.729d0*t913 + 9.82d0*t923)
       dbb  = -oneth*0.729d0*t9i23 + twoth*9.82d0*t9i13

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 7.3824387e-21
       dtermdt = 0.0d0
      end if

! rate
      fr    = ye * den * den * term
      dfrdt = ye * den * den * dtermdt * 1.0d-9
      dfrdd = ye * 2.0d0 * den * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_hep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! he3(p,e+nu)he4
      if (t9 .le. 3.0) then

       aa   = 8.78e-13 * t9i23 * exp(-6.141d0*t9i13)
       daa  = aa*(-twoth*t9i + oneth*6.141d0*t9i43)

       term    = aa
       dtermdt = daa

      else
       term    = 5.9733434e-15
       dtermdt = 0.0d0
      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end







      subroutine rate_png(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! p(n,g)d
! smith,kawano,malany 1992

       aa      = 1.0d0 - 0.8504*t912 + 0.4895*t9 &
                 - 0.09623*t932 + 8.471e-3*t92 &
                 - 2.80e-4*t952

       daa     =  -0.5d0*0.8504*t9i12 + 0.4895 &
                 - 1.5d0*0.09623*t912 + 2.0d0*8.471e-3*t9 &
                 - 2.5d0*2.80e-4*t932

       term    = 4.742e4 * aa
       dtermdt = 4.742e4 * daa


! wagoner,schramm 1977
!      aa      = 1.0d0 - 0.86*t912 + 0.429*t9
!      daa     =  -0.5d0*0.86*t9i12 + 0.429

!      term    = 4.4d4 * aa
!      dtermdt = 4.4d4 * daa



! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.71d+09 * t932 * exp(-25.82*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.82*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_dpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! d(p,g)he3
      aa      = 2.24d+03 * t9i23 * exp(-3.720*t9i13)
      daa     = aa*(-twoth*t9i + oneth*3.720*t9i43)

      bb      = 1.0d0 + 0.112*t913 + 3.38*t923 + 2.65*t9
      dbb     = oneth*0.112*t9i23 + twoth*3.38*t9i13 + 2.65

      term    = aa * bb
      dtermdt = daa * bb + aa * dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.63d+10 * t932 * exp(-63.750*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.750*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he3ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! he3(n,g)he4
      term    = 6.62 * (1.0d0 + 905.0*t9)
      dtermdt = 5.9911d3

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61d+10 * t932 * exp(-238.81*t9i)
      drevdt   = rev*(1.5d0*t9i + 238.81*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he3he3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he3(he3,2p)he4
      aa   = 6.04d+10 * t9i23 * exp(-12.276*t9i13)
      daa  = aa*(-twoth*t9i + oneth*12.276*t9i43)

      bb   = 1.0d0 + 0.034*t913 - 0.522*t923 - 0.124*t9 &
             + 0.353*t943 + 0.213*t953
      dbb  = oneth*0.034*t9i23 - twoth*0.522*t9i13 - 0.124 &
             + fourth*0.353*t913 + fiveth*0.213*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39e-10 * t9i32 * exp(-149.230*t9i)
      drevdt   = rev*(-1.5d0*t9i + 149.230*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end





      subroutine rate_he3he4(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(he4,g)be7
      aa      = 1.0d0 + 0.0495*t9
      daa     = 0.0495

      zz      = 1.0d0/aa
      t9a     = t9*zz
      dt9a    = (1.0d0 - t9a*daa)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 5.61d+6 * t9a56 * t9i32 * exp(-12.826/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 12.826/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e+10 * t932 * exp(-18.423*t9i)
      drevdt   = rev*(1.5d0*t9i + 18.423*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7em(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb

! be7(e-,nu+g)li7
      if (t9 .le. 3.0) then
       aa  = 0.0027 * t9i * exp(2.515e-3*t9i)
       daa = -aa*t9i - aa*2.515e-3*t9i2

       bb  = 1.0d0 - 0.537*t913 + 3.86*t923 + aa
       dbb = -oneth*0.537*t9i23 + twoth*3.86*t9i13 + daa

       term    = 1.34e-10 * t9i12 * bb
       dtermdt = -0.5d0*term*t9i + 1.34e-10*t9i12*dbb

      else
       term    = 0.0d0
       dtermdt = 0.0d0
      endif

! rates
      fr    = ye * den * term
      dfrdt = ye * den * dtermdt * 1.0d-9
      dfrdd = ye * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! be7(p,g)b8
      aa      = 3.11e+05 * t9i23 * exp(-10.262*t9i13)
      daa     = aa*(-twoth*t9i + oneth*10.262*t9i43)

      bb      = 2.53e+03 * t9i32 * exp(-7.306*t9i)
      dbb     = bb*(-1.5d0*t9i + 7.306*t9i2)

      term    = aa + bb
      dtermdt = daa + dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.30e+10 * t932 * exp(-1.595*t9i)
      drevdt   = rev*(1.5d0*t9i + 1.595*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_li7pag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz, &
                       term1,dterm1,term2,dterm2,rev1,drev1dt,rev2,drev2dt,q1
      parameter        (q1 = 1.0d0/2.876416d0)



! 7li(p,g)8be=>2a
      aa   = 1.56d5 * t9i23 * exp(-8.472d0 * t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*8.472*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.049d0*t913 + 2.498d0*t923 + 0.86d0*t9 &
             + 3.518d0*t943 + 3.08*t953
      dbb  = oneth*0.049*t9i23 + twoth*2.498*t9i13 + 0.86 &
             + fourth*3.518*t913 + fiveth*3.08*t923

      cc   = aa*bb
      dcc  = daa*bb + aa*dbb

      dd   =  1.55d6 * t9i32 * exp(-4.478d0 * t9i)
      ddd  = dd*(-1.5d0*t9i + 4.478*t9i2)

      term1  = cc + dd
      dterm1 = dcc + ddd

      rev1    = 6.55e+10 * t932 * exp(-200.225*t9i)
      drev1dt = rev1*(1.5d0*t9i + 200.225*t9i2)


! 7li(p,a)a
      aa     = 1.0d0 + 0.759*t9

      zz     = 1.0d0/aa
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*0.759)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa     = 1.096e+09 * t9i23 * exp(-8.472*t9i13)
      daa    = aa*(-twoth*t9i + oneth*8.472*t9i43)

      bb     = -4.830e+08 * t9a56 * t9i32 * exp(-8.472/t9a13)
      dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i + 8.472/t9a13**2*dt9a13)

      cc     = 1.06e+10 * t9i32 * exp(-30.442*t9i)
      dcc    = cc*(-1.5d0*t9i + 30.442*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      rev2    = 4.69 * exp(-201.291*t9i)
      drev2dt = rev2*201.291*t9i2


! sum the two forward rates (per f35 cf88)

      term    = term1 + term2
      dtermdt = dterm1 + dterm2


! per cf88 the reverse rate just the second one
      rev     = rev2
      drevdt  = drev2dt


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_b8ep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 0.77d0, &
                        con      = lntwo/halflife)


! b8(e+,nu)be8 => 2a

      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/2.25d0)


! c12(p,g)13n
      aa   = 2.04e+07 * t9i23 * exp(-13.69*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*13.69*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.03*t913 + 1.19*t923 + 0.254*t9 &
             + 2.06*t943 + 1.12*t953
      dbb  = oneth*0.03*t9i23 + twoth*1.19*t9i13 + 0.254 &
             + fourth*2.06*t913 + fiveth*1.12*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.08e+05 * t9i32 * exp(-4.925*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.925*t9i2)

      ee   = 2.15e+05 * t9i32 * exp(-18.179*t9i)
      dee  = ee*(-1.5d0*t9i + 18.179*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.84e+09 * t932 * exp(-22.553*t9i)
      drevdt   = rev*(1.5d0*t9i + 22.553*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n13em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 597.9d0, &
                        con      = lntwo/halflife)

! n13(e-nu)c13
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/4.0d0)


! c13(p,g)13n
      aa   = 8.01e+07 * t9i23 * exp(-13.717*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*13.717*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.030*t913 + 0.958*t923 + 0.204*t9 &
             + 1.39*t943 + 0.753*t953
      dbb  = oneth*0.030*t9i23 + twoth*0.958*t9i13 + 0.204 &
             + fourth*1.39*t913 + fiveth*0.753*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.21e+06 * t9i65 * exp(-5.701*t9i)
      ddd  = dd*(-sixfif*t9i + 5.701*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-87.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/10.850436d0)


! n14(p,g)o15
      aa  = 4.90e+07 * t9i23 * exp(-15.228*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.228*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.027*t913 - 0.778*t923 - 0.149*t9 &
             + 0.261*t943 + 0.127*t953
      dbb  = oneth*0.027*t9i23 - twoth*0.778*t9i13 - 0.149 &
             + fourth*0.261*t913 + fiveth*0.127*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.37e+03 * t9i32 * exp(-3.011*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.011*t9i2)

      ee   = 2.19e+04 * exp(-12.530*t9i)
      dee  = ee*12.530*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 2.70e+10 * t932 * exp(-84.678*t9i)
      drevdt = rev*(1.5d0*t9i + 84.678*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_o15em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 122.24d0, &
                        con      = lntwo/halflife)

! o15(e-nu)n15
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.776161d0)


! n14(a,g)f18
      aa  = 7.78d+09 * t9i23 * exp(-36.031*t9i13- t92*q1)
      daa = aa*(-twoth*t9i + oneth*36.031*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.012*t913 + 1.45*t923 + 0.117*t9 &
             + 1.97*t943 + 0.406*t953
      dbb  = oneth*0.012*t9i23 + twoth*1.45*t9i13 + 0.117 &
             + fourth*1.97*t913 + fiveth*0.406*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.36d-10 * t9i32 * exp(-2.798*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.798*t9i2)

      ee   = 2.03 * t9i32 * exp(-5.054*t9i)
      dee  = ee*(-1.5d0*t9i + 5.054*t9i2)

      ff   = 1.15d+04 * t9i23 * exp(-12.310*t9i)
      dff  = ff*(-twoth*t9i + 12.310*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-51.236*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.236*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.2025d0)


! n15(p,g)o16
      aa  = 9.78e+08 * t9i23 * exp(-15.251*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0  + 0.027*t913 + 0.219*t923 + 0.042*t9 &
             + 6.83*t943 + 3.32*t953
      dbb  = oneth*0.027*t9i23 + twoth*0.219*t9i13 + 0.042 &
             + fourth*6.83*t913 + fiveth*3.32*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.11e+04*t9i32*exp(-3.328*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.328*t9i2)

      ee   = 1.49e+04*t9i32*exp(-4.665*t9i)
      dee  = ee*(-1.5d0*t9i + 4.665*t9i2)

      ff   = 3.8e+06*t9i32*exp(-11.048*t9i)
      dff  = ff*(-1.5d0*t9i + 11.048*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.62e+10 * t932 * exp(-140.734*t9i)
      drevdt   = rev*(1.5d0*t9i + 140.734*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.272484d0)


! n15(p,a)c12
      aa  = 1.08d+12*t9i23*exp(-15.251*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.027*t913 + 2.62*t923 + 0.501*t9 &
             + 5.36*t943 + 2.60*t953
      dbb  = oneth*0.027*t9i23 + twoth*2.62*t9i13 + 0.501 &
             + fourth*5.36*t913 + fiveth*2.60*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.19d+08 * t9i32 * exp(-3.676*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.676*t9i2)

      ee   = 5.41d+08 * t9i12 * exp(-8.926*t9i)
      dee  = ee*(-0.5d0*t9i + 8.926*t9i2)

      ff   = theta * 4.72d+08 * t9i32 * exp(-7.721*t9i)
      dff  = ff*(-1.5d0*t9i + 7.721*t9i2)

      gg   = theta * 2.20d+09 * t9i32 * exp(-11.418*t9i)
      dgg  = gg*(-1.5d0*t9i + 11.418*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.06d-01*exp(-57.625*t9i)
      drevdt   = rev*57.625*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o16pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz


! o16(p,g)f17
      aa  = exp(-0.728*t923)
      daa = -twoth*aa*0.728*t9i13

      bb  = 1.0d0 + 2.13 * (1.0d0 - aa)
      dbb = -2.13*daa

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923*dbb

      dd   = exp(-16.692*t9i13)
      ddd  = oneth*dd*16.692*t9i43

      zz   = 1.0d0/cc
      ee   = dd*zz
      dee  = (ddd - ee*dcc)*zz

      term    = 1.50d+08 * ee
      dtermdt = 1.50d+08 * dee


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+09*t932*exp(-6.968*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.968*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,res1,dres1,res2,dres2,res3,dres3, &
                       res4,dres4,res5,dres5,res6,dres6,zz, &
                       theta,q1,q2
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.319225d0, &
                        q2    = 1.0d0/0.0016d0)



! o17(p,a)n14
! rate from jeff blackmons thesis, includes terms from fowler 75,
! landre 1990 (a&a 240, 85), and new results
! use rev factor from cf88 rate

      aa  = 1.53d+07 * t9i23 * exp(-16.712*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*16.712*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.025*t913 + 5.39*t923 + 0.940*t9 &
             + 13.5*t943 + 5.98*t953
      dbb  = oneth*0.025*t9i23 + twoth*5.39*t9i13 + 0.940 &
             + fourth*13.5*t913 + fiveth*5.98*t923

      res1  = aa * bb
      dres1 = daa*bb + aa*dbb

      res2  = 2.92d+06 * t9 * exp(-4.247*t9i)
      dres2 = res2*(t9i + 4.247*t9i2)


      aa    = 0.479 * t923 + 0.00312
      daa   = twoth*0.479*t9i13

      bb    = aa*aa
      dbb   = 2.0d0 * aa * daa

      cc    =  1.78d+05 * t9i23 * exp(-16.669*t9i13)
      dcc   = cc*(-twoth*t9i + oneth*16.669*t9i43)

      zz    = 1.0d0/bb
      res3  = cc*zz
      dres3 = (dcc - res3*dbb)*zz

      res4  = 8.68d+10 * t9 * exp(-16.667*t9i13 - t92*q2)
      dres4 = res4*(t9i + oneth*16.667*t9i43 - 2.0d0*t9*q2)

      res5  = 9.22d-04 * t9i32 * exp(-0.767*t9i)
      dres5 = res5*(-1.5d0*t9i + 0.767*t9i2)

      res6  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dres6 = res6*(-1.5d0*t9i + 2.077*t9i2)

      term    = res1 + res2 + res3 + res4 + res5 + res6
      dtermdt = dres1 + dres2 + dres3 + dres4 + dres5 + dres6

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.676 * exp(-13.825*t9i)
      drevdt   = rev*13.825*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_o17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       zz,theta
      parameter        (theta = 0.1d0)


! o17(p,g)18f
! from landre et al 1990 a&a 240, 85
      aa     = 1.0d0 + 2.69*t9
      zz     = 1.0d0/aa

      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*2.69)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa  = 7.97d+07 * t9a56 * t9i32 * exp(-16.712/t9a13)
      daa = aa*(dt9a56/t9a56 - 1.5d0*t9i + 16.712/t9a13**2*dt9a13)

      bb  = 1.0d0  + 0.025*t913 - 0.051*t923 - 8.82d-3*t9
      dbb = oneth*0.025*t9i23 - twoth*0.051*t9i13 - 8.82d-3
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = 1.51d+08 * t9i23 * exp(-16.712*t9i13)
      dcc = cc*(-twoth*t9i + oneth*16.712*t9i43)

      dd  = bb*cc
      ddd = dbb*cc + bb*dcc

      ee  = 1.56d+5 * t9i * exp(-6.272*t9i)
      dee = ee*(-t9i + 6.272*t9i2)

      ff  = 2.0d0 * theta * 3.16d-05 * t9i32 * exp(-0.767*t9i)
      dff = ff*(-1.5d0*t9i + 0.767*t9i2)

      gg  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dgg = gg*(-1.5d0*t9i + 2.077*t9i2)

      term    = aa + dd + ee + ff + gg
      dtermdt = daa + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.66d+10 * t932 * exp(-65.061*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.061*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_o18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.852321d0)


! o18(p,a)n15
      aa  = 3.63e+11 * t9i23 * exp(-16.729*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*16.729*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.025*t913 + 1.88*t923 + 0.327*t9 &
            + 4.66*t943 + 2.06*t953
      dbb = oneth*0.025*t9i23 + twoth*1.88*t9i13 + 0.327 &
            + fourth*4.66*t913 + fiveth*2.06*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 9.90e-14 * t9i32 * exp(-0.231*t9i)
      ddd = -1.5d0*dd*t9i + dd*0.231*t9i2

      ee  = 2.66e+04 * t9i32 * exp(-1.670*t9i)
      dee = -1.5d0*ee*t9i + ee*1.670*t9i2

      ff  = 2.41e+09 * t9i32 * exp(-7.638*t9i)
      dff = -1.5d0*ff*t9i + ff*7.638*t9i2

      gg  = 1.46e+09 * t9i * exp(-8.310*t9i)
      dgg = -gg*t9i + gg*8.310*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.66e-01 * exp(-46.191*t9i)
      drevdt   = rev*46.191*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.019321d0)


! o18(p,g)19f
      aa  = 3.45e+08 * t9i23 * exp(-16.729*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*16.729*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.025*t913 + 2.26*t923 + 0.394*t9 &
            + 30.56*t943 + 13.55*t953
      dbb = oneth*0.025*t9i23 + twoth*2.26*t9i13 + 0.394 &
            + fourth*30.56*t913 + fiveth*13.55*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 1.25e-15 * t9i32 * exp(-0.231*t9i)
      ddd = dd*(-1.5d0*t9i + 0.231*t9i2)

      ee  = 1.64e+02 * t9i32 * exp(-1.670*t9i)
      dee = ee*(-1.5d0*t9i + 1.670*t9i2)

      ff  = 1.28e+04 * t912 * exp(-5.098*t9i)
      dff = ff*(0.5d0*t9i + 5.098*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.20e+09 * t932 * exp(-92.769*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.769*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_f17em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 64.49d0, &
                        con      = lntwo/halflife)

! f17(e-nu)o17
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_f18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 6586.2d0, &
                        con      = lntwo/halflife)

! f18(e-nu)o18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,q1
      parameter        (q1 = 1.0d0/0.714025d0)


! f19(p,a)o16
      aa  = 3.55d+11 * t9i23 * exp(-18.113*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*18.113*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.023*t913 + 1.96*t923 + 0.316*t9 &
            + 2.86*t943 + 1.17*t953
      dbb = oneth*0.023*t9i23 + twoth*1.96*t9i13 + 0.316 &
            + fourth*2.86*t913 + fiveth*1.17*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.67d+06 * t9i32 * exp(-3.752*t9i)
      ddd = -1.5d0*dd*t9i + dd*3.752*t9i2

      ee  = 3.07d+08 * exp(-6.019*t9i)
      dee = ee*6.019*t9i2

      ff  = 4.0*exp(-2.090*t9i)
      dff = ff*2.090*t9i2

      gg  = 7.0*exp(-16.440*t9i)
      dgg = gg*16.440*t9i2

      hh  = 1.0d0 + ff + gg
      dhh = dff + dgg

      term    = (cc + dd + ee)/hh
      dtermdt = ((dcc + ddd + dee) - term*dhh)/hh


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.54e-01 * exp(-94.159*t9i)
      drevdt   = rev*94.159*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_n13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/0.69288976d0)


! n13(p,g)o14
! Keiner et al 1993 Nucl Phys A552, 66
      aa  = -1.727d+7 * t9i23 * exp(-15.168*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.168*t9i43 -2.0d0*t9*q1)

      bb  = 1.0d0 + 0.027*t913 - 17.54*t923 - 3.373*t9 &
            + 0.0176*t943 + 0.766d-2*t953
      dbb = oneth*0.027*t9i23 - twoth*17.54*t9i13 - 3.373 &
            + fourth*0.0176*t913 + fiveth*0.766d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 3.1d+05 * t9i32 * exp(-6.348*t9i)
      ddd = dd*(-1.5d0*t9i + 6.348*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! goes negative below about t7=1.5
! note cf88 rate stays positive
      if (term .lt. 0.0) then
       term    = 0.0d0
       dtermdt = 0.0d0
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.57d+10*t932*exp(-53.706*t9i)
      drevdt   = rev*(1.5d0*t9i + 53.706*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 70.606d0, &
                        con      = lntwo/halflife)

! o14(e-nu)n14
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.514089d0)


! o14(a,p)f17
      aa  = 1.68e+13 * t9i23 * exp(-39.388*t9i13- t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*39.388*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 + 13.117*t923 + 0.971*t9 &
            + 85.295*t943 + 16.061*t953
      dbb = oneth*0.011*t9i23 + twoth*13.117*t9i13 + 0.971 &
            + fourth*85.295*t913 + fiveth*16.061*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.31e+04 * t9i32 * exp(-11.733*t9i)
      ddd = -1.5d0*dd*t9i + dd*11.733*t9i2

      ee  = 1.79e+07 * t9i32 * exp(-22.609*t9i)
      dee = -1.5d0*ee*t9i + ee*22.609*t9i2

      ff  = 9.00e+03 * t9113 * exp(-12.517*t9i)
      dff = elvnth*ff*t9i + ff*12.517*t9i2

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01*exp(-13.820*t9i)
      drevdt   = rev*13.820*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh, &
                       q1,q2,q3
      parameter        (q1 = 1.0d0/9.0d0, &
                        q2 = 1.0d0/3.751969d0, &
                        q3 = 1.0d0/64.0d0)


! o15(a,g)ne19

      aa  = 3.57d+11 * t9i23 * exp(-39.584d+0*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*39.584d0*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 - 0.273*t923 - 0.020*t9
      dbb = oneth*0.011*t9i23 - twoth*0.273*t9i13 - 0.020

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 5.10d+10 * t9i23 * exp(-39.584d+0*t9i13 - t92*q2)
      ddd = dd*(-twoth*t9i + oneth*39.584*t9i43 - 2.0d0*t9*q2)

      ee  = 1.0d0 + 0.011*t913 + 1.59*t923 + 0.117*t9 &
            + 1.81*t943 + 0.338*t953
      dee = oneth*0.011*t9i23 + twoth*1.59*t9i13 + 0.117 &
            + fourth*1.81*t913 + fiveth*0.338*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee

      gg  = 3.95d-1 * t9i32 * exp(-5.849*t9i)
      dgg = gg*(-1.5d0*t9i + 5.849*t9i2)

      hh  = 1.90d+1 * t9**2.85 * exp(-7.356*t9i - t92*q3)
      dhh = hh*(2.85*t9i + 7.356*t9i2 - 2.0d0*t9*q3)


      term    = cc + ff + gg + hh
      dtermdt = dcc + dff + dgg + dhh


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54e+10 * t932 * exp(-40.957*t9i)
      drevdt   = rev*(1.5d0*t9i + 40.957*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_f17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee


! f17(p,g)ne18
! wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e+07 * t9i23 * exp(-18.03*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.03*t9i43)

      bb  = 2.194 + 0.050*t913 - 0.376*t923 - 0.061*t9 &
            + 0.026*t943 + 0.011*t953
      dbb = oneth*0.050*t9i23 - twoth*0.376*t9i13 - 0.061 &
            + fourth*0.026*t913 + fiveth*0.011*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 839.0 * t9i32 * exp(-6.93*t9i)
      ddd = dd*(-1.5d0*t9i + 6.93*t9i2)

      ee  = 33.56 * t9i32 * exp(-7.75*t9i)
      dee = ee*(-1.5d0*t9i + 7.75*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.087e+11 * t932 * exp(-45.501*t9i)
      drevdt   = rev*(1.5d0*t9i + 45.501*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 1.672d0, &
                        con      = lntwo/halflife)

! ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff


! f18(p,a)o15
! wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e-10 * t9i32 * exp(-0.302*t9i)
      daa = aa*(-1.5d0*t9i + 0.302*t9i2)

      bb  = 1.56e+05 * t9i32 * exp(-3.84*t9i)
      dbb = bb*(-1.5d0*t9i + 3.84*t9i2)

      cc  = 1.36e+06 * t9i32 * exp(-5.22*t9i)
      dcc = cc*(-1.5d0*t9i + 5.22*t9i2)

      dd  = 8.1e-05 * t9i32 * exp(-1.05*t9i)
      ddd = dd*(-1.5d0*t9i + 1.05*t9i2)

      ee  = 8.9e-04 * t9i32 * exp(-1.51*t9i)
      dee = ee*(-1.5d0*t9i + 1.51*t9i2)

      ff  = 3.0e+05 * t9i32 * exp(-4.29*t9i)
      dff = ff*(-1.5d0*t9i + 4.29*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01 * exp(-33.433*t9i)
      drevdt   = rev*33.433*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne18ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 10.0d0, &
                        a1   = 18.0d0, &
                        ztot = 2.0d0 * z1, &
                        ared = 4.0d0*a1/(4.0d0 + a1), &
                        r    = 5.1566081196876965d0, &
                        c1   = 4.9080044545315392d10, &
                        c2   = 4.9592784569936502d-2, &
                        c3   = 1.9288564401521285d1, &
                        c4   = 4.6477847042196437d1)

! note:
!      r    = 1.09 * a1**oneth + 2.3
!      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
!      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
!      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
!      c4   = 4.2487 * (ztot**2*ared)**oneth


! ne18ap(a,p)na21
! was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 0.0d0
      drevdt = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.304164d0)


! ne19(p,g)na20

      aa  = 1.71d+6 * t9i23 * exp(-19.431d0*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

      bb  = 1.0d0 + 0.021*t913 + 0.130*t923 + 1.95d-2*t9 &
            + 3.86d-2*t943 + 1.47d-02*t953
      dbb = oneth*0.021*t9i23 + twoth*0.130*t9i13 + 1.95d-2 &
            + fourth*3.86d-2*t913 + fiveth*1.47d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb


      dd  = 1.89d+5 * t9i23 * exp(-19.431d0*t9i13 - t92*q1)
      ddd = dd*(-twoth*t9i + oneth*19.431*t9i43 - 2.0d0*t9*q1)

      ee  = 1.0d0 + 0.021*t913 + 2.13*t923 + 0.320*t9 &
            + 2.80*t943 + 1.07*t953
      dee = oneth*0.021*t9i23 + twoth*2.13*t9i13 + 0.320 &
            + fourth*2.80*t913 + fiveth*1.07*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee

      gg  = 8.45d+3 * t9i54 * exp(-7.64d0*t9i)
      dgg = gg*(-fivfour*t9i + 7.64d0*t9i2)


      term    = cc + ff + gg
      dtermdt = dcc + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.39e+09 * t932 * exp(-25.519*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.519*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne19em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 17.3982d0, &
                        con      = lntwo/halflife)

! ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_si26ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 14.0d0, &
                        a1   = 26.0d0, &
                        ztot = 2.0d0 * z1, &
                        ared = 4.0d0*a1/(4.0d0 + a1), &
                        r    = 5.5291207145640335d0, &
                        c1   = 7.3266779970543091d10, &
                        c2   = 4.7895369289991982d-02, &
                        c3   = 2.4322657793918662d1, &
                        c4   = 5.9292366232997814d1)

! note:
!      r    = 1.09 * a1**oneth + 2.3
!      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
!      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
!      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
!      c4   = 4.2487 * (ztot**2*ared)**oneth



! si26ap(a,p)p29
! was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                       zz,q1
      parameter        (q1 = 1.0d0/12.222016d0)


! c12(a,g)o16
      aa   = 1.0d0 + 0.0489d0*t9i23
      daa  = -twoth*0.0489d0*t9i53

      bb   = t92*aa*aa
      dbb  = 2.0d0*(bb*t9i + t92*aa*daa)

      cc   = exp(-32.120d0*t9i13 - t92*q1)
      dcc  = cc * (oneth*32.120d0*t9i43 - 2.0d0*t9*q1)

      dd   = 1.0d0 + 0.2654d0*t9i23
      ddd  = -twoth*0.2654d0*t9i53

      ee   = t92*dd*dd
      dee  = 2.0d0*(ee*t9i + t92*dd*ddd)

      ff   = exp(-32.120d0*t9i13)
      dff  = ff * oneth*32.120d0*t9i43

      gg   = 1.25d3 * t9i32 * exp(-27.499*t9i)
      dgg  = gg*(-1.5d0*t9i + 27.499*t9i2)

      hh   = 1.43d-2 * t95 * exp(-15.541*t9i)
      dhh  = hh*(5.0d0*t9i + 15.541*t9i2)

      zz   = 1.0d0/bb
      f1   = cc*zz
      df1  = (dcc - f1*dbb)*zz

      zz   = 1.0d0/ee
      f2   = ff*zz
      df2  = (dff - f2*dee)*zz

      term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
      dtermdt = 1.04d8*df1 + 1.76d8*df2 + dgg + dhh


! 1.7 times cf88 value
      term     = 1.7d0 * term
      dtermdt  = 1.7d0 * dtermdt

      fr    = term * den
      dfrdt = dtermdt * den * 1.0d-9
      dfrdd = term

      rev    = 5.13d10 * t932 * exp(-83.111*t9i)
      drevdt = rev*(1.5d0*t9i + 83.111*t9i2)

      rr     = rev * term
      drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd  = 0.0d0

      return
      end






      subroutine rate_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                       drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1,rc28, &
                       q1,q2
      parameter        (rc28   = 0.1d0, &
                        q1     = 1.0d0/0.009604d0, &
                        q2     = 1.0d0/0.055225d0)



! triple alfa to c12
! this is a(a,g)be8
      aa    = 7.40d+05 * t9i32 * exp(-1.0663*t9i)
      daa   = aa*(-1.5d0*t9i  + 1.0663*t9i2)

      bb    = 4.164d+09 * t9i23 * exp(-13.49*t9i13 - t92*q1)
      dbb   = bb*(-twoth*t9i + oneth*13.49*t9i43 - 2.0d0*t9*q1)

      cc    = 1.0d0 + 0.031*t913 + 8.009*t923 + 1.732*t9 &
              + 49.883*t943 + 27.426*t953
      dcc   = oneth*0.031*t9i23 + twoth*8.009*t9i13 + 1.732 &
              + fourth*49.883*t913 + fiveth*27.426*t923

      r2abe    = aa + bb * cc
      dr2abedt = daa + dbb*cc + bb*dcc


! this is be8(a,g)c12
      dd    = 130.0d0 * t9i32 * exp(-3.3364*t9i)
      ddd   = dd*(-1.5d0*t9i + 3.3364*t9i2)

      ee    = 2.510d+07 * t9i23 * exp(-23.57*t9i13 - t92*q2)
      dee   = ee*(-twoth*t9i + oneth*23.57*t9i43 - 2.0d0*t9*q2)

      ff    = 1.0d0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + &
              19.176*t943 + 6.034*t953
      dff   = oneth*0.018*t9i23 + twoth*5.249*t9i13 + 0.650 &
              + fourth*19.176*t913 + fiveth*6.034*t923

      rbeac    = dd + ee * ff
      drbeacdt = ddd + dee * ff + ee * dff


! a factor
      xx    = rc28 * 1.35d-07 * t9i32 * exp(-24.811*t9i)
      dxx   = xx*(-1.5d0*t9i + 24.811*t9i2)


! high temperature rate
      if (t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac &
                 + 2.90d-16 * r2abe * drbeacdt &
                 + dxx


! low temperature rate
      else
       uu   = 0.8d0*exp(-(0.025*t9i)**3.263)
       yy   = 0.2d0 + uu
!       yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*t9i)**2.263 * (0.025*t9i2)
       vv   = 4.0d0*exp(-(t9/0.025)**9.227)
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = 0.01d0 + yy * aa
!       f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa
       term = 2.90d-16 * r2abe * rbeac * f1 +  xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac * f1 &
                 + 2.90d-16 * r2abe * drbeacdt * f1 &
                 + 2.90d-16 * r2abe * rbeac * df1 &
                 + dxx
      end if


! rates
!      term    = 1.2d0 * term
!      dtermdt = 1.2d0 * term

      fr    = term * den * den
      dfrdt = dtermdt * den * den * 1.0d-9
      dfrdd = 2.0d0 * term * den

      rev    = 2.00d+20*t93*exp(-84.424*t9i)
      drevdt = rev*(3.0d0*t9i + 84.424*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       aa,zz


! c12 + c12 reaction
      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9*zz
      dt9a    = (1.0d0 -  t9a*0.0396)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 4.27d+26 * t9a56 * t9i32 * &
                exp(-84.165/t9a13 - 2.12d-03*t93)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12c12npa(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

! locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,dd,ddd,zz, &
                       b24n,db24n,b24p,db24p,b24a,db24a


! c12(c12,n)mg23
! c12(c12,p)na23
! c12(c12,a)ne20


      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9 * zz
      dt9a    = (1.0d0 -  t9a*0.0396) * zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      aa = 4.27d+26 * t9a56 * t9i32 * &
           exp(-84.165/t9a13 - 2.12d-03*t93)

      daa = aa * (dt9a56/t9a56 - 1.5d0*t9i &
                + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)


! neutron branching from dayras switkowski and woosley 1976
      if (t9 .ge. 1.5) then

       bb    =  0.055 * exp(0.976 - 0.789*t9)
       dbb   = -bb*0.789

       b24n  = 0.055  - bb
       db24n = -dbb

      else

       bb    = 1.0d0 + 0.0789*t9 + 7.74*t92
       dbb   = 0.0789 + 2.0d0*7.74*t9

       cc    = 0.766*t9i3
       dcc   = -3.0d0*cc*t9i

       dd    = bb * cc
       ddd   = dbb*cc + bb*dcc

       b24n  = 0.859*exp(-dd)
       db24n = -b24n*ddd
      end if


! proton branching ratio
      if (t9.gt.3.) then
        b24p  = oneth*(1.0d0 - b24n)
        db24p = -oneth*db24n

        b24a  = 2.0d0 * b24p
        db24a = 2.0d0 * db24p

       else
        b24p  = 0.5d0*(1.0d0 - b24n)
        db24p = -0.5d0*db24n

        b24a  = b24p
        db24a = db24p

       end if

! rates

! c12(c12,n)mg23
      term    = aa * b24n
      dtermdt = daa*b24n + aa*db24n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(30.16100515d0*t9i)
       drevdt = -rev*30.16100515d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term


! c12(c12,p)na23
      term    = aa * b24p
      dtermdt = daa*b24p + aa*db24p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(-25.98325915d0*t9i)
       drevdt = rev*25.98325915d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! c12(c12,a)ne20
      term    = aa * b24a
      dtermdt = daa*b24a + aa*db24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.42 * exp(-53.576110995d0*t9i)
       drevdt = rev*53.576110995d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term

      return
      end






      subroutine rate_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz


! c12 + o16 reaction; see cf88 references 47-4

      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18*t9a*t9a)
       daa     = -aa * 0.36 * t9a * dt9a

       bb      = 1.06d-03*exp(2.562*t9a23)
       dbb     = bb * 2.562 * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0d0/cc
       term    = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                       + 106.594/t9a23*dt9a13 - zz*dcc)

      else
!       term    = 2.6288035d-29
       term    = 0.0d0
       dtermdt = 0.0d0
      endif


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c12o16npa(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

! locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13, &
                       t9a23,dt9a23,t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,b27n,b27p,b24a,zz


! c12(o16,n)si27
! c12(o16,p)al27
! c12(o16,a)mg24


      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa     = exp(-0.18*t9a*t9a)
       daa    = -aa * 0.36 * t9a * dt9a

       bb     = 1.06d-03*exp(2.562*t9a23)
       dbb    = bb * 2.562 * dt9a23

       cc     = aa + bb
       dcc    = daa + dbb

       zz     = 1.0d0/cc
       dd     = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) *zz
       ddd    = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 106.594/t9a23 * dt9a13 - dcc*zz)

      else
!       dd     = 2.6288035d-29
       dd     = 0.0d0
       ddd    = 0.0d0
      endif


! branching ratios from pwnsz data
        b27n = 0.1d0
        b27p = 0.5d0
        b24a = 0.4d0


! rates

! c12(o16,n)si27
      term    = dd * b27n
      dtermdt = ddd * b27n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(4.8972467d0*t9i)
       drevdt = -rev*4.8972467d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den*(drevdt*term + rev*dtermdt)*1.0d-9
      drr1dd = rev*term


! c12(o16,p)al27
      term    = dd * b27p
      dtermdt = ddd * b27p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(-59.9970745d0*t9i)
       drevdt = rev*59.9970745d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! c12(o16,a)mg24
      term    = dd * b24a
      dtermdt = ddd * b24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.83d0 * exp(-78.5648345d0*t9i)
       drevdt = rev*78.5648345d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


      return
      end








      subroutine rate_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt


! o16 + o16
      term  = 7.10d36 * t9i23 * &
              exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      dtermdt = -twoth*term*t9i &
                + term * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16o16npad(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd, &
                      fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd, &
                      fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa, &
                       b32n,db32n,b32p,db32p,b32a,db32a,b32d,db32d, &
                       ezro,dezro,dlt,ddlt,xxt,dxxt,thrs,dthrs


! o16(o16,n)s31
! o16(o16,p)p31
! o16(o16,a)si28
! o16(o16,d)p30

      aa  = 7.10d36 * t9i23 * exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      daa = -twoth*aa*t9i + aa * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! branching ratios highly uncertain;  guessed using fcz 1975
! deuteron channel is endoergic. apply error function cut-off.
       ezro = 3.9*t923
       dezro = twoth*ezro*t9i

       dlt  = 1.34*t9**fivsix
       ddlt = fivsix*dlt*t9i

       xxt  = 2.0d0*(2.406 - ezro)/dlt
       dxxt = -(2.0d0*dezro + xxt*ddlt)/dlt

       call fowthrsh(xxt,thrs,dthrs)

       b32d  = 0.05d0*thrs
       db32d = 0.05d0*dthrs*dxxt

! shut down deutrium channel
!       b32d  = 0.0d0
!       db32d = 0.0d0

       b32n  = 0.1d0
       db32n = 0.0d0

       b32a  = 0.25d0
       db32a = 0.0d0

       b32p  = 1.0d0 - b32d - b32a - b32n
       db32p = -db32d



! rates

! o16(o16,n)s31
      term    = aa * b32n
      dtermdt = daa*b32n + aa*db32n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92 * exp(-16.8038228d0*t9i)
       drevdt = rev*16.8038228d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den*(drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term

! o16(o16,p)p31
      term    = aa * b32p
      dtermdt = daa*b32p + aa*db32p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92*exp(-89.0788286d0*t9i)
       drevdt = rev*89.0788286d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! o16(o16,a)si28
      term    = aa * b32a
      dtermdt = daa*b32a + aa*db32a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.46*exp(-111.3137212d0*t9i)
       drevdt = rev*111.3137212d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


! o16(o16,d)p30
      term    = aa * b32d
      dtermdt = daa*b32d + aa*db32d
      fr4     = den * term
      dfr4dt  = den * dtermdt * 1.0d-9
      dfr4dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 0.984*exp(27.9908982d0*t9i)
       drevdt = -rev*27.9908982d0*t9i2
      end if
      rr4    = den * rev * term
      drr4dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr4dd = rev*term

      return
      end






      subroutine fowthrsh(x,thrs,dthrs)
      include 'implno.dek'

! fowler threshold fudge function.
! err func rational (abramowitz p.299)7.1.25 and its derivative

! declare
      double precision x,thrs,dthrs, &
                       ag,z,z2,t,t2,t3,tt,er,aa,daa,dt,dtt,der

      ag   = sign(1.0d0,x)
      z    = abs(x)
      z2   = z*z
      aa   = 1.0d0 + 0.47047d0*z

      t    = 1.0d0/aa
      dt   = -t*t*0.47047*ag

      t2   = t*t
      t3   = t2*t

      tt   = 0.3480242d0*t - 0.0958798d0*t2 + 0.7478556d0*t3
      dtt  = dt * (0.3480242d0 - 2.0d0*0.0958798d0*t &
                   + 3.0d0*0.7478556d0*t2)

      thrs  = 0.5d0
      dthrs = -5.6452433937900004d-1

      if (z .ne. 0) then
       aa   = exp(-z2)
       daa  = -2.0d0*aa*z*ag

       er   = 1.0d0 - tt * aa
       der  = -dtt*aa - tt*daa

       thrs  = 0.5d0 * (1.0d0 - ag*er)
       dthrs = -0.5d0*ag*der
      end if
      return
      end





      subroutine rate_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       cc,dcc,term2,dterm2,rev,drevdt,q1
      parameter        (q1 = 1.0d0/2.515396d0)



! o16(a,g)ne20
      term1   = 9.37d9 * t9i23 * exp(-39.757*t9i13 - t92*q1)
      dterm1  = term1*(-twoth*t9i + oneth*39.757*t9i43 - 2.0d0*t9*q1)

      aa      = 62.1 * t9i32 * exp(-10.297*t9i)
      daa     = aa*(-1.5d0*t9i + 10.297*t9i2)

      bb      = 538.0d0 * t9i32 * exp(-12.226*t9i)
      dbb     = bb*(-1.5d0*t9i + 12.226*t9i2)

      cc      = 13.0d0 * t92 * exp(-20.093*t9i)
      dcc     = cc*(2.0d0*t9i + 20.093*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      term    = term1 + term2
      dtermdt = dterm1 + dterm2


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.65d+10*t932*exp(-54.937*t9i)
      drevdt   = rev*(1.5d0*t9i + 54.937*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       term2,dterm2,term3,dterm3,rev,drevdt,zz,rc102,q1
      parameter        (rc102 = 0.1d0, &
                        q1    = 1.0d0/4.923961d0)



! ne20(a,g)mg24

      aa   = 4.11d+11 * t9i23 * exp(-46.766*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*46.766*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.009*t913 + 0.882*t923 + 0.055*t9 &
             + 0.749*t943 + 0.119*t953
      dbb  = oneth*0.009*t9i23 + twoth*0.882*t9i13 + 0.055 &
             + fourth*0.749*t913 + fiveth*0.119*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb


      aa   = 5.27d+03 * t9i32 * exp(-15.869*t9i)
      daa  = aa*(-1.5d0*t9i + 15.869*t9i2)

      bb   = 6.51d+03 * t912 * exp(-16.223*t9i)
      dbb  = bb*(0.5d0*t9i + 16.223*t9i2)

      term2  = aa + bb
      dterm2 = daa + dbb


      aa   = 42.1 * t9i32 * exp(-9.115*t9i)
      daa  = aa*(-1.5d0*t9i + 9.115*t9i2)

      bb   =  32.0 * t9i23 * exp(-9.383*t9i)
      dbb  = bb*(-twoth*t9i + 9.383*t9i2)

      term3  = rc102 * (aa + bb)
      dterm3 = rc102 * (daa + dbb)


      aa  = 5.0d0*exp(-18.960*t9i)
      daa = aa*18.960*t9i2

      bb  = 1.0d0 + aa
      dbb = daa

      zz      = 1.0d0/bb
      term    = (term1 + term2 + term3)*zz
      dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.01d+10 * t932 * exp(-108.059*t9i)
      drevdt   = rev*(1.5d0*t9i + 108.059*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,hhi,rev,drevdt,rc121
      parameter        (rc121 = 0.1d0)



! 24mg(a,g)28si

      aa    = 4.78d+01 * t9i32 * exp(-13.506*t9i)
      daa   = aa*(-1.5d0*t9i + 13.506*t9i2)

      bb    =  2.38d+03 * t9i32 * exp(-15.218*t9i)
      dbb   = bb*(-1.5d0*t9i + 15.218*t9i2)

      cc    = 2.47d+02 * t932 * exp(-15.147*t9i)
      dcc   = cc*(1.5d0*t9i + 15.147*t9i2)

      dd    = rc121 * 1.72d-09 * t9i32 * exp(-5.028*t9i)
      ddd   = dd*(-1.5d0*t9i + 5.028*t9i2)

      ee    = rc121* 1.25d-03 * t9i32 * exp(-7.929*t9i)
      dee   = ee*(-1.5d0*t9i + 7.929*t9i2)

      ff    = rc121 * 2.43d+01 * t9i * exp(-11.523*t9i)
      dff   = ff*(-t9i + 11.523*t9i2)

      gg    = 5.0d0*exp(-15.882*t9i)
      dgg   = gg*15.882*t9i2

      hh    = 1.0d0 + gg
      hhi   = 1.0d0/hh

      term    = (aa + bb + cc + dd + ee + ff) * hhi
      dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.27d+10 * t932 * exp(-115.862*t9i)
      drevdt   = rev*(1.5d0*t9i + 115.862*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                       rev,drevdt,rc148,q1
      parameter        (rc148 = 0.1d0, &
                        q1    = 1.0d0/0.024649d0)



! 24mg(a,p)al27
      aa     = 1.10d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa    = -twoth*aa*t9i + aa*(23.261*t9i43 - 2.0d0*t9*q1)

      bb     =  1.0d0 + 0.018*t913 + 12.85*t923 + 1.61*t9 &
               + 89.87*t943 + 28.66*t953
      dbb    = oneth*0.018*t9i23 + twoth*12.85*t9i13 + 1.61 &
                + fourth*89.87*t913 + fiveth*28.66*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb

      aa     = 129.0d0 * t9i32 * exp(-2.517*t9i)
      daa    = -1.5d0*aa*t9i + aa*2.517*t9i2

      bb     = 5660.0d0 * t972 * exp(-3.421*t9i)
      dbb    = 3.5d0*bb*t9i +  bb*3.421*t9i2

      cc     = rc148 * 3.89d-08 * t9i32 * exp(-0.853*t9i)
      dcc    = -1.5d0*cc*t9i + cc*0.853*t9i2

      dd     = rc148 * 8.18d-09 * t9i32 * exp(-1.001*t9i)
      ddd    = -1.5d0*dd*t9i + dd*1.001*t9i2

      term2  = aa + bb + cc + dd
      dterm2 = daa + dbb + dcc + ddd

      ee     = oneth*exp(-9.792*t9i)
      dee    = ee*9.792*t9i2

      ff     =  twoth * exp(-11.773*t9i)
      dff    = ff*11.773*t9i2

      gg     = 1.0d0 + ee + ff
      dgg    = dee + dff

      term    = (term1 + term2)/gg
      dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


! the rates
      rev      = 1.81 * exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! al27(p,g)si28
! champagne 1996

      aa  = 1.32d+09 * t9i23 * exp(-23.26*t9i13)
      daa = aa*(-twoth*t9i + oneth*23.26*t9i43)

      bb  = 3.22d-10 * t9i32 * exp(-0.836*t9i)*0.17
      dbb = bb*(-1.5d0*t9i + 0.836*t9i2)

      cc  = 1.74d+00 * t9i32 * exp(-2.269*t9i)
      dcc = cc*(-1.5d0*t9i + 2.269*t9i2)

      dd  = 9.92d+00 * t9i32 * exp(-2.492*t9i)
      ddd = dd*(-1.5d0*t9i + 2.492*t9i2)

      ee  = 4.29d+01 * t9i32 * exp(-3.273*t9i)
      dee = ee*(-1.5d0*t9i + 3.273*t9i2)

      ff  = 1.34d+02 * t9i32 * exp(-3.654*t9i)
      dff = ff*(-1.5d0*t9i + 3.654*t9i2)

      gg  = 1.77d+04 * (t9**0.53) * exp(-4.588*t9i)
      dgg = gg*(0.53*t9i + 4.588*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11 * t932 * exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp, &
                       rev,drevdt,rc147,q1
      parameter        (rc147 = 0.1d0, &
                        q1    = 1.0d0/0.024025d0)



! 27al(p,g)si28  cf88

      aa  = 1.67d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*23.261*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.018*t913 + 5.81*t923 + 0.728*t9 &
            + 27.31*t943 + 8.71*t953
      dbb = oneth*0.018*t9i23 + twoth*5.81*t9i13 + 0.728 &
            + fourth*27.31*t913 + fiveth*8.71*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 2.20d+00 * t9i32 * exp(-2.269*t9i)
      ddd = dd*(-1.5d0*t9i + 2.269*t9i2)

      ee  = 1.22d+01 * t9i32 * exp(-2.491*t9i)
      dee = ee*(-1.5d0*t9i + 2.491*t9i2)

      ff  =  1.50d+04 * t9 * exp(-4.112*t9i)
      dff = ff*(t9i + 4.112*t9i2)

      gg  = rc147 * 6.50d-10 * t9i32 * exp(-0.853*t9i)
      dgg = gg*(-1.5d0*t9i + 0.853*t9i2)

      hh  = rc147 * 1.63d-10 * t9i32 * exp(-1.001*t9i)
      dhh = hh*(-1.5d0*t9i + 1.001*t9i2)

      xx     = oneth*exp(-9.792*t9i)
      dxx    = xx*9.792*t9i2

      yy     =  twoth * exp(-11.773*t9i)
      dyy    = yy*11.773*t9i2

      zz     = 1.0d0 + xx + yy
      dzz    = dxx + dyy

      pp      = 1.0d0/zz
      term    = (cc + dd + ee + ff + gg + hh)*pp
      dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11*t932*exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
      if (z .eq. 10.0) then
       daa = 0
      else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*t9 - 3.0d0*2.900d-4*t92
      end if

      term    = 4.82d+22 * t9i23 * exp(-61.015 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 61.015*t9i13*(oneth*t9i*aa - daa))

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.461d+10 * t932 * exp(-80.643*t9i)
      drevdt   = rev*(1.5d0*t9i + 80.643*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,p)p31
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*t9 - 3.0d0*2.341d-4*t92
      end if

      term    = 4.16d+13 * t9i23 * exp(-25.631 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*25.631*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 0.5825d0 * exp(-22.224*t9i)
      drevdt   = rev*22.224*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! p31(p,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*t9 + 3.0d0*6.444d-4*t92
      end if

      term    = 1.08d+16 * t9i23 * exp(-27.042 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 27.042*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.764d+10 * t932 * exp(-102.865*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.865*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,g)ar36
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*t9 - 3.0d0*4.067d-4*t92
      end if

      term    = 1.16d+24 * t9i23 * exp(-66.690 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 66.690*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.616d+10 * t932 * exp(-77.080*t9i)
      drevdt   = rev*(1.5d0*t9i + 77.080*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,p)cl35
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*t9 + 3.0d0*6.969d-4*t92
      end if

      term    = 1.27d+16 * t9i23 * exp(-31.044 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*31.044*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.144 * exp(-21.643*t9i)
      drevdt   = rev*21.643*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt


! cl35(p,g)ar36
      aa    = 1.0d0 + 1.761d-1*t9 - 1.322d-2*t92 + 5.245d-4*t93
      daa   = 1.761d-1 - 2.0d0*1.322d-2*t9 + 3.0d0*5.245d-4*t92


      term    =  4.48d+16 * t9i23 * exp(-29.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 29.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.568d+10*t932*exp(-98.722*t9i)
      drevdt   = rev*(1.5d0*t9i + 98.722*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*t9 + 3.0d0*3.790d-4*t92
      end if

      term    = 2.81d+30 * t9i23 * exp(-78.271 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 78.271*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.740d+10 * t932 * exp(-81.711*t9i)
      drevdt   = rev*(1.5d0*t9i + 81.711*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,p)k39
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*t9 + 3.0d0*4.021d-4*t92
      end if

      term    = 2.76d+13 * t9i23 * exp(-34.922 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*34.922*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.128*exp(-14.959*t9i)
      drevdt   = rev*14.959*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! k39(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.622d-1 - 2.0d0*1.119d-2*t9 + 3.0d0*3.910d-4*t92
      end if

      term    = 4.09d+16 * t9i23 * exp(-31.727 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 31.727*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.600d+10 * t932 * exp(-96.657*t9i)
      drevdt   = rev*(1.5d0*t9i + 96.657*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,g)ti44
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*z - 3.0d0*3.889d-4*z2
      end if

      term    = 4.66d+24 * t9i23 * exp(-76.435 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 76.435*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.843d+10 * t932 * exp(-59.510*t9i)
      drevdt   = rev*(1.5d0*t9i + 59.510*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,p)sc43
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*t9 - 3.0d0*5.071d-4*t92
      end if

      term    = 4.54d+14 * t9i23 * exp(-32.177 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*32.177*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 2.229 * exp(-40.966*t9i)
      drevdt   = rev*40.966*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! sc43(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*t9 - 3.0d0*5.463d-5*t92
      end if

      term    = 3.85d+16 * t9i23 * exp(-33.234 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 33.234*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.525d+11 * t932 * exp(-100.475*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.475*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*t9 + 3.0d0*5.324d-4*t92
      end if

      term    = 1.37d+26 * t9i23 * exp(-81.227 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.227*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.928d+10*t932*exp(-89.289*t9i)
      drevdt   = rev*(1.5d0*t9i + 89.289*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ti44ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,p)v47
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.655d-2 - 2.0d0*3.947d-3*t9 + 3.0d0*2.522d-4*t92
      end if

      term    = 6.54d+20 * t9i23 * exp(-66.678 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*66.678*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.104 * exp(-4.723*t9i)
      drevdt   = rev*4.723*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_v47pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! v47(p,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*t9 - 3.0d0*6.662d-5*t92
      end if

      term    = 2.05d+17 * t9i23 * exp(-35.568 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 35.568*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.649d+10*t932*exp(-93.999*t9i)
      drevdt   = rev*(1.5d0*t9i + 93.999*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_cr48ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*t9 + 3.0d0*2.848d-4*t92
      end if

      term    = 1.04d+23 * t9i23 * exp(-81.420 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.420*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.001d+10 * t932 * exp(-92.177*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.177*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_cr48ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,p)mn51
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.384d-2 + 2.0d0*1.081d-3*t9 - 3.0d0*5.933d-5*t92
      end if

      term    = 1.83d+26 * t9i23 * exp(-86.741 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*86.741*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.6087*exp(-6.510*t9i)
      drevdt   = rev*6.510*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mn51pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! mn51(p,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*t9 - 3.0d0*9.453d-5*t92
      end if

      term    = 3.77d+17 * t9i23 * exp(-37.516 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 37.516*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.150d+11*t932*exp(-85.667*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.667*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_fe52ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*t9 + 3.0d0*3.723d-4*t92
      end if

      term    = 1.05d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 91.674*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.064d+10*t932*exp(-92.850*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.850*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,p)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.367d-2 + 2.0d0*7.428d-4*t9 - 3.0d0*3.050d-5*t92
      end if

      term    = 1.30d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*91.674*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.4597*exp(-9.470*t9i)
      drevdt   = rev*9.470*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_co55pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! co55(p,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*t9 - 3.0d0*2.160d-5*t92
      end if

      term    = 1.21d+18 * t9i23 * exp(-39.604 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 39.604*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.537d+11*t932*exp(-83.382*t9i)
      drevdt   = rev*(1.5d0*t9i + 83.382*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq2


! fe52(n,g)fe53
      tq2     = t9 - 0.348d0
      term    = 9.604d+05 * exp(-0.0626*tq2)
      dtermdt = -term*0.0626

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.43d+09 * t932 * exp(-123.951*t9i)
      drevdt   = rev*(1.5d0*t9i + 123.951*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_fe53ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2


! fe53(n,g)fe54
      tq1   = t9/0.348
      tq10  = tq1**0.10
      dtq10 = 0.1d0*tq10/(0.348*tq1)
      tq2   = t9 - 0.348d0

      term    = 1.817d+06 * tq10 * exp(-0.06319*tq2)
      dtermdt = term/tq10*dtq10 - term*0.06319

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.56d+11 * t932 * exp(-155.284*t9i)
      drevdt   = rev*(1.5d0*t9i + 155.284*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe54(n,g)fe55
      aa   =  2.307390d+01 - 7.931795d-02 * t9i + 7.535681d+00 * t9i13 &
            - 1.595025d+01 * t913 + 1.377715d+00 * t9 - 1.291479d-01 * t953 &
            + 6.707473d+00 * log(t9)

      daa  =  7.931795d-02 * t9i2 - oneth * 7.535681d+00 * t9i43 &
            - oneth * 1.595025d+01 *t9i23 + 1.377715d+00 - fiveth * 1.291479d-01 *t923 &
            + 6.707473d+00 * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 4.800293d+09 * t932 * exp(-1.078986d+02 * t9i)
      dbb = bb*(1.5d0*t9i + 1.078986d+02 * t9i2)

! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end



      subroutine rate_fe54pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,z,z2,z3


! fe54(p,g)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.593d-2 - 2.0d0*3.445d-3*t9 + 3.0d0*8.594d-5*t92
      end if

      term    = 4.51d+17 * t9i23 * exp(-38.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 38.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.400d+09 * t932 * exp(-58.605*t9i)
      drevdt   = rev*(1.5d0*t9i + 58.605*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe54(a,p)co57
      aa   =  3.97474900d+01 - 6.06543100d+00 * t9i + 1.63239600d+02 * t9i13 &
            - 2.20457700d+02 * t913 + 8.63980400d+00 * t9 - 3.45841300d-01 * t953 &
            + 1.31464200d+02 * log(t9)

      daa  =  6.06543100d+00 * t9i2 - oneth * 1.63239600d+02 * t9i43 &
            - oneth * 2.20457700d+02 * t9i23 + 8.63980400d+00 - fiveth * 3.45841300d-01 * t923 &
            + 1.31464200d+02  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 2.16896000d+00  * exp(-2.05631700d+01 * t9i)
      dbb = bb * 2.05631700d+01 * t9i2

! reverse rate
      drrdd = term
      rr    = term*den
      drrdt = dtermdt*den

! forward rate
      fr    = rr*bb
      dfrdt = drrdt*bb + rr*dbb*1.0d-9
      dfrdd = drrdd*bb

      return
      end







      subroutine rate_fe55ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe55(n,g)fe56
      aa   =  1.954115d+01 - 6.834029d-02 * t9i + 5.379859d+00 * t9i13 &
            - 8.758150d+00 * t913 + 5.285107d-01 * t9 - 4.973739d-02  * t953 &
            + 4.065564d+00  * log(t9)

      daa  =  6.834029d-02 * t9i2 - oneth * 5.379859d+00 * t9i43 &
            - oneth * 8.758150d+00 * t9i23 + 5.285107d-01 - fiveth * 4.973739d-02  *t923 &
            + 4.065564d+00  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 7.684279d+10  * t932 * exp(-1.299472d+02  * t9i)
      dbb = bb*(1.5d0*t9i + 1.299472d+02 * t9i2)

! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end




      subroutine rate_fe56pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe56(p,g)co57
   
      aa   =  1.755960d+02 - 7.018872d+00 * t9i + 2.800131d+02 * t9i13 &
            - 4.749343d+02 * t913 + 2.683860d+01 * t9 - 1.542324d+00  * t953 &
            + 2.315911d+02  * log(t9)

      daa  =  7.018872d+00 * t9i2 - oneth * 2.800131d+02 * t9i43 &
            - oneth * 4.749343d+02 * t9i23 + 2.683860d+01 - fiveth * 1.542324d+00  *t923 &
            + 2.315911d+02  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 2.402486d+09 * t932 * exp(-6.995192d+01 * t9i)
      dbb = bb*(1.5d0*t9i + 6.995192d+01 * t9i2)


! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end
!---------------------------------------------------------------------

      subroutine azbar_simple(ymass,aion,zion,ionbeg,ionend, &
                              abar,zbar,ye)
      include 'implno.dek'

! this routine calculates common composition variables

! input:
! molar abundances             = ymass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! start of isotopes            = ionbeg
! end of isotopes              = ionend
!
! output:
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! electron fraction       = ye                dimensionless


! declare the pass
      integer :: ionbeg,ionend
      real*8  :: ymass(ionbeg:ionend),aion(ionbeg:ionend),zion(ionbeg:ionend),abar,zbar,ye

! local variables
      real*8  :: sum1

! mean number of nucleons, numerator assumes sum xmass = 1.0d0
      abar  = 1.0d0/sum(ymass(ionbeg:ionend))

! mean electron fraction and mean charge
      ye    = sum(zion(ionbeg:ionend)*ymass(ionbeg:ionend))
      zbar  = abar * ye

      return
      end subroutine azbar_simple

      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! this routine computes the fermi-dirac integrals f(dk,eta,theta) of
! index dk, with degeneracy parameter eta and relativity parameter theta.

! input is dk the index of the fermi-dirac function,
! eta the degeneracy parameter, and theta the relativity parameter.

! output is fd is computed by applying three 10-point gauss-legendre
! and one 10-point gauss-laguerre rules over four appropriate subintervals.
! fdeta is the derivative of f with respect to eta.
! fdthetha is the derivative of f with theta.
! fdeta2 is the second derivative of f with respect to eta.
! fdtheta2 is the second derivative of f with respect to theta
! fdetadtheta is the mixed second derivative of f with respect to eta and theta.

! reference: j.m. aparicio, apjs 117, 632 1998


! declare
      external :: fdfunc1,fdfunc2
      real*8   :: dk,eta,theta,fd,fdeta,fdtheta, &
                  fdeta2,fdtheta2,fdetadtheta, &
                  d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3, &
                  eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3), &
                  res(4),drde(4),drdt(4),drde2(4),drdt2(4),drdet(4)


!   parameters defining the location of the breakpoints for the
!   subintervals of integration:
      data d   / 3.3609d0 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d0 /
      data b1  / 1.1418d0 /
      data c1  / 2.9826d0 /
      data a2  / 3.7601d0 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d1 /
      data e2  / 1.0056d0 /
      data a3  / 7.5669d0 /
      data b3  / 1.1695d0 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d0 /
      data e3  /-1.2819d-1 /


!   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


!   definition of xi:
      eta1=sg*(eta-d)
      if (eta1 .le. 5.0d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

!   definition of the x_i:
      x1=(a1 + b1*xi + c1*xi2)    / (1.d0 +c1*xi)
      x2=(a2 + b2*xi + c2*d2*xi2) / (1.d0 +e2*xi + c2*xi2)
      x3=(a3 + b3*xi + c3*d3*xi2) / (1.d0 +e3*xi + c3*xi2)

!   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=sqrt(s1)

!   quadrature integrations:

! 9 significant figure accuracy
!      call dqleg010(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
!                    drde2(1), drdt2(1), drdet(1), par,3)
!      call dqleg010(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
!                    drde2(2), drdt2(2), drdet(2), par,3)
!      call dqleg010(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
!                    drde2(3), drdt2(3), drdet(3), par,3)
!      call dqlag010(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
!                    drde2(4), drdt2(4), drdet(4), par,3)

! 14 significant figure accuracy
      call dqleg020(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
                    drde2(1), drdt2(1), drdet(1), par,3)
      call dqleg020(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
                    drde2(2), drdt2(2), drdet(2), par,3)
      call dqleg020(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
                    drde2(3), drdt2(3), drdet(3), par,3)
      call dqlag020(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
                    drde2(4), drdt2(4), drdet(4), par,3)

! 18 significant figure accuracy
!      call dqleg040(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
!                    drde2(1), drdt2(1), drdet(1), par,3)
!      call dqleg040(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
!                    drde2(2), drdt2(2), drdet(2), par,3)
!      call dqleg040(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
!                    drde2(3), drdt2(3), drdet(3), par,3)
!      call dqlag040(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
!                    drde2(4), drdt2(4), drdet(4), par,3)

! 32 significant figure accuracy
!      call dqleg080(fdfunc2, 0.d0,  s12, res(1), drde(1), drdt(1), &
!                    drde2(1), drdt2(1), drdet(1), par,3)
!      call dqleg080(fdfunc1,   s1,   s2, res(2), drde(2), drdt(2), &
!                    drde2(2), drdt2(2), drdet(1), par,3)
!      call dqleg080(fdfunc1,   s2,   s3, res(3), drde(3), drdt(3), &
!                    drde2(3), drdt2(3), drdet(1), par,3)
!      call dqlag080(fdfunc1,   s3, 1.d0, res(4), drde(4), drdt(4), &
!                    drde2(4), drdt2(4), drdet(1), par,3)


! sum the contributions
      fd          = res(1)   + res(2)   + res(3)   + res(4)
      fdeta       = drde(1)  + drde(2)  + drde(3)  + drde(4)
      fdtheta     = drdt(1)  + drdt(2)  + drdt(3)  + drdt(4)
      fdeta2      = drde2(1) + drde2(2) + drde2(3) + drde2(4)
      fdtheta2    = drdt2(1) + drdt2(2) + drdt2(3) + drdt2(4)
      fdetadtheta = drdet(1) + drdet(2) + drdet(3) + drdet(4)
      return
      end subroutine dfermi

      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta, &
                         fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! forms the fermi-dirac integrand and its first and second
! derivatives with eta and theta.

! input:
! x is the integration variable
! par(1) is the index
! par(2) is the degeneravy parameter
! par(3) is the relativity parameter

! output:
! fd is the integrand
! fdeta is the first derivative with eta
! fdeta2 is the second derivative with eta
! fdtheta is the first derivative with theta
! fdtheta2 is the second derivative with theta
! fdetadtheta is the mixed second derivative

! declare the pass
      integer     :: n
      real*8      :: x,par(n),fd, &
                     fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta

! local variables
      real*8      :: dk,eta,theta, &
                     factor,xst,dxst,denom,denomi,denom2,xdk


! initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xst   = 1.0d0 + 0.5d0*x*theta
      dxst  = sqrt(xst)

!   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor      = exp(x-eta)
       denom       = factor + 1.0d0
       denomi      = 1.0d0/denom
       fd          = xdk * dxst * denomi
       fdeta       = fd * factor * denomi
       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * x * denom2
       fdtheta2    = -fdtheta * x * denom2
       fdetadtheta = fdtheta * factor * denomi
      else
       factor      = exp(eta-x)
       fd          = xdk * dxst * factor
       fdeta       = fd
       fdeta2      = fd
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * x * denom2
       fdtheta2    = -fdtheta * x * denom2
       fdetadtheta = fdtheta
      endif

      return
      end subroutine fdfunc1

      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta, &
                         fdeta2,fdtheta2,fdetadtheta)
      include 'implno.dek'

! forms the fermi-dirac integrand and its first and second
! derivatives with eta and theta when the variable change z**2=x
! has been made.

! input:
! x is the integration variable
! par(1) is the index
! par(2) is the degeneravy parameter
! par(3) is the relativity parameter

! output:
! fd is the integrand
! fdeta is the first derivative with eta
! fdeta2 is the second derivative with eta
! fdtheta is the first derivative with theta
! fdtheta2 is the second derivative with theta
! fdetadtheta is the mixed second derivative

! declare the pass
      integer     ::  n
      real*8      ::  x,par(n),fd, &
                      fdeta,fdeta2,fdtheta,fdtheta2,fdetadtheta

! local variables
      real*8      ::  dk,eta,theta, &
                      factor,xst,dxst,denom,denomi,denom2,xdk,xsq


      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xst   = 1.0d0 + 0.5d0 * xsq * theta
      dxst  = sqrt(xst)

!   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.0d2) then
       factor      = exp(xsq - eta)
       denom       = factor + 1.0d0
       denomi      = 1.0d0/denom
       fd          = 2.0d0 * xdk * dxst * denomi
       fdeta       = fd * factor * denomi
       fdeta2      = (2.0d0 * factor * denomi - 1.0d0)*fdeta
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * xsq * denom2
       fdtheta2    = -fdtheta * xsq * denom2
       fdetadtheta = fdtheta * factor * denomi
      else
       factor      = exp(eta - xsq)
       fd          = 2.0d0 * xdk * dxst * factor
       fdeta       = fd
       fdeta2      = fd
       denom2      = 1.0d0/(4.0d0 * xst)
       fdtheta     = fd * xsq * denom2
       fdtheta2    = -fdtheta * xsq * denom2
       fdetadtheta = fdtheta
      endif

      return
      end subroutine fdfunc2

            subroutine dqleg020(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-legendre rule for the fermi-dirac function and
! its derivatives with respect to eta and theta.
! on input f is the name of the subroutine containing the integrand,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to subroutine f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-legendre rule,
! drdeta is the derivative with respect to eta, and drdtheta is the
! derivative with respect to theta.
!
! note: since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc1,absc2,center,hlfrun,wg(10),xg(10), &
                      fval(2),deta(2),dtheta(2), &
                      deta2(2),dtheta2(2),detadtheta(2)


! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 20-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 20-point rule
! wg     - weights of the 20-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.65265211334973337546404093988382110d-2 /
      data xg (  2) /   2.27785851141645078080496195368574624d-1 /
      data xg (  3) /   3.73706088715419560672548177024927237d-1 /
      data xg (  4) /   5.10867001950827098004364050955250998d-1 /
      data xg (  5) /   6.36053680726515025452836696226285936d-1 /
      data xg (  6) /   7.46331906460150792614305070355641590d-1 /
      data xg (  7) /   8.39116971822218823394529061701520685d-1 /
      data xg (  8) /   9.12234428251325905867752441203298113d-1 /
      data xg (  9) /   9.63971927277913791267666131197277221d-1 /
      data xg ( 10) /   9.93128599185094924786122388471320278d-1 /

      data wg (  1) /   1.52753387130725850698084331955097593d-1 /
      data wg (  2) /   1.49172986472603746787828737001969436d-1 /
      data wg (  3) /   1.42096109318382051329298325067164933d-1 /
      data wg (  4) /   1.31688638449176626898494499748163134d-1 /
      data wg (  5) /   1.18194531961518417312377377711382287d-1 /
      data wg (  6) /   1.01930119817240435036750135480349876d-1 /
      data wg (  7) /   8.32767415767047487247581432220462061d-2 /
      data wg (  8) /   6.26720483341090635695065351870416063d-2 /
      data wg (  9) /   4.06014298003869413310399522749321098d-2 /
      data wg ( 10) /   1.76140071391521183118619623518528163d-2 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      center       = 0.5d0 * (a+b)
      hlfrun       = 0.5d0 * (b-a)
      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)

        call f(absc1, par, n, fval(1), deta(1), dtheta(1), &
               deta2(1), dtheta2(1), detadtheta(1))
        call f(absc2, par, n, fval(2), deta(2), dtheta(2), &
               deta2(2), dtheta2(2), detadtheta(2))

        result       = result    + (fval(1)    + fval(2))*wg(j)
        drdeta       = drdeta    + (deta(1)    + deta(2))*wg(j)
        drdtheta     = drdtheta  + (dtheta(1)  + dtheta(2))*wg(j)
        drdeta2      = drdeta2   + (deta2(1)   + deta2(2))*wg(j)
        drdtheta2    = drdtheta2 + (dtheta2(1) + dtheta2(2))*wg(j)
        drdetadtheta = drdetadtheta+(detadtheta(1)+detadtheta(2))*wg(j)
      enddo

      result       = result * hlfrun
      drdeta       = drdeta * hlfrun
      drdtheta     = drdtheta * hlfrun
      drdeta2      = drdeta2 * hlfrun
      drdtheta2    = drdtheta2 * hlfrun
      drdetadtheta = drdetadtheta * hlfrun
      return
      end subroutine dqleg020

      subroutine dqlag020(f,a,b,result,drdeta,drdtheta, &
                          drdeta2,drdtheta2,drdetadtheta,par,n)
      include 'implno.dek'
!
! 20 point gauss-laguerre rule for the fermi-dirac function.
! on input f is the external function defining the integrand
! f(x)=g(x)*w(x), where w(x) is the gaussian weight
! w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
! a is the lower end point of the interval, b is the higher end point,
! par is an array of constant parameters to be passed to the function f,
! and n is the length of the par array. on output result is the
! approximation from applying the 20-point gauss-laguerre rule.
! since the number of nodes is even, zero is not an abscissa.
!
! declare
      external    ::  f
      integer     ::  j,n
      real*8      ::  a,b,result,drdeta,drdtheta, &
                      drdeta2,drdtheta2,drdetadtheta,par(n), &
                      absc,wg(20),xg(20),fval,deta,dtheta, &
                      deta2,dtheta2,detadtheta


! the abscissae and weights are given for the interval (0,+inf).
! xg     - abscissae of the 20-point gauss-laguerre rule
! wg     - weights of the 20-point gauss rule. since f yet
!          includes the weight function, the values in wg
!          are actually exp(xg) times the standard
!          gauss-laguerre weights
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.05398896919887533666890045842150958d-2 /
      data xg (  2) /   3.72126818001611443794241388761146636d-1 /
      data xg (  3) /   9.16582102483273564667716277074183187d-1 /
      data xg (  4) /   1.70730653102834388068768966741305070d0 /
      data xg (  5) /   2.74919925530943212964503046049481338d0 /
      data xg (  6) /   4.04892531385088692237495336913333219d0 /
      data xg (  7) /   5.61517497086161651410453988565189234d0 /
      data xg (  8) /   7.45901745367106330976886021837181759d0 /
      data xg (  9) /   9.59439286958109677247367273428279837d0 /
      data xg ( 10) /   1.20388025469643163096234092988655158d1 /
      data xg ( 11) /   1.48142934426307399785126797100479756d1 /
      data xg ( 12) /   1.79488955205193760173657909926125096d1 /
      data xg ( 13) /   2.14787882402850109757351703695946692d1 /
      data xg ( 14) /   2.54517027931869055035186774846415418d1 /
      data xg ( 15) /   2.99325546317006120067136561351658232d1 /
      data xg ( 16) /   3.50134342404790000062849359066881395d1 /
      data xg ( 17) /   4.08330570567285710620295677078075526d1 /
      data xg ( 18) /   4.76199940473465021399416271528511211d1 /
      data xg ( 19) /   5.58107957500638988907507734444972356d1 /
      data xg ( 20) /   6.65244165256157538186403187914606659d1 /

      data wg (  1) /   1.81080062418989255451675405913110644d-1 /
      data wg (  2) /   4.22556767878563974520344172566458197d-1 /
      data wg (  3) /   6.66909546701848150373482114992515927d-1 /
      data wg (  4) /   9.15352372783073672670604684771868067d-1 /
      data wg (  5) /   1.16953970719554597380147822239577476d0 /
      data wg (  6) /   1.43135498592820598636844994891514331d0 /
      data wg (  7) /   1.70298113798502272402533261633206720d0 /
      data wg (  8) /   1.98701589079274721410921839275129020d0 /
      data wg (  9) /   2.28663578125343078546222854681495651d0 /
      data wg ( 10) /   2.60583472755383333269498950954033323d0 /
      data wg ( 11) /   2.94978373421395086600235416827285951d0 /
      data wg ( 12) /   3.32539578200931955236951937421751118d0 /
      data wg ( 13) /   3.74225547058981092111707293265377811d0 /
      data wg ( 14) /   4.21423671025188041986808063782478746d0 /
      data wg ( 15) /   4.76251846149020929695292197839096371d0 /
      data wg ( 16) /   5.42172604424557430380308297989981779d0 /
      data wg ( 17) /   6.25401235693242129289518490300707542d0 /
      data wg ( 18) /   7.38731438905443455194030019196464791d0 /
      data wg ( 19) /   9.15132873098747960794348242552950528d0 /
      data wg ( 20) /   1.28933886459399966710262871287485278d1 /


!           list of major variables
!           -----------------------
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula

      result       = 0.0d0
      drdeta       = 0.0d0
      drdtheta     = 0.0d0
      drdeta2      = 0.0d0
      drdtheta2    = 0.0d0
      drdetadtheta = 0.0d0

      do j=1,20
       absc = a+b*xg(j)

       call f(absc, par, n, fval, deta, dtheta, &
              deta2,dtheta2,detadtheta)

       result       = result    + fval*wg(j)
       drdeta       = drdeta    + deta*wg(j)
       drdtheta     = drdtheta  + dtheta*wg(j)
       drdeta2      = drdeta2   + deta2*wg(j)
       drdtheta2    = drdtheta2 + dtheta2*wg(j)
       drdetadtheta = drdetadtheta + detadtheta*wg(j)
      enddo

      result       = result*b
      drdeta       = drdeta*b
      drdtheta     = drdtheta*b
      drdeta2      = drdeta2*b
      drdtheta2    = drdtheta2*b
      drdetadtheta = drdetadtheta*b
      return
      end subroutine dqlag020

      subroutine locate(xx,n,x,j)
      implicit none
      save


! given an array xx of length n, and a value of x, this routine returns
! a value j such that x is between xx(j) and xx(j+1). the array xx must be
! monotonic. j=0 or j=n indicates that x is out of range. bisection is used
! to find the entry

! declare
      integer           n,j,jl,ju,jm
      double precision  xx(n),x

! initialize
      jl = 0
      ju = n+1

! compute a midpoint, and replace either the upper or lower limit
 10   if (ju-jl .gt. 1) then
       jm = (ju+jl)/2
       if ( (xx(n) .ge. xx(1)) .eqv. (x .ge. xx(jm)) ) then
        jl = jm
       else
        ju = jm
       end if
       goto 10
      end if
      if (x .eq. xx(1))then
        j = 1
      else if(x .eq. xx(n))then
        j = n - 1
      else
        j = jl
      end if
      return
    end subroutine locate
!---------------------------------------------------------------------

      subroutine coords(j,k)
      include 'implno.dek'
      include 'network.dek'

! write the coordinates of the isotopes being linked - useful for network plots

      integer :: j,k
      write(44,10) int(nion(j)), int(zion(j)), int(nion(k)), int(zion(k)), & 
                   ionam(j),ionam(k)
10    format(1x,4i4,2a8)
      return
      end subroutine coords


