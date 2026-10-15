!timreact21
!timrat21
!init_aprox21
!tfactors(temp)
!aprox21rat
!screen_aprox21(y)
!screen5(temp,den,zbar,abar,z2bar
!rhs21(y,rate,dydt,deriva)
!rate_pp(btemp,bden, 
!rate_png(btemp,bden, 
!rate_dpg(btemp,bden, 
!rate_he3ng(btemp,bden, 
!rate_he3he3(btemp,bden, 
!rate_he3he4(btemp,bden, 
!rate_tripalf(btemp,bden, 
!rate_c12pg(btemp,bden, 
!rate_n14pg(btemp,bden, 
!rate_n15pg(btemp,bden, 
!rate_n15pa(btemp,bden, 
!rate_o16pg(btemp,bden, 
!rate_n14ag(btemp,bden, 
!rate_c12ag(btemp,bden, 
!rate_c12c12(btemp,bden, 
!rate_c12o16(btemp,bden, 
!rate_o16o16(btemp,bden, 
!rate_o16ag(btemp,bden, 
!rate_ne20ag(btemp,bden, 
!rate_mg24ag(btemp,bden, 
!rate_mg24ap(btemp,bden, 
!rate_al27pg(btemp,bden, 
!rate_si28ag(btemp,bden, 
!rate_si28ap(btemp,bden, 
!rate_p31pg(btemp,bden, 
!rate_s32ag(btemp,bden, 
!rate_s32ap(btemp,bden, 
!rate_cl35pg(btemp,bden, 
!rate_ar36ag(btemp,bden, 
!rate_ar36ap(btemp,bden, 
!rate_k39pg(btemp,bden, 
!rate_ca40ag(btemp,bden, 
!rate_ca40ap(btemp,bden, 
!rate_sc43pg(btemp,bden, 
!rate_ti44ag(btemp,bden, 
!rate_ti44ap(btemp,bden, 
!rate_v47pg(btemp,bden, 
!rate_cr48ag(btemp,bden, 
!rate_cr48ap(btemp,bden, 
!rate_mn51pg(btemp,bden, 
!rate_fe52ag(btemp,bden, 
!rate_fe52ap(btemp,bden, 
!rate_co55pg(btemp,bden, 
!rate_fe52ng(btemp,bden, 
!rate_fe53ng(btemp,bden, 
!rate_fe54ng(btemp,bden, 
!rate_fe54pg(btemp,bden, 
!rate_fe54ap(btemp,bden, 
!rate_fe55ng(btemp,bden, 
!rate_fe56pg(btemp,bden,
!azbar_simple

  !     Get Timmes's reaction rates from MG solution
      subroutine timreact21(nd, pa, ydot)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'burn_common.dek'
!     Input
      integer   :: nd            
      real*8    :: pa(*), ydot(50)
!     Local
      logical   :: deriva
      integer   :: i      
      real*8    :: y(50), dydt(50),temp, den, pres, enuc
      deriva = .false.
!     Set temperature, density and pressure
      temp = pa(25+nd)
      den = pa(1)
      pres = pa(nd+2)
      call tfactors(temp)
      bden = den
      btemp = temp
      bpres = pres
!     Set abundances
      do i=1,21
         y(i) = pa(2+nd+i)
      enddo
!     Get rates
      call aprox21rat
!     Do screening
      screen_on = 1
      call screen_aprox21(y)      
      call rhs21(y,ratdum,dydt,deriva)
      do i=1,21
         ydot(i) = dydt(i)
      enddo
      return
    end subroutine timreact21

      subroutine timrat21(nd, pa, RR)
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
      real*8    :: y(50), temp, den, pres, enuc
      deriva = .false.
!     Set temperature, density and pressure
      temp = pa(25+nd)
      den = pa(1)
      pres = pa(nd+2)      
      call tfactors(temp)
      bden = den
      btemp = temp
      bpres = pres
!     Set abundances
      do i=1,21
         y(i) = pa(2+nd+i)
      enddo
!     Get rates
      call aprox21rat
!     Do screening
      screen_on = 1
      call screen_aprox21(y)
      do i=1,67
         RR(i) = ratdum(i)
      enddo      
      return
    end subroutine timrat21

      subroutine aprox21_init(eaion, ezion, ebion)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this routine initializes stuff for the aprox13 network
!
!     Output
      real*8 :: eaion(100), ezion(100), ebion(100)          

! declare
      integer            :: i
      real*8, parameter  :: mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/(clight*clight)
      logical, parameter :: icoords = .true.


! set the size of the network and the number of rates
      idnet   = idaprox21
      ionmax  = 21
      ionbeg  = 1
      ionend  = ionmax
      iener   = ionmax + 1
      itemp   = ionmax + 2
      iden    = ionmax + 3
      ivelx   = ionmax + 4
      iposx   = ionmax + 5
      neqs    = iposx
      nrat    = 112
      netname = 'aprox21'


! set the id numbers of the elements
      ih1   = ionbeg
      ihe3  = ionbeg+ 1
      ihe4  = ionbeg + 2
      ic12  = ionbeg + 3
      in14  = ionbeg + 4
      io16  = ionbeg + 5
      ine20 = ionbeg + 6
      img24 = ionbeg + 7
      isi28 = ionbeg + 8
      is32  = ionbeg + 9
      iar36 = ionbeg + 10
      ica40 = ionbeg + 11
      iti44 = ionbeg + 12
      icr48 = ionbeg + 13
      icr56 = ionbeg + 14
      ife52 = ionbeg + 15
      ife54 = ionbeg + 16
      ife56 = ionbeg + 17
      ini56 = ionbeg + 18
      ineut = ionbeg + 19
      iprot = ionbeg + 20


! set the names of the elements
      ionam(ih1)   = 'h1  '
      ionam(ihe3)  = 'he3 '
      ionam(ihe4)  = 'he4 '
      ionam(ic12)  = 'c12 '
      ionam(in14)  = 'n14 '
      ionam(io16)  = 'o16 '
      ionam(ine20) = 'ne20'
      ionam(img24) = 'mg24'
      ionam(isi28) = 'si28'
      ionam(is32)  = 's32 '
      ionam(iar36) = 'ar36'
      ionam(ica40) = 'ca40'
      ionam(iti44) = 'ti44'
      ionam(icr48) = 'cr48'
      ionam(icr56) = 'cr56'
      ionam(ife52) = 'fe52'
      ionam(ife54) = 'fe54'
      ionam(ife56) = 'fe56'
      ionam(ini56) = 'ni56'
      ionam(ineut) = 'neut'
      ionam(iprot) = 'prot'

      ionam(iener) = 'ener '
      ionam(itemp) = 'temp '
      ionam(iden)  = 'den  '
      ionam(ivelx) = 'velx '
      ionam(iposx) = 'posx '


! set the number of nucleons in the element
      aion(ih1)   = 1.0d0
      aion(ihe3)  = 3.0d0
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(in14)  = 14.0d0
      aion(io16)  = 16.0d0
      aion(ine20) = 20.0d0
      aion(img24) = 24.0d0
      aion(isi28) = 28.0d0
      aion(is32)  = 32.0d0
      aion(iar36) = 36.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(icr56) = 56.0d0
      aion(ife52) = 52.0d0
      aion(ife54) = 54.0d0
      aion(ife56) = 56.0d0
      aion(ini56) = 56.0d0
      aion(ineut) = 1.0d0
      aion(iprot) = 1.0d0

! set the number of protons in the element
      zion(ih1)   = 1.0d0
      zion(ihe3)  = 2.0d0
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(in14)  = 7.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0
      zion(img24) = 12.0d0
      zion(isi28) = 14.0d0
      zion(is32)  = 16.0d0
      zion(iar36) = 18.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(icr56) = 24.0d0
      zion(ife52) = 26.0d0
      zion(ife54) = 26.0d0
      zion(ife56) = 26.0d0
      zion(ini56) = 28.0d0
      zion(ineut) = 0.0d0
      zion(iprot) = 1.0d0


! set the binding energy of the element
      bion(ih1)   = 0.0d0
      bion(ihe3)  = 7.71819d0
      bion(ihe4)  = 28.29603d0
      bion(ic12)  = 92.16294d0
      bion(in14)  = 104.65998d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0
      bion(img24) = 198.25790d0
      bion(isi28) = 236.53790d0
      bion(is32)  = 271.78250d0
      bion(iar36) = 306.72020d0
      bion(ica40) = 342.05680d0
      bion(iti44) = 375.47720d0
      bion(icr48) = 411.46900d0
      bion(icr56) = 488.4970d0
      bion(ife52) = 447.70800d0
      bion(ife54) = 471.7696d0
      bion(ife56) = 492.2450d0
      bion(ini56) = 484.00300d0
      bion(ineut) = 0.0d0
      bion(iprot) = 0.0d0


! set the number of neutrons and mass
       nion(ionbeg:ionend) = aion(ionbeg:ionend) - zion(ionbeg:ionend)

! mass of each isotope
       mion(ionbeg:ionend) = nion(ionbeg:ionend)*mn + zion(ionbeg:ionend)*mp - bion(ionbeg:ionend)*mev2gr

! molar mass
       wion(ionbeg:ionend) = avo * mion(ionbeg:ionend)

! a common approximation
      wion(ionbeg:ionend) = aion(ionbeg:ionend)

! set the partition functions - statistical weights, ground-state only here
      wpart(ionbeg:ionend) = 1.0d0
      wpart(ih1)   = 2.0d0
      wpart(ineut) = 2.0d0
      wpart(iprot) = 2.0d0


! set the id numbers of the reaction rates
      ir3a   = 1
      irg3a  = 2
      ircag  = 3
      ir1212 = 4
      ir1216 = 5
      ir1616 = 6
      iroga  = 7
      iroag  = 8
      irnega = 9
      irneag = 10
      irmgga = 11
      irmgag = 12
      irsiga = 13
      irmgap = 14
      iralpa = 15
      iralpg = 16
      irsigp = 17
      irsiag = 18
      irsga  = 19
      irsiap = 20
      irppa  = 21
      irppg  = 22
      irsgp  = 23
      irsag  = 24
      irarga = 25
      irsap  = 26
      irclpa = 27
      irclpg = 28
      irargp = 29
      irarag = 30
      ircaga = 31
      irarap = 32
      irkpa  = 33
      irkpg  = 34
      ircagp = 35
      ircaag = 36
      irtiga = 37
      ircaap = 38
      irscpa = 39
      irscpg = 40
      irtigp = 41
      irtiag = 42
      ircrga = 43
      irtiap = 44
      irvpa  = 45
      irvpg  = 46
      ircrgp = 47
      ircrag = 48
      irfega = 49
      ircrap = 50
      irmnpa = 51
      irmnpg = 52
      irfegp = 53
      irfeag = 54
      irniga = 55
      irfeap = 56
      ircopa = 57
      ircopg = 58
      irnigp = 59

! for fe54 photodisintegration
      ir52ng = 60
      ir53gn = 61
      ir53ng = 62
      ir54gn = 63
      irfepg = 64
      ircogp = 65

! for he4 photodisintegration
      irheng = 66
      irhegn = 67
      irhng  = 68
      irdgn  = 69
      irdpg  = 70
      irhegp = 71

! weak reactions
      irpen   = 72
      irnep   = 73
      irn56ec = 74

! ppchain
      irpp    = 75
      ir33    = 76
      irhe3ag = 77

! cno cycles
      ircpg  = 78
      irnpg  = 79
      ifa    = 80
      ifg    = 81
      iropg  = 82
      irnag  = 83

! for reactions to fe56 
      ir54ng   = 84
      ir55gn   = 85
      ir55ng   = 86
      ir56gn   = 87
      irfe54ap = 88
      irco57pa = 89
      irfe56pg = 90
      irco57gp = 91

! the equilibrium links
      irr1   = 92
      irs1   = 93
      irt1   = 94
      iru1   = 95
      irv1   = 96
      irw1   = 97
      irx1   = 98

      ir1f54 = 99
      ir2f54 = 100
      ir3f54 = 101
      ir4f54 = 102
      ir5f54 = 103
      ir6f54 = 104
      ir7f54 = 105
      ir8f54 = 106

      iralf1 = 107
      iralf2 = 108

      irfe56_aux1 = 109
      irfe56_aux2 = 110
      irfe56_aux3 = 111
      irfe56_aux4 = 112

      ! eaion, ezion, ebion
      eaion(1:ionmax) = aion(1:ionmax)
      ezion(1:ionmax) = zion(1:ionmax)
      ebion(1:ionmax) = bion(1:ionmax)      
      return
      end subroutine aprox21_init

    
      subroutine aprox21rat
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine generates nuclear reaction rates for the aprox21 network.

! declare
      integer  ::  i
      real*8   ::  rrate,drratedt,drratedd,ff1,dff1dt,dff1dd, &
                   ff2,dff2dt,dff2dd,tot,dtotdt,dtotdd,invtot


! zero the rates
      ratraw(1:nrat)    = 0.0d0
      dratrawdt(1:nrat) = 0.0d0
      dratrawdd(1:nrat) = 0.0d0

      if (btemp .lt. 1.0e6) return



! get the temperature factors
      call tfactors(btemp)


! p(p,e+nu)d
      call rate_pp(btemp,bden, &
                   ratraw(irpp),dratrawdt(irpp),dratrawdd(irpp), &
                   rrate,drratedt,drratedd)

! p(n,g)2d
      call rate_png(btemp,bden, &
                    ratraw(irhng),dratrawdt(irhng),dratrawdd(irhng), &
                    ratraw(irdgn),dratrawdt(irdgn),dratrawdd(irdgn))

! d(p,g)he3
      call rate_dpg(btemp,bden, &
                    ratraw(irdpg),dratrawdt(irdpg),dratrawdd(irdpg), &
                    ratraw(irhegp),dratrawdt(irhegp),dratrawdd(irhegp))

! he3(n,g)he4
      call rate_he3ng(btemp,bden, &
                      ratraw(irheng),dratrawdt(irheng),dratrawdd(irheng), &
                      ratraw(irhegn),dratrawdt(irhegn),dratrawdd(irhegn))

! he3(he3,2p)he4
      call rate_he3he3(btemp,bden, &
                       ratraw(ir33),dratrawdt(ir33),dratrawdd(ir33), &
                       rrate,drratedt,drratedd)

! he3(he4,g)be7
      call rate_he3he4(btemp,bden, &
                       ratraw(irhe3ag),dratrawdt(irhe3ag),dratrawdd(irhe3ag), &
                       rrate,drratedt,drratedd)

! triple alpha to c12
      call rate_tripalf(btemp,bden, &
                    ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
                    ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

! c12(p,g)13n
      call rate_c12pg(btemp,bden, &
                      ratraw(ircpg),dratrawdt(ircpg),dratrawdd(ircpg), &
                      rrate,drratedt,drratedd)

! n14(p,g)o15
      call rate_n14pg(btemp,bden, &
                      ratraw(irnpg),dratrawdt(irnpg),dratrawdd(irnpg), &
                      rrate,drratedt,drratedd)

! fraction fg of n15 goes (pg) to o16, fraction fa of n15 goes (pa) to c12
! result is  n14+2p  goes to o16    at a rate = rnpg*fg
!                    goes to c12+a  at a rate = rnpg*fa

      call rate_n15pg(btemp,bden, &
                      ff1,dff1dt,dff1dd, &
                      rrate,drratedt,drratedd)

      call rate_n15pa(btemp,bden, &
                      ff2,dff2dt,dff2dd, &
                      rrate,drratedt,drratedd)

      tot            = ff1 + ff2
      dtotdt         = dff1dt + dff2dt
      dtotdd         = dff1dd + dff2dd
      invtot         = 1.0d0/tot

      ratraw(ifa)    = ff2 * invtot
      dratrawdt(ifa) = dff2dt * invtot - ff2 * invtot*invtot * dtotdt
      dratrawdd(ifa) = dff2dd * invtot - ff2 * invtot*invtot * dtotdd

      ratraw(ifg)    = 1.0d0 - ratraw(ifa)
      dratrawdt(ifg) = -dratrawdt(ifa)
      dratrawdd(ifg) = -dratrawdd(ifa)


! o16(p,g)f17
      call rate_o16pg(btemp,bden, &
                      ratraw(iropg),dratrawdt(iropg),dratrawdd(iropg), &
                      rrate,drratedt,drratedd)

! n14(a,g)f18
      call rate_n14ag(btemp,bden, &
                      ratraw(irnag),dratrawdt(irnag),dratrawdd(irnag), &
                      rrate,drratedt,drratedd)

! c12(a,g)o16
      call rate_c12ag(btemp,bden, &
                      ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
                      ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

! c12 + c12
      call rate_c12c12(btemp,bden, &
                       ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
                       rrate,drratedt,drratedd)

! c12 + o16
      call rate_c12o16(btemp,bden, &
                       ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
                       rrate,drratedt,drratedd)

! o16 + o16
      call rate_o16o16(btemp,bden, &
                       ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
                       rrate,drratedt,drratedd)

! o16(a,g)ne20
      call rate_o16ag(btemp,bden, &
                      ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
                      ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

! ne20(a,g)mg24
      call rate_ne20ag(btemp,bden, &
                       ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
                       ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

! mg24(a,g)si28
      call rate_mg24ag(btemp,bden, &
                       ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
                       ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

! mg24(a,p)al27
      call rate_mg24ap(btemp,bden, &
                       ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap), &
                       ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))

! al27(p,g)si28
      call rate_al27pg(btemp,bden, &
                       ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg), &
                       ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))

! si28(a,g)s32
      call rate_si28ag(btemp,bden, &
                       ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag), &
                       ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))

! si28(a,p)p31
      call rate_si28ap(btemp,bden, &
                       ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap), &
                       ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))

! p31(p,g)s32
      call rate_p31pg(btemp,bden, &
                       ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg), &
                       ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))

! s32(a,g)ar36
      call rate_s32ag(btemp,bden, &
                       ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag), &
                       ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))

! s32(a,p)cl35
      call rate_s32ap(btemp,bden, &
                       ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap), &
                       ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))

! cl35(p,g)ar36
      call rate_cl35pg(btemp,bden, &
                       ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg), &
                       ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))

! ar36(a,g)ca40
      call rate_ar36ag(btemp,bden, &
                       ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag), &
                       ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))

! ar36(a,p)k39
      call rate_ar36ap(btemp,bden, &
                       ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap), &
                       ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))

! k39(p,g)ca40
      call rate_k39pg(btemp,bden, &
                       ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg), &
                       ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))

! ca40(a,g)ti44
      call rate_ca40ag(btemp,bden, &
                       ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
                       ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

! ca40(a,p)sc43
      call rate_ca40ap(btemp,bden, &
                       ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap), &
                       ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))

! sc43(p,g)ti44
      call rate_sc43pg(btemp,bden, &
                       ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg), &
                       ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))

! ti44(a,g)cr48
      call rate_ti44ag(btemp,bden, &
                       ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag), &
                       ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

! ti44(a,p)v47
      call rate_ti44ap(btemp,bden, &
                       ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
                       ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))

! v47(p,g)cr48
      call rate_v47pg(btemp,bden, &
                       ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg), &
                       ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))

! cr48(a,g)fe52
      call rate_cr48ag(btemp,bden, &
                       ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag), &
                       ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))

! cr48(a,p)mn51
      call rate_cr48ap(btemp,bden, &
                       ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap), &
                       ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))

! mn51(p,g)fe52
      call rate_mn51pg(btemp,bden, &
                       ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg), &
                       ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))

! fe52(a,g)ni56
      call rate_fe52ag(btemp,bden, &
                       ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag), &
                       ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))

! fe52(a,p)co55
      call rate_fe52ap(btemp,bden, &
                       ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap), &
                       ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))

! co55(p,g)ni56
      call rate_co55pg(btemp,bden, &
                       ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg), &
                       ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))

! fe52(n,g)fe53
      call rate_fe52ng(btemp,bden, &
                       ratraw(ir52ng),dratrawdt(ir52ng),dratrawdd(ir52ng), &
                       ratraw(ir53gn),dratrawdt(ir53gn),dratrawdd(ir53gn))

! fe53(n,g)fe54
      call rate_fe53ng(btemp,bden, &
                       ratraw(ir53ng),dratrawdt(ir53ng),dratrawdd(ir53ng), &
                       ratraw(ir54gn),dratrawdt(ir54gn),dratrawdd(ir54gn))

! fe54(n,g)fe55
      call rate_fe54ng(btemp,bden, &
                       ratraw(ir54ng),dratrawdt(ir54ng),dratrawdd(ir54ng), &
                       ratraw(ir55gn),dratrawdt(ir55gn),dratrawdd(ir55gn))


! fe54(p,g)co55
      call rate_fe54pg(btemp,bden, &
                       ratraw(irfepg),dratrawdt(irfepg),dratrawdd(irfepg), &
                       ratraw(ircogp),dratrawdt(ircogp),dratrawdd(ircopg))

! fe54(a,p)co57
      call rate_fe54ap(btemp,bden, &
                       ratraw(irfe54ap),dratrawdt(irfe54ap),dratrawdd(irfe54ap), &
                       ratraw(irco57pa),dratrawdt(irco57pa),dratrawdd(irco57pa))

! fe55(n,g)fe56
      call rate_fe55ng(btemp,bden, &
                       ratraw(ir55ng),dratrawdt(ir55ng),dratrawdd(ir55ng), &
                       ratraw(ir56gn),dratrawdt(ir56gn),dratrawdd(ir56gn))

! fe56(p,g)co57
      call rate_fe56pg(btemp,bden, &
                       ratraw(irfe56pg),dratrawdt(irfe56pg),dratrawdd(irfe56pg), &
                       ratraw(irco57gp),dratrawdt(irco57gp),dratrawdd(irco57gp))






! shut down
! for fe54 photodisintegration
!      ratraw(ir52ng)    = 0.0d0
!      dratrawdt(ir52ng) = 0.0d0
!      dratrawdd(ir52ng) = 0.0d0
!      ratraw(ir53gn)    = 0.0d0
!      dratrawdt(ir53gn) = 0.0d0
!      dratrawdd(ir53gn) = 0.0d0
!      ratraw(ir53ng)    = 0.0d0
!      dratrawdt(ir53ng) = 0.0d0
!      dratrawdd(ir53ng) = 0.0d0
!      ratraw(ir54gn)    = 0.0d0
!      dratrawdt(ir54gn) = 0.0d0
!      dratrawdd(ir54gn) = 0.0d0
!      ratraw(irfepg)    = 0.0d0
!      dratrawdt(irfepg) = 0.0d0
!      dratrawdd(irfepg) = 0.0d0
!      ratraw(ircogp)    = 0.0d0
!      dratrawdt(ircogp) = 0.0d0
!      dratrawdd(ircogp) = 0.0d0

! for he4 photodisintegration
!      ratraw(irheng)     = 0.0d0
!      dratrawdt(irheng)  = 0.0d0
!      dratrawdd(irheng)  = 0.0d0
!      ratraw(irhegn)     = 0.0d0
!      dratrawdt(irhegn)  = 0.0d0
!      dratrawdd(irhegn)  = 0.0d0
!      ratraw(irhng)      = 0.0d0
!      dratrawdt(irhng)   = 0.0d0
!      dratrawdd(irhng)   = 0.0d0
!      ratraw(irdgn)      = 0.0d0
!      dratrawdt(irdgn)   = 0.0d0
!      dratrawdd(irdgn)   = 0.0d0
!      ratraw(irdpg)      = 0.0d0
!      dratrawdt(irdpg)   = 0.0d0
!      dratrawdd(irdpg)   = 0.0d0
!      ratraw(irhegp)     = 0.0d0
!      dratrawdt(irhegp)  = 0.0d0
!      dratrawdd(irhegp)  = 0.0d0
!
! for weak reactions
!      ratraw(irpen)      = 0.0d0
!      dratrawdt(irpen)   = 0.0d0
!      dratrawdd(irpen)   = 0.0d0
!      ratraw(irnep)      = 0.0d0
!      dratrawdt(irnep)   = 0.0d0
!      dratrawdd(irnep)   = 0.0d0
!      ratraw(irn56ec)    = 0.0d0
!      dratrawdt(irn56ec) = 0.0d0
!      dratrawdd(irn56ec) = 0.0d0
!
! ppchain
!      ratraw(irpp)       = 0.0d0
!      dratrawdt(irpp)    = 0.0d0
!      dratrawdd(irpp)    = 0.0d0
!      ratraw(ir33)       = 0.0d0
!      dratrawdt(ir33)    = 0.0d0
!      dratrawdd(ir33)    = 0.0d0
!      ratraw(irhe3ag)    = 0.0d0
!      dratrawdt(irhe3ag) = 0.0d0
!      dratrawdd(irhe3ag) = 0.0d0
!
! cno cycles
!      ratraw(ircpg)      = 0.0d0
!      dratrawdt(ircpg)   = 0.0d0
!      dratrawdd(ircpg)   = 0.0d0
!      ratraw(irnpg)      = 0.0d0
!      dratrawdt(irnpg)   = 0.0d0
!      dratrawdd(irnpg)   = 0.0d0
!      ratraw(iropg)      = 0.0d0
!      dratrawdt(iropg)   = 0.0d0
!      dratrawdd(iropg)   = 0.0d0
!      ratraw(ifa)        = 0.0d0
!      dratrawdt(ifa)     = 0.0d0
!      dratrawdd(ifa)     = 0.0d0
!      ratraw(ifg)        = 0.0d0
!      dratrawdt(ifg)     = 0.0d0
!      dratrawdd(ifg)     = 0.0d0
!      ratraw(irnag)      = 0.0d0
!      dratrawdt(irnag)   = 0.0d0
!      dratrawdd(irnag)   = 0.0d0



! for a strict alpha chain with only (a,g) and (g,a) reactions
! shut down the (a,p) (p,a) (p,g) and (g,p) rates
!      ratraw(irmgap) = 1.0d-40
!      ratraw(iralpa) = 0.0d0
!      ratraw(iralpg) = 1.0d-40
!      ratraw(irsigp) = 1.0d-40

!      ratraw(irsiga)  = 0.0d0
!      ratraw(irsiap) = 1.0d-40
!      ratraw(irppa)  = 0.0d0
!      ratraw(irppg)  = 1.0d-40
!      ratraw(irsgp)  = 1.0d-40

!      ratraw(irsag)  = 0.0d0
!      ratraw(irsap)  = 0.0d0

!      ratraw(irclpa) = 0.0d0
!      ratraw(irclpg) = 1.0d-40
!      ratraw(irargp) = 1.0d-40

!      ratraw(irarap) = 1.0d-40
!      ratraw(irkpa)  = 0.0d0
!      ratraw(irkpg)  = 1.0d-40
!      ratraw(ircagp) = 1.0d-40

!      ratraw(ircaap) = 1.0d-40
!      ratraw(irscpa) = 0.0d0
!      ratraw(irscpg) = 1.0d-40
!      ratraw(irtigp) = 1.0d-40

!      ratraw(irtiap) = 1.0d-40
!      ratraw(irvpa)  = 0.0d0
!      ratraw(irvpg)  = 1.0d-40
!      ratraw(ircrgp) = 1.0d-40

!      ratraw(ircrap) = 1.0d-40
!      ratraw(irmnpa) = 0.0d0
!      ratraw(irmnpg) = 1.0d-40
!      ratraw(irfegp) = 1.0d-40

! shutting down this (a,p)(p,g) sequence
! will make aprox21 exactly like iso13
!      ratraw(irfeap) = 1.0d-40
!      ratraw(ircopa) = 0.0d0
!      ratraw(ircopg) = 1.0d-40
!      ratraw(irnigp) = 1.0d-40

      return
      end subroutine aprox21rat

      subroutine screen_aprox21(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine computes the screening factors
! and applies them to the raw reaction rates,
! producing the final reaction rates used by the
! right hand sides and jacobian matrix elements

! this routine assumes screen_on = 1 or = 0 and that
! screen_init = 1 or = 0 has been set at a higher level,
! presumably in the top level driver

! declare the pass
      real*8  :: y(*)

! declare
      integer :: i,jscr,init=1
      integer :: maxrate_loc(1),minrate_loc(1)
      real*8  :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                 sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar, &
                 denom,denomdt,denomdd,xx,zz

! initialize
       ratdum(1:nrat) = ratraw(1:nrat)  ; dratdumdt(1:nrat) = dratrawdt(1:nrat)  ; dratdumdd(1:nrat) = dratrawdd(1:nrat)
       scfac(1:nrat)  = 1.0d0  ; dscfacdt(1:nrat)  = 0.0d0  ; dscfacdd(1:nrat)  = 0.0d0
     

! if screening is on
      if (screen_on .ne. 0) then


! with the passed composition compute abar, zbar and z**2 bar

       call azbar_simple(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)
       z2bar  = abar * sum( zion(1:ionmax) * zion(1:ionmax) * y(1:ionmax) )


! first the always fun triple alpha and its inverse
      jscr = 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),4.0d0,8.0d0, &
                   jscr,init,sc2a,sc2adt,sc2add)

      sc3a   = sc1a * sc2a
      sc3adt = sc1adt*sc2a + sc1a*sc2adt
      sc3add = sc1add*sc2a + sc1a*sc2add

      ratdum(ir3a)    = ratraw(ir3a) * sc3a
      dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
      dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

      scfac(ir3a) = sc3a  ; dscfacdt(ir3a)  = sc3adt  ; dscfacdd(ir3a)  = sc3add

      ratdum(irg3a)    = ratraw(irg3a) * sc3a
      dratdumdt(irg3a) = dratrawdt(irg3a)*sc3a + ratraw(irg3a)*sc3adt
      dratdumdd(irg3a) = dratrawdd(irg3a)*sc3a + ratraw(irg3a)*sc3add

      scfac(irg3a)  = sc3a  ; dscfacdt(irg3a)  = sc3adt  ; dscfacdd(irg3a)  = sc3add



! c12 to o16
! c12(a,g)o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircag)     = ratraw(ircag) * sc1a
      dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
      dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

      scfac(ircag)  = sc1a  ; dscfacdt(ircag)   = sc1adt  ; dscfacdt(ircag)   = sc1add

      ratdum(iroga)     = ratraw(iroga) * sc1a
      dratdumdt(iroga)  = dratrawdt(iroga)*sc1a + ratraw(iroga)*sc1adt
      dratdumdd(iroga)  = dratrawdd(iroga)*sc1a + ratraw(iroga)*sc1add

      scfac(iroga)  = sc1a  ; dscfacdt(iroga)   = sc1adt  ; dscfacdt(iroga)   = sc1add



! c12 + c12
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ic12),aion(ic12), &
                   jscr,init,sc1a,sc1adt,sc1add)


      ratdum(ir1212)    = ratraw(ir1212) * sc1a
      dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
      dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

      scfac(ir1212)     = sc1a  ; dscfacdt(ir1212)  = sc1adt  ; dscfacdd(ir1212)  = sc1add



! c12 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1216)    = ratraw(ir1216) * sc1a
      dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
      dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

      scfac(ir1216)     = sc1a  ; dscfacdt(ir1216)  = sc1adt  ; dscfacdd(ir1216)  = sc1add


! o16 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)


      ratdum(ir1616)    = ratraw(ir1616) * sc1a
      dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
      dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

      scfac(ir1616)     = sc1a  ; dscfacdt(ir1616)  = sc1adt  ; dscfacdd(ir1616)  = sc1add



! o16 to ne20
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! o16(a,g)ne20
      ratdum(iroag)    = ratraw(iroag) * sc1a
      dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
      dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

      scfac(iroag)  = sc1a  ; dscfacdt(iroag)  = sc1adt  ; dscfacdd(iroag)  = sc1add

      ratdum(irnega)    = ratraw(irnega) * sc1a
      dratdumdt(irnega) = dratrawdt(irnega)*sc1a + ratraw(irnega)*sc1adt
      dratdumdd(irnega) = dratrawdd(irnega)*sc1a + ratraw(irnega)*sc1add

      scfac(irnega)  = sc1a  ; dscfacdt(irnega)  = sc1adt  ; dscfacdd(irnega)  = sc1add


! ne20 to mg24
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ine20),aion(ine20),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ne20(a,g)mg24
      ratdum(irneag)    = ratraw(irneag) * sc1a
      dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
      dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

      scfac(irneag) = sc1a  ; dscfacdt(irneag)  = sc1adt  ; dscfacdd(irneag)  = sc1add

      ratdum(irmgga)    = ratraw(irmgga) * sc1a
      dratdumdt(irmgga) = dratrawdt(irmgga)*sc1a + ratraw(irmgga)*sc1adt
      dratdumdd(irmgga) = dratrawdd(irmgga)*sc1a + ratraw(irmgga)*sc1add

      scfac(irmgga) = sc1a  ; dscfacdt(irmgga)  = sc1adt  ; dscfacdd(irmgga)  = sc1add



! mg24 to si28
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(img24),aion(img24),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! mg24(a,g)si28
      ratdum(irmgag)    = ratraw(irmgag) * sc1a
      dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
      dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

      scfac(irmgag) = sc1a ; dscfacdt(irmgag)  = sc1adt  ; dscfacdd(irmgag)  = sc1add

      ratdum(irsiga)    = ratraw(irsiga) * sc1a
      dratdumdt(irsiga) = dratrawdt(irsiga)*sc1a + ratraw(irsiga)*sc1adt
      dratdumdd(irsiga) = dratrawdd(irsiga)*sc1a + ratraw(irsiga)*sc1add

      scfac(irsiga) = sc1a ; dscfacdt(irsiga)  = sc1adt  ; dscfacdd(irsiga)  = sc1add


! mg24(a,p)al27
      ratdum(irmgap)    = ratraw(irmgap) * sc1a
      dratdumdt(irmgap) = dratrawdt(irmgap)*sc1a + ratraw(irmgap)*sc1adt
      dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + ratraw(irmgap)*sc1add

      scfac(irmgap)     = sc1a  ; dscfacdt(irmgap)  = sc1adt  ; dscfacdd(irmgap)  = sc1add

      ratdum(iralpa)    = ratraw(iralpa) * sc1a
      dratdumdt(iralpa) = dratrawdt(iralpa)*sc1a + ratraw(iralpa)*sc1adt
      dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + ratraw(iralpa)*sc1add

      scfac(iralpa)     = sc1a  ; dscfacdt(iralpa)  = sc1adt  ; dscfacdd(iralpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   13.0d0,27.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! al27(p,g)si28
      ratdum(iralpg)    = ratraw(iralpg) * sc1a
      dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
      dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

      scfac(iralpg)     = sc1a  ; dscfacdt(iralpg)  = sc1adt  ; dscfacdd(iralpg)  = sc1add

      ratdum(irsigp)    = ratraw(irsigp) * sc1a
      dratdumdt(irsigp) = dratrawdt(irsigp)*sc1a + ratraw(irsigp)*sc1adt
      dratdumdd(irsigp) = dratrawdd(irsigp)*sc1a + ratraw(irsigp)*sc1add

      scfac(irsigp)     = sc1a  ; dscfacdt(irsigp)  = sc1adt  ; dscfacdd(irsigp)  = sc1add



! si28 to s32
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(isi28),aion(isi28),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! si28(a,g)s32
      ratdum(irsiag)    = ratraw(irsiag) * sc1a
      dratdumdt(irsiag) = dratrawdt(irsiag)*sc1a + ratraw(irsiag)*sc1adt
      dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + ratraw(irsiag)*sc1add

      scfac(irsiag)     = sc1a  ; dscfacdt(irsiag)  = sc1adt  ; dscfacdd(irsiag)  = sc1add

      ratdum(irsga)    = ratraw(irsga) * sc1a
      dratdumdt(irsga) = dratrawdt(irsga)*sc1a + ratraw(irsga)*sc1adt
      dratdumdd(irsga) = dratrawdd(irsga)*sc1a + ratraw(irsga)*sc1add

      scfac(irsga)     = sc1a  ; dscfacdt(irsga)  = sc1adt  ; dscfacdd(irsga)  = sc1add


! si28(a,p)p31
      ratdum(irsiap)    = ratraw(irsiap) * sc1a
      dratdumdt(irsiap) = dratrawdt(irsiap)*sc1a + ratraw(irsiap)*sc1adt
      dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + ratraw(irsiap)*sc1add

      scfac(irsiap)     = sc1a  ; dscfacdt(irsiap)  = sc1adt  ; dscfacdd(irsiap)  = sc1add

      ratdum(irppa)     = ratraw(irppa) * sc1a
      dratdumdt(irppa)  = dratrawdt(irppa)*sc1a  + ratraw(irppa)*sc1adt
      dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + ratraw(irppa)*sc1add

      scfac(irppa)      = sc1a  ; dscfacdt(irppa)   = sc1adt  ; dscfacdd(irppa)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   15.0d0,31.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! p31(p,g)s32
      ratdum(irppg)     = ratraw(irppg) * sc1a
      dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
      dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

      scfac(irppg)      = sc1a  ; dscfacdt(irppg)   = sc1adt  ; dscfacdd(irppg)   = sc1add

      ratdum(irsgp)     = ratraw(irsgp) * sc1a
      dratdumdt(irsgp)  = dratrawdt(irsgp)*sc1a + ratraw(irsgp)*sc1adt
      dratdumdd(irsgp)  = dratrawdd(irsgp)*sc1a + ratraw(irsgp)*sc1add

      scfac(irsgp)      = sc1a  ; dscfacdt(irsgp)   = sc1adt  ; dscfacdd(irsgp)   = sc1add



! s32 to ar36
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(is32),aion(is32),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! s32(a,g)ar36
      ratdum(irsag)     = ratraw(irsag) * sc1a
      dratdumdt(irsag)  = dratrawdt(irsag)*sc1a + ratraw(irsag)*sc1adt
      dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + ratraw(irsag)*sc1add

      scfac(irsag)      = sc1a ; dscfacdt(irsag)   = sc1adt ; dscfacdd(irsag)   = sc1add

      ratdum(irarga)     = ratraw(irarga) * sc1a
      dratdumdt(irarga)  = dratrawdt(irarga)*sc1a + ratraw(irarga)*sc1adt
      dratdumdd(irarga)  = dratrawdd(irarga)*sc1a + ratraw(irarga)*sc1add

      scfac(irarga)      = sc1a ; dscfacdt(irarga)   = sc1adt ; dscfacdd(irarga)   = sc1add

! s32(a,p)cl35
      ratdum(irsap)     = ratraw(irsap) * sc1a
      dratdumdt(irsap)  = dratrawdt(irsap)*sc1a + ratraw(irsap)*sc1adt
      dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + ratraw(irsap)*sc1add

      scfac(irsap)      = sc1a  ; dscfacdt(irsap)   = sc1adt  ; dscfacdd(irsap)   = sc1add

      ratdum(irclpa)    = ratraw(irclpa) * sc1a
      dratdumdt(irclpa) = dratrawdt(irclpa)*sc1a + ratraw(irclpa)*sc1adt
      dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + ratraw(irclpa)*sc1add

      scfac(irclpa)     = sc1a  ; dscfacdt(irclpa)  = sc1adt  ; dscfacdt(irclpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   17.0d0,35.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! cl35(p,g)ar36
      ratdum(irclpg)    = ratraw(irclpg) * sc1a
      dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
      dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

      scfac(irclpg)     = sc1a  ; dscfacdt(irclpg)  = sc1adt  ; dscfacdd(irclpg)  = sc1add

      ratdum(irargp)    = ratraw(irargp) * sc1a
      dratdumdt(irargp) = dratrawdt(irargp)*sc1a + ratraw(irargp)*sc1adt
      dratdumdd(irargp) = dratrawdd(irargp)*sc1a + ratraw(irargp)*sc1add

      scfac(irargp)     = sc1a  ; dscfacdt(irargp)  = sc1adt  ; dscfacdd(irargp)  = sc1add



! ar36 to ca40
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iar36),aion(iar36),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ar36(a,g)ca40
      ratdum(irarag)    = ratraw(irarag) * sc1a
      dratdumdt(irarag) = dratrawdt(irarag)*sc1a + ratraw(irarag)*sc1adt
      dratdumdd(irarag) = dratrawdd(irarag)*sc1a + ratraw(irarag)*sc1add

      scfac(irarag)     = sc1a  ; dscfacdt(irarag)  = sc1adt  ; dscfacdd(irarag)  = sc1add

      ratdum(ircaga)    = ratraw(ircaga) * sc1a
      dratdumdt(ircaga) = dratrawdt(ircaga)*sc1a + ratraw(ircaga)*sc1adt
      dratdumdd(ircaga) = dratrawdd(ircaga)*sc1a + ratraw(ircaga)*sc1add

      scfac(ircaga)     = sc1a  ; dscfacdt(ircaga)  = sc1adt  ; dscfacdd(ircaga)  = sc1add


! ar36(a,p)k39
      ratdum(irarap)    = ratraw(irarap) * sc1a
      dratdumdt(irarap) = dratrawdt(irarap)*sc1a + ratraw(irarap)*sc1adt
      dratdumdd(irarap) = dratrawdd(irarap)*sc1a + ratraw(irarap)*sc1add

      scfac(irarap)     = sc1a  ; dscfacdt(irarap)  = sc1adt  ; dscfacdd(irarap)  = sc1add

      ratdum(irkpa)     = ratraw(irkpa) * sc1a
      dratdumdt(irkpa)  = dratrawdt(irkpa)*sc1a  + ratraw(irkpa)*sc1adt
      dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + ratraw(irkpa)*sc1add

      scfac(irkpa)      = sc1a  ; dscfacdt(irkpa)   = sc1adt  ; dscfacdd(irkpa)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   19.0d0,39.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! k39(p,g)ca40
      ratdum(irkpg)     = ratraw(irkpg) * sc1a
      dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
      dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

      scfac(irkpg)      = sc1a  ; dscfacdt(irkpg)   = sc1adt  ; dscfacdd(irkpg)   = sc1add

      ratdum(ircagp)     = ratraw(ircagp) * sc1a
      dratdumdt(ircagp)  = dratrawdt(ircagp)*sc1a  + ratraw(ircagp)*sc1adt
      dratdumdd(ircagp)  = dratrawdd(ircagp)*sc1a  + ratraw(ircagp)*sc1add

      scfac(ircagp)      = sc1a  ; dscfacdt(ircagp)   = sc1adt  ; dscfacdd(ircagp)   = sc1add



! ca40 to ti44
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ica40),aion(ica40),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ca40(a,g)ti44
      ratdum(ircaag)    = ratraw(ircaag) * sc1a
      dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
      dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

      scfac(ircaag)     = sc1a  ; dscfacdt(ircaag)  = sc1adt  ; dscfacdd(ircaag)  = sc1add

      ratdum(irtiga)    = ratraw(irtiga) * sc1a
      dratdumdt(irtiga) = dratrawdt(irtiga)*sc1a + ratraw(irtiga)*sc1adt
      dratdumdd(irtiga) = dratrawdd(irtiga)*sc1a + ratraw(irtiga)*sc1add

      scfac(irtiga)     = sc1a  ; dscfacdt(irtiga)  = sc1adt  ; dscfacdd(irtiga)  = sc1add


! ca40(a,p)sc43
      ratdum(ircaap)    = ratraw(ircaap) * sc1a
      dratdumdt(ircaap) = dratrawdt(ircaap)*sc1a + ratraw(ircaap)*sc1adt
      dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + ratraw(ircaap)*sc1add

      scfac(ircaap)     = sc1a  ; dscfacdt(ircaap)  = sc1adt  ; dscfacdd(ircaap)  = sc1add

      ratdum(irscpa)    = ratraw(irscpa) * sc1a
      dratdumdt(irscpa) = dratrawdt(irscpa)*sc1a + ratraw(irscpa)*sc1adt
      dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + ratraw(irscpa)*sc1add

      scfac(irscpa)     = sc1a  ; dscfacdt(irscpa)  = sc1adt  ; dscfacdd(irscpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   21.0d0,43.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! sc43(p,g)ti44
      ratdum(irscpg)    = ratraw(irscpg) * sc1a
      dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
      dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

      scfac(irscpg)     = sc1a  ; dscfacdt(irscpg)  = sc1adt  ; dscfacdd(irscpg)  = sc1add

      ratdum(irtigp)    = ratraw(irtigp) * sc1a
      dratdumdt(irtigp) = dratrawdt(irtigp)*sc1a + ratraw(irtigp)*sc1adt
      dratdumdd(irtigp) = dratrawdd(irtigp)*sc1a + ratraw(irtigp)*sc1add

      scfac(irtigp)     = sc1a  ; dscfacdt(irtigp)  = sc1adt  ; dscfacdd(irtigp)  = sc1add



! ti44 to cr48
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iti44),aion(iti44),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! ti44(a,g)cr48
      ratdum(irtiag)    = ratraw(irtiag) * sc1a
      dratdumdt(irtiag) = dratrawdt(irtiag)*sc1a + ratraw(irtiag)*sc1adt
      dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + ratraw(irtiag)*sc1add

      scfac(irtiag)     = sc1a  ; dscfacdt(irtiag)  = sc1adt  ; dscfacdd(irtiag)  = sc1add

      ratdum(ircrga)    = ratraw(ircrga) * sc1a
      dratdumdt(ircrga) = dratrawdt(ircrga)*sc1a + ratraw(ircrga)*sc1adt
      dratdumdd(ircrga) = dratrawdd(ircrga)*sc1a + ratraw(ircrga)*sc1add

      scfac(ircrga)     = sc1a  ; dscfacdt(ircrga)  = sc1adt  ; dscfacdd(ircrga)  = sc1add

! ti44(a,p)v47
      ratdum(irtiap)    = ratraw(irtiap) * sc1a
      dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
      dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + ratraw(irtiap)*sc1add

      scfac(irtiap)  = sc1a  ; dscfacdt(irtiap)  = sc1adt  ; dscfacdd(irtiap)  = sc1add

      ratdum(irvpa)     = ratraw(irvpa) * sc1a
      dratdumdt(irvpa)  = dratrawdt(irvpa)*sc1a  + ratraw(irvpa)*sc1adt
      dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + ratraw(irvpa)*sc1add

      scfac(irvpa)      = sc1a  ; dscfacdt(irvpa)   = sc1adt  ; dscfacdd(irvpa)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   23.0d0,47.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! v47(p,g)cr48
      ratdum(irvpg)     = ratraw(irvpg) * sc1a
      dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
      dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

      scfac(irvpg)      = sc1a  ; dscfacdt(irvpg)   = sc1adt  ; dscfacdd(irvpg)   = sc1add

      ratdum(ircrgp)     = ratraw(ircrgp) * sc1a
      dratdumdt(ircrgp)  = dratrawdt(ircrgp)*sc1a  + ratraw(ircrgp)*sc1adt
      dratdumdd(ircrgp)  = dratrawdd(ircrgp)*sc1a  + ratraw(ircrgp)*sc1add

      scfac(ircrgp)      = sc1a  ; dscfacdt(ircrgp)   = sc1adt  ; dscfacdd(ircrgp)   = sc1add



! cr48 to fe52
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(icr48),aion(icr48),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! cr48(a,g)fe52
      ratdum(ircrag)    = ratraw(ircrag) * sc1a
      dratdumdt(ircrag) = dratrawdt(ircrag)*sc1a + ratraw(ircrag)*sc1adt
      dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + ratraw(ircrag)*sc1add

      scfac(ircrag)     = sc1a  ; dscfacdt(ircrag)  = sc1adt  ; dscfacdd(ircrag)  = sc1add

      ratdum(irfega)    = ratraw(irfega) * sc1a
      dratdumdt(irfega) = dratrawdt(irfega)*sc1a + ratraw(irfega)*sc1adt
      dratdumdd(irfega) = dratrawdd(irfega)*sc1a + ratraw(irfega)*sc1add

      scfac(irfega)     = sc1a  ; dscfacdt(irfega)  = sc1adt  ; dscfacdd(irfega)  = sc1add


! cr48(a,p)mn51
      ratdum(ircrap)    = ratraw(ircrap) * sc1a
      dratdumdt(ircrap) = dratrawdt(ircrap)*sc1a + ratraw(ircrap)*sc1adt
      dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + ratraw(ircrap)*sc1add

      scfac(ircrap)     = sc1a  ; dscfacdt(ircrap)  = sc1adt  ; dscfacdd(ircrap)  = sc1add

      ratdum(irmnpa)    = ratraw(irmnpa) * sc1a
      dratdumdt(irmnpa) = dratrawdt(irmnpa)*sc1a + ratraw(irmnpa)*sc1adt
      dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + ratraw(irmnpa)*sc1add

      scfac(irmnpa)     = sc1a  ; dscfacdt(irmnpa)  = sc1adt  ; dscfacdd(irmnpa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   25.0d0,51.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! mn51(p,g)fe52
      ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
      dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
      dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

      scfac(irmnpg)     = sc1a  ; dscfacdt(irmnpg)  = sc1adt  ; dscfacdd(irmnpg)  = sc1add

      ratdum(irfegp)    = ratraw(irfegp) * sc1a
      dratdumdt(irfegp) = dratrawdt(irfegp)*sc1a + ratraw(irfegp)*sc1adt
      dratdumdd(irfegp) = dratrawdd(irfegp)*sc1a + ratraw(irfegp)*sc1add

      scfac(irfegp)     = sc1a  ; dscfacdt(irfegp)  = sc1adt  ; dscfacdd(irfegp)  = sc1add


! fe to ni 
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife52),aion(ife52),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! fe52(a,g)ni56
      ratdum(irfeag)    = ratraw(irfeag) * sc1a
      dratdumdt(irfeag) = dratrawdt(irfeag)*sc1a + ratraw(irfeag)*sc1adt
      dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + ratraw(irfeag)*sc1add

      scfac(irfeag)     = sc1a  ; dscfacdt(irfeag)  = sc1adt  ; dscfacdd(irfeag)  = sc1add

      ratdum(irniga)    = ratraw(irniga) * sc1a
      dratdumdt(irniga) = dratrawdt(irniga)*sc1a + ratraw(irniga)*sc1adt
      dratdumdd(irniga) = dratrawdd(irniga)*sc1a + ratraw(irniga)*sc1add

      scfac(irniga)     = sc1a  ; dscfacdt(irniga)  = sc1adt  ; dscfacdd(irniga)  = sc1add


! fe52(a,p)co55
      ratdum(irfeap) = ratraw(irfeap) * sc1a
      dratdumdt(irfeap) = dratrawdt(irfeap)*sc1a + ratraw(irfeap)*sc1adt
      dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + ratraw(irfeap)*sc1add

      scfac(irfeap)     = sc1a  ; dscfacdt(irfeap)  = sc1adt  ; dscfacdd(irfeap)  = sc1add

      ratdum(ircopa)    = ratraw(ircopa) * sc1a
      dratdumdt(ircopa) = dratrawdt(ircopa)*sc1a + ratraw(ircopa)*sc1adt
      dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + ratraw(ircopa)*sc1add

      scfac(ircopa)     = sc1a  ; dscfacdt(ircopa)  = sc1adt  ; dscfacdd(ircopa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   27.0d0,55.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)


! co55(p,g)ni56
      ratdum(ircopg)    = ratraw(ircopg) * sc1a
      dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
      dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

      scfac(ircopg)     = sc1a  ; dscfacdt(ircopg)  = sc1adt  ; dscfacdd(ircopg)  = sc1add

      ratdum(irnigp)    = ratraw(irnigp) * sc1a
      dratdumdt(irnigp) = dratrawdt(irnigp)*sc1a + ratraw(irnigp)*sc1adt
      dratdumdd(irnigp) = dratrawdd(irnigp)*sc1a + ratraw(irnigp)*sc1add

      scfac(irnigp)     = sc1a  ; dscfacdt(irnigp)  = sc1adt  ; dscfacdd(irniga)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife54),aion(ife54),1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! fe54(p,g)co55
      ratdum(irfepg)    = ratraw(irfepg) * sc1a
      dratdumdt(irfepg) = dratrawdt(irfepg)*sc1a + ratraw(irfepg)*sc1adt
      dratdumdd(irfepg) = dratrawdd(irfepg)*sc1a + ratraw(irfepg)*sc1add

      scfac(irfepg)     = sc1a  ; dscfacdt(irfepg)  = sc1adt  ; dscfacdd(irfepg)  = sc1add

      ratdum(ircogp)    = ratraw(ircogp) * sc1a
      dratdumdt(ircogp) = dratrawdt(ircogp)*sc1a + ratraw(ircogp)*sc1adt
      dratdumdd(ircogp) = dratrawdd(ircogp)*sc1a + ratraw(ircogp)*sc1add

      scfac(ircogp)     = sc1a  ; dscfacdt(ircogp)  = sc1adt  ; dscfacdd(ircogp)  = sc1add



! fe54(a,p)co57(g,p)fe56
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife54),aion(ife54),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)


! fe54(a,p)co57
      ratdum(irfe54ap)    = ratraw(irfe54ap) * sc1a
      dratdumdt(irfe54ap) = dratrawdt(irfe54ap)*sc1a + ratraw(irfe54ap)*sc1adt
      dratdumdd(irfe54ap) = dratrawdd(irfe54ap)*sc1a + ratraw(irfe54ap)*sc1add

      scfac(irfe54ap)     = sc1a  ; dscfacdt(irfe54ap)  = sc1adt  ; dscfacdd(irfe54ap)  = sc1add

      ratdum(irco57pa)    = ratraw(irco57pa) * sc1a
      dratdumdt(irco57pa) = dratrawdt(irco57pa)*sc1a + ratraw(irco57pa)*sc1adt
      dratdumdd(irco57pa) = dratrawdd(irco57pa)*sc1a + ratraw(irco57pa)*sc1add

      scfac(irco57pa)     = sc1a  ; dscfacdt(irco57pa)  = sc1adt  ; dscfacdd(irco57pa)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife56),aion(ife56),1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add) 

! fe56(p,g)co57
      ratdum(irfe56pg)    = ratraw(irfe56pg) * sc1a
      dratdumdt(irfe56pg) = dratrawdt(irfe56pg)*sc1a + ratraw(irfe56pg)*sc1adt
      dratdumdd(irfe56pg) = dratrawdd(irfe56pg)*sc1a + ratraw(irfe56pg)*sc1add

      scfac(irfe56pg)     = sc1a  ; dscfacdt(irfe56pg)  = sc1adt  ; dscfacdd(irfe56pg)  = sc1add

      ratdum(irco57gp)    = ratraw(irco57gp) * sc1a
      dratdumdt(irco57gp) = dratrawdt(irco57gp)*sc1a + ratraw(irco57gp)*sc1adt
      dratdumdd(irco57gp) = dratrawdd(irco57gp)*sc1a + ratraw(irco57gp)*sc1add

      scfac(irco57gp)     = sc1a  ; dscfacdt(irco57gp)  = sc1adt  ; dscfacdd(irco57gp)  = sc1add



! d(p,g)he4
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   1.0d0,2.0d0,zion(ih1),aion(ih1), &
                   jscr,init,sc1a,sc1adt,sc1add)


! d(p,g)he3
      ratdum(irdpg)    = ratraw(irdpg) * sc1a
      dratdumdt(irdpg) = dratrawdt(irdpg)*sc1a + ratraw(irdpg)*sc1adt
      dratdumdd(irdpg) = dratrawdd(irdpg)*sc1a + ratraw(irdpg)*sc1add

      scfac(irdpg)     = sc1a  ; dscfacdt(irdpg)  = sc1adt ; dscfacdd(irdpg) = sc1add

      ratdum(irhegp)    = ratraw(irhegp) * sc1a
      dratdumdt(irhegp) = dratrawdt(irhegp)*sc1a + ratraw(irhegp)*sc1adt
      dratdumdd(irhegp) = dratrawdd(irhegp)*sc1a + ratraw(irhegp)*sc1add

      scfac(irhegp)     = sc1a  ; dscfacdt(irhegp)  = sc1adt ; dscfacdd(irhegp) = sc1add


! pp
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ih1),aion(ih1),zion(ih1),aion(ih1), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irpp)    = ratraw(irpp) * sc1a
      dratdumdt(irpp) = dratrawdt(irpp)*sc1a + ratraw(irpp)*sc1adt
      dratdumdd(irpp) = dratrawdd(irpp)*sc1a + ratraw(irpp)*sc1add

      scfac(irpp)     = sc1a  ; dscfacdt(irpp)  = sc1adt  ; dscfacdd(irpp)  = sc1add


! he3 + he3
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3), &
                   jscr,init,sc1a,sc1adt,sc1add)

! he3(he3,2p)he4
      ratdum(ir33)    = ratraw(ir33) * sc1a
      dratdumdt(ir33) = dratrawdt(ir33)*sc1a + ratraw(ir33)*sc1adt
      dratdumdd(ir33) = dratrawdd(ir33)*sc1a + ratraw(ir33)*sc1add

      scfac(ir33)     = sc1a  ; dscfacdt(ir33)  = sc1adt  ; dscfacdd(ir33)  = sc1add


! he3 + he4
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

! he3(he4,g)be7
      ratdum(irhe3ag)    = ratraw(irhe3ag) * sc1a
      dratdumdt(irhe3ag) = dratrawdt(irhe3ag)*sc1a &
                           + ratraw(irhe3ag)*sc1adt
      dratdumdd(irhe3ag) = dratrawdd(irhe3ag)*sc1a &
                           + ratraw(irhe3ag)*sc1add

      scfac(irhe3ag)     = sc1a  ; dscfacdt(irhe3ag)  = sc1adt  ; dscfacdd(irhe3ag)  = sc1add



! cno cycles

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ih1),aion(ih1), &
                   jscr,init,sc1a,sc1adt,sc1add)

! c12(p,g)13n
      ratdum(ircpg)    = ratraw(ircpg) * sc1a
      dratdumdt(ircpg) = dratrawdt(ircpg)*sc1a + ratraw(ircpg)*sc1adt
      dratdumdd(ircpg) = dratrawdd(ircpg)*sc1a + ratraw(ircpg)*sc1add

      scfac(ircpg)     = sc1a  ; dscfacdt(ircpg)  = sc1adt  ; dscfacdd(ircpg)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(in14),aion(in14),zion(ih1),aion(ih1), &
                   jscr,init,sc1a,sc1adt,sc1add)

! n14(p,g)o15
      ratdum(irnpg)    = ratraw(irnpg) * sc1a
      dratdumdt(irnpg) = dratrawdt(irnpg)*sc1a + ratraw(irnpg)*sc1adt
      dratdumdd(irnpg) = dratrawdd(irnpg)*sc1a + ratraw(irnpg)*sc1add

      scfac(irnpg)     = sc1a ; dscfacdt(irnpg)  = sc1adt  ; dscfacdd(irnpg)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(ih1),aion(ih1), &
                   jscr,init,sc1a,sc1adt,sc1add)


! o16(p,g)f17
      ratdum(iropg)    = ratraw(iropg) * sc1a
      dratdumdt(iropg) = dratrawdt(iropg)*sc1a+ratraw(iropg)*sc1adt
      dratdumdd(iropg) = dratrawdd(iropg)*sc1a+ratraw(iropg)*sc1add

      scfac(iropg)     = sc1a  ; dscfacdt(iropg)  = sc1adt  ; dscfacdd(iropg)  = sc1add

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(in14),aion(in14),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

! n14(a,g)f18
      ratdum(irnag)    = ratraw(irnag) * sc1a
      dratdumdt(irnag) = dratrawdt(irnag)*sc1a + ratraw(irnag)*sc1adt
      dratdumdd(irnag) = dratrawdd(irnag)*sc1a + ratraw(irnag)*sc1add

      scfac(irnag)     = sc1a  ; dscfacdt(irnag)  = sc1adt  ; dscfacdd(irnag)  = sc1add


! reset the screen initialization flag
      init = 0

! end of screening if block
      end if




! now form those lovely equilibrium rates

! mg24(a,p)27al(p,g)28si
      ratdum(irr1)     = 0.0d0  ; dratdumdt(irr1)  = 0.0d0  ; dratdumdd(irr1)  = 0.0d0
      denom    = ratdum(iralpa) + ratdum(iralpg)
      denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
      denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
      if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irr1)    = ratdum(iralpa)*zz
       dratdumdt(irr1) = (dratdumdt(iralpa) - ratdum(irr1)*denomdt)*zz
       dratdumdd(irr1) = (dratdumdd(iralpa) - ratdum(irr1)*denomdd)*zz
      end if

! si28(a,p)p31(p,g)s32
      ratdum(irs1)     = 0.0d0  ; dratdumdt(irs1)  = 0.0d0  ; dratdumdd(irs1)  = 0.0d0
      denom    = ratdum(irppa) + ratdum(irppg)
      denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
      denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
      if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irs1)    = ratdum(irppa)*zz
       dratdumdt(irs1) = (dratdumdt(irppa) - ratdum(irs1)*denomdt)*zz
       dratdumdd(irs1) = (dratdumdd(irppa) - ratdum(irs1)*denomdd)*zz
      end if

! s32(a,p)cl35(p,g)ar36
      ratdum(irt1)     = 0.0d0  ; dratdumdt(irt1)  = 0.0d0  ; dratdumdd(irt1)  = 0.0d0
      denom    = ratdum(irclpa) + ratdum(irclpg)
      denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
      denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
      if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irt1)    = ratdum(irclpa)*zz
       dratdumdt(irt1) = (dratdumdt(irclpa) - ratdum(irt1)*denomdt)*zz
       dratdumdd(irt1) = (dratdumdd(irclpa) - ratdum(irt1)*denomdd)*zz
      end if

! ar36(a,p)k39(p,g)ca40
      ratdum(iru1)     = 0.0d0  ; dratdumdt(iru1)  = 0.0d0  ; dratdumdd(iru1)  = 0.0d0
      denom    = ratdum(irkpa) + ratdum(irkpg)
      denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
      denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
      if (denom .gt. 1.0d-50) then
       zz   = 1.0d0/denom
       ratdum(iru1)   = ratdum(irkpa)*zz
       dratdumdt(iru1) = (dratdumdt(irkpa) - ratdum(iru1)*denomdt)*zz
       dratdumdd(iru1) = (dratdumdd(irkpa) - ratdum(iru1)*denomdd)*zz
      end if

! ca40(a,p)sc43(p,g)ti44
      ratdum(irv1)     = 0.0d0  ; dratdumdt(irv1)  = 0.0d0  ; dratdumdd(irv1)  = 0.0d0
      denom    = ratdum(irscpa) + ratdum(irscpg)
      denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
      denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
      if (denom .gt. 1.0d-50) then
       zz  = 1.0d0/denom
       ratdum(irv1)    = ratdum(irscpa)*zz
       dratdumdt(irv1) = (dratdumdt(irscpa) - ratdum(irv1)*denomdt)*zz
       dratdumdd(irv1) = (dratdumdd(irscpa) - ratdum(irv1)*denomdd)*zz
      end if

! ti44(a,p)v47(p,g)cr48
      ratdum(irw1)    = 0.0d0  ; dratdumdt(irw1) = 0.0d0  ; dratdumdd(irw1) = 0.0d0
      denom    = ratdum(irvpa) + ratdum(irvpg)
      denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
      denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
      if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irw1)    = ratdum(irvpa)*zz
       dratdumdt(irw1) = (dratdumdt(irvpa) - ratdum(irw1)*denomdt)*zz
       dratdumdd(irw1) = (dratdumdd(irvpa) - ratdum(irw1)*denomdd)*zz
      end if

! cr48(a,p)mn51(p,g)fe52
      ratdum(irx1)    = 0.0d0  ; dratdumdt(irx1) = 0.0d0  ; dratdumdd(irx1) = 0.0d0
      denom    = ratdum(irmnpa) + ratdum(irmnpg)
      denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
      denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
      if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
      endif


! fe52(n,g)fe53(n,g)fe54 equilibrium links

      ratdum(ir1f54) = 0.0d0  ; dratdumdy1(ir1f54) = 0.0d0  ; dratdumdt(ir1f54)  = 0.0d0  ; dratdumdd(ir1f54)  = 0.0d0
      ratdum(ir2f54) = 0.0d0  ; dratdumdy1(ir2f54) = 0.0d0  ; dratdumdt(ir2f54)  = 0.0d0  ; dratdumdd(ir2f54)  = 0.0d0

      denom   = ratdum(ir53gn) + y(ineut)*ratdum(ir53ng)
      denomdt = dratdumdt(ir53gn) + y(ineut)*dratdumdt(ir53ng)
      denomdd = dratdumdd(ir53gn) + y(ineut)*dratdumdd(ir53ng)

      if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then
       zz      = 1.0d0/denom

       ratdum(ir1f54)     = ratdum(ir54gn)*ratdum(ir53gn)*zz
       dratdumdy1(ir1f54) = -ratdum(ir1f54)*zz * ratdum(ir53ng)
       dratdumdt(ir1f54)  = dratdumdt(ir54gn)*ratdum(ir53gn)*zz + ratdum(ir54gn)*dratdumdt(ir53gn)*zz - ratdum(ir1f54)*zz*denomdt
       dratdumdd(ir1f54)  = dratdumdd(ir54gn)*ratdum(ir53gn)*zz + ratdum(ir54gn)*dratdumdd(ir53gn)*zz - ratdum(ir1f54)*zz*denomdd

       ratdum(ir2f54)     = ratdum(ir52ng)*ratdum(ir53ng)*zz
       dratdumdy1(ir2f54) = -ratdum(ir2f54)*zz * ratdum(ir53ng)
       dratdumdt(ir2f54)  = dratdumdt(ir52ng)*ratdum(ir53ng)*zz + ratdum(ir52ng)*dratdumdt(ir53ng)*zz - ratdum(ir2f54)*zz*denomdt
       dratdumdd(ir2f54)  = dratdumdd(ir52ng)*ratdum(ir53ng)*zz + ratdum(ir52ng)*dratdumdd(ir53ng)*zz - ratdum(ir2f54)*zz*denomdd
      end if


! fe54(n,g)fe55(n,g)fe56 equilibrium links

      ratdum(irfe56_aux1) = 0.0d0  ; dratdumdy1(irfe56_aux1) = 0.0d0
      dratdumdt(irfe56_aux1)  = 0.0d0  ; dratdumdd(irfe56_aux1)  = 0.0d0
      ratdum(irfe56_aux2) = 0.0d0  ; dratdumdy1(irfe56_aux2) = 0.0d0
      dratdumdt(irfe56_aux2)  = 0.0d0  ; dratdumdd(irfe56_aux2)  = 0.0d0

      denom   = ratdum(ir55gn)    + y(ineut)*ratdum(ir55ng)
      denomdt = dratdumdt(ir55gn) + y(ineut)*dratdumdt(ir55ng)
      denomdd = dratdumdd(ir55gn) + y(ineut)*dratdumdd(ir55ng)

      if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then
       zz      = 1.0d0/denom

       ratdum(irfe56_aux1)     = ratdum(ir56gn)*ratdum(ir55gn)*zz
       dratdumdy1(irfe56_aux1) = -ratdum(irfe56_aux1)*zz * ratdum(ir55ng)
       dratdumdt(irfe56_aux1)  = dratdumdt(ir56gn)*ratdum(ir55gn)*zz &
        + ratdum(ir56gn)*dratdumdt(ir55gn)*zz - &
       ratdum(irfe56_aux1)*zz*denomdt
       dratdumdd(irfe56_aux1)  = dratdumdd(ir56gn)*ratdum(ir55gn)*zz &
       + ratdum(ir56gn)*dratdumdd(ir55gn)*zz - &
       ratdum(irfe56_aux1)*zz*denomdd

       ratdum(irfe56_aux2)     = ratdum(ir54ng)*ratdum(ir55ng)*zz
       dratdumdy1(irfe56_aux2) = -ratdum(irfe56_aux2)*zz * ratdum(ir55ng)
       dratdumdt(irfe56_aux2)  = dratdumdt(ir54ng)*ratdum(ir55ng)*zz &
       + ratdum(ir54ng)*dratdumdt(ir55ng)*zz - &
       ratdum(irfe56_aux2)*zz*denomdt
       dratdumdd(irfe56_aux2)  = dratdumdd(ir54ng)*ratdum(ir55ng)*zz &
       + ratdum(ir54ng)*dratdumdd(ir55ng)*zz  - &
       ratdum(irfe56_aux2)*zz*denomdd
      end if


! fe54(a,p)co57(g,p)fe56 equilibrium links 

      ratdum(irfe56_aux3)     = 0.0d0
      dratdumdy1(irfe56_aux3) = 0.0d0
      dratdumdt(irfe56_aux3)  = 0.0d0
      dratdumdd(irfe56_aux3)  = 0.0d0
      ratdum(irfe56_aux4)     = 0.0d0
      dratdumdy1(irfe56_aux4) = 0.0d0
      dratdumdt(irfe56_aux4)  = 0.0d0
      dratdumdd(irfe56_aux4)  = 0.0d0

      denom   = ratdum(irco57gp)    + y(iprot)*ratdum(irco57pa)
      denomdt = dratdumdt(irco57gp) + y(iprot)*dratdumdt(irco57pa)
      denomdd = dratdumdd(irco57gp) + y(iprot)*dratdumdd(irco57pa)

      if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then
       zz      = 1.0d0/denom

       ratdum(irfe56_aux3)     = ratdum(irfe56pg) * ratdum(irco57pa) * zz
       dratdumdy1(irfe56_aux3) = -ratdum(irfe56_aux3) * zz * ratdum(irco57pa)
       dratdumdt(irfe56_aux3)  = dratdumdt(irfe56pg) * ratdum(irco57pa) * zz &
       + ratdum(irfe56pg) * dratdumdt(irco57pa) * zz - &
       ratdum(irfe56_aux3) * zz * denomdt
       dratdumdd(irfe56_aux3)  = dratdumdd(irfe56pg) * ratdum(irco57pa) * zz &
       + ratdum(irfe56pg) * dratdumdd(irco57pa) * zz - &
       ratdum(irfe56_aux3) * zz * denomdd

       ratdum(irfe56_aux4)     = ratdum(irfe54ap) * ratdum(irco57gp) * zz
       dratdumdy1(irfe56_aux4) = -ratdum(irfe56_aux4) * zz * ratdum(irco57pa)
       dratdumdt(irfe56_aux4)  = dratdumdt(irfe54ap) * ratdum(irco57gp) * zz &
       + ratdum(irfe54ap) * dratdumdt(irco57gp) * zz - &
       ratdum(irfe56_aux4) * zz * denomdt
       dratdumdd(irfe56_aux4)  = dratdumdd(irfe54ap) * ratdum(irco57gp) * zz &
       + ratdum(irfe54ap) * dratdumdd(irco57gp) * zz  - &
       ratdum(irfe56_aux4) * zz * denomdd
      end if


! fe54(p,g)co55(p,g)ni56 equilibrium links r3f54 r4f54
! fe52(a,p)co55(g,p)fe54 equilibrium links r5f54 r6f54
! fe52(a,p)co55(p,g)ni56 equilibrium links r7f54 r8f54

      ratdum(ir3f54) = 0.0d0  ; dratdumdy1(ir3f54) = 0.0d0
      dratdumdt(ir3f54) = 0.0d0  ; dratdumdd(ir3f54) = 0.0d0
      ratdum(ir4f54) = 0.0d0  ; dratdumdy1(ir4f54) = 0.0d0
      dratdumdt(ir4f54) = 0.0d0  ; dratdumdd(ir4f54) = 0.0d0
      ratdum(ir5f54) = 0.0d0  ; dratdumdy1(ir5f54) = 0.0d0
      dratdumdt(ir5f54) = 0.0d0  ; dratdumdd(ir5f54) = 0.0d0
      ratdum(ir6f54) = 0.0d0  ; dratdumdy1(ir6f54) = 0.0d0
      dratdumdt(ir6f54) = 0.0d0  ; dratdumdd(ir6f54) = 0.0d0
      ratdum(ir7f54) = 0.0d0  ; dratdumdy1(ir7f54) = 0.0d0
      dratdumdt(ir7f54) = 0.0d0  ; dratdumdd(ir7f54) = 0.0d0
      ratdum(ir8f54) = 0.0d0  ; dratdumdy1(ir8f54) = 0.0d0
      dratdumdt(ir8f54) = 0.0d0  ; dratdumdd(ir8f54) = 0.0d0

      denom   = ratdum(ircogp) + y(iprot)*(ratdum(ircopg) + ratdum(ircopa))

      if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then

       denomdt = dratdumdt(ircogp) + y(iprot)*(dratdumdt(ircopg) + dratdumdt(ircopa))
       denomdd = dratdumdd(ircogp) + y(iprot)*(dratdumdd(ircopg) + dratdumdd(ircopa))

       zz      = 1.0d0/denom

       ratdum(ir3f54)     = ratdum(irfepg) * ratdum(ircopg) * zz
       dratdumdy1(ir3f54) = -ratdum(ir3f54) * zz * (ratdum(ircopg) + ratdum(ircopa))
       dratdumdt(ir3f54)  = dratdumdt(irfepg) * ratdum(ircopg) * zz &
       + ratdum(irfepg) * dratdumdt(ircopg) * zz - ratdum(ir3f54)*zz*denomdt
       dratdumdd(ir3f54)  = dratdumdd(irfepg) * ratdum(ircopg) * zz &
       + ratdum(irfepg) * dratdumdd(ircopg) * zz - ratdum(ir3f54)*zz*denomdd

       ratdum(ir4f54)     = ratdum(irnigp) * ratdum(ircogp) * zz
       dratdumdy1(ir4f54) = -ratdum(ir4f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir4f54)  =  dratdumdt(irnigp) * ratdum(ircogp) * zz &
       + ratdum(irnigp) * dratdumdt(ircogp) * zz - ratdum(ir4f54)*zz*denomdt
       dratdumdd(ir4f54)  = dratdumdd(irnigp) * ratdum(ircogp) * zz &
       + ratdum(irnigp) * dratdumdd(ircogp) * zz  - ratdum(ir4f54)*zz*denomdd

       ratdum(ir5f54)     = ratdum(irfepg) * ratdum(ircopa) * zz
       dratdumdy1(ir5f54) = -ratdum(ir5f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir5f54)  = dratdumdt(irfepg) * ratdum(ircopa) * zz &
       + ratdum(irfepg) * dratdumdt(ircopa) * zz - ratdum(ir5f54) * zz * denomdt
       dratdumdd(ir5f54)  = dratdumdd(irfepg) * ratdum(ircopa) * zz &
       + ratdum(irfepg) * dratdumdd(ircopa) * zz - ratdum(ir5f54) * zz * denomdd

       ratdum(ir6f54)     = ratdum(irfeap) * ratdum(ircogp) * zz
       dratdumdy1(ir6f54) = -ratdum(ir6f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir6f54)  = dratdumdt(irfeap) * ratdum(ircogp) * zz &
       + ratdum(irfeap) * dratdumdt(ircogp) * zz - ratdum(ir6f54) * zz * denomdt
       dratdumdd(ir6f54)  = dratdumdd(irfeap) * ratdum(ircogp) * zz &
       + ratdum(irfeap) * dratdumdd(ircogp) * zz - ratdum(ir6f54) * zz * denomdd

       ratdum(ir7f54)     = ratdum(irfeap) * ratdum(ircopg) * zz
       dratdumdy1(ir7f54) = -ratdum(ir7f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir7f54)  = dratdumdt(irfeap) * ratdum(ircopg) * zz &
       + ratdum(irfeap) * dratdumdt(ircopg) * zz - ratdum(ir7f54) * zz * denomdt
       dratdumdd(ir7f54)  = dratdumdd(irfeap) * ratdum(ircopg) * zz &
       + ratdum(irfeap) * dratdumdd(ircopg) * zz - ratdum(ir7f54) * zz * denomdd

       ratdum(ir8f54)     = ratdum(irnigp) * ratdum(ircopa) * zz
       dratdumdy1(ir8f54) = -ratdum(ir8f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir8f54)  = dratdumdt(irnigp) * ratdum(ircopa) * zz &
       + ratdum(irnigp) * dratdumdt(ircopa) * zz - ratdum(ir8f54) * zz * denomdt
       dratdumdd(ir8f54)  = dratdumdd(irnigp) * ratdum(ircopa) * zz &
       + ratdum(irnigp) * dratdumdd(ircopa) * zz - ratdum(ir8f54) * zz * denomdd
      end if


! p(n,g)h2(n,g)3h(p,g)he4   photodisintegrated n and p back to he4 equilibrium links
! p(n,g)h2(p,g)he3(n,g)he4

      ratdum(iralf1)     = 0.0d0  ; dratdumdy1(iralf1) = 0.0d0
      dratdumdy2(iralf1) = 0.0d0  ; dratdumdt(iralf1)  = 0.0d0
      dratdumdd(iralf1)  = 0.0d0  
      ratdum(iralf2)     = 0.0d0  ; dratdumdy1(iralf2) = 0.0d0
      dratdumdy2(iralf2) = 0.0d0  ; dratdumdt(iralf2)  = 0.0d0
      dratdumdd(iralf2)  = 0.0d0

      denom  = ratdum(irhegp)*ratdum(irdgn) &
      + y(ineut)*ratdum(irheng)*ratdum(irdgn) &
      + y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)

      if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then

       denomdt  = dratdumdt(irhegp)*ratdum(irdgn) + ratdum(irhegp)*dratdumdt(irdgn) & 
                +  y(ineut) * (dratdumdt(irheng)*ratdum(irdgn) + ratdum(irheng)*dratdumdt(irdgn)) &
                +  y(ineut)*y(iprot) * (dratdumdt(irheng)*ratdum(irdpg) + ratdum(irheng)*dratdumdt(irdpg))

       denomdd  = dratdumdd(irhegp)*ratdum(irdgn) + ratdum(irhegp)*dratdumdd(irdgn) &
                +  y(ineut) * (dratdumdd(irheng)*ratdum(irdgn) + ratdum(irheng)*dratdumdd(irdgn)) &
                +  y(ineut)*y(iprot) * (dratdumdd(irheng)*ratdum(irdpg) + ratdum(irheng)*dratdumdd(irdpg))

       zz = 1.0d0/denom

       ratdum(iralf1)     = ratdum(irhegn) * ratdum(irhegp) * ratdum(irdgn) * zz
       dratdumdy1(iralf1) = -ratdum(iralf1) * zz * (ratdum(irheng)*ratdum(irdgn) &
       + y(iprot)*ratdum(irheng)*ratdum(irdpg))
       dratdumdy2(iralf1) = -ratdum(iralf1) * zz * y(ineut) * ratdum(irheng) * ratdum(irdpg) 
       dratdumdt(iralf1)  = dratdumdt(irhegn)*ratdum(irhegp) * ratdum(irdgn) * zz &
                 + ratdum(irhegn)*dratdumdt(irhegp)*ratdum(irdgn)*zz &
                 + ratdum(irhegn)*ratdum(irhegp)*dratdumdt(irdgn)*zz &
                 - ratdum(iralf1)*zz*denomdt
       dratdumdd(iralf1)  = dratdumdd(irhegn) * ratdum(irhegp) * ratdum(irdgn) * zz &
                 + ratdum(irhegn)*dratdumdd(irhegp)*ratdum(irdgn)*zz &
                 + ratdum(irhegn)*ratdum(irhegp)*dratdumdd(irdgn)*zz &
                 - ratdum(iralf1)*zz*denomdt


       ratdum(iralf2)     = ratdum(irheng)*ratdum(irdpg) * ratdum(irhng)*zz
       dratdumdy1(iralf2) = -ratdum(iralf2) * zz * (ratdum(irheng)*ratdum(irdgn) &
       + y(iprot)*ratdum(irheng)*ratdum(irdpg))
       dratdumdy2(iralf2) = -ratdum(iralf2) * zz * y(ineut) * ratdum(irheng) * ratdum(irdpg)
       dratdumdt(iralf2)  = dratdumdt(irheng)*ratdum(irdpg) * ratdum(irhng) * zz &
                 + ratdum(irheng)*dratdumdt(irdpg)*ratdum(irhng)*zz &
                 + ratdum(irheng)*ratdum(irdpg)*dratdumdt(irhng)*zz &
                 - ratdum(iralf2)*zz*denomdt
       dratdumdd(iralf2)  = dratdumdd(irheng)*ratdum(irdpg) * ratdum(irhng)*zz &
                 + ratdum(irheng)*dratdumdd(irdpg)*ratdum(irhng)*zz &
                 + ratdum(irheng)*ratdum(irdpg)*dratdumdd(irhng)*zz &
                 - ratdum(iralf2)*zz*denomdd
      end if



! he3(a,g)be7(p,g)8b(e+nu)8be(2a)
! beta limit he3+he4 by the 8B decay half life

      if (y(ihe4) .gt. 1.0d-30) then
       xx            = 0.896d0/y(ihe4)
       ratdum(irhe3ag)  = min(ratdum(irhe3ag),xx)
       if (ratdum(irhe3ag) .eq. xx) then
        dratdumdy1(irhe3ag) = -xx/y(ihe4)
        dratdumdt(irhe3ag)  = 0.0d0
        dratdumdd(irhe3ag)  = 0.0d0
       else
        dratdumdy1(irhe3ag) = 0.0d0
       endif
      endif


! beta limit n14(p,g)o15(enu)o16  and o16(p,g)f17(e+nu)17o(p,a)n14

      if (y(ih1)  .gt. 1.0d-30) then
       xx = 5.68d-03/(y(ih1)*1.57d0)
       ratdum(irnpg) = min(ratdum(irnpg),xx)
       if (ratdum(irnpg) .eq. xx) then
        dratdumdy1(irnpg) = -xx/y(ih1)
        dratdumdt(irnpg)  = 0.0d0
        dratdumdd(irnpg)  = 0.0d0
       else
        dratdumdy1(irnpg) = 0.0d0
       end if

       xx = 0.0105d0/y(ih1)
       ratdum(iropg) = min(ratdum(iropg),xx)
       if (ratdum(iropg) .eq. xx) then
        dratdumdy1(iropg) = -xx/y(ih1)
        dratdumdt(iropg)  = 0.0d0
        dratdumdd(iropg)  = 0.0d0
       else
        dratdumdy1(iropg) = 0.0d0
       end if
      end if




! write the screened rates
      if (write_rates .eq. 1) then
       maxrate_loc = maxloc(ratraw(1:nrat))
       minrate_loc = minloc(ratraw(1:nrat))
       open(unit=22,file='rate_chek.dat',status='unknown')
       write(22,'(1x,1p2e14.6,i6)') btemp,bden,nrat
       write(22,'(1x,a,a,a,1pe10.2)') 'max rate: ',trim(ratnam(maxrate_loc(1))),' ',ratraw(maxrate_loc)
       write(22,'(1x,a,a,a,1pe10.2)') 'min rate: ',trim(ratnam(minrate_loc(1))),' ',ratraw(minrate_loc)
       do i=1,nrat
        write(22,'(1x,i4,a,a,a,t20,1p8e10.2)') i,' ',trim(ratnam(i)),' ',ratraw(i),scfac(i),ratdum(i)
       enddo
       close(unit=22)
       write(6,*) 
       write(6,*) 'wrote rates to rate_chek.dat'
       write(6,*) 
       write_rates = 0
      end if

      return
    end subroutine screen_aprox21

      subroutine screen5(temp,den,zbar,abar,z2bar, &
                         z1,a1,z2,a2,jscreen,init, &
                         scor,scordt,scordd)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'


! this subroutine calculates screening factors and their derivatives
! for nuclear reaction rates in the weak, intermediate and strong regimes.
! based on dewitt, graboske, cooper apj, 181, 439, 1973 and 
! graboske, dewit, grossman and cooper apj 181 457 1973 for weak screening. 
!
! based on alastuey and jancovici apj 226 1034 1978 
! and jancovici jsp 17, 5. 357 1977 with plasma parameters 
! from itoh et al apj 234 1079 1979, for strong screening.


! input:
! temp    = temperature
! den     = density
! zbar    = mean charge per nucleus
! abar    = mean number of nucleons per nucleus
! z2bar   = mean square charge per nucleus
! z1 a1   = charge and number in the entrance channel
! z2 a2   = charge and number in the exit channel
! jscreen = counter of which reaction is being calculated
! init    = flag to compute the more expensive functions just once

! output:
! scor    = screening correction
! scordt  = derivative of screening correction with temperature
! scordd  = derivative of screening correction with density


! declare the pass
      integer          jscreen,init
      double precision temp,den,zbar,abar,z2bar,z1,a1,z2,a2, &
                       scor,scordt,scordd


! local variables
      integer          iscreen_type
      double precision aa,daadt,daadd,bb,cc,dccdt,dccdd, &
                       pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd, &
                       ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd, &
                       vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni, &
                       qlam0z,qlam0zdt,qlam0zdd, &
                       h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd, &
                       h12x,dh12xdt,dh12xdd,alfa,beta, x1, x2, denom, &
                       taufac,taufacdt,gamp,gampdt,gampdd, &
                       gamef,gamefdt,gamefdd, &
                       tau12,tau12dt,alph12,alph12dt,alph12dd, &
                       xlgfac,dxlgfacdt,dxlgfacdd, &
                       gamp14,gamp14dt,gamp14dd, &
                       xni,dxnidd,ytot, &
                       temp_old,den_old,zbar_old,abar_old


! screening variables
! zs13    = (z1+z2)**(1./3.)
! zhat    = combination of z1 and z2 raised to the 5/3 power
! zhat2   = combination of z1 and z2 raised to the 5/12 power
! lzav    = log of effective charge
! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power


      integer          nscreen_max
      parameter        (nscreen_max = 2*abignet + 40)

      double precision zs13(nscreen_max),zhat(nscreen_max), &
                       zhat2(nscreen_max),lzav(nscreen_max), &
                       aznut(nscreen_max),zs13inv(nscreen_max)


! parameter fact is the cube root of 2
      double precision  x13,x14,x53,x532,x512,fact,co2,gamefx,gamefs, &
                        dgamma, blend_frac
      parameter        (x13    = 1.0d0/3.0d0, &
                        x14    = 1.0d0/4.0d0, &
                        x53    = 5.0d0/3.0d0, &
                        x532   = 5.0d0/32.0d0, &
                        x512   = 5.0d0/12.0d0, &
                        fact   = 1.25992104989487d0, &
                        co2    = x13 * 4.248719d3, &
                        gamefx = 0.1d0, &
                        gamefs = 0.4d0, &
                        dgamma     = 1.0d0/(gamefs - gamefx), &
                        blend_frac = 0.05d0)


      data     temp_old/-1.0d0/, den_old/-1.0d0/, &
               zbar_old/-1.0d0/, abar_old/-1.0d0/


! compute and store the more expensive screening factors
      if (init .eq. 1) then
       if (jscreen .gt. nscreen_max) &
       stop 'jscreen > nscreen_max in screen5'
       zs13(jscreen)    = (z1 + z2)**x13
       zs13inv(jscreen) = 1.0d0/zs13(jscreen)
       zhat(jscreen)    = (z1 + z2)**x53  - z1**x53 - z2**x53
       zhat2(jscreen)   = (z1 + z2)**x512 - z1**x512 -z2**x512
       lzav(jscreen)    = x53 * log(z1*z2/(z1 + z2))
       aznut(jscreen)   = (z1**2 * z2**2 * a1*a2 / (a1 + a2))**x13
      endif


! calculate average plasma, if need be
      if (temp_old .ne. temp .or. &
          den_old  .ne. den  .or. &
          zbar_old  .ne. zbar  .or. &
          abar_old  .ne. abar ) then

       temp_old = temp
       den_old  = den
       zbar_old  = zbar
       abar_old  = abar

       ytot     = 1.0d0/abar
       rr       = den * ytot
       tempi   = 1.0d0/temp
       dtempi  = -tempi*tempi
       deni    = 1.0d0/den

! equation 4 of dewitt et al becomes ztilde = z2bar + zbar
       pp       = sqrt(rr*tempi*(z2bar + zbar))
       qq       = 0.5d0/pp *(z2bar + zbar)
       dppdt    = qq*rr*dtempi
       dppdd    = qq*ytot*tempi

! equation 6 of dewitt et al 
       qlam0z   = 1.88d8 * tempi * pp
       qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
       qlam0zdd = 1.88d8 * tempi * dppdd

       taufac   = co2 * tempi**x13
       taufacdt = -x13*taufac*tempi

       qq      = rr*zbar
       xni     = qq**x13
       dxnidd  = x13 * xni * deni

       aa     = 2.27493d5 * tempi * xni
       daadt  = 2.27493d5 * dtempi * xni
       daadd  = 2.27493d5 * tempi * dxnidd
      end if


! calculate individual screening factors
      bb       = z1 * z2
      gamp     = aa
      gampdt   = daadt
      gampdd   = daadd

      qq       = fact * bb * zs13inv(jscreen)
      gamef    = qq * gamp
      gamefdt  = qq * gampdt
      gamefdd  = qq * gampdd


! equation 6 of itoh et al
      tau12    = taufac * aznut(jscreen)
      tau12dt  = taufacdt * aznut(jscreen)

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

       qq       = zs13(jscreen)/(fact * bb)
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

       gamp14   = gamp**x14
       rr       = 1.0d0/gamp
       qq       = 0.25d0*gamp14*rr
       gamp14dt = qq * gampdt
       gamp14dd = qq * gampdd

! equation 25 of alastuey
! coefficients from equation 16 of jancovici

       cc       =   0.896434d0 * gamp * zhat(jscreen) &
                  - 3.44740d0  * gamp14 * zhat2(jscreen) &
                  - 0.5551d0   * (log(gamp) + lzav(jscreen)) &
                  - 2.996d0

       dccdt    =   0.896434d0 * gampdt * zhat(jscreen) &
                  - 3.44740d0  * gamp14dt * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdt

       dccdd    =   0.896434d0 * gampdd * zhat(jscreen) &
                  - 3.44740d0  * gamp14dd * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdd

       a3     = alph12 * alph12 * alph12
       da3    = 3.0d0 * alph12 * alph12


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
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt
       dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd


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

! quantum correction of jancovici 1977
       h12    = log(xlgfac) + h12
       rr     = 1.0d0/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt
       dh12dd = rr*dxlgfacdd + dh12dd


! intermediate regime
       if (gamef .le. gamefs) then
        iscreen_type = intermediate_screening

        rr     =  dgamma*(gamefs - gamef)
        drrdt  = -dgamma*gamefdt
        drrdd  = -dgamma*gamefdd

        ss     = dgamma*(gamef - gamefx)
        dssdt  = dgamma*gamefdt
        dssdd  = dgamma*gamefdd


! store current values for blending with weak regime
        h12x    = h12
        dh12xdt = dh12dt
        dh12xdd = dh12dd

! real linear ; gives same result as below
!         x1    = gamefx
!         x2    = gamefs
!         denom = 1.0d0/(x1-x2)
!         alfa  = (gamef - x2)*denom
!         beta  = (x1 - gamef)*denom
!         h12   = alfa*h12w + beta*h12x

!         x1    = h12w*rr + h12x*ss

         h12    = h12w*rr + h12x*ss
         dh12dt = dh12wdt*rr + h12w*drrdt + dh12xdt*ss + h12x*dssdt
         dh12dd = dh12wdd*rr + h12w*drrdd + dh12xdd*ss + h12x*dssdd


! blend the transition region - from bill paxton
!        if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
!         write(6,'(1x,a,1p2e12.4)') 'blending weak and intemediate'
!         alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
!         alfa   = 0.5d0 * (1.0d0 - cos(pi*alfa))
!         beta   = 1.0d0 - alfa

!         h12    = alfa * h12    + beta * h12w
!         dh12dt = alfa * dh12dt + beta * dh12wdt
!         dh12dd = alfa * dh12dd + beta * dh12wdd

!         h12    = alfa * h12 + beta * h12x
!         dh12dt = alfa * dh12dt + beta * dh12xdt
!         dh12dd = alfa * dh12dd + beta * dh12xdd

!        end if

       else
        iscreen_type = strong_screening
       end if


! end of intermediate and strong screening if
       else
        iscreen_type = weak_screening
      end if



! machine limit the output
! further limit to avoid the pyconuclear regime
      h12    = max(min(h12,30.0d0),0.0d0)
      scor   = exp(h12)
      if (h12 .eq. 30.0d0) then
       scordt = 0.0d0
       scordd = 0.0d0
      else
       scordt = scor * dh12dt
       scordd = scor * dh12dd
      end if
      screen_regime(jscreen) = iscreen_type
      screen_gamef(jscreen) = gamef


!      write(6,*) jscreen, screen_regime(jscreen)
!         write(6,'(1x,1p10e12.4)') temp,den,gamef,h12w,h12,float(iscreen_type)
!         read(5,*)


!      write(6,111) 'weak =',h12w,' total =',h12,
!     1             ' 1-ratio =',1.0d0-h12w/h12,' correction',scor
! 111  format(1x,4(a,1pe13.6))
!      read(5,*)

      return
    end subroutine screen5

      subroutine rhs21(y,rate,dydt,deriva)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'


! evaluates the right hand side of the aprox21 odes
! input:
! y      = array of abundances
! rate   = array of rates
! deriva = logical deciding if the rhs (deriva = false) or its 
!          derivative wrt to temp or den (deriva = true) is being evalauted
!
!output:
! dydt = array of rate of change of y with time


! declare the pass
      logical   ::  deriva
      real*8    ::  y(*),rate(*),dydt(*)


! local variables
      integer           :: i
      real*8, parameter :: sixth = 1.0d0/6.0d0

! quad precision dydt sums
      real*8           :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
      real*8           :: qray(abignet)



! zero the odes
      dydt(1:neqs) = 0.0d0
      qray(1:neqs) = 0.0q0


! set up the system of odes:

! hydrogen reactions
      a1 = -1.5d0 * y(ih1) * y(ih1) * rate(irpp)
      a2 =  y(ihe3) * y(ihe3) * rate(ir33) 
      a3 = -y(ihe3) * y(ihe4) * rate(irhe3ag) 
      a4 = -2.0d0 * y(ic12) * y(ih1) * rate(ircpg) 
      a5 = -2.0d0 * y(in14) * y(ih1) * rate(irnpg) 
      a6 = -2.0d0 * y(io16) * y(ih1) * rate(iropg) 
      a7 = -3.0d0 * y(ih1) * rate(irpen)

      qray(ih1) = qray(ih1) + a1 + a2 + a3 + a4 + a5 + a6 + a7      

! he3 reactions
      a1  =  0.5d0 * y(ih1) * y(ih1) * rate(irpp) 
      a2  = -y(ihe3) * y(ihe3) * rate(ir33) 
      a3  = -y(ihe3) * y(ihe4) * rate(irhe3ag) 
      a4  =  y(ih1) * rate(irpen)

      qray(ihe3) = qray(ihe3) + a1 + a2 + a3 + a4

! he4 reactions
! heavy ion reactions
      a1  = 0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
      a2  = 0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
      a3  = 0.56d0 * 0.5d0 * y(io16) * y(io16) * rate(ir1616)
      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3


! (a,g) and (g,a) reactions
      a1  = -0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
      a2  =  3.0d0 * y(ic12) * rate(irg3a) 
      a3  = -y(ihe4)  * y(ic12) * rate(ircag) 
      a4  =  y(io16)  * rate(iroga) 
      a5  = -y(ihe4)  * y(io16) * rate(iroag) 
      a6  =  y(ine20) * rate(irnega) 
      a7  = -y(ihe4)  * y(ine20) * rate(irneag) 
      a8  =  y(img24) * rate(irmgga) 
      a9  = -y(ihe4)  * y(img24)* rate(irmgag) 
      a10 =  y(isi28) * rate(irsiga) 
      a11 = -y(ihe4)  * y(isi28)*rate(irsiag) 
      a12 =  y(is32)  * rate(irsga)

      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12
            
      a1  = -y(ihe4)  * y(is32) * rate(irsag) 
      a2  =  y(iar36) * rate(irarga) 
      a3  = -y(ihe4)  * y(iar36)*rate(irarag) 
      a4  =  y(ica40) * rate(ircaga) 
      a5  = -y(ihe4)  * y(ica40)*rate(ircaag) 
      a6  =  y(iti44) * rate(irtiga) 
      a7  = -y(ihe4)  * y(iti44)*rate(irtiag) 
      a8  =  y(icr48) * rate(ircrga) 
      a9  = -y(ihe4)  * y(icr48)*rate(ircrag) 
      a10 =  y(ife52) * rate(irfega) 
      a11 = -y(ihe4)  * y(ife52) * rate(irfeag) 
      a12 =  y(ini56) * rate(irniga)

      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12


! (a,p)(p,g) and (g,p)(p,a) reactions

      if (.not.deriva) then
       a1  =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) 
       a2  = -y(ihe4)  * y(img24) * rate(irmgap)*(1.0d0-rate(irr1))
       a3  =  y(isi28) * rate(irsigp) * rate(irr1) 
       a4  = -y(ihe4)  * y(isi28) * rate(irsiap)*(1.0d0-rate(irs1)) 
       a5  =  y(is32)  * rate(irsgp) * rate(irs1) 
       a6  = -y(ihe4)  * y(is32) * rate(irsap)*(1.0d0-rate(irt1)) 
       a7  =  y(iar36) * rate(irargp) * rate(irt1) 
       a8  = -y(ihe4)  * y(iar36) * rate(irarap)*(1.0d0-rate(iru1)) 
       a9  =  y(ica40) * rate(ircagp) * rate(iru1) 
       a10 = -y(ihe4)  * y(ica40) * rate(ircaap)*(1.0d0-rate(irv1)) 
       a11 =  y(iti44) * rate(irtigp) * rate(irv1)
       a12 = -y(ihe4)  * y(iti44) * rate(irtiap)*(1.0d0-rate(irw1)) 
       a13 =  y(icr48) * rate(ircrgp) * rate(irw1) 
       a14 = -y(ihe4)  * y(icr48) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a15 =  y(ife52) * rate(irfegp) * rate(irx1) 

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15


      else
       a1  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1) * rate(ir1616) 
       a2  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1) * ratdum(ir1616)
       a3  = -y(ihe4)*y(img24) * rate(irmgap)*(1.0d0 - ratdum(irr1)) 
       a4  =  y(ihe4)*y(img24) * ratdum(irmgap)*rate(irr1)
       a5  =  y(isi28) * ratdum(irsigp) * rate(irr1) 
       a6  =  y(isi28) * rate(irsigp) * ratdum(irr1)
       a7  = -y(ihe4)*y(isi28) * rate(irsiap)*(1.0d0 - ratdum(irs1)) 
       a8  =  y(ihe4)*y(isi28) * ratdum(irsiap) * rate(irs1)
       a9  =  y(is32)  * ratdum(irsgp) * rate(irs1)
       a10 =  y(is32)  * rate(irsgp) * ratdum(irs1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10

       a1  = -y(ihe4)*y(is32) * rate(irsap)*(1.0d0 - ratdum(irt1)) 
       a2  =  y(ihe4)*y(is32) * ratdum(irsap)*rate(irt1)
       a3  =  y(iar36) * ratdum(irargp) * rate(irt1) 
       a4  =  y(iar36) * rate(irargp) * ratdum(irt1)
       a5  = -y(ihe4)*y(iar36) * rate(irarap)*(1.0d0 - ratdum(iru1))
       a6  =  y(ihe4)*y(iar36) * ratdum(irarap)*rate(iru1)
       a7  =  y(ica40) * ratdum(ircagp) * rate(iru1)
       a8  =  y(ica40) * rate(ircagp) * ratdum(iru1)
       a9  = -y(ihe4)*y(ica40) * rate(ircaap)*(1.0d0-ratdum (irv1)) 
       a10 =  y(ihe4)*y(ica40) * ratdum(ircaap)*rate(irv1)
       a11 =  y(iti44) * ratdum(irtigp) * rate(irv1)
       a12 =  y(iti44) * rate(irtigp) * ratdum(irv1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12

       a1  = -y(ihe4)*y(iti44) * rate(irtiap)*(1.0d0 - ratdum(irw1))
       a2  =  y(ihe4)*y(iti44) * ratdum(irtiap)*rate(irw1)
       a3  =  y(icr48) * ratdum(ircrgp) * rate(irw1) 
       a4  =  y(icr48) * rate(ircrgp) * ratdum(irw1)
       a5  = -y(ihe4)*y(icr48) * rate(ircrap)*(1.0d0 - ratdum(irx1))
       a6  =  y(ihe4)*y(icr48) * ratdum(ircrap)*rate(irx1)
       a7  =  y(ife52) * ratdum(irfegp) * rate(irx1) 
       a8  =  y(ife52) * rate(irfegp) * ratdum(irx1)
       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8
      end if


! photodisintegration reactions
      a1 =  y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
      a2 = -y(ife52) * y(ihe4) * rate(ir6f54) 
      a3 = -y(ife52) * y(ihe4) * y(iprot) * rate(ir7f54) 
      a4 =  y(ini56) * y(iprot) * rate(ir8f54) 
      a5 = -y(ihe4) * rate(iralf1) 
      a6 =  y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2) 
      a7 =  y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
      a8 = -y(ife54) * y(ihe4) * rate(irfe56_aux4) 

      qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8


! ppchain
      a1 = 0.5d0 * y(ihe3) * y(ihe3) * rate(ir33)
      a2 = y(ihe3) * y(ihe4) * rate(irhe3ag)

      qray(ihe4) =  qray(ihe4) + a1 + a2 


! cno cycles
      a1 = y(io16) * y(ih1) * rate(iropg) 
      a2 = -y(ihe4) * y(in14) * rate(irnag) * 1.5d0

      qray(ihe4) =  qray(ihe4) + a1 + a2

      if (.not. deriva) then
       a1 = y(in14) * y(ih1) * rate(ifa) * rate(irnpg) 
       qray(ihe4) =  qray(ihe4) + a1
      else
       a1 = y(in14) * y(ih1) * rate(ifa) * ratdum(irnpg) 
       a2 = y(in14) * y(ih1) * ratdum(ifa) * rate(irnpg) 
       qray(ihe4) =  qray(ihe4) + a1 + a2
      end if



! c12 reactions
      a1 = -y(ic12) * y(ic12) * rate(ir1212) 
      a2 = -y(ic12) * y(io16) * rate(ir1216) 
      a3 =  sixth * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
      a4 = -y(ic12) * rate(irg3a) 
      a5 = -y(ic12) * y(ihe4) * rate(ircag) 
      a6 =  y(io16) * rate(iroga) 
      a7 = -y(ic12) * y(ih1) * rate(ircpg) 

      qray(ic12) =  qray(ic12) + a1 + a2 + a3 + a4 + a5 + a6 + a7

      if (.not. deriva) then
       a1 =  y(in14) * y(ih1) * rate(ifa) * rate(irnpg)
       qray(ic12) =  qray(ic12) + a1
      else
       a1 =  y(in14) * y(ih1) * rate(ifa) * ratdum(irnpg)
       a2 =  y(in14) * y(ih1) * ratdum(ifa) * rate(irnpg)
       qray(ic12) =  qray(ic12) + a1 + a2
      end if


! n14 reactions
      a1 =  y(ic12) * y(ih1) * rate(ircpg) 
      a2 = -y(in14) * y(ih1) * rate(irnpg) 
      a3 =  y(io16) * y(ih1) * rate(iropg) 
      a4 = -y(ihe4) * y(in14) * rate(irnag)

      qray(in14) =  qray(in14) + a1 + a2 + a3 + a4 


! o16 reactions
      a1 = -y(ic12) * y(io16) * rate(ir1216) 
      a2 = -y(io16) * y(io16) * rate(ir1616) 
      a3 =  y(ic12) * y(ihe4) * rate(ircag) 
      a4 = -y(io16) * y(ihe4) * rate(iroag) 
      a5 = -y(io16) * rate(iroga) 
      a6 =  y(ine20) * rate(irnega) 
      a7 = -y(io16) * y(ih1) * rate(iropg)

      qray(io16) =  qray(io16) + a1 + a2 + a3 + a4 + a5 + a6 + a7

      if (.not. deriva) then
       a1 =  y(in14) * y(ih1) * rate(ifg) * rate(irnpg) 
       qray(io16) =  qray(io16) + a1
      else
       a1 =  y(in14) * y(ih1) * rate(ifg) * ratdum(irnpg) 
       a2 =  y(in14) * y(ih1) * ratdum(ifg) * rate(irnpg) 
       qray(io16) =  qray(io16) + a1 + a2
      end if



! ne20 reactions
      a1 =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
      a2 =  y(io16) * y(ihe4) * rate(iroag) 
      a3 = -y(ine20) * y(ihe4) * rate(irneag) 
      a4 = -y(ine20) * rate(irnega) 
      a5 =  y(img24) * rate(irmgga) 
      a6 =  y(in14) * y(ihe4) * rate(irnag)

      qray(ine20) =  qray(ine20) + a1 + a2 + a3 + a4 + a5 + a6


! mg24 reactions
      a1 =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
      a2 =  y(ine20) * y(ihe4) * rate(irneag) 
      a3 = -y(img24) * y(ihe4) * rate(irmgag) 
      a4 = -y(img24) * rate(irmgga) 
      a5 =  y(isi28) * rate(irsiga)

      qray(img24) =  qray(img24) + a1 + a2 + a3 + a4 + a5 

      if (.not.deriva) then
       a1 = -y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) 
       a2 =  y(isi28) * rate(irr1) * rate(irsigp)

       qray(img24) =  qray(img24) + a1 + a2

      else
       a1 = -y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a2 =  y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
       a3 =  y(isi28) * ratdum(irr1) * rate(irsigp) 
       a4 =  y(isi28) * rate(irr1) * ratdum(irsigp)

       qray(img24) =  qray(img24) + a1 + a2 + a3 + a4
      end if



! si28 reactions
      a1 =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
      a2 =  0.56d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616) 
      a3 =  y(img24) * y(ihe4) * rate(irmgag) 
      a4 = -y(isi28) * y(ihe4) * rate(irsiag) 
      a5 = -y(isi28) * rate(irsiga) 
      a6 =  y(is32)  * rate(irsga)

      qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5 + a6

      if (.not.deriva) then
       a1 =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) 
       a2 =  y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) 
       a3 = -y(isi28) * rate(irr1) * rate(irsigp) 
       a4 = -y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) 
       a5 =  y(is32)  * rate(irs1) * rate(irsgp)

       qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5

      else
       a1  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1)*rate(ir1616) 
       a2  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1)*ratdum(ir1616)
       a3  =  y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a4  = -y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
       a5  = -y(isi28) * ratdum(irr1) * rate(irsigp) 
       a6  = -y(isi28) * rate(irr1) * ratdum(irsigp)
       a7  = -y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0 - ratdum(irs1))
       a8  =  y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1)
       a9  = y(is32) * ratdum(irs1) * rate(irsgp)
       a10 = y(is32) * rate(irs1) * ratdum(irsgp)

       qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10
      end if



! s32 reactions
      a1 =  0.1d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616) 
      a2 =  y(isi28) * y(ihe4) * rate(irsiag) 
      a3 = -y(is32) * y(ihe4) * rate(irsag) 
      a4 = -y(is32) * rate(irsga) 
      a5 =  y(iar36) * rate(irarga)

       qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5

      if (.not.deriva) then
       a1 =  0.34d0*0.5d0*y(io16)*y(io16)* rate(ir1616)*(1.0d0-rate(irs1)) 
       a2 =  y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) 
       a3 = -y(is32) * rate(irs1) * rate(irsgp) 
       a4 = -y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) 
       a5 =  y(iar36) * rate(irt1) * rate(irargp)

       qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5

      else
       a1  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(ir1616)*(1.0d0-ratdum(irs1))
       a2  = -0.34d0*0.5d0*y(io16)*y(io16) * ratdum(ir1616)*rate(irs1)
       a3  =  y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0-ratdum(irs1))
       a4  = -y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1) 
       a5  = -y(is32) * ratdum(irs1) * rate(irsgp) 
       a6  = -y(is32) * rate(irs1) * ratdum(irsgp)
       a7  = -y(is32)*y(ihe4) * rate(irsap)*(1.0d0-ratdum(irt1))
       a8  =  y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
       a9  =  y(iar36) * ratdum(irt1) * rate(irargp)
       a10 =  y(iar36) * rate(irt1) * ratdum(irargp)

       qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10
      end if



! ar36 reactions
      a1 =  y(is32)  * y(ihe4) * rate(irsag) 
      a2 = -y(iar36) * y(ihe4) * rate(irarag)
      a3 = -y(iar36) * rate(irarga) 
      a4 =  y(ica40) * rate(ircaga)

      qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 = y(is32)  * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) 
       a2 = -y(iar36) * rate(irt1) * rate(irargp) 
       a3 = -y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) 
       a4 =  y(ica40) * rate(ircagp) * rate(iru1)

       qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4

      else
       a1 =  y(is32)*y(ihe4) * rate(irsap)*(1.0d0 - ratdum(irt1)) 
       a2 = -y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
       a3 = -y(iar36) * ratdum(irt1) * rate(irargp) 
       a4 = -y(iar36) * rate(irt1) * ratdum(irargp)
       a5 = -y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
       a6 =  y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
       a7 =  y(ica40) * ratdum(ircagp) * rate(iru1) 
       a8 =  y(ica40) * rate(ircagp) * ratdum(iru1)

       qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if



! ca40 reactions
      a1 =  y(iar36) * y(ihe4) * rate(irarag)
      a2 = -y(ica40) * y(ihe4) * rate(ircaag)
      a3 = -y(ica40) * rate(ircaga)
      a4 =  y(iti44) * rate(irtiga)

      qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) 
       a2 = -y(ica40) * rate(ircagp) * rate(iru1) 
       a3 = -y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) 
       a4 =  y(iti44) * rate(irtigp) * rate(irv1)

      qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4

      else
       a1 =  y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
       a2 = -y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
       a3 = -y(ica40) * ratdum(ircagp) * rate(iru1) 
       a4 = -y(ica40) * rate(ircagp) * ratdum(iru1)
       a5 = -y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1)) 
       a6 =  y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
       a7 =  y(iti44) * ratdum(irtigp) * rate(irv1) 
       a8 =  y(iti44) * rate(irtigp) * ratdum(irv1)

       qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if



! ti44 reactions
      a1 =  y(ica40) * y(ihe4) * rate(ircaag) 
      a2 = -y(iti44) * y(ihe4) * rate(irtiag) 
      a3 = -y(iti44) * rate(irtiga) 
      a4 =  y(icr48) * rate(ircrga)

      qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) 
       a2 = -y(iti44) * rate(irv1) * rate(irtigp) 
       a3 = -y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) 
       a4 =  y(icr48) * rate(irw1) * rate(ircrgp)

      qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4

      else
       a1 =  y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1)) 
       a2 = -y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
       a3 = -y(iti44) * ratdum(irv1) * rate(irtigp) 
       a4 = -y(iti44) * rate(irv1) * ratdum(irtigp)
       a5 = -y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1))
       a6 =  y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
       a7 =  y(icr48) * ratdum(irw1) * rate(ircrgp)
       a8 =  y(icr48) * rate(irw1) * ratdum(ircrgp)

       qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if



! cr48 reactions
      a1 =  y(iti44) * y(ihe4) * rate(irtiag) 
      a2 = -y(icr48) * y(ihe4) * rate(ircrag) 
      a3 = -y(icr48) * rate(ircrga) 
      a4 =  y(ife52) * rate(irfega)

      qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) 
       a2 = -y(icr48) * rate(irw1) * rate(ircrgp) 
       a3 = -y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a4 =  y(ife52) * rate(irx1) * rate(irfegp)

       qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4

      else
        a1 =  y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1)) 
        a2 = -y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
        a3 = -y(icr48) * ratdum(irw1) * rate(ircrgp) 
        a4 = -y(icr48) * rate(irw1) * ratdum(ircrgp)
        a5 = -y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
        a6 =  y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
        a7 =  y(ife52) * ratdum(irx1) * rate(irfegp) 
        a8 =  y(ife52) * rate(irx1) * ratdum(irfegp)

       qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! cr56 reactions
      qray(icr56) = qray(icr56) + y(ife56) * 1.0d-04 * rate(irn56ec)



! fe52 reactions
      a1 =  y(icr48) * y(ihe4) * rate(ircrag) 
      a2 = -y(ife52) * y(ihe4) * rate(irfeag) 
      a3 = -y(ife52) * rate(irfega) 
      a4 =  y(ini56) * rate(irniga)

      qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a2 = -y(ife52) * rate(irx1) * rate(irfegp) 
       qray(ife52) =  qray(ife52) + a1 + a2

      else
       a1 =  y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
       a2 = -y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
       a3 = -y(ife52) * ratdum(irx1) * rate(irfegp) 
       a4 = -y(ife52) * rate(irx1) * ratdum(irfegp)
       qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4 
      end if


! photodisintegration reactions
      a1 =  y(ife54) * rate(ir1f54) 
      a2 = -y(ife52) * y(ineut) * y(ineut) * rate(ir2f54) 
      a3 =  y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
      a4 = -y(ife52) * y(ihe4) * rate(ir6f54) 
      a5 = -y(ife52) * y(ihe4) * y(iprot) * rate(ir7f54) 
      a6 =  y(ini56) * y(iprot) * rate(ir8f54) 

      qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4 + a5 + a6 



! fe54 reactions
      a1  = -y(ife54) * rate(ir1f54)
      a2  =  y(ife52) * y(ineut) * y(ineut) * rate(ir2f54) 
      a3  = -y(ife54) * y(iprot) * y(iprot) * rate(ir3f54) 
      a4  =  y(ini56) * rate(ir4f54) 
      a5  = -y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
      a6  =  y(ife52) * y(ihe4) * rate(ir6f54) 
      a7  =  y(ife56) * rate(irfe56_aux1) 
      a8  = -y(ife54) * y(ineut) * y(ineut) * rate(irfe56_aux2) 
      a9  =  y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
      a10 = -y(ife54) * y(ihe4) * rate(irfe56_aux4) 

      qray(ife54) =  qray(ife54) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10


! fe56 reactions
      a1 =  y(ini56) * rate(irn56ec) 
      a2 = -y(ife56) * 1.0d-04 * rate(irn56ec) 
      a3 = -y(ife56) * rate(irfe56_aux1) 
      a4 =  y(ife54) * y(ineut) * y(ineut) * rate(irfe56_aux2)  
      a5 = -y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
      a6 =  y(ife54) * y(ihe4) * rate(irfe56_aux4) 

      qray(ife56) =  qray(ife56) + a1 + a2 + a3 + a4 + a5 + a6 


! ni56 reactions
      a1 =  y(ife52) * y(ihe4) * rate(irfeag) 
      a2 = -y(ini56) * rate(irniga) 
      a3 = -y(ini56) * rate(irn56ec)

      qray(ini56) =  qray(ini56) + a1 + a2 + a3


! photodisintegration reactions
      a1 =  y(ife54) * y(iprot) * y(iprot) * rate(ir3f54) 
      a2 = -y(ini56) * rate(ir4f54) 
      a3 =  y(ife52) * y(ihe4)* y(iprot) * rate(ir7f54) 
      a4 = -y(ini56) * y(iprot) * rate(ir8f54)

      qray(ini56) =  qray(ini56) + a1 + a2 + a3 + a4



! photodisintegration neutrons
      a1 =  2.0d0 * y(ife54) * rate(ir1f54) 
      a2 = -2.0d0 * y(ife52) * y(ineut) * y(ineut) * rate(ir2f54) 
      a3 =  2.0d0 * y(ihe4) * rate(iralf1) 
      a4 = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2) 
      a5 =  y(iprot) * rate(irpen) 
      a6 = -y(ineut) * rate(irnep) 
      a7 =  2.0d0 * y(ife56) * rate(irfe56_aux1) 
      a8 = -2.0d0 * y(ife54) * y(ineut) * y(ineut) * rate(irfe56_aux2)

      qray(ineut) =  qray(ineut) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8


! photodisintegration protons
      a1  = -2.0d0 * y(ife54) * y(iprot) * y(iprot) * rate(ir3f54) 
      a2  =  2.0d0 * y(ini56) * rate(ir4f54) 
      a3  = -2.0d0 * y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
      a4  =  2.0d0 * y(ife52) * y(ihe4) * rate(ir6f54) 
      a5  =  2.0d0 * y(ihe4) * rate(iralf1) 
      a6  = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2) 
      a7  = -y(iprot) * rate(irpen) 
      a8  =  y(ineut) * rate(irnep)  
      a9  = -2.0d0 * y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
      a10 =  2.0d0 * y(ife54) * y(ihe4) * rate(irfe56_aux4)

      qray(iprot) =  qray(iprot) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10



! now set the double precision return argument dydt

      dydt(1:ionmax) = qray(1:ionmax)
 
      qwork1(1:ionmax) = qray(1:ionmax)

      return
      end subroutine rhs21


