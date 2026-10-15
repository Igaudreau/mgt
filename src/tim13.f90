!timreac13
!timrat13
!init_aprox13
!tfactors13(temp)
!aprox13rat
!screen_aprox13(y)
!screen_13_5(temp,den,zbar,abar,z2bar
!rhs13(y,rate,dydt,deriva)
!rate13_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
!rate13_ti44ap(temp,den, &
!rate13_v47pg(temp,den, &
!rate13_cr48ag(temp,den, &
!rate13_cr48ap(temp,den, &
!rate13_mn51pg(temp,den, &
!rate13_fe52ag(temp,den, &
!rate13_fe52ap(temp,den, &
!rate13_co55pg(temp,den, &
!azbar_simple

!     Get Timmes's reaction rates from MG solution
      subroutine timreact13(nd, pa, ydot)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'burn_common.dek'
!     Input
      integer   ::  nd
      real*8    ::  pa(*)
!     Output      
      real*8    ::  ydot(50)      
!     Local
      logical   :: deriva
      integer   :: i      
      real*8    :: y(50), dydt(50),temp, den, pres, enuc
      deriva = .false.
! Set temperature, density and pressure
      temp = pa(17+nd)
      den = pa(1)
      pres = pa(nd+2)
      call tfactors13(temp)
      bden = den
      btemp = temp
      bpres = pres
!     Set abundances
      do i= 1,13
         y(i) = pa(2+nd+i)
      enddo
!     Get rates
      call aprox13rat
!     Do screening
      screen_on = 1
      call screen_aprox13(y)      
      call rhs13(y,ratdum,dydt,deriva)
      do i= 1,13
         ydot(i) = dydt(i)
      enddo
      return
    end subroutine timreact13

      subroutine timrat13(nd, pa, RR)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'burn_common.dek'
!     Input
      integer   :: nd     
      real*8    :: pa(*), RR(*)
!     Local
      logical   :: deriva
      integer   :: i      
      real*8    :: y(50), temp, den, pres, enuc
      deriva = .false.
!     Set temperature, density and pressure
      temp = pa(17+nd)
      den = pa(1)
      pres = pa(nd+2)
      call tfactors13(temp)
      bden = den
      btemp = temp
      bpres = pres
!     Set abundances
      do i=1,13
         y(i) = pa(3+i)
      enddo
!     Get rates
      call aprox13rat
!     Do screening
      screen_on = 1
      call screen_aprox13(y)
      do i=1,67
         RR(i) = ratdum(i)
      enddo      
      return
    end subroutine timrat13
    

      subroutine aprox13_init(eaion, ezion, ebion)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

! this routine initializes stuff for the aprox13 network
!
!     Output
      real*8 :: eaion(100), ezion(100), ebion(100)                
! declare
      integer           :: i
      real*8, parameter :: mev2erg = ev2erg*1.0d6, mev2gr  = mev2erg/clight**2


! set the size of the network and the number of rates
      idnet   = idaprox13

      ionmax  = 13
      iener   = ionmax + 1
      itemp   = ionmax + 2
      iden    = ionmax + 3
      ivelx   = ionmax + 4
      iposx   = ionmax + 5
      neqs    = iposx


      nrat    = 67
      netname = 'aprox13'


! set the id numbers of the elements

      ionbeg = 1
      ionend = ionmax

      ihe4  = ionbeg
      ic12  = ionbeg + 1
      io16  = ionbeg + 2
      ine20 = ionbeg + 3
      img24 = ionbeg + 4
      isi28 = ionbeg + 5
      is32  = ionbeg + 6
      iar36 = ionbeg + 7
      ica40 = ionbeg + 8
      iti44 = ionbeg + 9
      icr48 = ionbeg + 10
      ife52 = ionbeg + 11
      ini56 = ionbeg + 12


! set the names of the elements
      ionam(ihe4)  = 'he4  '
      ionam(ic12)  = 'c12  '
      ionam(io16)  = 'o16  '
      ionam(ine20) = 'ne20 '
      ionam(img24) = 'mg24 '
      ionam(isi28) = 'si28 '
      ionam(is32)  = 's32  '
      ionam(iar36) = 'ar36 '
      ionam(ica40) = 'ca40 '
      ionam(iti44) = 'ti44 '
      ionam(icr48) = 'cr48 '
      ionam(ife52) = 'fe52 '
      ionam(ini56) = 'ni56 '

      ionam(iener) = 'ener '
      ionam(itemp) = 'temp '
      ionam(iden)  = 'den  '
      ionam(ivelx) = 'velx '
      ionam(iposx) = 'posx '


! set the number of nucleons in the element
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(io16)  = 16.0d0
      aion(ine20) = 20.0d0
      aion(img24) = 24.0d0
      aion(isi28) = 28.0d0
      aion(is32)  = 32.0d0
      aion(iar36) = 36.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(ife52) = 52.0d0
      aion(ini56) = 56.0d0

! set the number of protons in the element
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0
      zion(img24) = 12.0d0
      zion(isi28) = 14.0d0
      zion(is32)  = 16.0d0
      zion(iar36) = 18.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(ife52) = 26.0d0
      zion(ini56) = 28.0d0


! set the binding energy of the element
      bion(ihe4)  =  28.29603d0
      bion(ic12)  =  92.16294d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0
      bion(img24) = 198.25790d0
      bion(isi28) = 236.53790d0
      bion(is32)  = 271.78250d0
      bion(iar36) = 306.72020d0
      bion(ica40) = 342.05680d0
      bion(iti44) = 375.47720d0
      bion(icr48) = 411.46900d0
      bion(ife52) = 447.70800d0
      bion(ini56) = 484.00300d0


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


! set the id numbers of the reaction rates
      ir3a   = 1
      irg3a  = 2
      ircag  = 3
      iroga  = 4
      ir1212 = 5
      ir1216 = 6
      ir1616 = 7
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

      irr1   = 60
      irs1   = 61
      irt1   = 62
      iru1   = 63
      irv1   = 64
      irw1   = 65
      irx1   = 66
      iry1   = 67


! set the names of the reaction rates
      ratnam(ir3a)   = 'r3a  '
      ratnam(irg3a)  = 'rg3a '
      ratnam(ircag)  = 'rcag '
      ratnam(ir1212) = 'r1212'
      ratnam(ir1216) = 'r1216'
      ratnam(ir1616) = 'r1616'
      ratnam(iroga)  = 'roga '
      ratnam(iroag)  = 'roag '
      ratnam(irnega) = 'rnega'
      ratnam(irneag) = 'rneag'
      ratnam(irmgga) = 'rmgga'
      ratnam(irmgag) = 'rmgag'
      ratnam(irsiga) = 'rsiga'
      ratnam(irmgap) = 'rmgap'
      ratnam(iralpa) = 'ralpa'
      ratnam(iralpg) = 'ralpg'
      ratnam(irsigp) = 'rsigp'
      ratnam(irsiag) = 'rsiag'
      ratnam(irsga)  = 'rsga '
      ratnam(irsiap) = 'rsiap'
      ratnam(irppa)  = 'rppa '
      ratnam(irppg)  = 'rppg '
      ratnam(irsgp)  = 'rsgp '
      ratnam(irsag)  = 'rsag '
      ratnam(irarga) = 'rarga'
      ratnam(irsap)  = 'rsap '
      ratnam(irclpa) = 'rclpa'
      ratnam(irclpg) = 'rclpg'
      ratnam(irargp) = 'rargp'
      ratnam(irarag) = 'rarag'
      ratnam(ircaga) = 'rcaga'
      ratnam(irarap) = 'rarap'
      ratnam(irkpa)  = 'rkpa '
      ratnam(irkpg)  = 'rkpg '
      ratnam(ircagp) = 'rcagp'
      ratnam(ircaag) = 'rcaag'
      ratnam(irtiga) = 'rtiga'
      ratnam(ircaap) = 'rcaap'
      ratnam(irscpa) = 'rscpa'
      ratnam(irscpg) = 'rscpg'
      ratnam(irtigp) = 'rtigp'
      ratnam(irtiag) = 'rtiag'
      ratnam(ircrga) = 'rcrga'
      ratnam(irtiap) = 'rtiap'
      ratnam(irvpa)  = 'rvpa '
      ratnam(irvpg)  = 'rvpg '
      ratnam(ircrgp) = 'rcrgp'
      ratnam(ircrag) = 'rcrag'
      ratnam(irfega) = 'rfega'
      ratnam(ircrap) = 'rcrap'
      ratnam(irmnpa) = 'rmnpa'
      ratnam(irmnpg) = 'rmnpg'
      ratnam(irfegp) = 'rfegp'
      ratnam(irfeag) = 'rfeag'
      ratnam(irniga) = 'rniga'
      ratnam(irfeap) = 'rfeap'
      ratnam(ircopa) = 'rcopa'
      ratnam(ircopg) = 'rcopg'
      ratnam(irnigp) = 'rnigp'

      ratnam(irr1)   = 'irr1 '
      ratnam(irs1)   = 'irs1 '
      ratnam(irt1)   = 'irt1 '
      ratnam(iru1)   = 'iru1 '
      ratnam(irv1)   = 'irv1 '
      ratnam(irw1)   = 'irw1 '
      ratnam(irx1)   = 'irx1 '
      ratnam(iry1)   = 'iry1 '

      ! eaion, ezion, ebion
      eaion(1:ionmax) = aion(1:ionmax)
      ezion(1:ionmax) = zion(1:ionmax)
      ebion(1:ionmax) = bion(1:ionmax)
      return
    end subroutine aprox13_init

      subroutine tfactors13(temp)
      include 'implno.dek'
      include 'tfactors.dek'

! sets various popular temperature factors into common block
! this routine must be called before any of the rates are called

! declare the pass
      real*8 temp

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
    end subroutine tfactors13
    
      subroutine aprox13rat
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine generates unscreened
! nuclear reaction rates for the aprox13 network.

! declare
      integer          i
      double precision rrate,drratedt,drratedd


! zero the rates
      ratraw(1:nrat) = 0.0d0
      dratrawdt(1:nrat) = 0.0d0
      dratrawdd(1:nrat) = 0.0d0


      if (btemp .lt. 1.0e6) return


! get the temperature factors
      call tfactors13(btemp)


! triple alpha to c12
      call rate13_tripalf(btemp,bden, &
           ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
           ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

! c12(a,g)o16
      call rate13_c12ag(btemp,bden, &
           ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
           ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

! c12 + c12
      call rate13_c12c12(btemp,bden, &
           ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
           rrate,drratedt,drratedd)

! c12 + o16
      call rate13_c12o16(btemp,bden, &
           ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
           rrate,drratedt,drratedd)

! o16 + o16
      call rate13_o16o16(btemp,bden, &
           ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
           rrate,drratedt,drratedd)

! o16(a,g)ne20
      call rate13_o16ag(btemp,bden, &
           ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
           ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

! ne20(a,g)mg24
      call rate13_ne20ag(btemp,bden, &
           ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
           ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

! mg24(a,g)si28
      call rate13_mg24ag(btemp,bden, &
           ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
           ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

! mg24(a,p)al27
      call rate13_mg24ap(btemp,bden, &
           ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap), &
           ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))

! al27(p,g)si28
      call rate13_al27pg(btemp,bden, &
           ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg), &
           ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))

! si28(a,g)s32
      call rate13_si28ag(btemp,bden, &
           ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag), &
           ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))

! si28(a,p)p31
      call rate13_si28ap(btemp,bden, &
           ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap), &
           ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))

! p31(p,g)s32
      call rate13_p31pg(btemp,bden, &
           ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg), &
           ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))

! s32(a,g)ar36
      call rate13_s32ag(btemp,bden, &
           ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag), &
           ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))

! s32(a,p)cl35
      call rate13_s32ap(btemp,bden, &
           ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap), &
           ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))

! cl35(p,g)ar36
      call rate13_cl35pg(btemp,bden, &
           ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg), &
           ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))

! ar36(a,g)ca40
      call rate13_ar36ag(btemp,bden, &
           ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag), &
           ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))

! ar36(a,p)k39
      call rate13_ar36ap(btemp,bden, &
           ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap), &
           ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))

! k39(p,g)ca40
      call rate13_k39pg(btemp,bden, &
           ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg), &
           ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))

! ca40(a,g)ti44
      call rate13_ca40ag(btemp,bden, &
           ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
           ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

! ca40(a,p)sc43
      call rate13_ca40ap(btemp,bden, &
           ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap), &
           ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))

! sc43(p,g)ti44
      call rate13_sc43pg(btemp,bden, &
           ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg), &
           ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))

! ti44(a,g)cr48
      call rate13_ti44ag(btemp,bden, &
           ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag), &
           ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

! ti44(a,p)v47
      call rate13_ti44ap(btemp,bden, &
           ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
           ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))

! v47(p,g)cr48
      call rate13_v47pg(btemp,bden, &
           ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg), &
           ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))

! cr48(a,g)fe52
      call rate13_cr48ag(btemp,bden, &
           ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag), &
           ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))

! cr48(a,p)mn51
      call rate13_cr48ap(btemp,bden, &
           ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap), &
           ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))

! mn51(p,g)fe52
      call rate13_mn51pg(btemp,bden, &
           ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg), &
           ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))

! fe52(a,g)ni56
      call rate13_fe52ag(btemp,bden, &
           ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag), &
           ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))

! fe52(a,p)co55
      call rate13_fe52ap(btemp,bden, &
           ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap), &
           ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))

! co55(p,g)ni56
      call rate13_co55pg(btemp,bden, &
           ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg), &
           ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))


! write out the rates
!      write(6,133) btemp,bden
! 133  format(1x,1pe12.4)
!      do i=1,nrat
!       write(6,134) ratnam(i),ratraw(i)
! 134   format(1x,a,'  ',1pe14.6,1pe11.3,1pe14.6)
!      enddo
!      read(5,*)



! for a strict alpha chain with only (a,g) and (g,a) reactions
! shut down the (a,p) (p,g) and (g,p) (p,a) rates

! irr1
!           ratraw(irmgag) = 0.0d0
!           ratraw(irsiga) = 0.0d0
!      ratraw(irmgap)    = 0.0d0  ; dratrawdt(irmgap) = 0.0d0  ; dratrawdd(irmgap) = 0.0d0
!      ratraw(iralpa)    = 0.0d0  ; dratrawdt(iralpa) = 0.0d0  ; dratrawdd(iralpa) = 0.0d0
!      ratraw(iralpg)    = 0.0d0  ; dratrawdt(iralpg) = 0.0d0  ; dratrawdd(iralpg) = 0.0d0
!      ratraw(irsigp)    = 0.0d0  ; dratrawdt(irsigp) = 0.0d0  ; dratrawdd(irsigp) = 0.0d0

! irs1
!           ratraw(ir1616) = 0.0d0
!           ratraw(irsiag) = 0.0d0
!           ratraw(irsga) = 0.0d0
!      ratraw(irsiap) = 0.0d0  ; dratrawdt(irsiap) = 0.0d0  ; dratrawdd(irsiap) = 0.0d0
!      ratraw(irppa)  = 0.0d0  ; dratrawdt(irppa)  = 0.0d0  ; dratrawdd(irppa)  = 0.0d0
!      ratraw(irppg)  = 0.0d0  ; dratrawdt(irppg)  = 0.0d0  ; dratrawdd(irppg)  = 0.0d0
!      ratraw(irsgp)  = 0.0d0  ; dratrawdt(irsgp)  = 0.0d0  ; dratrawdd(irsgp)  = 0.0d0

! irt1 
!           ratraw(irsag) = 0.0d0
!           ratraw(irarga) = 0.0d0
!      ratraw(irsap)  = 0.0d0  ;  dratrawdt(irsap)  = 0.0d0  ; dratrawdd(irsap)  = 0.0d0
!      ratraw(irclpa) = 0.0d0  ;  dratrawdt(irclpa) = 0.0d0  ; dratrawdd(irclpa) = 0.0d0
!      ratraw(irclpg) = 0.0d0  ;  dratrawdt(irclpg) = 0.0d0  ; dratrawdd(irclpg) = 0.0d0
!      ratraw(irargp) = 0.0d0  ;  dratrawdt(irargp) = 0.0d0  ; dratrawdd(irargp) = 0.0d0

! iru1 
!           ratraw(irarag) = 0.0d0
!           ratraw(ircaga) = 0.0d0
!      ratraw(irarap) = 0.0d0  ; dratrawdt(irarap) = 0.0d0  ; dratrawdd(irarap) = 0.0d0
!      ratraw(irkpa)  = 0.0d0  ;  dratrawdt(irkpa)  = 0.0d0 ; dratrawdd(irkpa)  = 0.0d0
!      ratraw(irkpg)  = 0.0d0  ; dratrawdt(irkpg)  = 0.0d0  ; dratrawdd(irkpg)  = 0.0d0
!      ratraw(ircagp) = 0.0d0  ; dratrawdt(ircagp) = 0.0d0  ; dratrawdd(ircagp) = 0.0d0

! irv1 
!           ratraw(ircaag) = 0.0d0
!           ratraw(irtiga) = 0.0d0
!      ratraw(ircaap) = 0.0d0  ; dratrawdt(ircaap) = 0.0d0  ; dratrawdd(ircaap) = 0.0d0
!      ratraw(irscpa) = 0.0d0  ; dratrawdt(irscpa) = 0.0d0  ; dratrawdd(irscpa) = 0.0d0
!      ratraw(irscpg) = 0.0d0  ; dratrawdt(irscpg) = 0.0d0  ; dratrawdd(irscpg) = 0.0d0
!      ratraw(irtigp) = 0.0d0  ; dratrawdt(irtigp) = 0.0d0  ; dratrawdd(irtigp) = 0.0d0

! irw1 
!           ratraw(irtiag) = 0.0d0
!           ratraw(ircrga) = 0.0d0
!      ratraw(irtiap) = 0.0d0  ; dratrawdt(irtiap) = 0.0d0  ; dratrawdd(irtiap) = 0.0d0
!      ratraw(irvpa)  = 0.0d0  ; dratrawdt(irvpa)  = 0.0d0  ; dratrawdd(irvpa)  = 0.0d0
!      ratraw(irvpg)  = 0.0d0  ; dratrawdt(irvpg)  = 0.0d0  ; dratrawdd(irvpg)  = 0.0d0
!      ratraw(ircrgp) = 0.0d0  ; dratrawdt(ircrgp) = 0.0d0  ; dratrawdd(ircrgp) = 0.0d0


! irx1 
!           ratraw(ircrag) = 0.0d0
!           ratraw(irfega) = 0.0d0
!      ratraw(ircrap) = 0.0d0  ; dratrawdt(ircrap) = 0.0d0  ; dratrawdd(ircrap) = 0.0d0
!      ratraw(irmnpa) = 0.0d0  ; dratrawdt(irmnpa) = 0.0d0  ; dratrawdd(irmnpa) = 0.0d0
!      ratraw(irmnpg) = 0.0d0  ; dratrawdt(irmnpg) = 0.0d0  ; dratrawdd(irmnpg) = 0.0d0
!      ratraw(irfegp) = 0.0d0  ; dratrawdt(irfegp) = 0.0d0  ; dratrawdd(irfegp) = 0.0d0


! iry1
!           ratraw(irfeag) = 0.0d0
!           ratraw(irniga) = 0.0d0
!      ratraw(irfeap) = 0.0d0  ; dratrawdt(irfeap) = 0.0d0  ; dratrawdd(irfeap) = 0.0d0
!      ratraw(ircopa) = 0.0d0  ; dratrawdt(ircopa) = 0.0d0  ; dratrawdd(ircopa) = 0.0d0
!      ratraw(ircopg) = 0.0d0  ; dratrawdt(ircopg) = 0.0d0  ; dratrawdd(ircopg) = 0.0d0
!      ratraw(irnigp) = 0.0d0  ; dratrawdt(irnigp) = 0.0d0  ; dratrawdd(irnigp) = 0.0d0



      return
      end subroutine aprox13rat

      subroutine screen_aprox13(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine computes the screening factors
! and applies them to the raw reaction rates,
! producing the final reaction rates used by the
! right hand sides and jacobian matrix elements

! this routine assumes screen_on = 1 or = 0 has been set at a higher
! level presumably in the top level driver


! declare the pass
      real*8  :: y(*)

! declare
      integer :: i,jscr,init=1
      real*8  :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                 sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar, &
                 denom,denomdt,denomdd,xx,zz


! initialize
       ratdum(1:nrat) = ratraw(1:nrat)  ; dratdumdt(1:nrat) = dratrawdt(1:nrat)  ; dratdumdd(1:nrat) = dratrawdd(1:nrat)
       scfac(1:nrat)  = 1.0d0  ; dscfacdt(1:nrat)  = 0.0d0  ; dscfacdd(1:nrat)  = 0.0d0


! if screening is on
      if (screen_on .ne. 0) then


! with the passed composition compute abar, zbar and z**2 bar

       call azbar_simple13(y,aion,zion,ionbeg,ionmax,abar,zbar,ye)
       z2bar  = abar * sum( zion(1:ionmax) * zion(1:ionmax) * y(1:ionmax) )


! first the always fun triple alpha and its inverse
      jscr = 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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


! c12 to o16 jscr = 3
! c12(a,g)o16
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircag)     = ratraw(ircag) * sc1a
      dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
      dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

      scfac(ircag)  = sc1a  ; dscfacdt(ircag)   = sc1adt  ; dscfacdd(ircag)   = sc1add

      ratdum(iroga)     = ratraw(iroga) * sc1a
      dratdumdt(iroga)  = dratrawdt(iroga)*sc1a + ratraw(iroga)*sc1adt
      dratdumdd(iroga)  = dratrawdd(iroga)*sc1a + ratraw(iroga)*sc1add

      scfac(iroga)  = sc1a  ; dscfacdt(iroga)   = sc1adt  ; dscfacdd(iroga)   = sc1add


! c12 + c12 jscr - 4
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ic12),aion(ic12), &
                   jscr,init,sc1a,sc1adt,sc1add)


      ratdum(ir1212)    = ratraw(ir1212) * sc1a
      dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
      dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

      scfac(ir1212)     = sc1a  ; dscfacdt(ir1212)  = sc1adt  ; dscfacdd(ir1212)  = sc1add



! c12 + o16 jscr = 5
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1216)    = ratraw(ir1216) * sc1a
      dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
      dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

      scfac(ir1216)     = sc1a  ; dscfacdt(ir1216)  = sc1adt  ; dscfacdd(ir1216)  = sc1add


! o16 + o16 jscr = 6
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)


      ratdum(ir1616)    = ratraw(ir1616) * sc1a
      dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
      dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

      scfac(ir1616)     = sc1a  ; dscfacdt(ir1616)  = sc1adt  ; dscfacdd(ir1616)  = sc1add



! o16 to ne20 jscr = 7
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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


! ne20 to mg24 jscr = 8
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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



! mg24 to si28 jscr = 9
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   13.0d0,27.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! al27(p,g)si28 jscr = 10
      ratdum(iralpg)    = ratraw(iralpg) * sc1a
      dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
      dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

      scfac(iralpg)     = sc1a  ; dscfacdt(iralpg)  = sc1adt  ; dscfacdd(iralpg)  = sc1add

      ratdum(irsigp)    = ratraw(irsigp) * sc1a
      dratdumdt(irsigp) = dratrawdt(irsigp)*sc1a + ratraw(irsigp)*sc1adt
      dratdumdd(irsigp) = dratrawdd(irsigp)*sc1a + ratraw(irsigp)*sc1add

      scfac(irsigp)     = sc1a  ; dscfacdt(irsigp)  = sc1adt  ; dscfacdd(irsigp)  = sc1add


! si28 to s32 jscr = 11
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   15.0d0,31.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! p31(p,g)s32 jscr = 12
      ratdum(irppg)     = ratraw(irppg) * sc1a
      dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
      dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

      scfac(irppg)      = sc1a  ; dscfacdt(irppg)   = sc1adt  ; dscfacdd(irppg)   = sc1add

      ratdum(irsgp)     = ratraw(irsgp) * sc1a
      dratdumdt(irsgp)  = dratrawdt(irsgp)*sc1a + ratraw(irsgp)*sc1adt
      dratdumdd(irsgp)  = dratrawdd(irsgp)*sc1a + ratraw(irsgp)*sc1add

      scfac(irsgp)      = sc1a  ; dscfacdt(irsgp)   = sc1adt  ; dscfacdd(irsgp)   = sc1add



! s32 to ar36 jscr = 13
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   17.0d0,35.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! cl35(p,g)ar36 jscr = 14
      ratdum(irclpg)    = ratraw(irclpg) * sc1a
      dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
      dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

      scfac(irclpg)     = sc1a  ; dscfacdt(irclpg)  = sc1adt  ; dscfacdd(irclpg)  = sc1add

      ratdum(irargp)    = ratraw(irargp) * sc1a
      dratdumdt(irargp) = dratrawdt(irargp)*sc1a + ratraw(irargp)*sc1adt
      dratdumdd(irargp) = dratrawdd(irargp)*sc1a + ratraw(irargp)*sc1add

      scfac(irargp)     = sc1a  ; dscfacdt(irargp)  = sc1adt  ; dscfacdd(irargp)  = sc1add



! ar36 to ca40 jscr = 15
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   19.0d0,39.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! k39(p,g)ca40 jscr = 16
      ratdum(irkpg)     = ratraw(irkpg) * sc1a
      dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
      dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

      scfac(irkpg)      = sc1a  ; dscfacdt(irkpg)   = sc1adt  ; dscfacdd(irkpg)   = sc1add

      ratdum(ircagp)     = ratraw(ircagp) * sc1a
      dratdumdt(ircagp)  = dratrawdt(ircagp)*sc1a  + ratraw(ircagp)*sc1adt
      dratdumdd(ircagp)  = dratrawdd(ircagp)*sc1a  + ratraw(ircagp)*sc1add

      scfac(ircagp)      = sc1a  ; dscfacdt(ircagp)   = sc1adt  ; dscfacdd(ircagp)   = sc1add



! ca40 to ti44 jscr = 17
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   21.0d0,43.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! sc43(p,g)ti44 jscr = 18
      ratdum(irscpg)    = ratraw(irscpg) * sc1a
      dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
      dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

      scfac(irscpg)     = sc1a  ; dscfacdt(irscpg)  = sc1adt  ; dscfacdd(irscpg)  = sc1add

      ratdum(irtigp)    = ratraw(irtigp) * sc1a
      dratdumdt(irtigp) = dratrawdt(irtigp)*sc1a + ratraw(irtigp)*sc1adt
      dratdumdd(irtigp) = dratrawdd(irtigp)*sc1a + ratraw(irtigp)*sc1add

      scfac(irtigp)     = sc1a  ; dscfacdt(irtigp)  = sc1adt  ; dscfacdd(irtigp)  = sc1add



! ti44 to cr48 jscr = 19
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   23.0d0,47.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! v47(p,g)cr48 jscr = 20
      ratdum(irvpg)     = ratraw(irvpg) * sc1a
      dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
      dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

      scfac(irvpg)      = sc1a  ; dscfacdt(irvpg)   = sc1adt  ; dscfacdd(irvpg)   = sc1add

      ratdum(ircrgp)     = ratraw(ircrgp) * sc1a
      dratdumdt(ircrgp)  = dratrawdt(ircrgp)*sc1a  + ratraw(ircrgp)*sc1adt
      dratdumdd(ircrgp)  = dratrawdd(ircrgp)*sc1a  + ratraw(ircrgp)*sc1add

      scfac(ircrgp)      = sc1a  ; dscfacdt(ircrgp)   = sc1adt  ; dscfacdd(ircrgp)   = sc1add



! cr48 to fe52 jscr = 21
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   25.0d0,51.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

! mn51(p,g)fe52 jscr = 22
      ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
      dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
      dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

      scfac(irmnpg)     = sc1a  ; dscfacdt(irmnpg)  = sc1adt  ; dscfacdd(irmnpg)  = sc1add

      ratdum(irfegp)    = ratraw(irfegp) * sc1a
      dratdumdt(irfegp) = dratrawdt(irfegp)*sc1a + ratraw(irfegp)*sc1adt
      dratdumdd(irfegp) = dratrawdd(irfegp)*sc1a + ratraw(irfegp)*sc1add

      scfac(irfegp)     = sc1a  ; dscfacdt(irfegp)  = sc1adt  ; dscfacdd(irfegp)  = sc1add



! fe52 to ni56 jscr = 23
      jscr = jscr + 1
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
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
      call screen_13_5(btemp,bden,zbar,abar,z2bar, &
                   27.0d0,55.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)


! co55(p,g)ni56 jscr = 24
      ratdum(ircopg)    = ratraw(ircopg) * sc1a
      dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
      dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

      scfac(ircopg)     = sc1a  ; dscfacdt(ircopg)  = sc1adt  ; dscfacdd(ircopg)  = sc1add

      ratdum(irnigp)    = ratraw(irnigp) * sc1a
      dratdumdt(irnigp) = dratrawdt(irnigp)*sc1a + ratraw(irnigp)*sc1adt
      dratdumdd(irnigp) = dratrawdd(irnigp)*sc1a + ratraw(irnigp)*sc1add

      scfac(irnigp)     = sc1a  ; dscfacdt(irnigp)  = sc1adt  ; dscfacdd(irniga)  = sc1add



! reset the screen initialization flag
      init = 0

! end of screening block if
      end if


! now form those lovely dummy proton link rates

! mg24(a,p)27al(p,g)28si
      ratdum(irr1)     = 0.0d0  ; dratdumdt(irr1)  = 0.0d0  ; dratdumdd(irr1)  = 0.0d0
      denom    = ratdum(iralpa) + ratdum(iralpg)
      denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
      denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
      if (denom .gt. 1.0d-30) then
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
      if (denom .gt. 1.0d-30) then
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
      if (denom .gt. 1.0d-30) then
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
      if (denom .gt. 1.0d-30) then
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
      if (denom .gt. 1.0d-30) then
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
      if (denom .gt. 1.0d-30) then
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
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
      endif

! fe52(a,p)co55(p,g)ni56
      ratdum(iry1)    = 0.0d0  ; dratdumdt(iry1) = 0.0d0  ; dratdumdd(iry1) = 0.0d0
      denom    = ratdum(ircopa) + ratdum(ircopg)
      denomdt  = dratdumdt(ircopa) + dratdumdt(ircopg)
      denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
      if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(iry1)    = ratdum(ircopa)*zz
       dratdumdt(iry1) = (dratdumdt(ircopa) - ratdum(iry1)*denomdt)*zz
       dratdumdd(iry1) = (dratdumdd(ircopa) - ratdum(iry1)*denomdd)*zz
      end if



! debugs
!      do i=1,nrat
!       if (ratdum(i) .lt. 0.0) then
!        write(6,110) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
! 110    format(1x,i4,' ',a,' ',1p3e12.4)
!        stop 'negative rate'
!       end if
!      enddo
!      read(5,*)

!      write(6,*) btemp,bden
!      write(6,*) abar,zbar
!      do i=1,nrat
!       write(6,111) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
!       write(6,111) i,ratnam(i),ratdum(i),dratdumdt(i),dratdumdd(i)
! 111   format(1x,i4,' ',a,' ',1p3e14.6)
!      enddo
!      read(5,*)

      return
      end subroutine screen_aprox13

      subroutine screen_13_5(temp,den,zbar,abar,z2bar, &
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
       stop 'jscreen > nscreen_max in screen_13_5'
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
    end subroutine screen_13_5
    
      subroutine rhs13(y,rate,dydt,deriva)
      include 'implno.dek'
      include 'network.dek'

! evaluates the right hand side of the aprox13 odes

! declare the pass
      logical   ::  deriva
      real*8    ::  y(*),rate(*),dydt(*)



! local variables
      integer           :: i
      real*8, parameter :: sixth = 1.0d0/6.0d0

! quad precision dydt sums
      real*8           :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17
      real*8           :: qray(abignet)



! zero the odes
      dydt(1:neqs) = 0.0d0
      qray(1:neqs) = 0.0q0


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
       a16 = -y(ihe4)  * y(ife52) * rate(irfeap)*(1.0d0-rate(iry1))
       a17 =  y(ini56) * rate(irnigp) * rate(iry1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 & 
                                + a11 + a12 + a13 + a14 + a15 + a16 + a17

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
       a9  = -y(ihe4)*y(ife52) * rate(irfeap)*(1.0d0 - ratdum(iry1))
       a10 =  y(ihe4)*y(ife52) * ratdum(irfeap)*rate(iry1)
       a11 =  y(ini56) * ratdum(irnigp) * rate(iry1) 
       a12 =  y(ini56) * rate(irnigp) * ratdum(iry1)

       qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12
      end if


! c12 reactions
      a1 = -y(ic12) * y(ic12) * rate(ir1212) 
      a2 = -y(ic12) * y(io16) * rate(ir1216) 
      a3 =  sixth * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
      a4 = -y(ic12) * rate(irg3a) 
      a5 = -y(ic12) * y(ihe4) * rate(ircag) 
      a6 =  y(io16) * rate(iroga) 

      qray(ic12) =  qray(ic12) + a1 + a2 + a3 + a4 + a5 + a6


! o16 reactions
      a1 = -y(ic12) * y(io16) * rate(ir1216) 
      a2 = -y(io16) * y(io16) * rate(ir1616) 
      a3 =  y(ic12) * y(ihe4) * rate(ircag) 
      a4 = -y(io16) * y(ihe4) * rate(iroag) 
      a5 = -y(io16) * rate(iroga) 
      a6 =  y(ine20) * rate(irnega) 

      qray(io16) =  qray(io16) + a1 + a2 + a3 + a4 + a5 + a6 


! ne20 reactions
      a1 =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
      a2 =  y(io16) * y(ihe4) * rate(iroag) 
      a3 = -y(ine20) * y(ihe4) * rate(irneag) 
      a4 = -y(ine20) * rate(irnega) 
      a5 =  y(img24) * rate(irmgga) 

      qray(ine20) =  qray(ine20) + a1 + a2 + a3 + a4 + a5


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


! fe52 reactions
      a1 =  y(icr48) * y(ihe4) * rate(ircrag) 
      a2 = -y(ife52) * y(ihe4) * rate(irfeag) 
      a3 = -y(ife52) * rate(irfega) 
      a4 =  y(ini56) * rate(irniga)

      qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4

      if (.not.deriva) then
       a1 =  y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
       a2 = -y(ife52) * rate(irx1) * rate(irfegp) 
       a3 = -y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) 
       a4 =  y(ini56) * rate(iry1) * rate(irnigp)

       qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4

      else
       a1 =  y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
       a2 = -y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
       a3 = -y(ife52) * ratdum(irx1) * rate(irfegp) 
       a4 = -y(ife52) * rate(irx1) * ratdum(irfegp)
       a5 = -y(ife52)*y(ihe4) * rate(irfeap)*(1.0d0-ratdum(iry1))
       a6 =  y(ife52)*y(ihe4) * ratdum(irfeap)*rate(iry1)
       a7 =  y(ini56) * ratdum(iry1) * rate(irnigp)
       a8 =  y(ini56) * rate(iry1) * ratdum(irnigp)

       qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
      end if


! ni56 reactions
      a1 =  y(ife52) * y(ihe4) * rate(irfeag) 
      a2 = -y(ini56) * rate(irniga) 

      qray(ini56) =  qray(ini56) + a1 + a2 

      if (.not.deriva) then
       a1 =  y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) 
       a2 = -y(ini56) * rate(iry1) * rate(irnigp)

       qray(ini56) =  qray(ini56) + a1 + a2 

      else
       a1 =  y(ife52)*y(ihe4) * rate(irfeap)*(1.0d0-ratdum(iry1))
       a2 = -y(ife52)*y(ihe4) * ratdum(irfeap)*rate(iry1)
       a3 = -y(ini56) * ratdum(iry1) * rate(irnigp)
       a4 = -y(ini56) * rate(iry1) * ratdum(irnigp)

       qray(ini56) =  qray(ini56) + a1 + a2 + a3 + a4
      end if


! now set the double precision return argument dydt

      dydt(1:ionmax) = qray(1:ionmax)


      return
      end subroutine rhs13
      
      subroutine rate13_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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






      subroutine rate13_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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







      subroutine rate13_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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






      subroutine rate13_c12c12npa(temp,den, &
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






      subroutine rate13_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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




      subroutine rate13_c12o16npa(temp,den, &
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








      subroutine rate13_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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






      subroutine rate13_o16o16npad(temp,den, &
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

       call fowthrsh13(xxt,thrs,dthrs)

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






      subroutine fowthrsh13(x,thrs,dthrs)
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





      subroutine rate13_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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




      subroutine rate13_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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




      subroutine rate13_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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






      subroutine rate13_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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







      subroutine rate13_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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






      subroutine rate13_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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






      subroutine rate13_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! sc43(p,g)ti44
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






      subroutine rate13_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
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





      subroutine rate13_ti44ap(temp,den, &
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





      subroutine rate13_v47pg(temp,den, &
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






      subroutine rate13_cr48ag(temp,den, &
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





      subroutine rate13_cr48ap(temp,den, &
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





      subroutine rate13_mn51pg(temp,den, &
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






      subroutine rate13_fe52ag(temp,den, &
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





      subroutine rate13_fe52ap(temp,den, &
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





      subroutine rate13_co55pg(temp,den, &
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

      subroutine azbar_simple13(ymass,aion,zion,ionbeg,ionend, &
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
      end subroutine azbar_simple13
