! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team

!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************  
      subroutine mesareact45(nd, init, pa, ydot)
       use rates_def
       use chem_def
       use net_def
       use net_lib
       use chem_lib
!     Input
       integer   :: nd, init
       real*8    :: pa(*), ydot(50)
!      locals
       integer, SAVE :: screening_mode, i
       real(dp), SAVE :: xa(50), T, logT, Rho, logRho, eta, d_eta_dlnT, d_eta_dlnRho, &
            d_eps_nuc_dx(50), dabar_dx(50), dzbar_dx(50), dmc_dx(50), &
            weak_rate_factor, xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, xsum, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, &
            eps_nuc_categories(num_categories), eps_neu_total, &
            dxdt(50), d_dxdt_dRho(50), d_dxdt_dT(50), &
            d_dxdt_dx(50,50)
       real(dp), SAVE, target :: rate_factors_ary(370)
       real(dp), SAVE, pointer, dimension(:) :: rate_factors
       logical :: SAVE, skip_jacobian
       type (Net_Info), SAVE :: n
       character (len=80), SAVE :: string
       integer, SAVE :: ierr, handle, species 
       integer, SAVE, pointer :: chem_id(:), net_iso(:)
       character (len=100), SAVE :: net_file
!       character (len=64) :: mesa_dir       
       integer, SAVE :: num_reactions
!       integer  :: ierr                       
       if (init .eq. 1) then
!         mesa_dir = '../..'
        mesa_dir = ''
        ! choose the network to use
        net_file = 'mesa_45.net'
        ! initialize
        ierr = 0
        call initialize(mesa_dir, ierr)
!       set up the network         
!        call setup_net(net_file, handle, &
!          species, chem_id, net_iso, num_reactions, ierr)
        ierr = 0
        handle = alloc_net_handle(ierr)
        call net_start_def(handle, ierr)
        write(*,*) 'load ' // trim(net_file)
        call read_net_file(net_file, handle, ierr)
        call net_finish_def(handle, ierr)
        call net_setup_tables(handle, '', ierr)
        species = net_num_isos(handle, ierr)
        call get_chem_id_table_ptr(handle, chem_id, ierr)
        call get_net_iso_table_ptr(handle, net_iso, ierr)
        num_reactions = net_num_reactions(handle, ierr)
!       set some pointers and options
        ierr = 0
        rate_factors => rate_factors_ary
        eta = 0
        rate_factors(:) = 1
        weak_rate_factor = 1
        screening_mode = extended_screening
        skip_jacobian = .false.
        return
       endif
!      Set temperature, density and pressure
       temp = pa(50+nd)
       den = pa(1)
       pres = pa(nd+2)
!      Set abundances
       do i=1,species
          xa(i) = pa(2+nd+i)
       enddo
       logT   = log10(temp)
       logRho = log10(den)       
!      Get some composition variables
       call composition_info( &
            species, chem_id, xa, xh, xhe, z, abar, zbar, z2bar, z53bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
!      Get rates    
       call net_get( &
            handle, skip_jacobian, n, species, num_reactions, &
            xa, T, logT, Rho, logRho, & 
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, & 
            std_reaction_Qs, std_reaction_neuQs, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
            screening_mode, &     
            eps_nuc_categories, eps_neu_total, & 
            ierr)
      do i=1,species
         ydot(i) = dxdt(i)
      enddo
      ydot(46) = eps_nuc
      return
    end subroutine mesareact45
    
    subroutine initialize(mesa_dir, ierr)
      use const_lib, only: const_init
      use math_lib
      use chem_lib, only: chem_init
      use rates_lib, only: rates_init, rates_warning_init
      use net_lib, only : net_init
      character (len=*), intent(in) :: mesa_dir
      integer, intent(out) :: ierr
      ierr = 0
      call math_init()
      call const_init(mesa_dir,ierr)     
      call chem_init('isotopes.data', ierr)
      call rates_init('reactions.list', '', 'rate_tables', .false., .false.,&
           '', '', '',  ierr)
      call rates_warning_init(.true., 10d0)
      call net_init(ierr)

    end subroutine initialize

    subroutine setup_net(net_file, handle, &
         species, chem_id, net_iso, num_reactions, ierr)
      use net_lib
      use rates_def, only: rates_reaction_id_max
      character (len=100), intent(in) :: net_file
      integer, pointer :: chem_id(:), net_iso(:) ! set, but not allocated
      integer, intent(out) :: handle, species, num_reactions, ierr
      ierr = 0
      handle = alloc_net_handle(ierr)
      call net_start_def(handle, ierr)
      call read_net_file(net_file, handle, ierr)
      call net_finish_def(handle, ierr)
      call net_setup_tables(handle, '', ierr)
      species = net_num_isos(handle, ierr)
      call get_chem_id_table_ptr(handle, chem_id, ierr)   
      call get_net_iso_table_ptr(handle, net_iso, ierr)
      num_reactions = net_num_reactions(handle, ierr)
    end subroutine setup_net







