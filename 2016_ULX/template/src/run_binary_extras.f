! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use crlibm_lib
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pinters to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: beta
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
         extras_binary_startup = keep_going
      end function  extras_binary_startup
      
      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: r_isco, Z1, Z2
         call binary_ptr(binary_id, b, ierr)

         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (b% s1% model_number > 1000 .and. b% ignore_rlof) then
               extras_binary_finish_step = terminate
               write(*,*) "Terminate due to pre-MS evolution taking too long"
               return
         end if

         if (b% ignore_rlof .and. &
            abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) < 0.005 .and. &
            b% s1% star_age > 1d2) then
            ! if here, primary reached thermal equilibrium (reached ZAMS), so activate RLOF
            ! this is the amount of overflow of a q=1 system at L2, anything more than this
            ! is too much
            if (b% rl_relative_gap(1) > 0.320819224) then
               extras_binary_finish_step = terminate
               write(*,*) "Terminate due to overflow of L2 at ZAMS"
               return
            else if (b% rl_relative_gap(1) > 0d0) then
               write(*,*) "model is overflowing at ZAMS"
               ! initially overflowing system will evolve rapidly to q=1, limit mdot
               ! for contact system to work and remove hard limit on mdot change
               b% max_implicit_abs_mdot = 1d-3
               b% s1% delta_mdot_limit = -1
               b% s2% delta_mdot_limit = -1
               b% s1% delta_mdot_hard_limit = -1
               b% s2% delta_mdot_hard_limit = -1
            else
               write(*,*) "model is not overflowing at ZAMS"
            end if
            b% ignore_rlof = .false.
            b% terminate_if_L2_overflow = .true.
            b% s1% max_mdot_redo_cnt = 100
            b% s2% max_mdot_redo_cnt = 100
            write(*,*) "Engage RLOF!"
         else if (b% ignore_rlof .and. &
            (abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) > 0.005 .or. &
            b% s1% star_age < 1d2)) then
            ! if here, still not in ZAMS, keep period fixed
            b% period = b% initial_period_in_days*(24d0*60d0*60d0)
            b% separation = &
               pow_cr((b% s1% cgrav(1)*(b% m(1)+b% m(2)))*(b% period/(2*pi))**2,1d0/3d0)
            b% angular_momentum_j = b% m(1) * b% m(2) * sqrt( b% s1% cgrav(1) *&
               b% separation / (b% m(1) + b% m(2)) )
            b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
            b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
            b% rl_relative_gap(1) = (b% r(1) - b% rl(1)) / b% rl(1) ! gap < 0 means out of contact 
            b% rl_relative_gap(2) = (b% r(2) - b% rl(2)) / b% rl(2) ! gap < 0 means out of contact
         end if

         !relax fr for deep contact
         if (.not. b% ignore_rlof) then
            if (b% rl_relative_gap(1) > 0.1) then
               b% fr = 1d0
            else
               b% fr = 0.02d0
            end if
         else
            b% fr = 1d0
         end if

         ! check if stars are evolving homogeneously
         if (b% s1% center_h1 > 1d-3) then
            if (b% s1% center_he4 - b% s1% surface_he4 > 0.2) then
               extras_binary_finish_step = terminate
               write(*,*) "Terminate due to primary not evolving homogeneously"
               return
            end if
         end if
         !if (b% s2% center_h1 > 1d-3) then
         !   if (b% s2% center_he4 - b% s2% surface_he4 > 0.2) then
         !      extras_binary_finish_step = terminate
         !      write(*,*) "Terminate due to secondary not evolving homogeneously"
         !      return
         !   end if
         !end if

         ! turn a star into point mass once it depletes helium
         ! terminate once both deplete helium
         if (b% point_mass_i == 0) then
            ! if modelling q=1 system then just terminate at He depletion
            if (b% s1% center_he4 < 1d-5) then
               ! terminate when the result is a PISNe
               if (b% s1% star_mass > 60 .and. b% s1% star_mass < 130) then
                  extras_binary_finish_step = terminate
                  write(*,*) "Terminate due to primary resulting in PISN"
                  return
               end if
               b% point_mass_i = 1
               b% d_i = 2
               b% a_i = 1
               b% s_donor => b% s2
               b% s_accretor => b% s1
               b% limit_retention_by_mdot_edd = .true.
               b% mdot_scheme = "roche_lobe"
               b% initial_bh_spin = b% s1% total_angular_momentum * clight &
                  / (b% m(1)**2 * standard_cgrav)
               write(*,*) "initial black hole spin is", b% initial_bh_spin
               b% initial_bh_spin = min(b% initial_bh_spin,1d0)
               Z1 = 1d0 + pow_cr(1d0 - b% initial_bh_spin**2,one_third) &
                  * (pow_cr(1d0 + b% initial_bh_spin,one_third) + pow_cr(1d0 - b% initial_bh_spin,one_third))
               Z2 = sqrt(3d0*b% initial_bh_spin**2 + Z1**2)
               r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
               b% eq_initial_bh_mass = b% m(b% point_mass_i) * sqrt(r_isco/6d0)
            else if (b% s2% center_he4 < 1d-5 .and. b% s2% center_c12 < 1d-3) then
               b% point_mass_i = 2
               b% d_i = 1
               b% a_i = 2
               b% s_donor => b% s1
               b% s_accretor => b% s2
               b% limit_retention_by_mdot_edd = .true.
               b% mdot_scheme = "roche_lobe"
            end if
         end if
         ! terminate when the other star is close to carbon depletion
         ! avoid going to zero carbon to save myself random convergence issues
         if (b% point_mass_i == 1) then
            if (b% s2% center_he4 < 1d-5 .and. b% s2% center_c12 < 1d-3) then
               extras_binary_finish_step = terminate
               write(*,*) "Terminate due to carbon depletion"
               return
            end if
         else if (b% point_mass_i == 2) then
            if (b% s1% center_he4 < 1d-5) then
               extras_binary_finish_step = terminate
               write(*,*) "Terminate due to helium depletion"
               return
            end if
         end if

         ! stop at different phases of mass transfer
         if (.not. b% ignore_rlof .and. &
            (b% rl_relative_gap(1) >-1d-3 .or. b% rl_relative_gap(2) >-1d-3)) then
            if(b% point_mass_i == 0 .and. b% d_i == 1) then
               extras_binary_finish_step = terminate
               write(*,*) "Mass transfer from primary to secondary"
               return
            end if
            if(b% point_mass_i == 0 .and. b% d_i == 2) then
               if (b% s1% center_h1 < 1d-6) then
                  extras_binary_finish_step = terminate
                  write(*,*) "Mass transfer from secondary to primary, primary off MS"
                  return
               else
                  extras_binary_finish_step = terminate
                  write(*,*) "Mass transfer from secondary to primary, primary on the MS"
                  return
               end if
            end if
            if(b% point_mass_i == 1) then
               if (b% s2% center_he4 < 1d-6) then
                  write(*,*) "Mass transfer from secondary to point mass, case BB!"
                  if (b% s2% surface_h1 < 1d-6) then
                     extras_binary_finish_step = terminate
                     write(*,*) "Terminate due to deep case BB!"
                     return
                  end if
               else if (b% s2% center_h1 < 1d-6) then
                  write(*,*) "Mass transfer from secondary to point mass, secondary off MS"
               else
                  write(*,*) "Mass transfer from secondary to point mass, secondary on the MS"
               end if
            end if
            if(b% point_mass_i == 2) then
               if (b% s2% center_h1 < 1d-6) then
                  write(*,*) "Mass transfer from primary to point mass, primary off MS"
               else
                  write(*,*) "Mass transfer from primary to point mass, primary on the MS"
               end if
            end if
         end if

         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow_cr(m1/m2,one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
      end function eval_rlobe
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
