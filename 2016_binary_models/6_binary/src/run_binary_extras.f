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
      use utils_lib
      
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
         
         b% use_other_jdot_ls = .true.
         b% other_jdot_ls => correct_jdot_ls

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls

      subroutine correct_jdot_ls(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: delta_J
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_ls = 0
         ! ignore in first step, or if not doing rotation
         if (b% doing_first_model_of_run) &
             return
         ! bulk change in spin angular momentum takes tides into account
         delta_J = b% s_donor% total_angular_momentum_old - &
             b% s_donor% total_angular_momentum
         ! ignore angular momentum lost through winds
         if (b% s_donor% mstar_dot < 0) &
             delta_J = delta_J - b% s_donor% angular_momentum_removed * &
                 abs(b% mdot_system_wind(b% d_i) / b% s_donor% mstar_dot)
         ! Repeat for accretor
         if (b% point_mass_i == 0) then
             delta_J = delta_J + b% s_accretor% total_angular_momentum_old - &
                 b% s_accretor% total_angular_momentum
             if (b% s_accretor% mstar_dot < 0) then
                 ! all AM lost from the accretor is lost from the system
                 delta_J = delta_J - b% s_accretor% angular_momentum_removed
             end if
         end if
         b% jdot_ls = delta_J / b% s_donor% dt
      end subroutine correct_jdot_ls

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 3
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

         names(1) = "mdot_limit_low"
         names(2) = "mdot_limit_high"
         names(3) = "ignore_rlof"

         if (b% point_mass_i == 0) then
            vals(1) = -7.196 &
               + safe_log10_cr(b% s2% photosphere_r) &
               + safe_log10_cr(b% s2% photosphere_L + b% s1% photosphere_L) &
               - safe_log10_cr(b% s2% mstar/Msun)

            vals(2) = -7.196 &
               + safe_log10_cr(b% rl(2)/Rsun) &
               + safe_log10_cr(b% s2% photosphere_L + b% s1% photosphere_L) &
               - safe_log10_cr(b% s2% mstar/Msun)
         else
            vals(1) = 1d99
            vals(2) = 1d99
         end if

         if (b% ignore_rlof) then
            vals(3) = 1
         else
            vals(3) = 0
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
         
         if (.not. restart) then
            ! ixtra values are used as counters for mass transfer events (cases A,B,C)
            b% s1% ixtra1 = 0
            b% s1% ixtra2 = 0
            b% s1% ixtra3 = 0
            ! lxtra are used to register events
            b% s1% lxtra1 = .false. ! OC Ne ignition for primary.
            b% s1% lxtra16 = .false. ! OC Ne ignition for secondary.
            b% s1% lxtra2 = .false. ! Primary depleted carbon.
            b% s1% lxtra19 = .false. ! Primary depleted helium with mass > 14msun.
            b% s1% lxtra3 = .false. ! Secondary depleted carbon.
            b% s1% lxtra4 = .false. ! Primary terminated last.
            b% s1% lxtra5 = .false. ! Secondary terminated last.
            b% s1% lxtra6 = .false. ! ZAMS L2 overflow.
            b% s1% lxtra7 = .false. ! ZAMS RLOF.
            b% s1% lxtra8 = .false. ! L2OF after ZAMS.
            b% s1% lxtra9 = .false. ! Contact.
            b% s1% lxtra10 = .false. ! Contact with inverse mass transfer.
            b% s1% lxtra11 = .false. ! Reached lower mdot limit.
            b% s1% lxtra12 = .false. ! Reached upper mdot limit.
            b% s1% lxtra13 = .false. ! Inverse mass transfer past secondary H depletion, primary in MS.
            b% s1% lxtra18 = .false. ! Inverse mass transfer past secondary H depletion, primary off MS.
            b% s1% lxtra14 = .false. ! Inverse mass transfer before secondary H depletion, primary off MS.
            b% s1% lxtra17 = .false. ! Terminated due to mass transfer in q<=0.25 system
            b% s1% lxtra20 = .false. ! Reached maximum mass transfer rate, assume merger
            b% s1% lxtra21 = .false. ! System is Darwin unstable
            b% s1% lxtra15 = .false. ! Normal termination
         end if

         if (restart) then
            if (b% s1% center_h1 < 1d-2) then
               b% s1% fp_min = 1d0-0.25*b% s1% center_h1/1d-2
               b% s1% ft_min = 1d0-0.05*b% s1% center_h1/1d-2
            else
               b% s1% fp_min = 0.75d0
               b% s1% ft_min = 0.95d0
            end if
         end if
         
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
         integer :: ierr, star_id, i
         real(dp) :: q, mdot_limit_low, mdot_limit_high, &
            center_h1, center_h1_old, center_he4, center_he4_old
         logical :: is_ne_biggest

         extras_binary_finish_step = keep_going

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (b% s1% model_number > 1000 .and. b% ignore_rlof) then
               write(*,*) "Terminate due to pre-MS evolution taking too long"
               extras_binary_finish_step = terminate
               write(*,*) "TERMINATING: overflow of L2 at ZAMS"
               b% s1% lxtra6 = .true.
               call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
               call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
         end if

         if (b% ignore_rlof .and. &
            abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) < 0.005 .and. &
            b% s1% star_age > 1d2) then
            ! if here, primary reached thermal equilibrium (reached ZAMS), so activate RLOF
            ! this is the amount of overflow of a q=1 system at L2, anything more than this
            ! is too much
            if (b% m(1) > b% m(2)) then
               q = b% m(2) / b% m(1)
               star_id = 2
            else
               q = b% m(1) / b% m(2)
               star_id = 1
            end if
            if (b% rl_relative_gap(star_id) > 0.29858997d0*atan_cr(1.83530121d0*pow_cr(q,0.39661426d0))) then
               extras_binary_finish_step = terminate
               write(*,*) "TERMINATING: overflow of L2 at ZAMS"
               b% s1% lxtra6 = .true.
               call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
               call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
               return
            else if (b% separation < &
                  sqrt(3*(b% s1% total_angular_momentum / b% s1% omega_avg_surf &
                   + b% s2% total_angular_momentum / b% s2% omega_avg_surf) &
                   / (b% m(1) * b% m(2) /(b% m(1) + b% m(2))))) then
               extras_binary_finish_step = terminate
               write(*,*) "TERMINATING: system is Darwin unstable"
               b% s1% lxtra21 = .true.
               call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
               call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
               return
            else if (b% rl_relative_gap(1) > 0d0) then
               write(*,*) "model is overflowing at ZAMS"
               b% s1% lxtra7 = .true.
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
            b% s1% max_mdot_redo_cnt = 100
            b% s2% max_mdot_redo_cnt = 100
            write(*,*) "Engage RLOF!"
            ! save ZAMS profiles
            call star_write_profile_info(b% s1% id, "LOGS1/prof_1ZAMS.data", b% s1% id, ierr)
            if (ierr /= 0) return ! failure in profile
            call star_write_profile_info(b% s2% id, "LOGS2/prof_1ZAMS.data", b% s2% id, ierr)
            if (ierr /= 0) return ! failure in profile
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
            if (b% rl_relative_gap(1) > 0.005) then
               b% fr = 0.05d0
            else
               b% fr = 0.01d0
            end if
            if (b% rl_relative_gap(1) > -1d-2 .and. b% m2 / b% m1 < 0.251) then
               extras_binary_finish_step = terminate
               write(*,*) "Terminated due to mass transfer in q_i<=0.25 system"
               b% s1% lxtra17 = .true.
               if (b% s1% center_h1 > 1d-6) then
                  b% s1% ixtra1 = b% s1% ixtra1 + 1
               else if (center_he4 > 1d-6) then
                  b% s1% ixtra2 = b% s1% ixtra2 + 1
               else
                  b% s1% ixtra3 = b% s1% ixtra3 + 1
               end if
               call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
               call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
               return
            end if
         else
            b% fr = 1d0
         end if

         if (b% point_mass_i == 0) then
            ! check for limits on mass loss
            mdot_limit_low = -7.196 &
               + safe_log10_cr(b% s2% photosphere_r) &
               + safe_log10_cr(b% s2% photosphere_L + b% s1% photosphere_L) &
               - safe_log10_cr(b% s2% mstar/Msun)

            mdot_limit_high = -7.196 &
               + safe_log10_cr(b% rl(2)/Rsun) &
               + safe_log10_cr(b% s2% photosphere_L + b% s1% photosphere_L) &
               - safe_log10_cr(b% s2% mstar/Msun)

            if (safe_log10_cr(abs(b% mdot_system_wind(1)+b% mdot_system_wind(2))/Msun*secyer) &
               > mdot_limit_low) then
               write(*,*) "mass loss above first limit"
               b% s1% lxtra11 = .true.
            end if

            if (safe_log10_cr(abs(b% mdot_system_wind(1)+b% mdot_system_wind(2))/Msun*secyer) &
               > mdot_limit_high) then
               extras_binary_finish_step = terminate
               write(*,*) "TERMINATING: mass loss above second limit"
               b% s1% lxtra12 = .true.
            end if

            ! check for inverse mass transfer at late stages
            if (b% s2% center_h1 < 1d-10 .and. b% d_i == 2 .and. b% rl_relative_gap(2) > -1d-2) then
               if (b% s1% center_h1 > 1d-10) then
                  extras_binary_finish_step = terminate
                  write(*,*) "TERMINATING: Inverse mass transfer past secondary H depletion, primary in MS"
                  b% s1% lxtra13 = .true.
               else
                  extras_binary_finish_step = terminate
                  write(*,*) "TERMINATING: Inverse mass transfer past secondary H depletion, primary off MS"
                  b% s1% lxtra18 = .true.
               end if
            else if (b% s1% center_h1 < 1d-10 .and. b% d_i == 2 .and. b% rl_relative_gap(2) > -1d-2) then
               extras_binary_finish_step = terminate
               write(*,*) "TERMINATING: Inverse mass transfer before secondary H depletion, primary off MS"
               b% s1% lxtra14 = .true.
            end if
         end if

         !remove gradL_composition term after MS, it can cause the convective helium core to recede
         if (b% s1% center_h1 < 1d-6) then
            b% s1% num_cells_for_smooth_gradL_composition_term = 0
         end if
         if (b% s2% center_h1 < 1d-6) then
            b% s2% num_cells_for_smooth_gradL_composition_term = 0
         end if

         !check if mass transfer rate reached maximun, assume merger if it happens
         if(abs(b% mtransfer_rate) >= b% max_implicit_abs_mdot*Msun/secyer) then
            extras_binary_finish_step = terminate
            write(*,*) "TERMINATING: Reached maximum mass transfer rate"
            b% s1% lxtra20 = .true.
         end if

         ! check for termination due to carbon depletion or off center neon ignition for primary
         if (b% point_mass_i /= 1) then
            ! ignore centrifugal foreces past the MS in the donor
            ! as these can cause problems during late inflated phases
            if (b% s1% center_h1 < 1d-2) then
               b% s1% fp_min = 1d0-0.25*b% s1% center_h1/1d-2
               b% s1% ft_min = 1d0-0.05*b% s1% center_h1/1d-2
            else
               b% s1% fp_min = 0.75d0
               b% s1% ft_min = 0.95d0
            end if
            !turn on MLT++ in the primary for the more massive systems.
            !For these I mainly care about getting their lifetimes right, so theres not much problem
            !in doing this.
            if(.not. b% s1% okay_to_reduce_gradT_excess .and. b% s1% he_core_mass > 13d0) then
               if(b% s1% power_he_burn > b% s1% power_h_burn) then
                  write(*,*) "Turn on MLT++ for the primary"
                  b% s1% okay_to_reduce_gradT_excess = .true.
                  b% s1% gradT_excess_age_fraction = 0.99
                  b% s1% gradT_excess_max_change = 0.001
               end if
            end if
            !write(*,*) "check", b% s1% fp_min, b% s1% ft_min
            if (b% s1% center_c12 < 5d-3 .and. b% s1% center_he4 < 1d-6) then
               b% s1% lxtra2 = .true.
               if (b% point_mass_i == 0) then
                  write(*,*) "Primary has depleted central carbon"
               else
                  extras_binary_finish_step = terminate
                  b% s1% lxtra4 = .true.
                  write(*,*) "Terminate due to primary depleting carbon (inverse sn?)"
               end if
            else if (b% s1% star_mass > 13d0 .and. b% s1% center_he4 < 1d-6) then
               ! if star is too massive when helium is depleted, terminate as well to avoid convergence issues
               b% s1% lxtra19 = .true.
               if (b% point_mass_i == 0) then
                  write(*,*) "Primary has depleted central helium (M > 13 Msun)"
               else
                  extras_binary_finish_step = terminate
                  b% s1% lxtra4 = .true.
                  write(*,*) "Terminate due to primary depleting helium (M > 13 Msun) (inverse sn?)"
               end if
            else
               ! check if neon is by far greatest source of energy
               is_ne_biggest = .true.
               do i=1, num_categories
                  if(i /= i_burn_ne .and. b% s1% L_by_category(i_burn_ne) < 10*b% s1% L_by_category(i)) then
                     is_ne_biggest = .false.
                     exit
                  end if
               end do
               if (is_ne_biggest .and. b% s1% max_eps_z_m/b% s1% xmstar > 0.01) then
                  b% s1% lxtra1 = .true.
                  if (b% point_mass_i == 0) then
                     write(*,*) "offcenter neon ignition for primary at q=",  b% s1% max_eps_z_m/b% s1% xmstar, &
                        b% s1% max_eps_z_m
                  else
                     write(*,*) "Terminate due to offcenter neon ignition for primary (inverse sn?) at q=", &
                     b% s1% max_eps_z_m/b% s1% xmstar, b% s1% max_eps_z_m
                     extras_binary_finish_step = terminate
                     b% s1% lxtra4 = .true.
                  end if
               end if
            end if
            if ((b% s1% lxtra2 .or. b% s1% lxtra1 .or. b% s1% lxtra19) .and. b% point_mass_i == 0) then
               b% point_mass_i = 1
               b% d_i = 2
               b% a_i = 1
               b% s_donor => b% s2
               b% s_accretor => b% s1
               b% period = 1d99
               b% mtransfer_rate = 0d0
               b% separation = &
                  pow_cr((b% s1% cgrav(1)*(b% m(1)+b% m(2)))*(b% period/(2*pi))**2,1d0/3d0)
               b% angular_momentum_j = b% m(1) * b% m(2) * sqrt( b% s1% cgrav(1) *&
                  b% separation / (b% m(1) + b% m(2)) )
               b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
               b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
               b% rl_relative_gap(1) = (b% r(1) - b% rl(1)) / b% rl(1) ! gap < 0 means out of contact 
               b% rl_relative_gap(2) = (b% r(2) - b% rl(2)) / b% rl(2) ! gap < 0 means out of contact
               b% do_tidal_sync = .false.
               b% mdot_scheme = "roche_lobe"
               ! turn off outer_xa smooth and premix_omega in secondary
               ! se these can cause trouble at late phases
               b% s2% premix_omega = .false.
               b% s2% smooth_outer_xa_big = -1d0
               b% s2% smooth_outer_xa_small = -1d0
            end if
         end if

         ! check for termination due to carbon depletion or off center neon ignition for secondary
         if (b% point_mass_i /= 2) then
            if (b% s2% center_c12 < 5d-3 .and. b% s2% center_he4 < 1d-6) then
               b% s1% lxtra3 = .true.
               if (b% point_mass_i == 0) then
                  write(*,*) "Secondary has depleted central carbon"
               else
                  extras_binary_finish_step = terminate
                  b% s1% lxtra5 = .true.
                  write(*,*) "Terminate due to secondary depleting carbon"
               end if
            else
               ! check if neon is by far greatest source of energy
               is_ne_biggest = .true.
               do i=1, num_categories
                  if(i /= i_burn_ne .and. b% s2% L_by_category(i_burn_ne) < 10*b% s2% L_by_category(i)) then
                     is_ne_biggest = .false.
                     exit
                  end if
               end do
               if (is_ne_biggest .and. b% s2% max_eps_z_m/b% s2% xmstar > 0.01) then
                  b% s1% lxtra16 = .true.
                  if (b% point_mass_i == 0) then
                     write(*,*) "offcenter neon ignition for secondary at q=",  b% s2% max_eps_z_m/b% s2% xmstar, &
                        b% s2% max_eps_z_m
                  else
                     write(*,*) "Terminate due to offcenter neon ignition for secondary at q=",  b% s2% max_eps_z_m/b% s2% xmstar, &
                        b% s2% max_eps_z_m
                     extras_binary_finish_step = terminate
                     b% s1% lxtra5 = .true.
                  end if
               end if
            end if
            if ((b% s1% lxtra3 .or. b% s1% lxtra16) .and. b% point_mass_i == 0) then
               b% point_mass_i = 2
               b% d_i = 1
               b% a_i = 2
               b% s_donor => b% s1
               b% s_accretor => b% s2
               b% period = 1d99
               b% mtransfer_rate = 0d0
               b% separation = &
                  pow_cr((b% s1% cgrav(1)*(b% m(1)+b% m(2)))*(b% period/(2*pi))**2,1d0/3d0)
               b% angular_momentum_j = b% m(1) * b% m(2) * sqrt( b% s1% cgrav(1) *&
                  b% separation / (b% m(1) + b% m(2)) )
               b% rl(1) = eval_rlobe(b% m(1), b% m(2), b% separation)
               b% rl(2) = eval_rlobe(b% m(2), b% m(1), b% separation)
               b% rl_relative_gap(1) = (b% r(1) - b% rl(1)) / b% rl(1) ! gap < 0 means out of contact 
               b% rl_relative_gap(2) = (b% r(2) - b% rl(2)) / b% rl(2) ! gap < 0 means out of contact
               b% do_tidal_sync = .false.
               b% mdot_scheme = "roche_lobe"
            end if
         end if

         ! check for L2 overflow after ZAMS
         if(.not. b% ignore_rlof .and. extras_binary_finish_step /= terminate) then
            if (b% m(1) > b% m(2)) then
               q = b% m(2) / b% m(1)
               star_id = 2
            else
               q = b% m(1) / b% m(2)
               star_id = 1
            end if
            if (b% rl_relative_gap(star_id) > 0.29858997d0*atan_cr(1.83530121d0*pow_cr(q,0.39661426d0))) then
               write(*,*) "Terminate due to L2 overflow"
               b% s1% lxtra8 = .true.
               extras_binary_finish_step = terminate
            end if
         end if

         if (extras_binary_finish_step == terminate) then
            call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
            if (ierr /= 0) return ! failure in profile
            call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
            if (ierr /= 0) return ! failure in profile
            b% s1% lxtra15 = .true.
         else
            !additional profiles to be saved
            center_h1 = b% s1% xa(b% s1% net_iso(ih1),b% s1% nz)
            center_h1_old = b% s1% xa_old(b% s1% net_iso(ih1),b% s1% nz_old)
            center_he4 = b% s1% xa(b% s1% net_iso(ihe4),b% s1% nz)
            center_he4_old = b% s1% xa_old(b% s1% net_iso(ihe4), b% s1% nz_old)
            if (center_h1 < 0.5 .and. center_h1_old > 0.5) then
               call star_write_profile_info(b% s1% id, "LOGS1/prof_2H50.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 0.25 .and. center_h1_old > 0.25) then
               call star_write_profile_info(b% s1% id, "LOGS1/prof_3H25.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 1d-6 .and. center_h1_old > 1d-6) then
               call star_write_profile_info(b% s1% id, "LOGS1/prof_4H00.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 1d-6) then
               if (center_he4 < 0.75 .and. center_he4_old > 0.75) then
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_5He75.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 0.5 .and. center_he4_old > 0.5) then
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_6He50.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 0.25 .and. center_he4_old > 0.25) then
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_7He25.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 1d-5 .and. center_he4_old > 1d-5) then
                  !create the profile at Y<1d-5, to avoid profiles to be created every single
                  !timestep if we only evolve until helium depletion, defined at Y=1d-6
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_8He00.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               end if
            end if

            ! check for beginning of mass transfer phases
            if (b% point_mass_i == 0 .and. abs(b% mtransfer_rate) > abs(b% mtransfer_rate_older) &
               .and. log10_cr(abs(b% mtransfer_rate_older)) < -90) then
               ! a mass transfer phase started
               if (center_h1 > 1d-6) then
                  b% s1% ixtra1 = b% s1% ixtra1 + 1
               else if (center_he4 > 1d-6) then
                  b% s1% ixtra2 = b% s1% ixtra2 + 1
               else
                  b% s1% ixtra3 = b% s1% ixtra3 + 1
               end if
            end if

            ! check for contact phases
            if (b% point_mass_i == 0 .and. log10_cr(abs(b% mtransfer_rate_older)) > -90 &
               .and. b% rl_relative_gap(1) > 0d0 .and. b% rl_relative_gap(2) > 0d0) then
               b% s1% lxtra9 = .true.
               if (b% d_i == 2) then
                  b% s1% lxtra10 = .true.
               end if
            end if

            center_h1 = b% s2% xa(b% s2% net_iso(ih1),b% s2% nz)
            center_h1_old = b% s2% xa_old(b% s2% net_iso(ih1),b% s2% nz_old)
            center_he4 = b% s2% xa(b% s2% net_iso(ihe4), b% s2% nz)
            center_he4_old = b% s2% xa_old(b% s2% net_iso(ihe4),b% s2% nz_old)
            if (center_h1 < 0.5 .and. center_h1_old > 0.5) then
               call star_write_profile_info(b% s2% id, "LOGS2/prof_2H50.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 0.25 .and. center_h1_old > 0.25) then
               call star_write_profile_info(b% s2% id, "LOGS2/prof_3H25.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 1d-6 .and. center_h1_old > 1d-6) then
               call star_write_profile_info(b% s2% id, "LOGS2/prof_4H00.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
            else if (center_h1 < 1d-6) then
               if (center_he4 < 0.75 .and. center_he4_old > 0.75) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_5He75.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 0.5 .and. center_he4_old > 0.5) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_6He50.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 0.25 .and. center_he4_old > 0.25) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_7He25.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_he4 < 1d-6 .and. center_he4_old > 1d-6) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_8He00.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               end if
            end if
         end if
         
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
         integer :: iounit
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         write(*,*) "Saving summary data"
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            write(*,*) "Error allocating iounit"
            return
         end if      

         open(unit=iounit, file=trim('summary.txt'), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open file summary.txt to write'
            stop 1
         end if

         write(iounit,*) "Case A:", b% s1% ixtra1 
         write(iounit,*) "Case B:", b% s1% ixtra2 
         write(iounit,*) "Case C:", b% s1% ixtra3 

         if (b% s1% lxtra1) &
            write(iounit,*) "Off-center Neon ignition for primary"
         if (b% s1% lxtra16) &
            write(iounit,*) "Off-center Neon ignition for secondary"
         if (b% s1% lxtra2) &
            write(iounit,*) "Primary depleted carbon"
         if (b% s1% lxtra19) &
            write(iounit,*) "Primary depleted helium with M>13 Msun"
         if (b% s1% lxtra3) &
            write(iounit,*) "Secondary depleted carbon"
         if (b% s1% lxtra4) &
            write(iounit,*) "Primary terminated last"
         if (b% s1% lxtra5) &
            write(iounit,*) "Secondary terminated last"
         if (b% s1% lxtra6) &
            write(iounit,*) "ZAMS L2 overflow"
         if (b% s1% lxtra7) &
            write(iounit,*) "ZAMS RLOF"
         if (b% s1% lxtra8) &
            write(iounit,*) "L2 overflow after ZAMS"
         if (b% s1% lxtra9 .and. .not. b% s1% lxtra10) &
            write(iounit,*) "Contact without inverse mass transfer"
         if (b% s1% lxtra10) &
            write(iounit,*) "Contact with inverse mass transfer"
         if (b% s1% lxtra11) &
            write(iounit,*) "System reached lower mdot limit"
         if (b% s1% lxtra12) &
            write(iounit,*) "System reached upper mdot limit"
         if (b% s1% lxtra13) &
            write(iounit,*) "Inverse mass transfer past secondary H depletion, primary in MS."
         if (b% s1% lxtra18) &
            write(iounit,*) "Inverse mass transfer past secondary H depletion, primary off MS."
         if (b% s1% lxtra14) &
            write(iounit,*) "Inverse mass transfer before secondary H depletion, primary off MS."
         if (b% s1% lxtra17) &
            write(iounit,*) "Terminated due to mass transfer in q_i<=0.25 system"
         if (b% s1% lxtra20) &
            write(iounit,*) "Terminated due to mass transfer rate reaching maximum"
         if (b% s1% lxtra21) &
            write(iounit,*) "Initial system is darwin unstable"
         if (b% s1% lxtra15 .or. b% s1% lxtra17) then
            write(iounit,*) "Normal termination"
         else
            write(iounit,*) "Abnormal termination"
         end if

         close(iounit)
 
      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
