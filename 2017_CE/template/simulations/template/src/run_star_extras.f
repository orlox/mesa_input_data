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
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use crlibm_lib
      
      implicit none
      
      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 15 ! minutes

      integer, parameter :: restart_info_alloc = 1
      integer, parameter :: restart_info_get = 2
      integer, parameter :: restart_info_put = 3
      real(dp), pointer :: m(:), entropy(:), xa_initial(:,:), helper(:)
      real(dp), pointer :: U_in(:), U_out(:), Omega_in(:), Omega_out(:)
      real(dp) :: initial_radius
      integer :: num_points, num_species
      
      
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% other_wind => brott_wind
         s% other_energy => my_other_energy
         s% accrete_given_mass_fractions = .true.
         s% accrete_same_as_surface = .false.
      end subroutine extras_controls
      
      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! Z=0.017 is Z from Grevesse et al. 1996
         Z_div_Z_solar = s% Zbase/0.017d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10_cr(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         if (T1 < Teff_jump) then
            ! low T wind
            w = max(vink_wind, nieu_wind)
         else
            ! high T wind
            alfa = 0d0
            if (Xs > 0.7d0) then
               alfa = 1d0
            else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
               alfa = (Xs - 0.4d0)/0.3d0
            end if
            w = alfa * vink_wind + (1d0-alfa) * hamann_wind
         end if

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if
            
            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow_cr(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10_cr(L1/Lsun/1d5) &
                  - 1.313d0*log10_cr(M1/Msun/30) &
                  - 1.226d0*log10_cr(vinf_div_vesc/2d0) &
                  + 0.933d0*log10_cr(T1/4d4) &
                  - 10.92d0*pow2(log10_cr(T1/4d4)) &
                  + 0.85d0*log10_cr(Z_div_Z_solar)
               w1 = exp10_cr(logMdot)
            else
               w1 = 0
            end if
            
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow_cr(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10_cr(L1/Lsun/1d5) &
                  - 1.339d0*log10_cr(M1/Msun/30) &
                  - 1.601d0*log10_cr(vinf_div_vesc/2d0) &
                  + 1.07d0*log10_cr(T1/2d4) &
                  + 0.85d0*log10_cr(Z_div_Z_solar)
               w2 = exp10_cr(logMdot)
            else
               w2 = 0
            end if
            
            w = alfa*w1 + (1 - alfa)*w2
            
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10_cr(L1/Lsun) &
                     +0.16d0*log10_cr(M1/Msun) &
                     +0.81d0*log10_cr(R1/Rsun) &
                     +0.85d0*log10_cr(Z_div_Z_solar)
            w = exp10_cr(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10_cr(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10_cr(Z_div_Z_solar)
            w = exp10_cr(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      integer function extras_startup(id, restart, ierr)
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interp_pm
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         integer :: op_err
         type (star_info), pointer :: s
         integer :: restart_time, prev_time_used, k, i
         real(dp), pointer :: interp_work(:), work(:), p(:)
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. restart) then
            call system_clock(time0,clock_rate)
            call alloc_restart_info(s)
         else
            call unpack_restart_info(s)
            call system_clock(restart_time,clock_rate)
            prev_time_used = time1 - time0
            time1 = restart_time
            time0 = time1 - prev_time_used
         end if
         extras_startup = keep_going
         allocate(m(s% nz), entropy(4*s% nz), xa_initial(s% species, 4*s% nz), helper(4*s% nz))
         allocate(U_in(4*s% nz), U_out(4*s% nz), Omega_in(4*s% nz), Omega_out(4*s% nz))
         m(:) = s% m(:s% nz)
         do k=1, s% nz
            entropy(4*k-3) = exp_cr(s% lnS(k))
         end do
         do i=1, s% species
            do k=1, s% nz
               xa_initial(i,4*k-3) = s% xa(i,k)
            end do
         end do
         U_out(1) = s% energy(1)*s% dm(1)
         Omega_out(1) = - standard_cgrav*s% m(1)*s% dm_bar(1)/s% r(1)
         do k=2, s% nz
            U_out(4*k-3) = U_out(4*(k-1)-3) + s% energy(k)*s% dm(k)
            Omega_out(4*k-3) = Omega_out(4*(k-1)-3) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
         end do
         U_in(4*s% nz-3) = s% energy(s% nz)*s% dm(s% nz)
         Omega_in(4*s% nz-3) = - standard_cgrav*s% m(s% nz)*s% dm_bar(s% nz)/s% r(s% nz)
         do k=s% nz-1, 1, -1
            U_in(4*k-3) = U_in(4*(k+1)-3) + s% energy(k)*s% dm(k)
            Omega_in(4*k-3) = Omega_in(4*(k+1)-3) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
         end do

         initial_radius = s% r(1)
         num_points = s% nz

         allocate(interp_work(s% nz*pm_work_size), stat=ierr)
         call interp_pm(m, s% nz, entropy, pm_work_size, interp_work, 'entropy interpolant', op_err)
         call interp_pm(m, s% nz, U_in, pm_work_size, interp_work, 'U_in interpolant', op_err)
         call interp_pm(m, s% nz, U_out, pm_work_size, interp_work, 'U_out interpolant', op_err)
         call interp_pm(m, s% nz, Omega_in, pm_work_size, interp_work, 'Omega_in interpolant', op_err)
         call interp_pm(m, s% nz, Omega_out, pm_work_size, interp_work, 'Omega_out interpolant', op_err)
         do i=1, s% species
            helper(:) = xa_initial(i,:)
            call interp_pm(m, s% nz, helper, pm_work_size, interp_work, 'constant_xa', op_err)
            xa_initial(i,:) = helper(:)
         end do
         deallocate(interp_work)
         write(*,*) "created interpolant", num_points

      end function extras_startup
      
      subroutine my_other_energy(id, ierr)
         use const_def, only: Rsun
         use interp_1d_lib, only: interp_values
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, op_err
         real(dp), pointer :: vals(:)
         real(dp) :: timescale

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         allocate(vals(s% nz))
         timescale = abs(s% mstar/s% mstar_dot/1000d0)
         call interp_values(m, num_points, entropy, s% nz, s% m, vals, op_err)
         s% extra_heat(:) = 0
         do k = 1, s% nz
            s% extra_heat(k) = ( 1d0 - exp_cr(s% lnS(k))/vals(k) ) * exp_cr(s%lnE(k)) \
                / timescale
            !write(*,*) k, s% m(k)/Msun, s% extra_heat(k), exp_cr(s%lnS(k)), vals(k), exp_cr(s%lnE(k))
         end do


         deallocate(vals)
      end subroutine my_other_energy

      integer function extras_check_model(id, id_extra)
         use chem_def, only: chem_isos
         use interp_1d_lib, only: interp_value
         integer, intent(in) :: id, id_extra
         type (star_info), pointer :: s
         integer :: ierr, j, op_err
         real(dp) :: xaccrete
         real(dp) :: sum_xa

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% max_timestep = abs(s% mstar / s% mstar_dot/1000)
         !write(*,*) s% max_timestep/secyer

         ! set xaccrete
         !write(*,*) "check comp", s% m(1)/Msun
         sum_xa = 0d0
         s% num_accretion_species = s% species
         do j = 1, s% species
            s% accretion_species_id(j) = chem_isos% name(s% chem_id(j))
            helper(:) = xa_initial(j,:)
            call interp_value(m, num_points, helper, s% m(1), xaccrete, op_err)
            s% accretion_species_xa(j) = xaccrete
            !write(*,*) xaccrete
            sum_xa = sum_xa + xaccrete
         end do
         !write(*,*) "sum_xa", sum_xa
         s% accretion_species_xa(:) = s% accretion_species_xa(:)/sum_xa

         extras_check_model = keep_going

         !if(s% mass_change > 0d0 .and. s% r(1)>initial_radius/3d0) then
         !   write(*,*) "Terminating because R>R_i/3"
         !   extras_check_model = terminate
         !end if
         if(s% X(1) < 1e-4 .and. s% mass_change < 0d0) then
            !s% mass_change = 1e-10
            !s% use_other_energy = .false.
            write(*,*) "Terminating because He core is reached"
            extras_check_model = terminate
         end if
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 7
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         use interp_1d_lib, only: interp_value
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt
         type (star_info), pointer :: s
         real(dp) :: U_removed, Omega_removed, U_inold, Omega_inold, Ecore, Ecore_i, Ebind
         integer :: k, op_err
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dt = dble(time1 - time0) / clock_rate / 60
         names(1) = 'runtime_minutes'
         vals(1) = dt

         names(2) = "Ebind"
         names(3) = "delta_Ecore"
         names(4) = "Uout"
         names(5) = "lambda"
         names(6) = "lambda_plus_Ecore"
         names(7) = "lambda_sub_Uout"

         Ecore = 0
         do k=1, s% nz
            Ecore = Ecore + s% energy(k)*s% dm(k) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
         end do

         call interp_value(m, num_points, U_out, s% m(1), U_removed, op_err)
         call interp_value(m, num_points, Omega_out, s% m(1), Omega_removed, op_err)
         call interp_value(m, num_points, U_in, s% m(1), U_inold, op_err)
         call interp_value(m, num_points, Omega_in, s% m(1), Omega_inold, op_err)

         Ecore_i = U_inold + Omega_inold
         Ebind = U_removed + Omega_removed

         !write(*,*) Ecore, Ecore_i, Ecore-Ecore_i

         vals(2) = U_removed + Omega_removed
         vals(3) = Ecore - Ecore_i
         vals(4) = U_removed
         vals(5) = -standard_cgrav*m(1)*(m(1)-s% m(1))/initial_radius/Ebind
         vals(6) = -standard_cgrav*m(1)*(m(1)-s% m(1))/initial_radius/(Ebind-Ecore+Ecore_i)
         vals(7) = -standard_cgrav*m(1)*(m(1)-s% m(1))/initial_radius/(Ebind-U_removed)

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 6
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = "U_out"
         names(2) = "U_in"
         names(3) = "Omega_out"
         names(4) = "Omega_in"
         names(5) = "Ecore"
         names(6) = "Ebind"

         vals(1,2) = s% energy(1)*s% dm(1)
         vals(1,4) = - standard_cgrav*s% m(1)*s% dm_bar(1)/s% r(1)
         do k=2, s% nz
            vals(k,2) = vals(k-1,2) + s% energy(k)*s% dm(k)
            vals(k,4) = vals(k-1,4) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
         end do
         vals(s% nz,1) = s% energy(s% nz)*s% dm(s% nz)
         vals(s% nz,3) = - standard_cgrav*s% m(s% nz)*s% dm_bar(s% nz)/s% r(s% nz)
         do k=s% nz-1, 1, -1
            vals(k,1) = vals(k+1,1) + s% energy(k)*s% dm(k)
            vals(k,3) = vals(k+1,3) - standard_cgrav*s% m(k)*s% dm_bar(k)/s% r(k)
         end do
         do k=1, s% nz
            vals(k,5) = vals(k,1) + vals(k,3)
            vals(k,6) = vals(k,2) + vals(k,4)
         end do

      end subroutine data_for_extra_profile_columns
      

      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if(abs(s% xtra25_old - s% center_h1) > 0.0005) then
             s% dt_next = min(s% dt_next, s% dt * s% min_timestep_factor)
             write(*,*) "reducing dt due to large change in central hydrogen"
          else if(abs(s% xtra26_old - s% center_he4) > 0.001) then
             s% dt_next = min(s% dt_next, s% dt * s% min_timestep_factor)
             write(*,*) "reducing dt due to large change in central helium"
          else if(abs(s% xtra27_old - s% center_c12) > 0.001) then
             s% dt_next = min(s% dt_next, s% dt * s% min_timestep_factor)
             write(*,*) "reducing dt due to large change in central carbon"
         end if
         s% xtra25 = s% center_h1
         s% xtra26 = s% center_he4
         s% xtra27 = s% center_c12
         extras_finish_step = keep_going
         call system_clock(time1,clock_rate)
         call store_restart_info(s)
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         dt = dble(time1 - time0) / clock_rate / 60
         if (dt > 10*expected_runtime) then
            write(*,'(/,a30,2f18.6,a,/)') '>>>>>>> EXCESSIVE runtime', &
               dt, expected_runtime, '   <<<<<<<<<  ERROR'
         else
            write(*,'(/,a50,2f18.6,99i10/)') 'runtime, retries, backups, steps', &
               dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
         end if
         write(*,*) "deallocate arrays"
         deallocate(m,entropy,xa_initial,helper,U_in,U_out,Omega_in,Omega_out)
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring data so can do restarts

      
      subroutine alloc_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_alloc)
      end subroutine alloc_restart_info
      
      
      subroutine unpack_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_get)
      end subroutine unpack_restart_info
      
      
      subroutine store_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_put)
      end subroutine store_restart_info
      
      
      subroutine move_restart_info(s,op)
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg 
         call move_int(time0)
         call move_int(time1)
         
         num_ints = i
         
         i = 0
         ! call move_dbl 
         
         num_dbls = i
         
         if (op /= restart_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (restart_info_get)
               dbl = s% extra_work(i)
            case (restart_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            include 'formats'
            i = i+1
            select case (op)
            case (restart_info_get)
               !write(*,3) 'restore int', i, s% extra_iwork(i)
               int = s% extra_iwork(i)
            case (restart_info_put)
               !write(*,3) 'save int', i, int
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (restart_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (restart_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_restart_info
      
      


      end module run_star_extras
      
