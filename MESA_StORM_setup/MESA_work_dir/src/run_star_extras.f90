! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
      use math_lib
      use chem_def

      implicit none

!#######################################
! Initialise a bunch of parameters
!#######################################

      ! declarations for xtra_coeff_os (extra resolution of overshooting zones must now be defined directly in run_star_extras)
      real(dp) :: xtra_coef_os_full_on, xtra_coef_os_full_off, xtra_coef_os_above_nonburn, xtra_coef_os_below_nonburn, &
                  xtra_coef_os_above_burn_h, xtra_coef_os_below_burn_h, xtra_coef_os_above_burn_he, xtra_coef_os_below_burn_he, &
                  xtra_coef_os_above_burn_z, xtra_coef_os_below_burn_z, xtra_dist_os_above_nonburn, xtra_dist_os_below_nonburn, &
                  xtra_dist_os_above_burn_h, xtra_dist_os_below_burn_h, xtra_dist_os_above_burn_he, xtra_dist_os_below_burn_he, &
                  xtra_dist_os_above_burn_z, xtra_dist_os_below_burn_z
      namelist /xtra_coeff_os/ xtra_coef_os_full_on, xtra_coef_os_full_off, xtra_coef_os_above_nonburn, xtra_coef_os_below_nonburn, &
                               xtra_coef_os_above_burn_h, xtra_coef_os_below_burn_h, xtra_coef_os_above_burn_he, xtra_coef_os_below_burn_he, &
                               xtra_coef_os_above_burn_z, xtra_coef_os_below_burn_z, xtra_dist_os_above_nonburn, xtra_dist_os_below_nonburn, &
                               xtra_dist_os_above_burn_h, xtra_dist_os_below_burn_h, xtra_dist_os_above_burn_he, xtra_dist_os_below_burn_he, &
                               xtra_dist_os_above_burn_z, xtra_dist_os_below_burn_z

      ! Saving parameters during MS stage
      real(dp) :: Xc_save, Xc_save_step1, Xc_save_step2, Xc_save_step3, Xc_save_log_iterator, Xc_save_log_exponent, Xc_precision !--> saving parameters\

      ! Control minimum change in central H and He burning
      real(dp), parameter :: min_Xc_chb = 1d-12 !--> end of MS
      real(dp), parameter :: min_Yc_cheb = 1d-12 !--> end of core He burning
      real(dp), parameter :: lowest_delta_lg_XH_cntr_limit = 1d-5

      logical :: ask_reached_zams, reached_zams

      ! M_ini, X_ini, Y_ini, Z_ini --> initial parameters
      real(dp) :: M_ini, X_ini, Y_ini, Z_ini, logD_mix

      ! Galactic enrichment law parameters
      real(dp) :: Yp, dY_dZ

      ! output_dir --> where to save history, profile, mod files
      ! hist_name --> name for history files
      ! evol_stage --> current global evolutionary stage
      character(len=256)  :: output_dir, preMS_dir, hist_name, profile_name, evol_stage, ZAMS_file_name

      logical :: need_to_save
      ! these routines are called by the standard run_star check_model

      contains

 !********************************************************************************************
 ! *** Extra controls to be set, i.e. setting overshooting, mixing, save folders, etc.
 !********************************************************************************************
      subroutine extras_controls(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        character(len=256) :: end_model_name
        real(dp) :: a_ov, f_ov, f0
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

        ! Meshing is now defined in run_star_extras through hook
        call read_inlist_xtra_coeff_os(ierr)
        if (ierr /= 0) return
        s% other_mesh_delta_coeff_factor => other_mesh_delta_coeff_factor

        s% job% warn_run_star_extras = .true.

        ! Read input from bash script
        print *
        print "('Insert output_directory, M_ini, Z_ini, logD_mix, a_ov, f_ov')"
        read(*,*) output_dir, Z_ini, M_ini, logD_mix, a_ov, f_ov
        write(*,*)'Zini = ',Z_ini
        write(*,*)'Mini = ',M_ini
        write(*,*)'logD_mix = ',logD_mix
        write(*,*)'a_ov = ',a_ov
        write(*,*)'f_ov = ',f_ov

        f0 = 0.005d0

        ! Set initial fractions according to Galactic enrichtment law
        Yp = 0.2453d0 ! Primordial helium abundance (Aver et al. 2021)
        dY_dZ = 2.193d0 ! Galactic enrichment ratio. Calibrated by requiring that for solar (Xini, Yini, Zini) = (0.710, 0.276, 0.014)
        Y_ini = Yp + dY_dZ*Z_ini
        X_ini = 1 - (Z_ini + Y_ini)

        write(*,*) 'Xini= ',X_ini
        write(*,*) 'Yini= ',Y_ini
        write(*,*) 'Zini= ',Z_ini

        ! MS saving properties (see subroutine save_at_Xc, which is called during every step)
        s% x_logical_ctrl(3) = .false. ! Whether we're in a bisection looking for Xc near Xc_save
        s% x_integer_ctrl(1) = 0       ! Which iteration of the bisection we're in

        Xc_save_step1 = s% x_ctrl(3) !* X_ini ! Defined in inlist
        Xc_save_step2 = s% x_ctrl(4) !* X_ini ! Defined in inlist

        Xc_save_log_exponent = s% x_ctrl(5) ! Defines the steepness of the log saving on the hook
        Xc_save_log_iterator = s% x_ctrl(6) ! Determines the amount of models (in steps of 0.1) one wants to save on hook (1.0 --> ten models, 1.5 --> 15 models, 0.5 --> 5 models etc.)
        Xc_precision = 1d-4 ! Precision to which we save at a certain Xc

        if (s% job% create_pre_main_sequence_model) then
          s% job% save_model_when_terminate = .true.
          ! s% scale_max_correction = 0.1       ! to help pre-MS convergence
          ! Set initial composition and mass
          call set_XYZ_ini_preMS(id, X_ini, Y_ini, Z_ini, ierr)
          call set_initial_mass(id, M_ini, ierr)
          ! Set mixing parameters
          write(*,*)'logD is started as:', logD_mix
          call set_logD_mix(id, 1d0, ierr)            ! Set output directory and the name for the history file

          if (M_ini<10) then
              write (hist_name, '(A,F5.3,A,F4.2,A)') 'Z', Z_ini, '_M', M_ini, '_preMS'
              write (ZAMS_file_name, '(A,F5.3,A,F4.2,A)') 'Z', Z_ini, '_M', M_ini, '_ZAMS.mod'
          else
              write (hist_name, '(A,F5.3,A,F5.2,A)') 'Z', Z_ini, '_M', M_ini, '_preMS'
              write (ZAMS_file_name, '(A,F5.3,A,F5.2,A)') 'Z', Z_ini, '_M', M_ini, '_ZAMS.mod'
          endif
          output_dir = 'preMS/' ! CAREFUL: I overwrite the output_dir provided in the
          s% job% save_model_filename = trim(output_dir) // trim(ZAMS_file_name)
          s% log_directory= trim(output_dir)
        else
          s% kap_rq% Zbase = Z_ini
          ! Extra meshing around convective core boundary
          s% job% save_model_when_terminate = .true.
          if (ierr /= 0) return
          s% other_D_mix => IGW_D_mix_rho
          ! s% job% set_initial_dt = .true.
          ! s% job% years_for_initial_dt = 1d3
          s% job% load_saved_model = .true.
          ! Set output directory and the name for the history file
          preMS_dir = 'preMS/'
          s% log_directory= trim(output_dir)
          s% photo_directory=trim(preMS_dir)
          if (M_ini<10) then
            write (hist_name, '(A,F5.3,A,F4.2,A,F3.1,A,F5.3)') 'Z', Z_ini, '_M', M_ini, '_logD', logD_mix, '_fov', f_ov
            write (ZAMS_file_name, '(A,F5.3,A,F4.2,A)') 'Z', Z_ini, '_M', M_ini,'_ZAMS.mod'
          else
            write (hist_name, '(A,F5.3,A,F5.2,A,F3.1,A,F5.3)') 'Z', Z_ini, '_M', M_ini, '_logD', logD_mix, '_fov', f_ov
            write (ZAMS_file_name, '(A,F5.3,A,F5.2,A)') 'Z', Z_ini, '_M', M_ini, '_ZAMS.mod'
          endif
          s% job% load_model_filename = trim(preMS_dir) // '/' // trim(ZAMS_file_name)
          s% job% save_model_filename = trim(output_dir) // '/' // trim(hist_name) // '_TAMS.mod'
          ! Set mixing parameters
          write(*,*)'logD is started as:', logD_mix
          call set_logD_mix(id, logD_mix, ierr)

          ! As you vary Z, X_ini changes and thus it is possible that the initial Xc_save has already been passed, which will cause endless retries.
          ! To prevent this, the initial Xc_save is lowered by the initial Xc_save_step until Xc_save drops below X_ini.
          ! This effectively skips a couple of steps of output for metal-rich stars, though the remaining output is still at the same Xc as the solar grid.
          Xc_save = 0.701d0
          do while (X_ini < Xc_save)
            Xc_save = Xc_save - s% x_ctrl(3)
          end do
          ! For consistency, also raise the initial Xc_save if it is too far below the X_ini, to keep a fairly even coverage of the MS evolution
          ! However, here we use a wider tolerance of twice the initial Xc_save_step as Xc can already decrease a bit during the preMS
          do while(X_ini > Xc_save + 2d0*s% x_ctrl(3))
            Xc_save = Xc_save + s% x_ctrl(3)
          end do
        endif

        ! Naming and specifying history output. Final model is named in if-else above
        s% star_history_name = trim(hist_name) // '.hist'

        ! Sets core overshooting scheme based on user input.
        if (f_ov > 0.00d0 .and. a_ov == 0.00d0) then
           ! Exponential overshooting
           s% overshoot_scheme(1) = 'exponential'
           s% overshoot_zone_type(1) = 'burn_H'
           s% overshoot_zone_loc(1) = 'core'
           s% overshoot_bdy_loc(1) = 'any'
           s% overshoot_f(1) = f_ov
           s% overshoot_f0(1) = f0
           s% overshoot_D_min = 1d-2
           write(*,*) 'Diffusive exponential overshoot enabled'
        else if (f_ov == 0.00 .and. a_ov > 0.00) then
           ! Step overshooting
           s% overshoot_scheme(1) = 'step'
           s% overshoot_zone_type(1) = 'burn_H'
           s% overshoot_zone_loc(1) = 'core'
           s% overshoot_bdy_loc(1) = 'any'
           s% overshoot_f(1) = a_ov
           s% overshoot_f0(1) = f0
           s% overshoot_D_min = 1d-2
           write(*,*) 'Step overshoot enabled'
        else
            write(*,*) 'No overshoot enabled'
        end if

        ! Set some overshooting for:
        ! - the convection zones near the surface,
        ! - H-burning shell developing end-MS/post-MS.
        ! Turned on in inlist_project.
        ! Diffusive exponential used
        if (s% x_logical_ctrl(2) .and. s% x_ctrl(2) > f0)  then
            ! Envelope non-burn convection zones
            s% overshoot_scheme(2) = 'exponential'
            s% overshoot_zone_type(2) = 'nonburn'
            s% overshoot_zone_loc(2) = 'any'
            s% overshoot_bdy_loc(2) = 'any'
            s% overshoot_f(2) = s% x_ctrl(2)
            s% overshoot_f0(2) = f0
            ! H-shell
            s% overshoot_scheme(3) = 'exponential'
            s% overshoot_zone_type(3) = 'burn_H'
            s% overshoot_zone_loc(3) = 'shell'
            s% overshoot_bdy_loc(3) = 'any'
            s% overshoot_f(3) = s% x_ctrl(2)
            s% overshoot_f0(3) = f0
            write(*,*) 'Overshooting in nonburn convection zones set at: ', s% x_ctrl(2)
            write(*,*) 'Overshooting in shell burning layers set at: ', s% x_ctrl(2)
        else if (s% x_logical_ctrl(1) .and. s% x_ctrl(2) <= f0) then
            write(*,*) 'Set undershooting/overshooting must be larger than f0'
        end if

        ! Different mixing is set here
        if (s% x_logical_ctrl(1)) then
            s% other_D_mix  => IGW_D_mix_rho
            write(*,*) "IGW mixing profile (rho^-1) set with logD at base of envelope: ",logD_mix
        endif

      end subroutine extras_controls

!**********************************************************************************************
! Code for xtra_coeff_os (copied from $MESA_dir/star/test_suite/agb/src/run_star_extras.f)      !
! Must be defined manually now instead of as an inlist control (see v.12778 code update notes)!
! See inlist_xtra_coeff_os for changes                                                        !
!**********************************************************************************************

    ! Reads inlist_xtra_coeff_os
      subroutine read_inlist_xtra_coeff_os(ierr)
         use utils_lib
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit

         filename = 'inlist_xtra_coeff_os'

         write(*,*) 'read_inlist_xtra_coeff_os'

         xtra_coef_os_full_on = 1d-4
         xtra_coef_os_full_off = 0.1d0
         xtra_coef_os_above_nonburn = 0.2d0
         xtra_coef_os_below_nonburn = 0.2d0
         xtra_coef_os_above_burn_h = 0.2d0
         xtra_coef_os_below_burn_h = 0.2d0
         xtra_coef_os_above_burn_he = 0.2d0
         xtra_coef_os_below_burn_he = 0.2d0
         xtra_coef_os_above_burn_z = 0.2d0
         xtra_coef_os_below_burn_z = 0.2d0
         xtra_dist_os_above_nonburn = 2.0d0
         xtra_dist_os_below_nonburn = 2.0d0
         xtra_dist_os_above_burn_h = 2.0d0
         xtra_dist_os_below_burn_h = 2.0d0
         xtra_dist_os_above_burn_he = 2.0d0
         xtra_dist_os_below_burn_he = 2.0d0
         xtra_dist_os_above_burn_z = 2.0d0
         xtra_dist_os_below_burn_z = 2.0d0

         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=xtra_coeff_os, iostat=ierr)
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *)
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=xtra_coeff_os)
               close(unit)
            end if
         end if
      end subroutine read_inlist_xtra_coeff_os

      ! Other mesh routine (copied from agb testsuite included with MESA v.12778)
      subroutine other_mesh_delta_coeff_factor(id, ierr)
         use const_def
         use chem_def
         integer, intent(in) :: id
         real(dp), allocatable, dimension(:) :: eps_h, eps_he, eps_z
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: he_cntr, full_off, full_on, alfa_os
         integer :: k, kk, nz, max_eps_loc
         real(dp) :: xtra_coef, xtra_dist, coef, Hp, r_extra, max_eps, eps
         logical :: in_convective_region
         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         if (xtra_coef_os_above_nonburn == 1d0 .and. &
             xtra_coef_os_below_nonburn == 1d0 .and. &
             xtra_coef_os_above_burn_h == 1d0 .and. &
             xtra_coef_os_below_burn_h == 1d0 .and. &
             xtra_coef_os_above_burn_he == 1d0 .and. &
             xtra_coef_os_below_burn_he == 1d0 .and. &
             xtra_coef_os_above_burn_z == 1d0 .and. &
             xtra_coef_os_below_burn_z == 1d0) return

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Only turn on once reached ZAMS
         if (.not.reached_zams) return

         write(*,*) 'Turning additional meshing on, i.e. xtra_coef_os_above_burn_h:', xtra_coef_os_above_burn_h
         nz = s% nz
         he_cntr = s% xa(s% net_iso(ihe4),nz)
         full_off = xtra_coef_os_full_off
         full_on = xtra_coef_os_full_on
         if (he_cntr >= full_off) then
            alfa_os = 0
         else if (he_cntr <= full_on) then
            alfa_os = 1
         else
            alfa_os = (full_off - he_cntr)/(full_off - full_on)
         end if
         !write(*,1) 'alfa_os', alfa_os
         if (alfa_os == 0) return

         ! Manually defining eps_h, eps_he and eps_z
         allocate(eps_h(s% nz), eps_he(s% nz), eps_z(s% nz))
         do k=1,s% nz
            eps_h(k) = s% eps_nuc_categories(ipp,k) + s% eps_nuc_categories(icno,k)
            eps_he(k) = s% eps_nuc_categories(i3alf,k)
            eps_z(k) = s% eps_nuc(k) - (eps_h(k) + eps_he(k))
         end do

         ! first go from surface to center doing below convective boundaries
         in_convective_region = (s% mixing_type(1) == convective_mixing)
         k = 2
         max_eps = -1d99
         max_eps_loc = -1
         do while (k <= nz)
            eps = eps_h(k) + eps_he(k) + eps_z(k)
            if (in_convective_region) then
               if (s% mixing_type(k) == convective_mixing) then
                  if (eps > max_eps) then
                     max_eps = eps
                     max_eps_loc = k
                  end if
               else
                  in_convective_region = .false.
                  if (max_eps < 1d0) then
                     xtra_coef = xtra_coef_os_below_nonburn
                     xtra_dist = xtra_dist_os_below_nonburn
                  else if (eps_h(max_eps_loc) > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_below_burn_h
                     xtra_dist = xtra_dist_os_below_burn_h
                  else if (eps_he(max_eps_loc) > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_below_burn_he
                     xtra_dist = xtra_dist_os_below_burn_he
                  else
                     xtra_coef = xtra_coef_os_below_burn_z
                     xtra_dist = xtra_dist_os_below_burn_z
                  end if
                  xtra_coef = xtra_coef*alfa_os + (1-alfa_os)
                  if (xtra_coef > 0 .and. xtra_coef /= 1) then
                     coef = xtra_coef
                     do
                        if (s% mixing_type(k) /= overshoot_mixing) exit
                        if (coef < s% mesh_delta_coeff_factor(k)) then
                            s% mesh_delta_coeff_factor(k) = coef
                        end if
                        if (k == nz) exit
                        k = k+1
                     end do
                     if (xtra_dist > 0) then
                        Hp = s% Peos(k)/(s% rho(k)*s% grav(k))
                        r_extra = max(0d0, s% r(k) - xtra_dist*Hp)
                        if (dbg) write(*,2) 'extra below overshoot region', &
                           k, s% r(k)/Rsun, Hp/Rsun, r_extra/Rsun
                        do
                           if (s% r(k) < r_extra) exit
                           if (coef < s% mesh_delta_coeff_factor(k)) then
                              s% mesh_delta_coeff_factor(k) = coef
                           end if
                           if (k == nz) exit
                           k = k+1
                        end do
                     end if
                  end if
                  if (dbg) write(*,2) 'done with extra below overshoot region', k
                  if (dbg) write(*,*)
               end if
            else if (s% mixing_type(k) == convective_mixing) then
               in_convective_region = .true.
               max_eps = eps
               max_eps_loc = k
            end if
            k = k+1
         end do

         ! now go from center to surface doing above convective boundaries
         in_convective_region = (s% mixing_type(nz) == convective_mixing)
         k = nz-1
         max_eps = -1d99
         max_eps_loc = -1
         do while (k >= 1)
            eps = eps_h(k) + eps_he(k) + eps_z(k)
            if (in_convective_region) then
               if (s% mixing_type(k) == convective_mixing) then
                  if (eps > max_eps) then
                     max_eps = eps
                     max_eps_loc = k
                  end if
               else
                  in_convective_region = .false.
                  if (max_eps < 1d0) then
                     xtra_coef = xtra_coef_os_above_nonburn
                     xtra_dist = xtra_dist_os_above_nonburn
                  else if (eps_h(max_eps_loc) > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_above_burn_h
                     xtra_dist = xtra_dist_os_above_burn_h
                  else if (eps_he(max_eps_loc) > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_above_burn_he
                     xtra_dist = xtra_dist_os_above_burn_he
                  else
                     xtra_coef = xtra_coef_os_above_burn_z
                     xtra_dist = xtra_dist_os_above_burn_z
                  end if
                  xtra_coef = xtra_coef*alfa_os + (1-alfa_os)
                  if (dbg) write(*,2) 'xtra_coeff to surf', s% model_number, xtra_coef

                  if (xtra_coef > 0 .and. xtra_coef /= 1) then
                     coef = xtra_coef
                     do
                        if (s% mixing_type(k) /= overshoot_mixing) exit
                        if (coef < s% mesh_delta_coeff_factor(k)) then
                            s% mesh_delta_coeff_factor(k) = coef
                        end if
                        if (k == 1) exit
                        k = k-1
                     end do
                     if (xtra_dist > 0) then
                        Hp = s% Peos(k)/(s% rho(k)*s% grav(k))
                        r_extra = min(s% r(1), s% r(k) + xtra_dist*Hp)
                        if (dbg) write(*,2) 'extra above overshoot region', &
                           k, s% r(k)/Rsun, Hp/Rsun, r_extra/Rsun
                        do
                           if (s% r(k) > r_extra) exit
                           if (coef < s% mesh_delta_coeff_factor(k)) then
                               s% mesh_delta_coeff_factor(k) = coef
                           end if
                           if (k == 1) exit
                           k = k-1
                        end do
                     end if
                  end if
                  if (dbg) write(*,2) 'done with extra above overshoot region', k
                  if (dbg) write(*,*)
               end if
            else if (s% mixing_type(k) == convective_mixing) then
               in_convective_region = .true.
               max_eps = eps
               max_eps_loc = k
            end if
            k = k-1
         end do

      end subroutine other_mesh_delta_coeff_factor

!********************************************************************************************
!****     end of code for xtra_coeff_os
!********************************************************************************************

!********************************************************************************************
! Set initial model settings
!********************************************************************************************
subroutine extras_startup(id, restart, ierr)
   integer, intent(in) :: id
   logical, intent(in) :: restart
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   write(*,*)'logD is started as:', logD_mix
   call set_logD_mix(id, logD_mix, ierr)

   ! Keep track of evolutionary stage
   ! Assume initial models are run from birth
   reached_zams = .false.
   ask_reached_zams = .true.
 end subroutine extras_startup

!********************************************************************************************
! Check model, before we end the step (i.e. which evolutionary stage)
!********************************************************************************************
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         logical :: do_retry
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going

        ! Assume saving or retrying is not necessary unless specific conditions are triggered
        do_retry = .false.
        need_to_save = .false.

        ! Checks current evolutionary stage
       if (ask_reached_zams) then
            call check_reached_zams(id, reached_zams, ierr)
            if (reached_zams) ask_reached_zams = .false.
        endif

        ! Once ZAMS reached save at intervals of certain intervals when Xc>0.10*Xini
        ! Save at two percent intervals when Xc<0.10*Xini
        ! Once we finish in MS Xc_save is set to -1 in order to stop the loop.
        if (.not. do_retry .and. .not. need_to_save .and. Xc_save .ne. -1d0) then
        !    call save_at_Xc(id, Xc_save, Xc_precision, need_to_save, ierr)
            call save_at_Xc_retry_method(id, Xc_save, Xc_precision, need_to_save, do_retry, ierr)
            if (do_retry) extras_check_model = retry
            if (need_to_save .and. .not. do_retry) then
                ! Depending on our Xc value we save more/less
                ! Set to values which worked for model grid for bcep stars.
                ! >0.25 (first half), < 0.25 (second half), <0.013 (special sampling of hook)
                if (Xc_save .gt. 0.25d0) then
                    Xc_save = Xc_save-Xc_save_step1
                elseif (Xc_save .le. 0.251d0 .and. Xc_save .gt. 0.014d0) then
                    Xc_save = Xc_save - Xc_save_step2
                 ! In this last little bit of the MS-evolution we save ten times (between Xc=0.013 and s% xa_central_lower_limit(1))
                 ! This generally allows us to sample the blue hook after the REMS (at least 10 times), where Xc not 'linear' with changes in Teff, logL
                 ! We start with Xc_save_log_iterator = %s x_ctrl(6), and take off 0.1 every time we save until we hit -0.1.
                 ! The log iterator is used as the base of an exponent and makes us save more often closer to our lower limit
                 ! With higher values making the saving more steep. From tests 2 seems to be a good value.
                 ! After which Xc_save is set to -1. And Xc is not used as a save indicator anymore.
                elseif (Xc_save .le. 0.014d0) then
                    Xc_save_log_iterator = Xc_save_log_iterator - 0.1
                    if (Xc_save_log_iterator .gt. -0.05d0 .and. s% xa_central_lower_limit(1) .lt. 0.013d0) then
                        Xc_save = pow(Xc_save_log_iterator,Xc_save_log_exponent)*(0.014d0 - s% xa_central_lower_limit(1)) + s% xa_central_lower_limit(1)
                    elseif (Xc_save_log_iterator .le. -0.05d0 .or. s% xa_central_lower_limit(1) .gt. 0.013d0) then
                        Xc_save = -1d0
                    endif
                end if
            endif
        end if

         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model

      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

        if (need_to_save) then
            print*, 'MESA save_profile_model_number', s% model_number
            call save_profile_model_number(id, prefix=hist_name, ierr=ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: run_star_extras: extras_finish_step: save_profile_model_number failed'
                return
            endif
        endif

        if (need_to_save) then
            print*, 'MESA save_gyre', s% model_number
            call save_gyre_model_number(id, prefix=hist_name, add_atm = s% add_atmosphere_to_pulse_data, keep_surf = s% keep_surface_point_for_pulse_data, ierr=ierr)
            if (ierr /= 0) then
               write(*,*) 'Error: run_star_extras: extras_finish_step: save_gyre failed'
               return
            endif
        endif

         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step

        integer function how_many_extra_history_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_columns = 1
        end function how_many_extra_history_columns

        subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            character (len=4) :: element_array(32)  ! Array with names of surface abundances to track
            integer :: k
            integer :: vals_nr  ! The number assigned to the values
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            vals_nr = 1 ! Keeps track of added values.

            ! Routines by M. Michielsen
            call hist_asymptotic_dP(id, names, vals, vals_nr, ierr)
        end subroutine data_for_extra_history_columns

        integer function how_many_extra_profile_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_columns = 0
        end function how_many_extra_profile_columns

        subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
            integer, intent(in) :: id, n, nz
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
            integer :: vals_nr  ! The number assigned to the values
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
        end subroutine data_for_extra_profile_columns


        subroutine extras_after_evolve(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
        end subroutine extras_after_evolve

        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Additional subroutines defined (from github/mathiasm or by Siemen))
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Set the initial mass of the model (M. Michielsen)
        subroutine set_initial_mass(id, M_ini, ierr)
          integer, intent(in) :: id
          real(dp), intent(in) :: M_ini
          integer, intent(out) :: ierr

          type(star_info), pointer :: s

          call star_ptr(id, s, ierr)
          if (ierr /= 0) then
              write(*,*) 'Error: runstarex_input_physics: set_initial_mass: star_ptr failed'
              return
          endif
          if (M_ini<0.01d0 .or. M_ini>150d0) then
              write(*,*) 'Error: runstarex_input_physics: set_initial_mass: M_ini out of range of 0.01 to 150 Msun: M_ini= ',M_ini
              write(*,*) 'Make sure the unit of M_ini is Solar Mass, Msun.'
              ierr = -1
              return
          endif
          s% initial_mass = M_ini
        end subroutine set_initial_mass

    ! Set the initial chemical composition, for pre-main sequence models (M. Michielsen)
        subroutine set_XYZ_ini_preMS(id, X_ini, Y_ini, Z_ini, ierr)
          integer, intent(in) :: id
          real(dp), intent(in) :: X_ini, Y_ini, Z_ini
          integer, intent(out) :: ierr
          type(star_info), pointer :: s
          ! h2/h1 and he3/he4 mass fraction ratio based on Lodders09 in data/chem_data/lodders09.data
          real(dp), parameter :: ratio_h2_to_h1 = 2.7556428657d-5/7.0945477357d-1, &
                               ratio_he3_to_he4 = 8.4641515456d-5/2.7501644504d-1

          call star_ptr(id, s, ierr)
          if (ierr /= 0) then
              write(*,*) 'Error: runstarex_input_physics: set_XYZ_ini: star_ptr failed'
              return
          endif
          s% job% set_uniform_initial_composition = .true.

          s% initial_y = Y_ini
          s% initial_z = Z_ini
          s% kap_rq% Zbase = Z_ini  ! set Zbase
          s% initial_he3 = Y_ini * ratio_he3_to_he4

          ! What about these? MESA crashes if they are not set. But job and controls both have initial_he3.
          s% job% initial_h1  = X_ini * (1d0 - ratio_h2_to_h1)
          s% job% initial_h2  = X_ini * ratio_h2_to_h1
          s% job% initial_he3 = Y_ini * ratio_he3_to_he4
          s% job% initial_he4 = Y_ini * (1d0 - ratio_he3_to_he4)
          s% job% dump_missing_metals_into_heaviest = .false.
        end subroutine set_XYZ_ini_preMS

        ! set_logD_mix() should be called only from extras_startup(), otherwise, the mass coordinate of the convective core is messed up (M. Michielsen)
        subroutine set_logD_mix(id, logD, ierr)
            integer, intent(in) :: id
            real(dp), intent(in) :: logD
            integer, intent(out) :: ierr

            type(star_info), pointer :: s

            call star_ptr(id, s, ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: grid_io: set_logD_mix: star_ptr failed'
                return
            endif

            if (logD .le. 1d-6) then
                s% set_min_D_mix = .false.
                s% min_D_mix = -1d99
            else
                s% set_min_D_mix = .true.
                s% min_D_mix = 10d0**logD
                write(*,*)'Mixing set:', logD
            endif
        end subroutine set_logD_mix

        ! Routine from github/mathiasm (contribution by M. Pedersen)
        ! Replace the constant min_D_mix by a profile going according to a power of the density
        ! Based on the internal gravity wave (IGW) mixing profile of Rogers & McElwaine (2017),
        ! see Pedersen et al. (2018) for more details
        subroutine IGW_D_mix_rho(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            real(dp) :: rho_switch, D_ext, new_Dmix, igw_exponent
            integer :: k, k_DmixMin
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! If these statements are not true, then no extra mixing is applied in the envelope.
            ! IGW profile will not be implemented if MESA is doing the startup of a run
            if ((.not. s% doing_first_model_of_run) .and. (.not. s% doing_relax ) .and. (s% mixing_type(s%nz) .eq. convective_mixing)) then
                !------- Determining position where to switch to IGW profile -------
                ! Going from core towards the surface of the star/model, check when Dmix <= s% min_D_mix
                ! and store the index at this position to be used for rescaling of the x-axis below
                D_ext = s% min_D_mix
                do k = s% nz, 1, -1
                    if (s% D_mix(k) <= D_ext) then
                        k_DmixMin = k
                        exit
                    end if
                end do

                igw_exponent = s% x_ctrl(1)
                if (k_DmixMin > 0) then
                    ! Rescale the y-axis of the loaded Dmix profile according to the min_D_mix
                    rho_switch = s% Rho(k_DmixMin)
                    do k = 1, k_DmixMin
                        new_Dmix = D_ext * pow(s% Rho(k) / rho_switch, igw_exponent)
                        if (new_Dmix > s% D_mix(k)) then
                            s% D_mix(k) = new_Dmix
                            s% mixing_type(k) = anonymous_mixing
                        end if
                    end do
                end if
            end if
        end subroutine IGW_D_mix_rho

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Saving routines
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Trigger a retry if we've gone over our goal Xc with a manually adapted timestep depending on if we overshot or undershot Xc.
        subroutine save_at_Xc_retry_method(id, Xc_save, Xc_precision, make_save, do_retry, ierr)
            integer, intent(in) :: id
            real(dp), intent(in) :: Xc_save, Xc_precision
            logical, intent(out) :: make_save, do_retry
            integer, intent(out) :: ierr
            type(star_info), pointer :: s

            call star_ptr(id, s, ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: grid_io: set_logD_mix: star_ptr failed'
                return
            endif

            make_save = .false.
            do_retry  = .false.
            if (ABS(s% center_h1 - Xc_save) < Xc_precision) then
                make_save = .true.
                s% x_ctrl(7) = Xc_save
                s% x_logical_ctrl(3) = .false.
                s% x_integer_ctrl(1) = 0
                s% x_ctrl(8) = -1d0
                s% timestep_factor_for_retries = 0.5d0 ! relinquish control of timestep
            endif

            ! retry if target Xc was missed
            ! Check if we're in the process of converging on an Xc_save. If so, use the bisection method
            if (s% x_logical_ctrl(3) .and. .not. make_save) then
                do_retry = .true.
                s% timestep_factor_for_retries = 1.0d0 ! Take control of timestep
                s% x_integer_ctrl(1) = s% x_integer_ctrl(1) + 1
                if (s% center_h1 < Xc_save) then
                    write(*,*) 'overshot Xc_save, reducing dt'
                    s% dt = s% dt - s% x_ctrl(8) * pow(2d0, - s% x_integer_ctrl(1))
                else if (s% center_h1 > Xc_save) then
                    write(*,*) 'undershot Xc_save, increasing dt'
                    s% dt = s% dt + s% x_ctrl(8) * pow(2d0, - s% x_integer_ctrl(1))
                endif
                write(*,*) 'Retry: bisection iteration ', s% x_integer_ctrl(1)
                return
            else if (s% center_h1 < Xc_save .and. .not. make_save) then
                s% x_ctrl(8) = s% dt
                s% timestep_factor_for_retries = 1.0d0 ! Take control of timestep
                s% x_logical_ctrl(3) = .true.
                s% x_integer_ctrl(1) = 1
                s%dt = 0.5d0*s%dt ! If you overshoot, we always halve the timestep in the first iteration
                do_retry = .true.
                write(*,*) 'Retry: overshot Xc_save, bisection iteration', s% x_integer_ctrl(1)
                return
            endif
        end subroutine save_at_Xc_retry_method

        ! save a profile with the model number as identifier (by Siemen and Mathijs)
        subroutine save_profile_model_number(id, prefix, ierr)
            integer, intent(in) :: id
            character(len=*), intent(in) :: prefix
            integer, intent(out) :: ierr
            character(len=256) :: full_path, suffix

            type(star_info), pointer :: s
            call star_ptr(id, s, ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: runstarex_numerics: save_profile_model_number: star_ptr failed'
                return
            endif

            write(suffix, '(A, F7.5, A, i0)') 'Xc', s% x_ctrl(7), '_mn', s% model_number

            full_path = trim(s% log_directory) // '/profiles/' // trim(prefix) // '_' // trim(suffix) // '.prof'

            call star_write_profile_info(id, full_path, ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: runstarex_numerics: save_profile_model_number: star_write_profile_info failed'
                return
            endif
        end subroutine save_profile_model_number


        ! save a pulsation file with model number as identifier (by Siemen Burssens, adapted by Mathijs Vanrespaille)
        subroutine save_gyre_model_number(id, prefix, keep_surf, add_atm, ierr)
            integer, intent(in) :: id
            character(len=*), intent(in) :: prefix
            logical, intent(in), optional, value :: keep_surf, add_atm
            integer, intent(out) :: ierr

            logical :: keep_surface_point, add_atmosphere, add_center_point
            character(len=256) :: full_path, suffix

            type(star_info), pointer :: s
            call get_star_ptr(id, s, ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: runstarex_numerics: save_gyre: get_star_ptr failed'
                return
            endif

            write(suffix, '(A, F7.5, A, i0)') 'Xc', s% x_ctrl(7), '_mn', s% model_number

            full_path = trim(s% log_directory) // '/gyre/input_models/' // trim(prefix) // '_' // trim(suffix) // '.GYRE'

            if (keep_surf) then
                keep_surface_point = .true.
            else
                keep_surface_point = .false.  ! default
            endif

            if (add_atm) then
                add_atmosphere = .true. ! atm_option = 'table' doesn't allow to add the atmosphere (does it?)
            else
                add_atmosphere = .false.  ! default
            endif

          add_center_point = .true.       ! default

          call star_export_pulse_data(id, 'GYRE', full_path, add_center_point, keep_surface_point, add_atmosphere, ierr)

            if (ierr /= 0) then
                write(*,*) 'Error: runstarex_numerics: save_gyre: star_write_gyre failed'
                return
            endif
        end subroutine save_gyre_model_number

        subroutine hist_asymptotic_dP(id, names, vals, vals_nr, ierr)
            ! author: M. Michielsen
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            character (len=maxlen_history_column_name), intent(inout) :: names(:)
            real(dp), intent(inout) :: vals(:)
            double precision, allocatable :: brunt_N(:)
            double precision :: integral
            integer :: k
            integer, intent(inout) :: vals_nr

            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! Preparing the Brunt-Vaisala frequency N in a separate variable
            ! with units "1/day", and selecting the regions with the positive
            ! values of N^2, i.e. the gravity mode propagation region
            allocate(brunt_N(s%nz))

            do k = 1, s%nz
              brunt_N(k) = (86400.0d0/(2.0d0*pi))*sqrt(max(0.0d0,s% brunt_N2(k)))
            end do

            ! Compute the integral int{(N/r)*dr}
            integral = 0.0d0
            do k = 1, s%nz - 1
              integral = integral + 0.5d0*(brunt_N(k) + brunt_N(k+1))*abs(s% r(k+1) - s% r(k)) / s% r(k)
            end do

            ! Compute the asymptotic period spacing without the spherical degree information in units of per second
            names(vals_nr) = "Asymptotic_dP"
            vals(vals_nr) = pi / integral*86400.0d0
            vals_nr = vals_nr+1
            deallocate(brunt_N)
        end subroutine hist_asymptotic_dp


        subroutine check_reached_zams(id, reached_zams, ierr)
            integer, intent(in) :: id
            logical, intent(out) :: reached_zams
            integer, intent(out) :: ierr
            type(star_info), pointer :: s

            call star_ptr(id, s, ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: check_reached_zams: star_ptr failed'
                return
            endif

            ! Reached zams if nuclear burning (L_nuc_burn_total) accounts for 99% of the total luminosity (log_L_surf)
            if  (abs(safe_log10(s% L_nuc_burn_total) - s% log_surface_luminosity) <=  2d0 .and. (s%mixing_type(s%nz) .eq. convective_mixing)) then
                reached_zams = .true.
                return
            endif

            reached_zams = .false.
            return
        end subroutine check_reached_zams

!********************************************************************************************
! routines for saving and restoring extra data so can do restarts
!********************************************************************************************

        subroutine store_extra_info(s)
           integer, parameter :: extra_info_put = 3
           type (star_info), pointer :: s
           call move_extra_info(s,extra_info_put)
        end subroutine store_extra_info

        subroutine move_extra_info(s,op)
           integer, parameter :: extra_info_alloc = 1
           integer, parameter :: extra_info_get = 2
           integer, parameter :: extra_info_put = 3
           type (star_info), pointer :: s
           integer, intent(in) :: op
           integer :: i, j, num_ints, num_dbls, ierr

           i = 0
           num_ints = i
           num_dbls = i

           if (op /= extra_info_alloc) return
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
        end subroutine move_extra_info

end module run_star_extras
