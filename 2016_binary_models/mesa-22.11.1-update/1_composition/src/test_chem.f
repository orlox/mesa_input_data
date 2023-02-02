      module mod_test_chem
      use chem_lib
      use chem_def
      use const_def, only: dp
      
      implicit none
      
      contains
      
      
      subroutine do_get_composition
         use const_lib, only: const_init
         integer, parameter :: iounit_LMC = 33
         integer, parameter :: iounit_SMC = 34
         integer, parameter :: iounit_GAL = 35
         integer, parameter :: num_elements = 6
         integer :: ierr, i, j
         character (len=64) :: my_mesa_dir 
         real(dp) :: X_LMC, X_SMC, X_GAL, Y_LMC, Y_SMC, Y_GAL, Z_LMC, Z_SMC, Z_GAL
         real(dp) :: Zn_LMC, Zn_SMC, Zn_GAL
         real(dp) :: Zfill_LMC, Zfill_SMC, Zfill_GAL
         ! elements that will be properly computed
         integer :: elements(num_elements)
         ! extra "fillup" element. One that does not react (not included in net), but is
         ! used to account all the extra non reacting metals.
         integer :: fillup_element
         real(dp) :: AGS05(num_chem_elements)
         real(dp) :: LMC_INES(num_chem_elements)
         real(dp) :: SMC_INES(num_chem_elements)
         real(dp) :: GAL_INES(num_chem_elements)

         include 'formats'

         call GET_ENVIRONMENT_VARIABLE( "MESA_DIR", my_mesa_dir, status=ierr, trim_name=.true.)
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            stop 1
         end if        
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: failed in chem_init'
            stop 1
         end if

         ! set elements that will go into the output file, together with fillup element
         elements(1) = e_c
         elements(2) = e_n
         elements(3) = e_o
         elements(4) = e_ne
         elements(5) = e_mg
         elements(6) = e_fe
         fillup_element = e_ca

         AGS05(:) = 0
         LMC_INES(:) = 0
         SMC_INES(:) = 0
         GAL_INES(:) = 0

         ! first store log abundances from the paper (photosphere unless otherwise noted)
         ! relative to log abundance of H = 12.00
         AGS05(e_li) = 3.25 !meteor
         AGS05(e_be) = 1.38
         AGS05(e_b ) = 2.70
         AGS05(e_c ) = 8.39
         AGS05(e_n ) = 7.78
         AGS05(e_o ) = 8.66
         AGS05(e_f ) = 4.56
         AGS05(e_ne) = 7.84 !indirect
         AGS05(e_na) = 6.17
         AGS05(e_mg) = 7.53
         AGS05(e_al) = 6.37
         AGS05(e_si) = 7.51
         AGS05(e_p ) = 5.36
         AGS05(e_s ) = 7.14
         AGS05(e_cl) = 5.50
         AGS05(e_ar) = 6.18 !indirect
         AGS05(e_k ) = 5.08
         AGS05(e_ca) = 6.31
         AGS05(e_sc) = 3.05
         AGS05(e_ti) = 4.90
         AGS05(e_v ) = 4.00
         AGS05(e_cr) = 5.64
         AGS05(e_mn) = 5.39
         AGS05(e_fe) = 7.45
         AGS05(e_co) = 4.92
         AGS05(e_ni) = 6.23
         AGS05(e_cu) = 4.21
         AGS05(e_zn) = 4.60
         AGS05(e_ga) = 2.88
         AGS05(e_ge) = 3.58
         AGS05(e_as) = 2.29 !meteor
         AGS05(e_se) = 3.33 !meteor
         AGS05(e_br) = 2.56 !meteor
         AGS05(e_kr) = 3.28 !indirect
         AGS05(e_rb) = 2.60
         AGS05(e_sr) = 2.92
         AGS05(e_y ) = 2.21
         AGS05(e_zr) = 2.59
         AGS05(e_nb) = 1.42
         AGS05(e_mo) = 1.92
         AGS05(e_Ru) = 1.84
         AGS05(e_Rh) = 1.12
         AGS05(e_Pd) = 1.69
         AGS05(e_Ag) = 0.94
         AGS05(e_Cd) = 1.77
         AGS05(e_In) = 1.60
         AGS05(e_Sn) = 2.00
         AGS05(e_Sb) = 1.00
         AGS05(e_Te) = 2.19 !meteor
         AGS05(e_I ) = 1.51 !meteor
         AGS05(e_Xe) = 2.27 !indirect
         AGS05(e_Cs) = 1.07 !meteor
         AGS05(e_Ba) = 2.17
         AGS05(e_La) = 1.13
         AGS05(e_Ce) = 1.58
         AGS05(e_Pr) = 0.71
         AGS05(e_Nd) = 1.45
         AGS05(e_Sm) = 1.01
         AGS05(e_Eu) = 0.52
         AGS05(e_Gd) = 1.12
         AGS05(e_Tb) = 0.28
         AGS05(e_Dy) = 1.14
         AGS05(e_Ho) = 0.51
         AGS05(e_Er) = 0.93
         AGS05(e_Tm) = 0.00
         AGS05(e_Yb) = 1.08
         AGS05(e_Lu) = 0.06
         AGS05(e_Hf) = 0.88
         AGS05(e_Ta) = -0.17 !meteor
         AGS05(e_W ) = 1.11
         AGS05(e_Re) = 0.23 !meteor
         AGS05(e_Os) = 1.45
         AGS05(e_Ir) = 1.38
         AGS05(e_Pt) = 1.64 !meteor
         AGS05(e_Au) = 1.01
         AGS05(e_Hg) = 1.13 !meteor
         AGS05(e_Tl) = 0.90
         AGS05(e_Pb) = 2.00
         AGS05(e_Bi) = 0.65 !meteor
         AGS05(e_Th) = 0.06 !meteor
         AGS05(e_U) = -0.52

         ! get modified abundances
         do i = e_li, e_u
            LMC_INES(i) = AGS05(i) - 0.4d0
            SMC_INES(i) = AGS05(i) - 0.7d0
            GAL_INES(i) = AGS05(i)
         end do
         !C
         LMC_INES(e_c) = 7.75
         SMC_INES(e_c) = 7.37
         GAL_INES(e_c) = 8.13
         !N
         LMC_INES(e_n) = 6.90
         SMC_INES(e_n) = 6.50
         GAL_INES(e_n) = 7.64
         !O
         LMC_INES(e_o) = 8.35
         SMC_INES(e_o) = 7.98
         GAL_INES(e_o) = 8.55
         !Mg
         LMC_INES(e_mg) = 7.05
         SMC_INES(e_mg) = 6.72
         GAL_INES(e_mg) = 7.32
         !Si
         LMC_INES(e_si) = 7.20
         SMC_INES(e_si) = 6.80
         GAL_INES(e_si) = 7.41
         !Fe
         LMC_INES(e_fe) = 7.05
         SMC_INES(e_fe) = 6.78
         GAL_INES(e_fe) = 7.40

         !turn into mass fractions
         ! we must satisfy three equations:
         ! 1 = X+Y+Z
         ! Y = 0.2477 + (0.28-0.2477)Z/Zsun = 0.2477 + 1.9 Z
         ! Z = X * alpha
         ! where alpha is given by
         ! alpha = sum_(i=e_li,e_u) m_i/m_H * 10^(abundance_i-12)
         ! from which Z is
         ! Z = 0.7523 / (alpha^(-1)+2.9)
         ! and X is
         ! X = 0.7523 / (1 + 2.9*alpha)
         Z_LMC = 0
         Z_SMC = 0
         Z_GAL = 0
         ! this computes alpha and stores it in Z_LMC, Z_SMC, Z_GAL
         do i = e_li, e_u
            if (element_atomic_weight(i) == 0.0) cycle
            LMC_INES(i) = element_atomic_weight(i) / element_atomic_weight(e_h) * 10**(LMC_INES(i) - 12d0)
            Z_LMC = Z_LMC + LMC_INES(i)
            SMC_INES(i) = element_atomic_weight(i) / element_atomic_weight(e_h) * 10**(SMC_INES(i) - 12d0)
            Z_SMC = Z_SMC + SMC_INES(i)
            GAL_INES(i) = element_atomic_weight(i) / element_atomic_weight(e_h) * 10**(GAL_INES(i) - 12d0)
            Z_GAL = Z_GAL + GAL_INES(i)
         end do

         !update X using X = 0.7523 / (1+2.9*alpha)
         X_LMC = 0.7523 / (1.0+2.9*Z_LMC)
         X_SMC = 0.7523 / (1.0+2.9*Z_SMC)
         X_GAL = 0.7523 / (1.0+2.9*Z_GAL)

         !update Z using Z = 0.7523 / (alpha^(-1)+2.9)
         Z_LMC = 0.7523 / (1.0/Z_LMC+2.9)
         Z_SMC = 0.7523 / (1.0/Z_SMC+2.9)
         Z_GAL = 0.7523 / (1.0/Z_GAL+2.9)

         !update Y using Y = 0.2477 + 1.9 Z 
         Y_LMC = 0.2477 + 1.9*Z_LMC
         Y_SMC = 0.2477 + 1.9*Z_SMC
         Y_GAL = 0.2477 + 1.9*Z_GAL

         write(*,*) 'LMC X,Y,Z,sum', X_LMC, Y_LMC, Z_LMC, X_LMC+Y_LMC+Z_LMC
         write(*,*) 'SMC X,Y,Z,sum', X_SMC, Y_SMC, Z_SMC, X_SMC+Y_SMC+Z_SMC
         write(*,*) 'GAL X,Y,Z,sum', X_GAL, Y_GAL, Z_GAL, X_GAL+Y_GAL+Z_GAL

         ! update individual metal abundances, recompute Z from here so sum of all adds to 1
         Z_LMC = 0
         Z_SMC = 0
         Z_GAL = 0
         do i = e_li, e_u
            if (element_atomic_weight(i) == 0.0) cycle
            LMC_INES(i) = X_LMC*LMC_INES(i)
            Z_LMC = Z_LMC + LMC_INES(i)
            SMC_INES(i) = X_SMC*SMC_INES(i)
            Z_SMC = Z_SMC + SMC_INES(i)
            GAL_INES(i) = X_GAL*GAL_INES(i)
            Z_GAL = Z_GAL + GAL_INES(i)
         end do

         ! just to be sure compute Y so abundances add up to 1 within MP
         Y_LMC = 1.0 - X_LMC - Z_LMC
         Y_SMC = 1.0 - X_SMC - Z_SMC
         Y_GAL = 1.0 - X_GAL - Z_GAL

         write(*,*) 'final LMC X,Y,Z,sum', X_LMC, Y_LMC, Z_LMC, X_LMC+Y_LMC+Z_LMC
         write(*,*) 'final SMC X,Y,Z,sum', X_SMC, Y_SMC, Z_SMC, X_SMC+Y_SMC+Z_SMC
         write(*,*) 'final GAL X,Y,Z,sum', X_GAL, Y_GAL, Z_GAL, X_GAL+Y_GAL+Z_GAL

         !add up elements that go into fillup
         Zfill_LMC = 0
         Zfill_SMC = 0
         Zfill_GAL = 0
         element_loop: do i = e_li, e_u
            if (element_atomic_weight(i) == 0.0) cycle
            do j = 1, num_elements
               if (i==elements(j)) cycle element_loop
            end do
            Zfill_LMC = Zfill_LMC + LMC_INES(i)
            Zfill_SMC = Zfill_SMC + SMC_INES(i)
            Zfill_GAL = Zfill_GAL + GAL_INES(i)
         end do element_loop

         ! add up relative number abundances for metals
         Zn_LMC = 0
         Zn_SMC = 0
         Zn_GAL = 0
         do i = e_li, e_u
            if (element_atomic_weight(i) == 0.0) cycle
            Zn_LMC = Zn_LMC + LMC_INES(i) / element_atomic_weight(i)
            Zn_SMC = Zn_SMC + SMC_INES(i) / element_atomic_weight(i)
            Zn_GAL = Zn_GAL + GAL_INES(i) / element_atomic_weight(i)
         end do

         !write down xa files for MESA
         ierr = 0
         open(unit=iounit_LMC, file=trim('xa_LMC.data'), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open file xa_LMC.data to write'
            stop 1
         end if
         open(unit=iounit_SMC, file=trim('xa_SMC.data'), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open file xa_SMC.data to write'
            stop 1
         end if
         open(unit=iounit_GAL, file=trim('xa_GAL.data'), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open file xa_GAL.data to write'
            stop 1
         end if
         write(iounit_LMC,*) "!File contains abundances per mass of each element"
         write(iounit_LMC,*) "!commented out are the abundances in 12+[X/H] format (marked 'abund=')"
         write(iounit_LMC,*) "!and also relative metal abundances per number (marked 'by number='), which can be inputted"
         write(iounit_LMC,*) "!in the opal website to produce custom tables."
         write(iounit_LMC,*) "!Relative metal abundances by mass are marked 'by mass='."
         write(iounit_LMC,*) "!As not all elements are included in the network,"
         write(iounit_LMC,*) "!we add up all non-included abundances in a single inert one."
         write(iounit_SMC,*) "!File contains abundances per mass of each element"
         write(iounit_SMC,*) "!commented out are the abundances in 12+[X/H] format (marked 'abund=')"
         write(iounit_SMC,*) "!and also relative metal abundances per number (marked 'by number=', which can be inputted"
         write(iounit_SMC,*) "!in the opal website to produce custom tables."
         write(iounit_SMC,*) "!Relative metal abundances by mass are marked 'by mass='."
         write(iounit_SMC,*) "!As not all elements are included in the network,"
         write(iounit_SMC,*) "!we add up all non-included abundances in a single inert one."
         write(iounit_GAL,*) "!File contains abundances per mass of each element"
         write(iounit_GAL,*) "!commented out are the abundances in 12+[X/H] format (marked 'abund=')"
         write(iounit_GAL,*) "!and alsorelative metal abundances per number (marked 'by number=', which can be inputted"
         write(iounit_GAL,*) "!in the opal website to produce custom tables."
         write(iounit_GAL,*) "!Relative metal abundances by mass are marked 'by mass='."
         write(iounit_GAL,*) "!As not all elements are included in the network,"
         write(iounit_GAL,*) "!we add up all non-included abundances in a single inert one."

         write(iounit_LMC,*) chem_element_main_iso_name(e_h), X_LMC
         write(iounit_SMC,*) chem_element_main_iso_name(e_h), X_SMC
         write(iounit_GAL,*) chem_element_main_iso_name(e_h), X_GAL
         write(iounit_LMC,*) chem_element_main_iso_name(e_he), Y_LMC
         write(iounit_SMC,*) chem_element_main_iso_name(e_he), Y_SMC
         write(iounit_GAL,*) chem_element_main_iso_name(e_he), Y_GAL
         element_loop2: do i = e_li, e_u
            if (element_atomic_weight(i) == 0.0) cycle
            do j = 1, num_elements
               if (i==elements(j)) then
                  write(iounit_LMC,*) chem_element_main_iso_name(i), LMC_INES(i)
                  write(iounit_LMC,*) "!abund=",12.0+log10(LMC_INES(i)/X_LMC*element_atomic_weight(e_h)/element_atomic_weight(i))
                  write(iounit_LMC,*) "!by number=", chem_element_main_iso_name(i), LMC_INES(i) / element_atomic_weight(i)/Zn_LMC
                  write(iounit_LMC,*) "!by mass=", chem_element_main_iso_name(i), LMC_INES(i) / Z_LMC
                  write(iounit_SMC,*) chem_element_main_iso_name(i), SMC_INES(i)
                  write(iounit_SMC,*) "!abund=",12.0+log10(SMC_INES(i)/X_SMC*element_atomic_weight(e_h)/element_atomic_weight(i))
                  write(iounit_SMC,*) "!by number=", chem_element_main_iso_name(i), SMC_INES(i) / element_atomic_weight(i)/Zn_SMC
                  write(iounit_SMC,*) "!by mass=", chem_element_main_iso_name(i), SMC_INES(i) / Z_SMC
                  write(iounit_GAL,*) chem_element_main_iso_name(i), GAL_INES(i)
                  write(iounit_GAL,*) "!abund=",12.0+log10(GAL_INES(i)/X_GAL*element_atomic_weight(e_h)/element_atomic_weight(i))
                  write(iounit_GAL,*) "!by number=", chem_element_main_iso_name(i), GAL_INES(i) / element_atomic_weight(i)/Zn_GAL
                  write(iounit_GAL,*) "!by mass=", chem_element_main_iso_name(i), GAL_INES(i) / Z_GAL
                  cycle element_loop2
               end if
            end do
            if (i == fillup_element) then
               write(iounit_LMC,*) chem_element_main_iso_name(i), Zfill_LMC
               write(iounit_LMC,*) "!this element is fillup, original value is", LMC_INES(i)
               write(iounit_LMC,*) "!abund=",12.0+log10(LMC_INES(i)/X_LMC*element_atomic_weight(e_h)/element_atomic_weight(i))
               write(iounit_LMC,*) "!by number=", chem_element_main_iso_name(i), LMC_INES(i) / element_atomic_weight(i)/Zn_LMC
               write(iounit_LMC,*) "!by mass=", chem_element_main_iso_name(i), LMC_INES(i) / Z_LMC
               write(iounit_SMC,*) chem_element_main_iso_name(i), Zfill_SMC
               write(iounit_SMC,*) "!this element is fillup, original value is", SMC_INES(i)
               write(iounit_SMC,*) "!abund=",12.0+log10(SMC_INES(i)/X_SMC*element_atomic_weight(e_h)/element_atomic_weight(i))
               write(iounit_SMC,*) "!by number=", chem_element_main_iso_name(i), SMC_INES(i) / element_atomic_weight(i)/Zn_SMC
               write(iounit_SMC,*) "!by mass=", chem_element_main_iso_name(i), SMC_INES(i) / Z_SMC
               write(iounit_GAL,*) chem_element_main_iso_name(i), Zfill_GAL
               write(iounit_GAL,*) "!this element is fillup, original value is", GAL_INES(i)
               write(iounit_GAL,*) "!abund=",12.0+log10(GAL_INES(i)/X_GAL*element_atomic_weight(e_h)/element_atomic_weight(i))
               write(iounit_GAL,*) "!by number=", chem_element_main_iso_name(i), GAL_INES(i) / element_atomic_weight(i)/Zn_GAL
               write(iounit_GAL,*) "!by mass=", chem_element_main_iso_name(i), GAL_INES(i) / Z_GAL
            else
               write(iounit_LMC,*) "!", chem_element_main_iso_name(i), LMC_INES(i) 
               write(iounit_LMC,*) "!abund=",12.0+log10(LMC_INES(i)/X_LMC*element_atomic_weight(e_h)/element_atomic_weight(i))
               write(iounit_LMC,*) "!by number=", chem_element_main_iso_name(i), LMC_INES(i) / element_atomic_weight(i)/Zn_LMC
               write(iounit_LMC,*) "!by mass=", chem_element_main_iso_name(i), LMC_INES(i) / Z_LMC
               write(iounit_SMC,*) "!", chem_element_main_iso_name(i), SMC_INES(i)
               write(iounit_SMC,*) "!abund=",12.0+log10(SMC_INES(i)/X_SMC*element_atomic_weight(e_h)/element_atomic_weight(i))
               write(iounit_SMC,*) "!by number=", chem_element_main_iso_name(i), SMC_INES(i) / element_atomic_weight(i)/Zn_SMC
               write(iounit_SMC,*) "!by mass=", chem_element_main_iso_name(i), SMC_INES(i) / Z_SMC
               write(iounit_GAL,*) "!", chem_element_main_iso_name(i), GAL_INES(i)
               write(iounit_GAL,*) "!abund=",12.0+log10(GAL_INES(i)/X_GAL*element_atomic_weight(e_h)/element_atomic_weight(i))
               write(iounit_GAL,*) "!by number=", chem_element_main_iso_name(i), GAL_INES(i) / element_atomic_weight(i)/Zn_GAL
               write(iounit_GAL,*) "!by mass=", chem_element_main_iso_name(i), GAL_INES(i) / Z_GAL
            end if
         end do element_loop2
         close(iounit_LMC)
         close(iounit_SMC)
         close(iounit_GAL)
         
      end subroutine do_get_composition
      
      end module mod_test_chem
      
      program test_chem
      use mod_test_chem
      call do_get_composition
      end program




