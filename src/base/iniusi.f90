!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2021 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine iniusi

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE APPELANT LES ROUTINES UTILISATEUR POUR L'ENTREE DES
!   PARAMETRES DE CALCUL : ON PASSE ICI POUR TOUT CALCUL

! CETTE ROUTINE PERMET DE CACHER A L'UTILISATEUR LES APPELS
!   A VARPOS ET AU LECTEUR XML DE L'IHM

! LE DECOUPAGE DE L'ANCIEN USINI1 PERMET EGALEMENT DE MIEUX
!   CONTROLER LES ZONES OU SONT INITIALISES LES VARIABLES (PAR
!   LE BIAIS DE PARAMETRES PASSES EN ARGUMENT)


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use albase
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use radiat
use cs_coal_incl
use cs_c_bindings
use cs_cf_bindings
use field
use cdomod

!===============================================================================

implicit none

! Arguments

! Local variables

integer          nmodpp
integer          nscmax
integer          l_size, f_id
double precision relaxp, l_cp(1), l_xmasm(1), l_cv(1)

type(var_cal_opt) :: vcopt

!===============================================================================

interface

  subroutine cs_gui_physical_properties()  &
       bind(C, name='cs_gui_physical_properties')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_physical_properties

  subroutine cs_gui_porous_model()  &
       bind(C, name='cs_gui_porous_model')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_porous_model

  subroutine cs_gui_radiative_transfer_parameters()  &
       bind(C, name='cs_gui_radiative_transfer_parameters')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_radiative_transfer_parameters

  subroutine cs_velocity_pressure_set_solid()  &
       bind(C, name='cs_velocity_pressure_set_solid')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_velocity_pressure_set_solid

  subroutine cs_runaway_check_define_field_max(f_id, value)  &
       bind(C, name='cs_runaway_check_define_field_max')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: f_id
    real(c_double), value :: value
  end subroutine cs_runaway_check_define_field_max

end interface

!===============================================================================

! Check for restart and read matching time steps

call parameters_read_restart_info

!===============================================================================
! 1. Initialize model settings
!===============================================================================

! Flow model selection through GUI

call cs_gui_physical_model_select

! Flow model selection through user Fortran subroutine

call usppmo(1)

! Other models selection through GUI

! ALE parameters
call uialin (nalinf, nalimx, epalim)

! thermal model
call csther

! turbulence model choice
call cs_gui_turb_model

! constant or variable specific heat
call cscpva

! Other models selection through user Fortran subroutine

call usipph(1, iturb, itherm, iale)

! Flow and other models selection through user C function
call cs_user_model

! Activate CDO for ALE
if (iale.eq.2) then
  call cs_ale_activate
endif

! Other model parameters, including user-defined scalars

call cs_gui_user_variables
call cs_gui_user_arrays

! Solid zones

call cs_velocity_pressure_set_solid

!===============================================================================
! 2. Initialize parameters for specific physics
!===============================================================================

call cfnmtd(ficfpp, len(ficfpp))

! --- Activation du module transferts radiatifs

!     Il est necessaire de connaitre l'activation du module transferts
!     radiatifs tres tot de maniere a pouvoir reserver les variables
!     necessaires dans certaines physiques particuliere

!   - Interface Code_Saturne
!     ======================

call cs_gui_radiative_transfer_parameters

! Define fields for variables, check and build iscapp
! and computes the number of user scalars (nscaus)
if (icdo.lt.2) then
  call fldvar(nmodpp)
endif

if (iale.ge.1) then
  call uialvm
endif
call csivis

nscmax = nscamx

! ---> Physique particuliere : darcy

if (ippmod(idarcy).ge.0) then
  call daini1
endif

call field_get_id_try('velocity', f_id)
if (f_id .ge. 0) then
  if (ippmod(icompf).ge.0) then
    call cs_runaway_check_define_field_max(f_id, 1.0d5)
  else
    call cs_runaway_check_define_field_max(f_id, 1.0d4)
  endif
endif

!===============================================================================
! 3. INITIALISATION DE PARAMETRES "GLOBAUX"
!===============================================================================

! --- Parametres globaux

!     Pas de temps
!     Couplage vitesse/pression
!     Prise en compte de la pression hydrostatique
!     Estimateurs (pas encore dans l'IHM)


!   - Interface Code_Saturne
!     ======================

call csidtv()

call csiphy()

! Postprocessing

call cspstb(ipstdv)

! Define main properties (pointers, checks, ipp) if not in CDO mode only
if (icdo.lt.2) then
   call fldprp
endif

!===============================================================================
! 4. INITIALISATION DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
!===============================================================================

! --- Format des fichiers aval (entsor.h)
! --- Options du calcul (optcal.h)
! --- Constantes physiques (cstphy.h)


!   - Interface Code_Saturne
!     ======================

! Restart, read auxiliary file, frozen velocity field

call csisui(ntsuit, iccvfg)

! Time step (only ntmabs, dtref)
call cstime()

! Local numerical options

call uinum1(cdtvar)

! If CDO mode only, no pressure is defined at this stage
if (icdo.lt.2) then

  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

  ! Global numeric options
  relaxp = -1.d0
  call csnum2 (relaxp, imrgra)
  if (idtvar.ge.0) vcopt%relaxv = relaxp

  call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

endif

! Gravity, physical properties
call cs_gui_physical_properties

! Turbulence reference values (uref, almax)
call cs_gui_turb_ref_values

! Set turbulence constants according to model choices.
! This can be overwritten by the user in cs_user_parameters()
call cs_f_turb_complete_constants

! Scamin, scamax, turbulent flux model
call cssca2()

! Diffusivities
call cssca3()

! Porosity model
call cs_gui_porous_model()

! Init fan
call uifans()

! Init error estimator
call uieres(iescal, iespre, iesder, iescor, iestot)

!   - User functions
!     ==============

call usipsu(nmodpp)
call user_parameters

! If time step is local or variable, pass information to C layer, as it
! may be needed for some field (or moment) definitions.
if (idtvar.ne.0) then
  call time_step_define_variable(1)
endif
if (idtvar.eq.2.or.idtvar.eq.-1) then
  call time_step_define_local(1)
endif

call indsui(isuite)

! Default values of physical properties for the compressible model
if (ippmod(icompf).ge.0) then
  ! ieos has been set above in uippmo with the GUI or in usppmo without the GUI.
  ! The variability of the thermal conductivity
  ! (diffusivity_id for itempk) and the volume viscosity (iviscv) has
  ! been set in fldprp.

  ! Here call to uscfx2 to get visls_0(itempk), viscv0, xmasmr, ivivar and
  ! psginf, gammasg, cv0 in stiffened gas thermodynamic.
  ! With GUI, visls_0(itempk), viscv0, xmasmr and ivivar have already been read
  ! above in the call to csphys.
  call uscfx2

  ! Compute cv0 according to chosen EOS.
  l_size = 1
  l_cp(1) = cp0 ! dummy argument in stiffened gas
  l_xmasm(1) = xmasmr ! dummy argument in stiffened gas
  call cs_cf_thermo_cv(l_cp, l_xmasm, l_cv, l_size)
  cv0 = l_cv(1)
endif

! Choose the porous model
call cs_porous_model_set_model(iporos)

! --- Varpos
! If CDO mode only, skip this stage
if (icdo.lt.2) then
  call varpos
endif
! --- Internal coupling
call cs_gui_internal_coupling
call cs_user_internal_coupling

call cs_internal_coupling_setup

!----
! Formats
!----

return
end subroutine
