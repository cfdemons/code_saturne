/*============================================================================
 * Code couplings definition with SYRTHES and Code_Saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of Code_Saturne
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/



 /*----------------------------------------------------------------------------*/
/*!
 * Define new LUMA coupling.
 *
 * \param[in] luma_name         matching LUMA application name
 * \param[in] boundary_criteria surface selection criteria, or NULL  // NOTE: I think this is the name of the boundary in CS. 
 * \param[in] volume_criteria   volume selection criteria, or NULL   // NOTE: I think this is hte name of the region containing the desired cells in CS.
 * \param[in] projection_axis   x', 'y', or 'y' for 2D projection axis (case
 *                              independent), or ' ' for standard 3D coupling // NOTE: I don't think I need this for LUMA
 * \param[in] allow_nonmatching allow nearest-neighbor mapping where matching
 *                              within tolerance is not available (useful
 *                              when meshes have a different level of detail)
 * \param[in] tolerance         addition to local extents of each element
 *                              extent = base_extent * (1 + tolerance)
 * \param[in] verbosity         verbosity level
 * \param[in] visualization     visualization output level (0 or 1)
 *
 * In the case of a single Code_Saturne and single LUMA instance, the
 * 'luma_name' argument is ignored, as there is only one matching
 * possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * LUMA instances based on the 'luma_name' argument.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_luma_coupling(void);


