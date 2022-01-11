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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "cs_user_coupling.h"
#include "cs_luma_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_coupling.c
 *
 * \brief Code couplings definition with SYRTHES and Code_Saturne.
 *
 * See \ref user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with other instances of Code_Saturne.
 *
 * This is done by calling the \ref cs_sat_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_saturne_coupling
void
cs_user_saturne_coupling(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with LUMA code.
 *
 * This is done by calling the \ref cs_luma_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

 /*----------------------------------------------------------------------------*/
/*!
 * \brief Define new LUMA coupling.
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
 * \param[in] vars_in           initial of the variables to write into CS: 
                                "v" velocity, "t" temperature. So "vt" is velocity and temperature. 
 * \param[in] vars_out          initial of the variables to read from CS and transfer to LUMA:
                                "v" velocity, "t" temperature. So "vt" is velocity and temperature. 
 *
 * In the case of a single Code_Saturne and single LUMA instance, the
 * 'luma_name' argument is ignored, as there is only one matching
 * possibility.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * LUMA instances based on the 'luma_name' argument.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_luma_coupling
void
cs_user_luma_coupling(void)
{
	// variables to couple
	// v = velocity
	// t = temperature

	int verbosity = 1;
	cs_luma_coupling_define("LEFT", "Left",  "box[0.54, 0., 0., 0.55, 1.0, 0.2]", 'L', true, 1.0, verbosity, 0, "v", "v");

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with SYRTHES code.
 *
 * This is done by calling the \ref cs_syr_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_syrthes_coupling
void
cs_user_syrthes_coupling(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a volume exchange coefficient for SYRTHES couplings.
 *
 * \param[in]   coupling_id   Syrthes coupling id
 * \param[in]   syrthes_name  name of associated Syrthes instance
 * \param[in]   n_elts        number of associated cells
 * \param[in]   elt_ids       associated cell ids
 * \param[out]  h_vol         associated exchange coefficient (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_syrthes_coupling_volume_h
void
cs_user_syrthes_coupling_volume_h(int               coupling_id,
                                  const char       *syrthes_name,
                                  cs_lnum_t         n_elts,
                                  const cs_lnum_t   elt_ids[],
                                  cs_real_t         h_vol[])
{
  CS_UNUSED(coupling_id);
  CS_UNUSED(syrthes_name);
  CS_UNUSED(n_elts);
  CS_UNUSED(elt_ids);
  CS_UNUSED(h_vol);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
