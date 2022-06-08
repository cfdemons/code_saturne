#ifndef __CS_LUMA_COUPLING_H__
#define __CS_LUMA_COUPLING_H__

/*============================================================================
 * LUMA coupling
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

/* Structure associated to LUMA coupling */

typedef struct _cs_luma_coupling_t  cs_luma_coupling_t;

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the implicit/explicit contribution for source terms in a SYRTHES
 * volume coupling
 *
 * Fortran Interface:
 *
 * SUBROUTINE CTBVSY (NUMSYR, TFLUID, CTBIMP, CTBEXP)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION CTBIMP      : --> : Implicit contribution
 * DOUBLE PRECISION CTBEXP      : --> : Explicit contribution
 *----------------------------------------------------------------------------*/

void CS_PROCF (ctbvsy, CTBVSY)
(
 int        *numsyr,
 cs_real_t  *tfluid,
 cs_real_t  *ctbimp,
 cs_real_t  *ctbexp
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/
 
 /*----------------------------------------------------------------------------
 * Send coupling variables to LUMA
 *
 * parameters:
 *   luma_coupling     <-- LUMA coupling
 *   t_luma       <-- CS temperature, NULL if no coupled temperature
 *   v_luma       <-- CS velocity, NULL if no coupled velocity
 *   mode         <-- 0: surface coupling; 1: volume coupling  // I probably don't need this anymore
 *----------------------------------------------------------------------------*/
//TODO: Add a new argument with the id of the coupled entity in luma_coupling we want to access to. 
 
void
cs_luma_coupling_send_data(cs_luma_coupling_t      *luma_coupling,
							 cs_real_t            t_cs[],
							 cs_real_t            v_cs[],
                             int                  mode);
						  // int  i);

 /*----------------------------------------------------------------------------*/
/*!
 * \brief  Exchange cell values relative to a LUMA coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_luma_coupling_send_volume(void);
 
  /*----------------------------------------------------------------------------*/
/*!
 * \brief  Read boundary field/variable values relative to a LUMA coupling.
 *
 * \param[in]       nvar     number of variables
 * \param[in]       bc_type  boundary condition type
 * \param[in, out]  icodcl   boundary condition codes
 * \param[in, out]  rcodcl   boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_luma_coupling_recv_boundary(int        nvar,
                              int        bc_type[],
                              int        icodcl[],
                              cs_real_t  rcodcl[]);
 
  /*----------------------------------------------------------------------------
 * Receive coupling variables from LUMA
 *
 * parameters:
 *   luma_ent     <-- LUMA coupling entity
 *   t_luma       --> LUMA temperature, NULL if no coupled temperature
 *   v_luma       --> LUMA velocity, NULL if no coupled velocity
 *   mode         <-- 0: surface coupling; 1: volume coupling  // I probably don't need this anymore
 *----------------------------------------------------------------------------*/
//TODO: Add a new argument with the id of the coupled entity in luma_coupling we want to access to. 
 
void
cs_luma_coupling_recv_data(cs_luma_coupling_t   *luma_coupling,
							 cs_real_t            t_luma[],
							 cs_real_t            v_luma[],
                             int                  mode);

  /*----------------------------------------------------------------------------
 * Get local numbering of coupled elements
 *
 * parameters:
 *   luma_coupling  <-- LUMA coupling structure
 *   cpl_elt_ids   --> Ids of coupled elements (0 to n-1)
 *   mode          <-- 0 (surface); 1 (volume)
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_get_elt_ids(const cs_luma_coupling_t     *luma_coupling,
                             cs_lnum_t                  cpl_elt_ids[],
                             int                        mode);
 
  /*----------------------------------------------------------------------------
 * Get number of associated coupled elements in main mesh
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *   mode          <-- 0 (surface); 1 (volume)
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_luma_coupling_get_n_elts(const cs_luma_coupling_t *luma_coupling,
                            int                       mode);
							
/*----------------------------------------------------------------------------
 * Return 1 if this coupling is a volume coupling else 0
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *
 * returns:
 *   1 or 0
 *----------------------------------------------------------------------------*/

int
cs_luma_coupling_is_vol(const cs_luma_coupling_t  *luma_coupling);
 
  /*----------------------------------------------------------------------------
 * Return 1 if this coupling is a surface coupling else 0
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *
 * returns:
 *   1 or 0
 *----------------------------------------------------------------------------*/

int
cs_luma_coupling_is_surf(const cs_luma_coupling_t  *luma_coupling);
 
  /*----------------------------------------------------------------------------*/
/*
 * Create coupled meshes and setup PLE locator for LUMA couplings.
 */
/*----------------------------------------------------------------------------*/

void
cs_luma_coupling_init_meshes(void);
 
  /*----------------------------------------------------------------------------
 * Define coupled mesh and send it to LUMA
 *
 * Optional post-processing output is also built at this stage.
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_init_mesh(cs_luma_coupling_t  *luma_coupling);
 
  /*----------------------------------------------------------------------------*/
/*
 * Create coupled meshes and setup PLE locator for LUMA couplings.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_init_meshes(void);
 
 
  /*----------------------------------------------------------------------------
 * Get name of LUMA coupling.
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *
 * returns:
 *   pointer to LUMA coupling name
 *----------------------------------------------------------------------------*/

const char *
cs_luma_coupling_get_name(cs_luma_coupling_t  *luma_coupling);
 
  /*----------------------------------------------------------------------------
 * Initialize communicator for LUMA coupling
 *
 * parameters:
 *   luma_coupling  <-> LUMA coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *   luma_root_rank <-- LUMA root rank
 *   n_luma_ranks   <-- Number of ranks associated with LUMA
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_init_comm(cs_luma_coupling_t *luma_coupling,
                           int                 coupling_id,
                           int                 luma_root_rank,
                           int                 n_luma_ranks);
						   
 /*----------------------------------------------------------------------------
 * Get pointer to LUMA coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of LUMA coupling
 *
 * returns:
 *   pointer to LUMA coupling structure
 *----------------------------------------------------------------------------*/

cs_luma_coupling_t *
cs_luma_coupling_by_id(int  coupling_id);
 
 /*----------------------------------------------------------------------------
 * Get number of LUMA couplings.
 *
 * returns:
 *   number of LUMA couplings
 *----------------------------------------------------------------------------*/

int
cs_luma_coupling_n_couplings(void);
 
 /*----------------------------------------------------------------------------
 * Initialize LUMA couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_all_init(void);
 
 /*----------------------------------------------------------------------------*/
/*
 * Define new LUMA coupling.
 *
 * \param[in] luma_name         	 matching LUMA application name
 * \param[in] boundary_criteria 	 surface selection criteria, or NULL  // NOTE: I think this is the name of the boundary in CS. 
 * \param[in] volume_criteria   	 volume selection criteria, or NULL   // NOTE: I think this is hte name of the region containing the desired cells in CS.
 * \param[in] projection_axis   	 x', 'y', or 'y' for 2D projection axis (case
 *                              	 independent), or ' ' for standard 3D coupling // NOTE: I don't think I need this for LUMA
 * \param[in] allow_nonmatching 	 allow nearest-neighbor mapping where matching
 *                              	 within tolerance is not available (useful
 *                              	 when meshes have a different level of detail)
 * \param[in] tolerance         	 addition to local extents of each element
 *                              	 extent = base_extent * (1 + tolerance)
 * \param[in] verbosity         	 verbosity level
 * \param[in] visualization     	 visualization output level (0 or 1)
  * \param[in] vars_in               initial of the variables to write into CS: 
                                     "v" velocity, "t" temperature. So "vt" is velocity and temperature. 
 * \param[in] vars_out               initial of the variables to read from CS and transfer to LUMA:
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

void
cs_luma_coupling_define(const char  *luma_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       bool         allow_nonmatching,
                       float        tolerance,
                       int          verbosity,
                       int          visualization,
					   const char  *vars_in,
					   const char  *vars_out);

/*----------------------------------------------------------------------------
 * Create or redefine a luma_coupling_t structure.
 *
 * If a structure is redefined, associated locations are reset.
 *
 * parameters:
 *   dim                <-- spatial mesh dimension
 *   ref_axis           <-- reference axis
 *   luma_name          <-- LUMA application name
 *   allow_nonmatching  <-- nearest-neighbor search for non-matching faces flag
 *   tolerance          <-- addition to local extents of each element
 *                          extent = base_extent * (1 + tolerance)
 *   verbosity          <-- verbosity level
 *   visualization      <-- visualization output flag
 *	 vars_in            <-- first letter of the variables to write into CS
 *	 vars_out           <-- first letter of the variables to read from CS into LUMA
 *----------------------------------------------------------------------------*/

cs_luma_coupling_t *
cs_luma_coupling_create(int          dim,
                        int          ref_axis,
                        const char  *luma_name,
                        bool         allow_nonmatching,
                        float        tolerance,
                        int          verbosity,
                        int          visualization,
						const char  *vars_in,
						const char  *vars_out);

/*----------------------------------------------------------------------------
 * Add a mesh location to a luma_coupling_t structure.
 *
 * parameters:
 *   luma_coupling  <-- LUMA coupling structure
 *   location_id   <-- id of mesh location to add (boundary faces or cells)
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_add_location(cs_luma_coupling_t  *luma_coupling,
                              int                  location_id);

/*----------------------------------------------------------------------------
 * Destroy cs_luma_coupling_t structures
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_all_destroy(void);

/*----------------------------------------------------------------------------
 * Get name of LUMA coupling.
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *
 * returns:
 *   pointer to LUMA coupling name
 *----------------------------------------------------------------------------*/

//const char *
//cs_luma_coupling_get_name(cs_luma_coupling_t  *luma_coupling);

/*----------------------------------------------------------------------------
 * Set conservativity forcing flag to True (1) or False (0) for all defined
 * LUMA couplings
 *
 * parameter:
 *   flag     <--  Conservativity forcing flag to set
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_set_conservativity(int  flag);


/*----------------------------------------------------------------------------
 * Initialize communicator for LUMA coupling
 *
 * parameters:
 *   luma_coupling  <-> LUMA coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *   luma_root_rank <-- LUMA root rank
 *   n_luma_ranks   <-- Number of ranks associated with LUMA
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_init_comm(cs_luma_coupling_t *luma_coupling,
//                           int                 coupling_id,
//                           int                 luma_root_rank,
//                           int                 n_luma_ranks);

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to LUMA
 *
 * Optional post-processing output is also built at this stage.
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_init_mesh(cs_luma_coupling_t  *luma_coupling);


/*----------------------------------------------------------------------------
 * Get number of associated coupled elements in main mesh
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *   mode          <-- 0 (surface); 1 (volume)
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

//cs_lnum_t
//cs_luma_coupling_get_n_elts(const cs_luma_coupling_t *luma_coupling,
//                            int                       mode);

/*----------------------------------------------------------------------------
 * Get local numbering of coupled elements
 *
 * parameters:
 *   luma_coupling  <-- LUMA coupling structure
 *   cpl_elt_ids   --> List of coupled elements (0 to n-1)
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_get_elt_ids(const cs_luma_coupling_t  *luma_coupling,
//                             cs_lnum_t                  cpl_elt_ids[]);

/*----------------------------------------------------------------------------
 * Receive coupling variables from LUMA
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *   tsolid       --> solid temperature   // This needs to be changed to the name of the variable and the value of the variable. 
 *   mode         <-- 0: surface coupling; 1: volume coupling
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_recv_tsolid(cs_syr4_coupling_t  *syr_coupling,
//                             cs_real_t            tsolid[],
//                             int                  mode);

/*----------------------------------------------------------------------------
 * Send coupling variables to LUMA
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *   cpl_elt_ids  <-- ids of coupled elements
 *   tf           <-- fluid temperature         // SAME: This needs to be modified so that CS can chose from a group of variables to send to LUMA.
 *   hf           <-- fluid heat exchange coef. (numerical or user-defined)
 *   mode          <-- 0: surface coupling; 1: volume coupling
 *----------------------------------------------------------------------------*/

//void
//cs_luma_coupling_send_tf_hf(cs_syr4_coupling_t  *syr_coupling,
//                            const cs_lnum_t      cpl_elt_ids[],
//                            cs_real_t            tf[],
//                            cs_real_t            hf[],
//                           int                  mode);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LUMA_COUPLING_H__ */
