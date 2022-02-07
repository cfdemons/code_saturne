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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>
#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal.h"
#include "fvm_nodal_extract.h"
#include "fvm_nodal_project.h"
#include "fvm_selector.h"

#include "cs_coupling.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_selector.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"

#include "cs_field.h"
#include "cs_physical_model.h"
#include "cs_parameters.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_luma_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

const int  cs_luma_coupling_tag = 'C'+'S'+'_'+'C'+'O'+'U'+'P'+'L'+'I'+'N'+'G';

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Structure associated with LUMA coupling. Equivalent to a coupling interface */

typedef struct _cs_luma_coupling_ent_t {

  /* Mesh-related members */

  ple_locator_t     *locator;        /* Associated locator */

  int                elt_dim;        /* Element dimension */
  cs_lnum_t          n_elts;         /* Number of coupled elements */

  fvm_nodal_t       *elts;           /* Coupled elements */

  /* Saved arrays for post processing (float for reduced memory use) */

  int                post_mesh_id;   /* 0 if post-processing is not active,
                                        or post-processing mesh_id (< 0) */
										
  bool				is_temp_in;		/* Boolean to indicate which variables are coupled for this entity. Note: I guess this could be a bitmask.  */
  bool				is_temp_out;
  bool				is_vel_in;
  bool				is_vel_out;
  
  // For the moment only velocity and temperature can be coupled. Add more if needed. They will only contain data if indicated by the variable_names char.
  cs_real_t         *temp_in;     		/* Temperature received from LUMA */      
  cs_real_t         *vel_in;           /* Velocity received from LUMA. Note: check that velocity in CS is stored in a cs_real_t variable */
  
  double            *temp_out;     		/* Temperature sent to LUMA */      
  double            *vel_out;           /* Velocity sent to LUMA. Note: check that velocity in CS is stored in a cs_real_t variable */

 /* Saved array for volume coupling. I'm not sure I'll need it */
 
  double			*hvol;         /* Volumetric data?. */ 
  
  
  //float             *flux;           /* Flux (calculated) */
  //float             *tfluid_tmp;     /* Fluid temperature (points to flux in
  //                                      transient stage where solid_temp and
  //                                      fluid_temp are known, NULL once
  //                                      flux is calculated) */

} cs_luma_coupling_ent_t;

/* Structure associated with LUMA. Equivalent to the coupling adapter */

struct _cs_luma_coupling_t {

  /* Mesh-related members */

  int                      dim;             /* Coupled mesh dimension */
  int                      ref_axis;        /* Axis for edge extraction */

  char                    *luma_name;        /* Application name */
  char 					  *vars_in;		     /* First letter of the variables to be written into CS: v = velocity, t = temperature*/
  char					  *vars_out;         /* First letter of the variables to be read from CS and transferred to LUMA: v = velocity, t = temperature */ 

  int                      n_b_locations;   /* Number of boundary locations */
  int                      n_v_locations;   /* Number of volume locations */
  int                     *b_location_ids;  /* Boundary location ids */
  int                     *v_location_ids;  /* Volume location ids */

  cs_luma_coupling_ent_t  *faces;           /* Wall coupling structure */  // NOTE: This is an array of cs_luma_coupling_ent_t, so an array of "interfaces"
  cs_luma_coupling_ent_t  *cells;           /* Volume coupling structure */

  bool                     allow_nearest;   /* Allow nearest-neighbor
                                               mapping beyond basic matching
                                               tolerance */
  float                    tolerance;       /* Tolerance */
  int                      verbosity;       /* Verbosity level */
  int                      visualization;   /* Visualization output flag */

  /* Communication-related members */

#if defined(HAVE_MPI)

  MPI_Comm           comm;           /* Associated MPI communicator */

  int                n_luma_ranks;    /* Number of associated LUMA ranks */
  int                luma_root_rank;  /* First associated LUMA rank */

#endif

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static int                   cs_glob_luma_n_couplings = 0;
static cs_luma_coupling_t  **cs_glob_luma_couplings = NULL;

/* Start and end (negative) numbers associated with
   dedicated post processing meshes */

static int  cs_glob_luma_post_mesh_ext[2] = {0, 1};

//static int  cs_luma_coupling_implicit = 0;

/*============================================================================
 *  Private functions definitions
 *============================================================================*/

 /*----------------------------------------------------------------------------
 * Post process variables associated with LUMA couplings
 *
 * parameters:
 *   coupling        <--  Void pointer to LUMA coupling structure
 *   ts              <--  time step status structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_luma_coupling_post_function(void                  *coupling,
                                const cs_time_step_t  *ts)
{
  int type_id;

  const cs_luma_coupling_t  *luma_coupling = coupling;
  cs_luma_coupling_ent_t *coupling_ent = NULL;

  for (type_id = 0; type_id < 2; type_id++) {

    if (type_id == 0)
      coupling_ent = luma_coupling->faces;
    else
      coupling_ent = luma_coupling->cells;

    if (coupling_ent != NULL) {

      if (coupling_ent->post_mesh_id != 0) {

        const cs_real_t *cell_temp = NULL;
        const cs_real_t *face_temp = NULL;
		const cs_real_t *cell_vel = NULL;
        const cs_real_t *face_vel = NULL;

        if (type_id == 0)
		{
			if (coupling_ent->is_temp_in)
				face_temp = coupling_ent->temp_in;
			if (coupling_ent->is_vel_in)
				face_vel = coupling_ent->vel_in;
		}
        else
        {
			if (coupling_ent->is_temp_in)
				cell_temp = coupling_ent->temp_in;
			if (coupling_ent->is_vel_in)
				cell_vel = coupling_ent->vel_in;
		}

        cs_post_write_var(coupling_ent->post_mesh_id,
                          CS_POST_WRITER_ALL_ASSOCIATED,
                          _("LUMA T"),
                          1,
                          false,
                          false,
                          CS_POST_TYPE_cs_real_t,
                          cell_temp,
                          NULL,
                          face_temp,
                          ts);
						  
		cs_post_write_var(coupling_ent->post_mesh_id,
                          CS_POST_WRITER_ALL_ASSOCIATED,
                          _("LUMA velocity"),
                          3,
                          true,
                          false,
                          CS_POST_TYPE_cs_real_t,
                          cell_vel,
                          NULL,
                          face_vel,
                          ts);

        /*if (type_id == 1)
          cs_post_write_var(coupling_ent->post_mesh_id,
                            CS_POST_WRITER_ALL_ASSOCIATED,
                            _("Solid heat flux"),
                            1,
                            false,
                            false,
                            CS_POST_TYPE_float,
                            coupling_ent->flux,
                            NULL,
                            NULL,
                            ts);*/

      }
    }
  }
}
 
 /*----------------------------------------------------------------------------
 * Initialize post-processing of a LUMA coupling
 *
 * parameters:
 *   luma_coupling <-- partially initialized LUMA coupling structure
 *   coupling_ent <-- associated coupling mesh entity
 *----------------------------------------------------------------------------*/

static void
_post_init(cs_luma_coupling_t      *luma_coupling,
           cs_luma_coupling_ent_t  *coupling_ent)
{
  int dim_shift = 0;
  int coupling_id = -1;

  const int writer_id = -1;
  const int writer_ids[] = {writer_id};

  /* Determine coupling id */

  for (coupling_id = 0;
       (   coupling_id < cs_glob_luma_n_couplings
        && cs_glob_luma_couplings[coupling_id] != luma_coupling);
       coupling_id++);

  /* Exit silently if associated writer is not available */

  if (cs_post_writer_exists(writer_id) != true)
    return;

  int t_top_id
    = cs_timer_stats_switch(cs_timer_stats_id_by_name("postprocessing_stage"));

  /* Initialize post processing flag */

  coupling_ent->post_mesh_id = cs_post_get_free_mesh_id();

  /* Allocate arrays if not already present. NOTE: The arrays should be initialised already. 
     The ones that are not initialised are not going to be used.   */

/*if (coupling_ent->n_elts > 0) {
    if (coupling_ent->solid_temp == NULL) // surface coupling 
      BFT_MALLOC(coupling_ent->solid_temp, coupling_ent->n_elts, cs_real_t);
    if (coupling_ent->elt_dim == syr_coupling->dim) { // volume coupling
      if (coupling_ent->flux == NULL)
        BFT_MALLOC(coupling_ent->flux, coupling_ent->n_elts, float);
    }
  }
  coupling_ent->tfluid_tmp = NULL;*/

  /* Associate external mesh description with post processing subsystem */

  if (luma_coupling->dim == 2)
    dim_shift = 1;

  cs_post_define_existing_mesh(coupling_ent->post_mesh_id,
                               coupling_ent->elts,
                               dim_shift,
                               false,
                               false,
                               1,
                               writer_ids);

  /* Register post processing function */

  cs_post_add_time_dep_output(_cs_luma_coupling_post_function,
                              (void *)luma_coupling);

  /* Update start and end (negative) numbers associated with
     dedicated post processing meshes */

  if (cs_glob_luma_post_mesh_ext[0] == 0)
    cs_glob_luma_post_mesh_ext[0] = coupling_ent->post_mesh_id;

  cs_glob_luma_post_mesh_ext[1] = coupling_ent->post_mesh_id;

  cs_timer_stats_switch(t_top_id);
}
 
  /*----------------------------------------------------------------------------
 * Exchange synchronization messages between Code_Saturne and LUMA.
 *
 * parameters:
 *   luma_coupling  <--  LUMA coupling structure
 *   op_name_send  <--  operation name to send, or NULL. Only the 32
 *                      first characters are sent if the nae is longer.
 *   op_name_recv  <--  operation name to receive, or NULL (size: 33)
 *----------------------------------------------------------------------------*/

static void
_exchange_sync(cs_luma_coupling_t  *luma_coupling,
               const char          *op_name_send,
               char                *op_name_recv)
{
#if defined(HAVE_MPI)

  if (cs_glob_rank_id < 1) {

    MPI_Status status;

    if (op_name_send != NULL) {

      char _op_name_send[33];
      strncpy(_op_name_send, op_name_send, 32);
      _op_name_send[32] = '\0';

    /* Exchange command messages */
      if (op_name_recv != NULL) {
        MPI_Sendrecv(_op_name_send, 32, MPI_CHAR,
                     luma_coupling->luma_root_rank, cs_luma_coupling_tag,
                     op_name_recv, 32, MPI_CHAR,
                     luma_coupling->luma_root_rank, cs_luma_coupling_tag,
                     luma_coupling->comm, &status);
      }

      else
        MPI_Send(_op_name_send, 32, MPI_CHAR,
                 luma_coupling->luma_root_rank, cs_luma_coupling_tag,
                 luma_coupling->comm);

    }
    else if (op_name_recv != NULL) {
      MPI_Recv(op_name_recv, 32, MPI_CHAR,
               luma_coupling->luma_root_rank, cs_luma_coupling_tag,
               luma_coupling->comm, &status);
    }

  }

  if (op_name_recv != NULL && cs_glob_rank_id > -1) {
	printf("CS: Boardcasting message %s .\n", op_name_recv);
    MPI_Bcast(op_name_recv, 32, MPI_CHAR, 0, cs_glob_mpi_comm);
    op_name_recv[32] = '\0';
  }

#endif
}
 
 /*----------------------------------------------------------------------------
 * Exchange location synchronization status
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_sync_after_location(cs_luma_coupling_t  *luma_coupling)
{
  int retval = 1;

  char op_name_send[32 + 1];
  char op_name_recv[32 + 1];

  /* Communication with LUMA */
  /*----------------------------*/

  /* Ready to start time iterations */

  strcpy(op_name_send, "coupling:start");

  _exchange_sync(luma_coupling, op_name_send, op_name_recv);

  if (!strcmp(op_name_recv, "coupling:error:location")) {

    cs_coupling_set_sync_flag(PLE_COUPLING_STOP);

    cs_base_warn(__FILE__, __LINE__);

    bft_printf(_(" Message received from LUMA: \"%s\"\n"
                 " indicates meshes have not been matched correctly.\n\n"
                 " The calculation will not run.\n\n"),
               op_name_recv);

  }
  else if (strcmp(op_name_recv, "coupling:start"))
    bft_error(__FILE__, __LINE__, 0,
              _(" Message received from LUMA: \"%s\"\n"
                " indicates an error or is unexpected."),
              op_name_recv);

  else
    retval = 0;

  return retval;
}

/*----------------------------------------------------------------------------
 * Check if coupling location is complete
 *
 * parameters:
 *   luma_coupling <-- partially initialized LUMA coupling structure
 *   coupling_ent <-- coupling entity
 *   n_ext        --> number of unlocated points
 *   ext_luma      --> 1 if LUMA has some unlocted elements, 0 otherwise
 *
 * returns:
 *   true if location is complete, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_is_location_complete(cs_luma_coupling_t      *luma_coupling,
                      cs_luma_coupling_ent_t  *coupling_ent,
                      cs_gnum_t               *n_ext,
                      bool                    *ext_luma)
{
  bool location_complete = true;

  /* Check that all points are effectively located */

  ple_lnum_t n_exterior = ple_locator_get_n_exterior(coupling_ent->locator);
  *n_ext = n_exterior;

  char  op_name_send[32 + 1];
  char  op_name_recv[32 + 1];

  cs_parall_counter(n_ext, 1);

  if (*n_ext > 0) {
    strcpy(op_name_send, "coupling:location:incomplete");
    location_complete = false;
  }
  else
    strcpy(op_name_send, "coupling:location:ok");

  _exchange_sync(luma_coupling, op_name_send, op_name_recv);
  if (!strcmp(op_name_recv, "coupling:location:incomplete")) {
    location_complete = false;
    *ext_luma = true;
	printf("CS: LUMA sent location incomplete message. \n");
  }
  else
    *ext_luma = false;

  return location_complete;
} 

 /*----------------------------------------------------------------------------
 * Define nodal mesh for LUMA coupling from selection criteria.
 *
 * parameters:
 *   luma_coupling    <-- partially initialized LUMA coupling structure
 *   n_locations     <-- number of associated locations
 *   location_ids    <-- associated location ids
 *   elt_dim         <-- element dimension
 *
 * returns:
 *   pointer to created LUMA coupling entity helper structure
 *----------------------------------------------------------------------------*/

static cs_luma_coupling_ent_t *
_create_coupled_ent(cs_luma_coupling_t  *luma_coupling,
                    int                  n_locations,
                    int                  location_ids[],
                    int                  elt_dim)
{
  char *coupled_mesh_name = NULL;
  bool      ext_luma = false;
  cs_lnum_t n_exterior = 0;
  cs_gnum_t n_ext = 0;
  cs_coord_t *elt_centers = NULL;
  fvm_nodal_t *location_elts = NULL;
  float *cs_to_luma_dist = NULL;
  float *luma_to_cs_dist = NULL;

  bool location_complete = false;
  cs_luma_coupling_ent_t *coupling_ent = NULL;

  int locator_options[PLE_LOCATOR_N_OPTIONS];
  locator_options[PLE_LOCATOR_NUMBERING] = 1;

  assert(luma_coupling != NULL);

  /* Initialization */

  BFT_MALLOC(coupling_ent, 1, cs_luma_coupling_ent_t);

  coupling_ent->locator = NULL;
  coupling_ent->elt_dim = elt_dim;

  coupling_ent->n_elts = 0;
  coupling_ent->elts = NULL;

  coupling_ent->post_mesh_id = 0;
  coupling_ent->temp_in = NULL;
  coupling_ent->vel_in = NULL;
  coupling_ent->temp_out = NULL;
  coupling_ent->vel_out = NULL;
  
  coupling_ent->is_temp_in = false;
  coupling_ent->is_temp_out = false;
  coupling_ent->is_vel_in = false;
  coupling_ent->is_vel_out = false;
  
  // Again, I'm not sure I'll keep hvol
  coupling_ent->hvol = NULL;

  if (luma_coupling->verbosity > 0) {
    bft_printf(_("\nExtracting coupled mesh             ..."));
    bft_printf_flush();
  }

  /* Select elements */

  cs_lnum_t  n_elts = 0;
  cs_lnum_t *elt_list = NULL;

  for (int l_i = 0; l_i < n_locations; l_i++)
    n_elts += cs_mesh_location_get_n_elts(location_ids[l_i])[0];

  BFT_MALLOC(elt_list, n_elts, cs_lnum_t);

  n_elts = 0;
  for (int l_i = 0; l_i < n_locations; l_i++) {
    int loc_id = location_ids[l_i];
    const cs_lnum_t n = cs_mesh_location_get_n_elts(loc_id)[0];
    const cs_lnum_t *ids = cs_mesh_location_get_elt_list(loc_id);
    if (ids != NULL) {
      for (cs_lnum_t i = 0; i < n; i++)
        elt_list[n_elts++] = ids[i] + 1;
    }
    else {
      for (cs_lnum_t i = 0; i < n; i++)
        elt_list[n_elts++] = i + 1;
    }
  }
  
  /* Check which variables need coupling and allocate the corresponding buffers*/
  // Only the "faces" coupled entity will have vars_in, no vars_out
  if (elt_dim == luma_coupling->dim -1)
  {
    for (size_t i = 0; i < strlen(luma_coupling->vars_in); i++)
    {
      printf("CS: vars_in = %s , vars_in[%lu] = %c \n", luma_coupling->vars_in, i, luma_coupling->vars_in[i]);
    
      if ((luma_coupling->vars_in[i] == 'v') ||
        (luma_coupling->vars_in[i] == 'V'))
      {
        printf("CS: Allocating LUMA velocity, number of elements: %d \n", n_elts);
        BFT_MALLOC(coupling_ent->vel_in, 3 * n_elts, cs_real_t);
        coupling_ent->is_vel_in = true;
      }
      else if ((luma_coupling->vars_in[i] == 't') ||
        (luma_coupling->vars_in[i] == 'T'))
      {
        printf("CS: Allocating LUMA temperature?\n");
        BFT_MALLOC(coupling_ent->temp_in, n_elts, cs_real_t);
        coupling_ent->is_temp_in = true;
      }
      else
      {
       bft_printf(_("Coupled variable %c not recognised."),luma_coupling->vars_in[i]);
        bft_printf_flush();
      }
    }
  }
  
  // Only the "cells" coupled entity will have vars out, no vars_in
  if (elt_dim == luma_coupling->dim)
  {
    printf("CS: lenght of out variables %lu \n", strlen(luma_coupling->vars_out));
    for (size_t i = 0; i < strlen(luma_coupling->vars_out); i++)
    {
      printf("CS: vars_out = %s , vars_out[%lu] = %c \n", luma_coupling->vars_out, i, luma_coupling->vars_out[i]);
    
      if ((luma_coupling->vars_out[i] == 'v') ||
        (luma_coupling->vars_out[i] == 'V'))
      {
        printf("CS: Allocating LUMA velocity, number of elements: %d \n", n_elts);
        BFT_MALLOC(coupling_ent->vel_out, 3 * n_elts, cs_real_t);
        coupling_ent->is_vel_out = true;
      }
      else if ((luma_coupling->vars_out[i] == 'T') ||
        (luma_coupling->vars_out[i] == 'T'))
      {
        printf("CS: Allocating LUMA temperature?\n");
        BFT_MALLOC(coupling_ent->temp_out, n_elts, cs_real_t);
        coupling_ent->is_temp_out = true;
      }
      else
      {
        bft_printf(_("Coupled variable %c not recognised."),luma_coupling->vars_in[i]);
        bft_printf_flush();
      }
    }
  }
  
  
  /* Creation of a new nodal mesh from selected cells */

  if (elt_dim == luma_coupling->dim) {

    BFT_MALLOC(coupled_mesh_name,
                 strlen(_("LUMA %s cells"))
               + strlen(luma_coupling->luma_name) + 1, char);
    sprintf(coupled_mesh_name, _("LUMA %s cells"), luma_coupling->luma_name);

    coupling_ent->n_elts = n_elts;

    coupling_ent->elts
      = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                       coupled_mesh_name,
                                       false,
                                       coupling_ent->n_elts,
                                       elt_list);

    /* Allocate additional buffers */

    BFT_MALLOC(coupling_ent->hvol, coupling_ent->n_elts, double);

  }

  /* Creation of a new nodal mesh from selected border faces */

  else if (elt_dim == luma_coupling->dim - 1) {

    BFT_MALLOC(coupled_mesh_name,
               strlen("LUMA  faces") + strlen(luma_coupling->luma_name) + 1,
               char);
    sprintf(coupled_mesh_name, _("LUMA %s faces"), luma_coupling->luma_name);

    coupling_ent->n_elts = n_elts;

    coupling_ent->elts
      = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                       coupled_mesh_name,
                                       false,
                                       0,
                                       coupling_ent->n_elts,
                                       NULL,
                                       elt_list);

  }

  BFT_FREE(elt_list);

  BFT_FREE(coupled_mesh_name);

  if (luma_coupling->verbosity > 0) {
    bft_printf(" [ok]\n");
    bft_printf_flush();
  }

  if (fvm_nodal_get_n_g_vertices(coupling_ent->elts) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Selected mesh locations:\n"
                " lead to an empty mesh for LUMA coupling .\n"
                " \"%s\"\n"),
              luma_coupling->luma_name);

  /* In case of 2D coupling, project coupled elements to 2D */

  location_elts = coupling_ent->elts;

  if (luma_coupling->dim == 2) {

    double  a[6];
    cs_lnum_t  n_errors = 0;

    if (luma_coupling->verbosity > 0) {
      bft_printf(_("Projecting the extracted mesh to 2D ..."));
      bft_printf_flush();
    }

    fvm_nodal_project(coupling_ent->elts, luma_coupling->ref_axis, &n_errors);

    if (n_errors > 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error projecting the extracted mesh."));

    if (luma_coupling->verbosity > 0) {
      bft_printf(_(" [ok]\n"));
      bft_printf_flush();
    }

    location_elts = fvm_nodal_copy(coupling_ent->elts);

    if (luma_coupling->ref_axis == 0) {
      a[0] = 0.; a[1] = 1.; a[2] = 0.; a[3] = 0.; a[4] = 0.; a[5] = 1.;
    }
    else if (luma_coupling->ref_axis == 1) {
      a[0] = 1.; a[1] = 0.; a[2] = 0.; a[3] = 0.; a[4] = 0.; a[5] = 1.;
    }
    else if (luma_coupling->ref_axis == 2) {
      a[0] = 1.; a[1] = 0.; a[2] = 0.; a[3] = 0.; a[4] = 1.; a[5] = 0.;
    }

    fvm_nodal_project_coords(location_elts, a);
  }

  /* Element information */

  if (luma_coupling->verbosity > 0) {
    cs_gnum_t n_g_elts = coupling_ent->n_elts;
    cs_parall_counter(&n_g_elts, 1);
    bft_printf(_("\nExtracted mesh built of %llu elements.\n"),
               (unsigned long long)n_g_elts);
    bft_printf_flush();
  }

  /* Initialize post-processing */

  /* Precaution: deactivate visualization for time-dependent meshes,
     as this would require regenerating visualization at each time step */

  if (cs_post_get_writer_time_dep(-1) != FVM_WRITER_FIXED_MESH)
    luma_coupling->visualization = 0;

  if (luma_coupling->visualization != 0)
    _post_init(luma_coupling, coupling_ent);

  /* Build and initialize associated locator */

  if (luma_coupling->verbosity > 0) {
    bft_printf(_("\nLocator structure and mesh creation ..."));
    bft_printf_flush();
  }

  /* Retrieve coordinates using FVM functions rather than previous list and
     coordinates, in case the extracted nodal mesh contains elements in a
     different order (multiple element types) or coordinates are projected
     in 2D. */

  if (coupling_ent->n_elts > 0) {

    if (luma_coupling->visualization != 0)
      BFT_MALLOC(cs_to_luma_dist, coupling_ent->n_elts, float);

    BFT_MALLOC(elt_centers,
               coupling_ent->n_elts*luma_coupling->dim,
               cs_coord_t);
    fvm_nodal_get_element_centers(location_elts,
                                  CS_INTERLACE,
                                  coupling_ent->elt_dim,
                                  elt_centers);
  }

  /* Locate entities */

#if defined(PLE_HAVE_MPI)
  coupling_ent->locator = ple_locator_create(luma_coupling->comm,
                                             luma_coupling->n_luma_ranks,
                                             luma_coupling->luma_root_rank);
#else
  coupling_ent->locator = ple_locator_create();
#endif

  printf("CS: Hi before ple_locator_set_mesh \n");

  ple_locator_set_mesh(coupling_ent->locator,
                       location_elts,
                       locator_options,
                       0.,
                       luma_coupling->tolerance,
                       luma_coupling->dim,
                       coupling_ent->n_elts,
                       NULL,
                       NULL,
                       elt_centers,
                       cs_to_luma_dist,
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh);
					   
	printf("CS: Hi after ple_locator_set_mesh \n");

  /* Check that all points are effectively located */
  
  // Print cs_to_luma_dist ////
   //printf( "CS: number of coupled elements %d \n", coupling_ent->n_elts);

   /*if(luma_coupling->visualization != 0) // cs_to_luma_dist doesn't exist if visualisation is 0.
   {
		for (int num = 0; num < 200; num++)
			printf("CS: Distance from mesh to LUMA element %d is %f \n", num, cs_to_luma_dist[num]);
   }*/
   

  location_complete = _is_location_complete(luma_coupling,
                                            coupling_ent,
                                            &n_ext,
                                            &ext_luma);
											
  
  printf("location complete? %d. Number of points not located: %d \n", location_complete, ple_locator_get_n_exterior(coupling_ent->locator));

  if (luma_coupling->allow_nearest) {

    float tolerance = luma_coupling->tolerance;

    while (location_complete == false) {

      tolerance *= 4;

      if (luma_coupling->verbosity > 0) {
        bft_printf(_(" [failed]\n"));
        if (n_ext > 0)  // NOTE: I don't understand this. Because as far as I know, the locator only has information about CS. How would it receive coordinates from Syrthes?
          bft_printf(_(" %llu fluid mesh elements not located on solid mesh\n"),
                     (unsigned long long) n_ext);
        if (ext_luma)
          bft_printf(_(" Some solid mesh elements not located on fluid mesh\n"));
        bft_printf(_("\n   Extending search with tolerance factor %f..."),
                   tolerance);
        bft_printf_flush();
      }

      ple_locator_extend_search(coupling_ent->locator,
                                location_elts,
                                locator_options,
                                0.,
                                tolerance,
                                coupling_ent->n_elts,
                                NULL,
                                NULL,
                                elt_centers,
                                cs_to_luma_dist,
                                cs_coupling_mesh_extents,
                                cs_coupling_point_in_mesh);

      location_complete = _is_location_complete(luma_coupling,
                                                coupling_ent,
                                                &n_ext,
                                                &ext_luma);

    }

  }

  if (luma_coupling->verbosity > 0) {
    bft_printf(_(" [ok]\n"));
    bft_printf_flush();
  }

  if (location_elts != coupling_ent->elts)
    fvm_nodal_destroy(location_elts);

  if (elt_centers != NULL)
    BFT_FREE(elt_centers);

  if (luma_coupling->visualization != 0) {

    cs_post_activate_writer(-1, 1);
    cs_post_write_meshes(cs_glob_time_step);

    const cs_real_t *b_dist = NULL, *v_dist = NULL;

    if (coupling_ent->elt_dim == luma_coupling->dim - 1)
      b_dist = (cs_real_t*)cs_to_luma_dist;
    else if (coupling_ent->elt_dim == luma_coupling->dim)
      v_dist = (cs_real_t*)cs_to_luma_dist;

    cs_post_write_var(coupling_ent->post_mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      _("distance_to_LUMA"),
                      1,
                      false,
                      false, // use_parent, 
                      CS_POST_TYPE_float,
                      v_dist,
                      NULL,
                      b_dist,
                      NULL);  // time-independent variable 

    BFT_FREE(cs_to_luma_dist);

  }

  /* Post-process distances from LUMA points to Code_Saturne faces */

  if (elt_dim == luma_coupling->dim - 1) {

    printf("CS: Element dimensions %u, luma dimensions %u\n", elt_dim, luma_coupling->dim);

    cs_lnum_t n_dist_elts = ple_locator_get_n_dist_points(coupling_ent->locator);

    BFT_MALLOC(luma_to_cs_dist, n_dist_elts, float);

    ple_locator_exchange_point_var(coupling_ent->locator,
                                   luma_to_cs_dist,
                                   NULL,
                                   NULL,
                                   sizeof(float),
                                   1,
                                   1);

    if (   luma_coupling->visualization != 0
        && luma_coupling->allow_nearest == false) {

      cs_lnum_t i;
      int writer_ids[] = {-1};
      int mesh_id = coupling_ent->post_mesh_id - 1;
      cs_lnum_t *p_vtx_num = NULL;
      fvm_io_num_t *vtx_io_num = NULL;
      fvm_nodal_t *luma_points = fvm_nodal_create("LUMA face centers",
                                                 luma_coupling->dim);

      BFT_MALLOC(p_vtx_num, n_dist_elts, cs_lnum_t);

      for (i = 0; i < (cs_lnum_t)n_dist_elts; i++)
        p_vtx_num[i] = i+1;

      fvm_nodal_define_vertex_list(luma_points, n_dist_elts, p_vtx_num);
      fvm_nodal_set_shared_vertices
        (luma_points,
         ple_locator_get_dist_coords(coupling_ent->locator));

      if (cs_glob_n_ranks > 1) {

        vtx_io_num = fvm_io_num_create_from_scan(n_dist_elts);

        fvm_nodal_init_io_num(luma_points,
                              fvm_io_num_get_global_num(vtx_io_num),
                              0);

      }

      cs_post_define_existing_mesh(mesh_id,
                                   luma_points,
                                   0,
                                   true,
                                   false,
                                   1,
                                   writer_ids);

      cs_post_activate_writer(-1, 1);
      cs_post_write_meshes(cs_glob_time_step);

      cs_post_write_vertex_var(mesh_id,
                               CS_POST_WRITER_ALL_ASSOCIATED,
                               _("distance_to_Code_Saturne"),
                               1,
                               false,
                               false, // use parent
                               CS_POST_TYPE_float,
                               luma_to_cs_dist,
                               NULL); // time-independent variable 

      cs_post_free_mesh(mesh_id);

      if (cs_glob_n_ranks > 1)
        fvm_io_num_destroy(vtx_io_num);

    } /* Do post-processing */

    BFT_FREE(luma_to_cs_dist);

  }

  if (n_ext) {

    int i;
    int writer_ids[] = {-1};
    int mesh_id = cs_post_get_free_mesh_id();
    cs_lnum_t *post_vtx_num = NULL;
    cs_coord_t *exterior_coords = NULL;
    cs_coord_t *el_list = NULL;
    fvm_io_num_t *vtx_io_num = NULL;
    fvm_nodal_t *ulck_points = fvm_nodal_create("unlocated elements (centers)",
                                                3);
    n_exterior = ple_locator_get_n_exterior(coupling_ent->locator);
    const ple_lnum_t *exterior_list
      = ple_locator_get_exterior_list(coupling_ent->locator);

    BFT_MALLOC(post_vtx_num, n_exterior, cs_lnum_t);
    BFT_MALLOC(exterior_coords, 3*n_exterior, cs_coord_t);
    BFT_MALLOC(el_list,
               coupling_ent->n_elts*3,
               cs_coord_t);

    fvm_nodal_get_element_centers(coupling_ent->elts,
                                  CS_INTERLACE,
                                  coupling_ent->elt_dim,
                                  el_list);

    for (i = 0; i < (cs_lnum_t)n_exterior; i++) {
      post_vtx_num[i] = i+1;
      if (exterior_list[i] > coupling_ent->n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: invalid exterior elements selection."));
      exterior_coords[3*i   ] = el_list[3*(exterior_list[i]-1)   ];
      exterior_coords[3*i +1] = el_list[3*(exterior_list[i]-1) +1];
      exterior_coords[3*i +2] = el_list[3*(exterior_list[i]-1) +2];
    }

    fvm_nodal_define_vertex_list(ulck_points,
                                 (cs_lnum_t)n_exterior,
                                 post_vtx_num);
    fvm_nodal_set_shared_vertices
      (ulck_points,
       exterior_coords);

    if (cs_glob_n_ranks > 1) {
      vtx_io_num = fvm_io_num_create_from_scan(n_exterior);
      fvm_nodal_init_io_num(ulck_points,
                            fvm_io_num_get_global_num(vtx_io_num),
                            0);
    }

    cs_post_define_existing_mesh(mesh_id,
                                 ulck_points,
                                 0,
                                 true,
                                 false,
                                 1,
                                 writer_ids);

    cs_post_activate_writer(writer_ids[0], 1);
    cs_post_write_meshes(cs_glob_time_step);
    cs_post_free_mesh(mesh_id);

    if (cs_glob_n_ranks > 1)
      fvm_io_num_destroy(vtx_io_num);

    BFT_FREE(el_list);
    BFT_FREE(exterior_coords);

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Coupling with LUMA impossible:\n"
                 "%llu element centers from mesh \"%s\"\n"
                 "not located on LUMA mesh."),
               (unsigned long long)n_ext,
               fvm_nodal_get_name(coupling_ent->elts));

  }

  /* Ensure clean stop */

  if (location_complete == false)
    cs_coupling_set_sync_flag(PLE_COUPLING_STOP);

  return coupling_ent;
}
 
 /*----------------------------------------------------------------------------
 * Print information on yet unmatched LUMA couplings.
 *
 * parameters:
 *   n_unmatched    <--  number of unmatched couplings
 *   unmatched_ids  <--  array of unmatched couplings
 *----------------------------------------------------------------------------*/

static void
_print_all_unmatched_luma(int        n_unmatched,
                          const int  unmatched_ids[])
{
  /* Loop on defined LUMA instances */

  for (int i = 0; i < n_unmatched; i++) {

    cs_luma_coupling_t *luma_coupling
      = cs_luma_coupling_by_id(unmatched_ids[i]);
    const char *local_name = cs_luma_coupling_get_name(luma_coupling);

    bft_printf(_(" LUMA coupling:\n"
                 "   coupling id:              %d\n"
                 "   local name:               \"%s\"\n\n"),
               i, local_name);

  }

  bft_printf_flush();
}
 

 /*----------------------------------------------------------------------------
 * Initialize communicator for LUMA coupling
 *
 * parameters:
 *   luma_coupling  <-> LUMA coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *----------------------------------------------------------------------------*/

static void
_init_comm(cs_luma_coupling_t *luma_coupling,
           int                 coupling_id)

{
#if defined(HAVE_MPI)

  int  mpi_flag = 0;
  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  MPI_Initialized(&mpi_flag);

  if (mpi_flag == 0)
    return;

  bft_printf(_(" LUMA coupling %d: initializing MPI communication ... "),
             coupling_id);
  bft_printf_flush();

  // NOTE: I'm not sure this will work with LUMA, since LUMA is using the ple_communicator
  ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                    cs_glob_mpi_comm,
                                    luma_coupling->luma_root_rank,
                                    &(luma_coupling->comm),
                                    local_range,
                                    distant_range);

  bft_printf(_("[ok]\n"));
  bft_printf(_("  Local ranks = [%d..%d], distant ranks = [%d..%d].\n\n"),
             local_range[0], local_range[1] - 1,
             distant_range[0], distant_range[1] - 1);
  bft_printf_flush();

  luma_coupling->n_luma_ranks = distant_range[1] - distant_range[0];
  luma_coupling->luma_root_rank = distant_range[0];

#endif
}
 
 #if defined(HAVE_MPI)
 
 /*----------------------------------------------------------------------------
 * Initialize MPI LUMA couplings using MPI.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   n_unmatched    <->  pointer to number of unmatched couplings
 *   unmatched_ids  <->  pointer to array of unmatched couplings
 *----------------------------------------------------------------------------*/

static void
_init_all_mpi_luma(int  *n_unmatched,
                   int  **unmatched_ids)
{
  int _n_unmatched = *n_unmatched;
  int *_unmatched_ids = *unmatched_ids;

  const int n_couplings = cs_luma_coupling_n_couplings();

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return;

  const int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);
  
  printf("CS: Initialising MPI LUMA coupling. Number of MPI apps = %d \n", n_apps);

  /* Loop on applications */

  for (int i = 0; i < n_apps; i++) {

    ple_coupling_mpi_set_info_t ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
	
	printf("CS: App number %d, app_type %s \n", i, ai.app_type);

    if (strncmp(ai.app_type, "LUMA", 4) == 0) {

      int  match_queue_id = -1;
      int  coupling_id = -1;

      if (n_apps == 2 && n_couplings == 1 && _n_unmatched == 1) {
        match_queue_id = 0;
        coupling_id = 0;
      }
      else if (ai.app_name != NULL) {
        for (int j = 0; j < _n_unmatched; j++) {
          int k = _unmatched_ids[j];
          cs_luma_coupling_t *lcpl = cs_luma_coupling_by_id(k);
          if (strcmp(ai.app_name, cs_luma_coupling_get_name(lcpl)) == 0) {
            coupling_id = k;
            match_queue_id = j;
            break;
          }
        }
      }

      if (coupling_id > -1) {

        /* Remove from unmatched queue */
        _n_unmatched -= 1;
        for (int l = match_queue_id; l < _n_unmatched; l++)
          _unmatched_ids[l] = _unmatched_ids[l+1];
        if (_n_unmatched == 0)
          BFT_FREE(_unmatched_ids);

		printf("CS: Coupling ID %d \n", coupling_id);
	  
        /* Set communicator */
        cs_luma_coupling_init_comm(cs_luma_coupling_by_id(coupling_id),
                                   coupling_id,
                                   ai.root_rank,
                                   ai.n_ranks);

        /* Print matching info */

        const char *luma_version = cs_empty_string;
        const char *local_name = cs_empty_string;
        const char *distant_name = cs_empty_string;

        if (ai.app_name != NULL)
          local_name = ai.app_name;
        if (ai.app_type != NULL)
          luma_version = ai.app_type;
        if (ai.app_name != NULL)
          distant_name = ai.app_name;

        bft_printf(_(" LUMA coupling:\n"
                     "   coupling id:              %d\n"
                     "   version:                  \"%s\"\n"
                     "   local name:               \"%s\"\n"
                     "   distant application name: \"%s\"\n"
                     "   MPI application id:       %d\n"
                     "   MPI root rank:            %d\n"
                     "   number of MPI ranks:      %d\n\n"),
                   coupling_id, luma_version, local_name, distant_name,
                   i, ai.root_rank, ai.n_ranks);
      }

      /* Note that a LUMA app may be present in the coupling set, but
         not coupled to the current code_saturne instance, so
         coupling_id < 0 here should not be reported as an error or
         complained about here. In case of missing matches, only the
         codes having defined and missing couplings should complain. */

    }

  } /* End of loop on applications */

  bft_printf_flush();

  /* Set return values */

  *n_unmatched = _n_unmatched;
  *unmatched_ids = _unmatched_ids;
}
 
/*----------------------------------------------------------------------------
 * Find name of single LUMA coupling using MPI.
 *
 * If no coupling or multiple couplings are present, the default cannot be
 * determined, so NULL is returned.
 *----------------------------------------------------------------------------*/

static const char *
_mpi_luma_default_name(void)
{
  const char *retval = NULL;

  int n_luma_apps = 0;

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps == NULL)
    return NULL;

  int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);

  /* First pass to count available LUMA couplings */

  for (int i = 0; i < n_apps; i++) {
    const ple_coupling_mpi_set_info_t
      ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
	// DEBUG!!!
	printf("CS: Coupled app type %s \n", ai.app_type);
	
    if (strncmp(ai.app_type, "LUMA", 4) == 0) {
      if (n_luma_apps == 0)
        retval = ai.app_name;
      else
        retval = NULL;
      n_luma_apps += 1;
    }
  }

  return retval;
}

#endif /* defined(HAVE_MPI) */

/*============================================================================
 *  Public functions definitions
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
                             int                  mode)
						  // int  i)	 
{
	//printf("CS: Hi I'm in the sending data function. \n");

  cs_luma_coupling_ent_t  *coupling_ent = NULL;

  assert(mode == 0 || mode == 1);

  if (mode == 0)
    coupling_ent = luma_coupling->faces; //faces[i]
  else
    coupling_ent = luma_coupling->cells; //cells[i]

  if (coupling_ent == NULL)
    return;


  // Get the number of coupled points and the distances 
  cs_lnum_t n_dist = ple_locator_get_n_dist_points(coupling_ent->locator);
  const cs_lnum_t* dist_loc = ple_locator_get_dist_locations(coupling_ent->locator);
  
  /* "Interpolate" (nearest neighbour) and send data */

  if(coupling_ent->is_temp_out)
  {
	double* send_var = NULL;
	BFT_MALLOC(send_var, n_dist, double);
	
	for (int i = 0; i< n_dist; i++)
		send_var[i] = t_cs[dist_loc[i]];  // -1 has to do with Fortran and C indexing but I don't understand why 
	                                          // this is not done when the data is received from LUMA. 
		
	ple_locator_exchange_point_var(coupling_ent->locator,
                                 send_var,
                                 NULL,
                                 NULL,
                                 sizeof(double),
                                 1,
                                 0);
								 
	if (coupling_ent->n_elts > 0) 
	{
		cs_lnum_t i;
		for (i = 0; i < coupling_ent->n_elts; i++)
			coupling_ent->temp_out[i] = t_cs[i];
	}
	
	BFT_FREE(send_var);
  }
  
  if(coupling_ent->is_vel_out)
  {	
	double* send_var = NULL;
	BFT_MALLOC(send_var, 3 * n_dist, double);

  double flowrate1 = 0.;
  for (int i = 0; i < n_dist; i++)
  {
		for (int ic = 0; ic < 3; ic++)
      send_var[3 * i + ic] = v_cs[3 * i + ic];//v_cs[3 * dist_loc[i] + ic];  // -1 has to do with Fortran and C indexing but I don't understand why 
	                                          // this is not done when the data is received from LUMA.
/*     printf("---------I am here check the velocity ux \n---------");*/
//    printf("%e ", send_var[3 * i + 0]); 

    //flowrate1 += send_var[3 * i] * 0.01*0.01;
    
/*      printf("---------I am checking what are dist_loc \n---------");
     printf("%d \n", dist_loc[i]);  */

  }
  //printf("++++++CS: send the flowrate1 to LUMA: %e++++++ \n", flowrate1);

  //printf("CS: Data ready to send to LUMA. \n");
		
	ple_locator_exchange_point_var(coupling_ent->locator,
                                 send_var,
                                 NULL,
                                 NULL,
                                 sizeof(double),
                                 3,
                                 0);
								 
	//printf("CS: Data sent to LUMA. \n");
								 
	if (coupling_ent->n_elts > 0) 
	{
		cs_lnum_t i;
    //double flowrate2 = 0.;
    for (i = 0; i < coupling_ent->n_elts; i++)
    {
			for (int j = 0; j < 3; j++)
				coupling_ent->vel_out[3*i+j] = v_cs[3*i+j];
      //flowrate2 += v_cs[3 * i]*0.01*0.01;
    }
   //printf("++++++CS: send the flowrate2 to LUMA: %f++++++ \n", flowrate2);
  }

  BFT_FREE(send_var);
  //printf("CS: just freed send_var \n");
  }
  
}
 
 
 /*----------------------------------------------------------------------------*/
/*!
 * \brief  Send cell values to a LUMA coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_luma_coupling_send_volume(void)
{
  //printf("CS: I'm in the send volume to LUMA function \n");

  const int kcpluma = cs_field_key_id("luma_coupling");

  /* Get number of coupling cases */

  int n_cpl = cs_luma_coupling_n_couplings();
  // This is Marta's work send, the data back to LUMA
  /* Loop on couplings: get boundary temperature and/or velocity array for each coupling
    and apply matching boundary condition. */

   for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++)
    {

      cs_luma_coupling_t *luma_coupling = cs_luma_coupling_by_id(cpl_id);

      //printf("CS: Number of couplings with LUMA: %d . \n", n_cpl);

      if (!cs_luma_coupling_is_vol(luma_coupling)) // ignore if surface-only
        continue;

      //printf("CS: The coupling is volume, right? \n");

      // TODO: This will have to be adapted to having multiple coupled meshes. This will have to be a double loop.
      // So from here on it will loop through the different coupled meshes. luma_coupling->cells[m].

      cs_lnum_t n_cpl_cells = cs_luma_coupling_get_n_elts(luma_coupling, 1);

      // Get list of coupled cells

      cs_lnum_t *c_ids = NULL;
      BFT_MALLOC(c_ids, n_cpl_cells, cs_lnum_t);
      cs_luma_coupling_get_elt_ids(luma_coupling, c_ids, 1);

      //printf("CS: cellIDS 30: %u \n", c_ids[30]);

      // Obtain velocity / temperature to luma 

      cs_real_t *t_cs = NULL;
      cs_real_t *v_cs = NULL;

      // Only allocate the data that will be sent to LUMA
      if (luma_coupling->cells->is_temp_out)
        BFT_MALLOC(t_cs, n_cpl_cells, cs_real_t);
      if (luma_coupling->cells->is_vel_out)
      {
        BFT_MALLOC(v_cs, 3 * n_cpl_cells, cs_real_t);
        //printf("CS: v_cs allocated. \n");
      }

      // Scalar field coupling (i e. temperature) //
      // If we need to send temperature from LUMA
      if (luma_coupling->cells->is_temp_out)
      {
        //printf("CS: Hi! I need to send temperature data to LUMA! \n");

        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0)
        {
          bft_error(__FILE__, __LINE__, 0,
                    _("Coupling with LUMA is not possible for compressible flow."));
        }

        // Get CS temperature field
        // NOTE: Fragment extracted from Dirichlet for velocity
        // in cs_gui_boundary_conditions.c, lines 1953 to 1968
        const cs_field_t *ft = cs_field_by_name_try("temperature");
        //const int var_key_id = cs_field_key_id("variable_id");
        //int ivarv = cs_field_get_key_int(fv, var_key_id) - 1;

        const cs_real_t *cvar_t = (const cs_real_t *)ft->val;

        // NOTE: Does it have to first create an array with only the cells in the coupled region
        // because the point location functions in CS give the indices related to the coupled cells only?
        // Not indices related to the whole mesh?
        for (cs_lnum_t i = 0; i < n_cpl_cells; i++)
          t_cs[i] = cvar_t[c_ids[i]];

      } // End temperature coupling

      // If we need to send velocity to LUMA
      if (luma_coupling->cells->is_vel_out)
      {
        //printf("CS: Hi! I need to send velocity data to LUMA! \n");

        // Get CS velocity field
        const cs_field_t *fv = cs_field_by_name_try("velocity");
        const cs_real_t *cvar_v = (const cs_real_t *)fv->val;

        // NOTE: Does it have to first create an array with only the cells in the coupled region
        // because the point location functions in CS give the indices related to the coupled cells only?
        // Not indices related to the whole mesh?
        for (cs_lnum_t i = 0; i < n_cpl_cells; i++)
        {
          cs_lnum_t cid = c_ids[i];

          for (cs_lnum_t ic = 0; ic < 3; ic++)
            v_cs[3 * i + ic] = cvar_v[3 * cid + ic]; // NOTE: I'm really not sure this is how the velocity data is stored...
        
          // I am testing the data extract from CS is correct
          //printf("%e ", v_cs[3 * i]);
        }

      } // End velocity coupling

      //printf("CS: Hello before sending data to LUMA \n");

      cs_luma_coupling_send_data(luma_coupling, t_cs, v_cs, 1);

      //printf("CS: Hello after sending data to LUMA. \n");

      //NOTE: I don't understand why freeing c_ids gives a free(): invalid pointer error.
      BFT_FREE(c_ids);
      //printf("CS: just freed c_ids\n");
      if (luma_coupling->cells->is_temp_out)
      {
        BFT_FREE(t_cs);
        //printf("CS: just freed t_cs\n");
      }

      if (luma_coupling->cells->is_vel_out)
      {
        BFT_FREE(v_cs);
        //printf("CS: just freed v_cs\n");
      }

    }  // End loop on couplings
	
}
 
 /*----------------------------------------------------------------------------
 * Receive coupling variables from LUMA
 *
 * parameters:
 *   luma_coupling     <-- LUMA coupling
 *   t_luma       --> LUMA temperature, NULL if no coupled temperature
 *   v_luma       --> LUMA velocity, NULL if no coupled velocity
 *   mode         <-- 0: surface coupling; 1: volume coupling  // I probably don't need this anymore
 *----------------------------------------------------------------------------*/
//TODO: Add a new argument with the id of the coupled entity in luma_coupling we want to access to. 
 
void
cs_luma_coupling_recv_data(cs_luma_coupling_t      *luma_coupling,
							 cs_real_t            t_luma[],
							 cs_real_t            v_luma[],
                             int                  mode) 
						  // int  i)	 
{
	//printf("CS: Hi I'm in the receiving data function. \n");

  cs_luma_coupling_ent_t  *coupling_ent = NULL;

  assert(mode == 0 || mode == 1);

  if (mode == 0)
    coupling_ent = luma_coupling->faces; //faces[i]
  else
    coupling_ent = luma_coupling->cells; //cells[i]

  if (coupling_ent == NULL)
    return;

  //printf("CS: Hey! The coupling entity is not null. Is temp in? %d \n", coupling_ent->is_temp_in);

  /* Receive data */

  if(coupling_ent->is_temp_in)
  {
	ple_locator_exchange_point_var(coupling_ent->locator,
                                 NULL,
                                 t_luma,
                                 NULL,
                                 sizeof(cs_real_t),
                                 1,
                                 0);
								 
	if (coupling_ent->n_elts > 0) 
	{
		cs_lnum_t i;
		for (i = 0; i < coupling_ent->n_elts; i++)
			coupling_ent->temp_in[i] = t_luma[i];
	}
  }
  
  if(coupling_ent->is_vel_in)
  {	
	v_luma[0] = -1.0;
	//printf("CS: Receiving velocity data... %f \n", v_luma[0]);
	  
	ple_locator_exchange_point_var(coupling_ent->locator,
                                 NULL,
                                 v_luma,
                                 NULL,
                                 sizeof(cs_real_t),
                                 3,
                                 0);
								 
	//printf("CS: Velocity data received... %f, %f, %f \n", v_luma[0], v_luma[1], v_luma[2]);
	
	if (coupling_ent->n_elts > 0) 
	{
		cs_lnum_t i;
		for (i = 0; i < coupling_ent->n_elts; i++)
		{
			for (int j = 0; j < 3; j++)
				coupling_ent->vel_in[3*i+j] = v_luma[3*i+j];
		}
			
	}
  }
  
}
 
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
                             int                        mode)
{
  cs_luma_coupling_ent_t  *coupling_ent = NULL;

  /* Sanity checks */

  assert(luma_coupling != NULL);
  assert(mode == 0 || mode == 1);

  if (mode == 0)
    coupling_ent = luma_coupling->faces;
  else
    coupling_ent = luma_coupling->cells;

  if (coupling_ent != NULL)
    fvm_nodal_get_parent_id(coupling_ent->elts,
                            coupling_ent->elt_dim,
                            cpl_elt_ids);
}

 
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
                            int                       mode)
{
  cs_luma_coupling_ent_t  *coupling_ent = NULL;
  cs_lnum_t retval = 0;

  /* Sanity checks */

  assert(luma_coupling != NULL);
  assert(mode == 0 || mode == 1);

  if (mode == 0)
    coupling_ent = luma_coupling->faces;
  else
    coupling_ent = luma_coupling->cells;

  if (coupling_ent != NULL)
    retval = coupling_ent->n_elts;

  return retval;
}

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
cs_luma_coupling_is_vol(const cs_luma_coupling_t  *luma_coupling)
{
  int retval = 0;

  assert(luma_coupling != NULL);

  if (luma_coupling->n_v_locations > 0)
    retval = 1;

  return retval;
}
 
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
cs_luma_coupling_is_surf(const cs_luma_coupling_t  *luma_coupling)
{
  int retval = 0;

  assert(luma_coupling != NULL);

  if (luma_coupling->n_b_locations > 0)
    retval = 1;

  return retval;
}
 
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
                              cs_real_t  rcodcl[])
{
	/* LUMA coupling: get variables from LUMA: temperature and/or velocity. 
	 The values received from LUMA will be incorporated into the coupled CS boundary
	 ====================================== */
	 
	//printf("CS: Hi! CS is calling the LUMA coupling function!\n");

	const int kcpluma  = cs_field_key_id("luma_coupling");

	const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
	const cs_lnum_t n_vars = nvar;  /* cast to cs_lnum_t because this
									 is used in address computations here */

	/* Get number of coupling cases */

	int n_cpl = cs_luma_coupling_n_couplings();

	/* Loop on couplings: get boundary temperature and/or velocity array for each coupling
	   and apply matching boundary condition. */

	for (int cpl_id = 0; cpl_id < n_cpl; cpl_id++) {

		cs_luma_coupling_t *luma_coupling = cs_luma_coupling_by_id(cpl_id);
		
		//printf("CS: Number of couplings with LUMA(reveice layer): %d . \n", n_cpl);

		if (! cs_luma_coupling_is_surf(luma_coupling))  /* ignore if volume-only */
			continue;
			
		//printf("CS: The coupling is surface, right? \n");

		// TODO: This will have to be adapted to having multiple coupled meshes. This will have to be a double loop. 
		// So from here on it will loop through the different coupled meshes. luma_coupling->faces[m]. 

		cs_lnum_t n_cpl_faces = cs_luma_coupling_get_n_elts(luma_coupling, 0);

		/* Get list of coupled faces */

		cs_lnum_t  *f_ids;
		BFT_MALLOC(f_ids, n_cpl_faces, cs_lnum_t);
		cs_luma_coupling_get_elt_ids(luma_coupling, f_ids, 0);

		/* Read boundary velocity and / or temperature and interpolate if necessary */

		cs_real_t *t_luma = NULL;
		cs_real_t *v_luma = NULL;

		// Only allocate the data that will be received from LUMA
		if(luma_coupling->faces->is_temp_in)
			BFT_MALLOC(t_luma, n_cpl_faces, cs_real_t);
		if(luma_coupling->faces->is_vel_in)
			BFT_MALLOC(v_luma, 3 * n_cpl_faces, cs_real_t);

		//printf("CS: Hello before receiving data from LUMA \n");
	
		cs_luma_coupling_recv_data(luma_coupling, t_luma, v_luma, 0);
		
		//printf("CS: Hello after receiving data from LUMA \n");

		// Scalar field coupling (i e. temperature) //
		// If we received temperature from LUMA
		if(luma_coupling->faces->is_temp_in)
		{
			//printf("CS: Hi! I've received temperature data from LUMA! \n");
            
      // Dirichlet for temperature
      // Check in cs_gui_boundary_conditions.s file the setting

      // Get CS temperature field
      const cs_field_t *f = cs_field_by_name_try("temperature");
      const int var_key_id = cs_field_key_id("variable_id");
      //! YANG: Why this place minus 1
      int ivarv = cs_field_get_key_int(f, var_key_id) - 1;
      
      // Access the boundary elements that we need to modify
      double temp_recv = 0;
      for (cs_lnum_t i = 0; i < n_cpl_faces; i++)
      {
        cs_lnum_t face_id = f_ids[i];
        
        // Set the boundary condition type of the temperature to Dirichlet
        icodcl[ivarv * n_b_faces + face_id] = 1;

        // Set the value to the CS boundary
        rcodcl[ivarv * n_b_faces +face_id] = t_luma[i];

        temp_recv = t_luma[i];
      }
    } // End temperature coupling

			/*  For scalars coupled with LUMA, prescribe a Dirichlet
			  condition at coupled faces.
			  For the time being, pass here only once, as only one scalar is
			  coupled with LUMA.
			  For the compressible module, solve in energy, but save the
			  temperature separately, for BC's to be clearer (?).
			  TODO: Dirichlet boundary condition needs to be implemented for the velocity
			  vector too. How do I do that?*/

			// Loop through all the variables to see which ones are coupled with LUMA
      
/* 			int n_fields = cs_field_n_fields();
			for (int field_id = 0 ; field_id <  n_fields; field_id++) 
			{
				cs_field_t  *f = cs_field_by_id(field_id);

				int icpluma = 0;
				if (f->type & CS_FIELD_VARIABLE)
				  icpluma = cs_field_get_key_int(f, kcpluma);

				if (icpluma < 1)
				  continue;

				const int k_var_id = cs_field_key_id("variable_id");
				int var_id = cs_field_get_key_int(f, k_var_id) - 1;

				if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) 
				{
				 bft_error
					(__FILE__, __LINE__, 0,
					 _("Coupling with LUMA is not possible for compressible flow."));
				}

				//NOTE: I think that this gives us the three components of var_id. But I'm not sure, maybe I should message Yvan about it. 
				// Yes, because I'll need to know how to do velocity. And how to do an inlet. 
				int  *_icodcl = icodcl + (var_id*n_b_faces);
				cs_real_t  *_rcodcl1 = rcodcl + (var_id*n_b_faces);
				cs_real_t  *_rcodcl2 = rcodcl + (n_b_faces*n_vars + var_id*n_b_faces);
				cs_real_t  *_rcodcl3 = rcodcl + (2*n_b_faces*n_vars + var_id*n_b_faces);

				for (cs_lnum_t i = 0; i < n_cpl_faces; i++) 
				{
					cs_lnum_t face_id = f_ids[i];

					if (   _icodcl[face_id] != CS_INDEF
						&& _icodcl[face_id] != CS_INLET) 
					{
					  if (bc_type[face_id] == CS_INLET)
						_icodcl[face_id] = CS_INLET;
					}

					_rcodcl1[face_id] = t_luma[i];
					_rcodcl2[face_id] = cs_math_infinite_r;
					_rcodcl3[face_id] = 0.;

				}
				  
				// Require temperature -> enthalpy conversion

				if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY)
				{
					if (f == cs_thermal_model_field()) 
					{
					  for (cs_lnum_t i = 0; i < n_cpl_faces; i++) {
						cs_lnum_t face_id = f_ids[i];
						_icodcl[face_id] *= -1;
					  }
					}

				} //End case for enthalpy 
				
			} // End loop on fields

		} // End temperature coupling */
	  
		// If we received velocity from LUMA
		if(luma_coupling->faces->is_vel_in)
		{
			//printf("CS: Hi! I've received velocity data from LUMA! \n"); 
			
			// NOTE: Fragment extracted from /*Dirichlet for velocity*/ 
			// in cs_gui_boundary_conditions.c, lines 1953 to 1968
			
			// Get CS velocity field
			const cs_field_t *fv = cs_field_by_name_try("velocity");
			const int var_key_id = cs_field_key_id("variable_id");
			int ivarv = cs_field_get_key_int(fv, var_key_id) - 1;
			
			// Access the boundary elements that we need to modify
			double flowRate = 0.;
			for (cs_lnum_t i = 0; i < n_cpl_faces; i++) 
			{
				//printf("CS: %d face %d\n", i, f_ids[i]);
				cs_lnum_t face_id = f_ids[i];

				for (cs_lnum_t ic = 0; ic < 3; ic++)
				{
					/*int  *_icodcl = icodcl + ((ivarv + ic)*n_b_faces);
					
					if (   _icodcl[face_id] != CS_INDEF
					&& _icodcl[face_id] != CS_INLET) 
					{
						if (bc_type[face_id] == CS_INLET)
							_icodcl[face_id] = CS_INLET;
					}*/
					
					// Set the boundary condition type for each component of the velocity to Dirichlet
					icodcl[(ivarv + ic) * n_b_faces + face_id] = 1;
					
					// Set the value of each velocity component
					rcodcl[(ivarv + ic) * n_b_faces + face_id] = v_luma[3 * i + ic];
        }
        //printf("%d ", (ivarv + 0) * n_b_faces + face_id);
        flowRate += v_luma[3 * i] * 0.005 * 0.005;
			}
			//double totalFlowRate = 0;
			//MPI_Allreduce(flowRate, totalFlowRate,1,MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
			printf("CS: Received flow rate from LUMA %f \n", flowRate);
			
			
	  
		} // End velocity coupling

		BFT_FREE(f_ids);
		if (luma_coupling->faces->is_temp_in)
			BFT_FREE(t_luma);
		if (luma_coupling->faces->is_vel_in)
			BFT_FREE(v_luma);
		
	} /* End loop on couplings */
  
}
 
 /*----------------------------------------------------------------------------
 * Define coupled mesh and send it to LUMA
 *
 * Optional post-processing output is also built at this stage.
 *
 * parameters:
 *   luma_coupling <-- LUMA coupling structure
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_init_mesh(cs_luma_coupling_t  *luma_coupling)
{
  const cs_lnum_t verbosity = luma_coupling->verbosity;

  if (verbosity > 0)
    bft_printf(_("\n ** Processing the mesh for LUMA coupling "
                 "\"%s\"\n\n"),
                 luma_coupling->luma_name);

  /* Define coupled mesh */

 // luma_couplings gets the dimensions from the projection_axis input in cs_luma_coupling_define
 assert(luma_coupling->dim == 3 || luma_coupling->dim == 2);

  int match_flag = 0;

  // DEBUG!: Check if the input to the coupling and the CS configuration works. 
  // "inlet" is a boundary, then n_b_location = 1
  for (int i = 0; i < luma_coupling->n_b_locations; i++)
	printf("CS: Number of boundary locations %d and location id %d \n", luma_coupling->n_b_locations, luma_coupling->b_location_ids[i]);
  for (int i = 0; i < luma_coupling->n_v_locations; i++)
	printf("CS: Number of volume locations %d and location id %d \n", luma_coupling->n_v_locations, luma_coupling->v_location_ids[i]);
  
 
  if (luma_coupling->n_b_locations > 0) {
	// NOTE: does _create_coupling_ent creates one element of faces for each n_b_location? No, it only creates one cs_luma_coupling_ent_t
    luma_coupling->faces = _create_coupled_ent(luma_coupling,
                                              luma_coupling->n_b_locations,
                                              luma_coupling->b_location_ids,
                                              luma_coupling->dim - 1);
    match_flag += _sync_after_location(luma_coupling);
  }

  if (luma_coupling->n_v_locations > 0) {
    luma_coupling->cells = _create_coupled_ent(luma_coupling,
                                              luma_coupling->n_v_locations,
                                              luma_coupling->v_location_ids,
                                              luma_coupling->dim);
	printf("CS: Finished creating the volume coupled entity \n");
    match_flag += _sync_after_location(luma_coupling);
	printf("CS: Syncronisation with LUMA complete \n");
  }

  /* Communication with LUMA */
  /*----------------------------*/

  if (match_flag == 0 && verbosity > 0) {
    bft_printf(_("\n ** Mesh located for LUMA coupling \"%s\".\n\n"),
               luma_coupling->luma_name);
    bft_printf_flush();
  }
}
 
 
 /*----------------------------------------------------------------------------*/
/*
 * Create coupled meshes and setup PLE locator for LUMA couplings.
 */
/*----------------------------------------------------------------------------*/

void
cs_luma_coupling_init_meshes(void)
{
  int n_coupl = cs_luma_coupling_n_couplings();

  for (int coupl_id = 0; coupl_id < n_coupl; coupl_id++) {
    cs_luma_coupling_t *luma_coupling = cs_luma_coupling_by_id(coupl_id);
    cs_luma_coupling_init_mesh(luma_coupling);
  }
}
 
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
cs_luma_coupling_get_name(cs_luma_coupling_t  *luma_coupling)
{
  const char *retval = cs_empty_string;

  if (luma_coupling->luma_name != NULL)
    retval = luma_coupling->luma_name;

  return retval;
}
 
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
                           int                 n_luma_ranks)
{
#if defined(HAVE_MPI)

  char  volume_flag = ' ';
  char  boundary_flag = ' ';
  char  conservativity_flag = '1';
  char  allow_nearest_flag = '1';
  char  op_name_send[32 + 1];
  char  op_name_recv[32 + 1];

  luma_coupling->n_luma_ranks = n_luma_ranks;
  luma_coupling->luma_root_rank = luma_root_rank;

  _init_comm(luma_coupling, coupling_id);

  /* Exchange coupling options */

  if (luma_coupling->n_b_locations > 0)
    boundary_flag = 'b';
  if (luma_coupling->n_v_locations > 0)
    volume_flag = 'v';
 // if (cs_luma_coupling_conservativity == 0)
  //  conservativity_flag = '0';
  if (luma_coupling->allow_nearest == false)
    allow_nearest_flag = '0';

   // I'm not checking that the coupling options for LUMA from CS are the same as the ones from LUMA
   // I can just get LUMA to get the info from CS, right? 
  //snprintf(op_name_send, 32, "coupling:type:%c%c%c \2\2%c(%6.2g)",
  //         boundary_flag, volume_flag, conservativity_flag,
  //         allow_nearest_flag, (double)luma_coupling->tolerance);

  /*_exchange_sync(luma_coupling, op_name_send, op_name_recv);

  if (strncmp(op_name_recv, op_name_send, 16))
    bft_error
      (__FILE__, __LINE__, 0,
       _("========================================================\n"
         "   ** Incompatible LUMA coupling options:\n"
         "      ------------------------------------------------\n"
         "      Code_Saturne options: \"%s\"\n"
         "      LUMA options:      \"%s\"\n"
         "========================================================\n"),
       op_name_send, op_name_recv);*/

#endif
}

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
cs_luma_coupling_by_id(int  coupling_id)
{
  cs_luma_coupling_t  *retval = NULL;

  if (   coupling_id > -1
      && coupling_id < cs_glob_luma_n_couplings)
    retval = cs_glob_luma_couplings[coupling_id];

  return retval;
}
 
 /*----------------------------------------------------------------------------
 * Get number of LUMA couplings.
 *
 * returns:
 *   number of LUMA couplings
 *----------------------------------------------------------------------------*/

int
cs_luma_coupling_n_couplings(void)
{
  return cs_glob_luma_n_couplings;
}
 
 /*----------------------------------------------------------------------------
 * Initialize LUMA couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_all_init(void)
{
  int n_unmatched = cs_luma_coupling_n_couplings();

  int *unmatched_ids;
  BFT_MALLOC(unmatched_ids, n_unmatched, int);

  for (int i = 0; i < n_unmatched; i++)
    unmatched_ids[i] = i;

  /* First try using MPI */

#if defined(HAVE_MPI)

  if (n_unmatched > 0)
    _init_all_mpi_luma(&n_unmatched, &unmatched_ids);

#endif

  if (n_unmatched > 0) {

    bft_printf("Unmatched LUMA couplings:\n"
               "----------------------------\n\n");

    _print_all_unmatched_luma(n_unmatched, unmatched_ids);

    BFT_FREE(unmatched_ids);

    bft_error(__FILE__, __LINE__, 0,
              _("At least 1 LUMA coupling was defined for which\n"
                "no communication with a LUMA instance is possible."));
  }
}
 
 
 /*----------------------------------------------------------------------------
 * Add a mesh location to a luma_coupling_t structure.
 *
 * parameters:
 *   luma_coupling  <-- LUMA coupling structure
 *   location_id   <-- id of mesh location to add (boundary faces or cells)
 *----------------------------------------------------------------------------*/

void
cs_luma_coupling_add_location(cs_luma_coupling_t  *luma_coupling,
                              int                  location_id)
{
  cs_mesh_location_type_t l_type = cs_mesh_location_get_type(location_id);

  if (l_type == CS_MESH_LOCATION_BOUNDARY_FACES) {
	
	printf("CS: Boundary mesh located\n");
	
    int i = luma_coupling->n_b_locations;
    luma_coupling->n_b_locations += 1;
    BFT_REALLOC(luma_coupling->b_location_ids, luma_coupling->n_b_locations, int);

    luma_coupling->b_location_ids[i] = location_id;
  }

  else if (l_type == CS_MESH_LOCATION_CELLS) {
	  
	printf("CS: Cell mesh located\n");  
	  
    int i = luma_coupling->n_v_locations;
    luma_coupling->n_v_locations += 1;
    BFT_REALLOC(luma_coupling->v_location_ids, luma_coupling->n_v_locations, int);

    luma_coupling->v_location_ids[i] = location_id;
  }
}

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
					   const char  *vars_out)
{
  int dim = 3;
  int ref_axis = -1;

  switch (projection_axis) {
  case 'x':
  case 'X':
    dim = 2;
    ref_axis = 0;
    break;
  case 'y':
  case 'Y':
    dim = 2;
    ref_axis = 1;
    break;
  case 'z':
  case 'Z':
    dim = 2;
    ref_axis = 2;
    break;
  default:
    break;
  }

  /* Ensure name is available */

#if defined(HAVE_MPI)
  if (luma_name == NULL)
    luma_name = _mpi_luma_default_name();
#endif

  if (luma_name == NULL)
    luma_name = cs_empty_string;

  /* Define additional coupling */

  cs_luma_coupling_t  *luma_coupling = cs_luma_coupling_create(dim,
                                                              ref_axis,
                                                              luma_name,
                                                              allow_nonmatching,
                                                              tolerance,
                                                              verbosity,
                                                              visualization,
															  vars_in,
															  vars_out);

  /* Add locations if done at that stage (deprecated) */

  int n_locations = cs_mesh_location_n_locations();

  const char *sel_criteria[2] = {boundary_criteria, volume_criteria};
  const char *type_name[2] = {"faces", "cells"};
  cs_mesh_location_type_t type_filter[2] = {CS_MESH_LOCATION_BOUNDARY_FACES,
                                            CS_MESH_LOCATION_CELLS};

  for (int i = 0; i < 2; i++) {

    if (sel_criteria[i] != NULL) {
      for (int j = 0; j < n_locations && sel_criteria[i] != NULL; j++) {
        cs_mesh_location_type_t l_type = cs_mesh_location_get_type(j);
        if (l_type & type_filter[i]) {
          const char *c = cs_mesh_location_get_selection_string(j);
          if (c != NULL) {
            if (strcmp(c, sel_criteria[i]) == 0) {
              cs_luma_coupling_add_location(luma_coupling, j);
              sel_criteria[i] = NULL;
            }
          }
        }
      }
    }

    if (sel_criteria[i] != NULL) {

      char *name;
      size_t l = strlen(luma_name) + strlen(type_name[i]) + 2;
      BFT_MALLOC(name, l, char);
      snprintf(name, l, "%s_%s", luma_name, type_name[i]);

      int j = cs_mesh_location_add(name, type_filter[i], sel_criteria[i]);

      BFT_FREE(name);

      cs_luma_coupling_add_location(luma_coupling, j);

    }

  }

}
 
 
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
						const char  *vars_out)
{
  cs_luma_coupling_t *luma_coupling = NULL;

  /* Search in existing couplings */
  
  for (int i = 0;
       i < cs_glob_luma_n_couplings;
       i++) {

    if (strcmp(cs_glob_luma_couplings[i]->luma_name, luma_name) == 0) {
      luma_coupling = cs_glob_luma_couplings[i];

      BFT_FREE(luma_coupling->luma_name);
      BFT_FREE(luma_coupling->b_location_ids);
      BFT_FREE(luma_coupling->v_location_ids);
	  BFT_FREE(luma_coupling->vars_in);
	  BFT_FREE(luma_coupling->vars_out);

      assert(luma_coupling->faces == NULL);  /* not built yet at this stage */
      assert(luma_coupling->cells == NULL);
    }
  }

  /* Allocate _cs_luma_coupling_t structure */

  if (luma_coupling == NULL) {
    BFT_REALLOC(cs_glob_luma_couplings,
                cs_glob_luma_n_couplings + 1, cs_luma_coupling_t *);
    BFT_MALLOC(luma_coupling, 1, cs_luma_coupling_t);

    cs_glob_luma_couplings[cs_glob_luma_n_couplings] = luma_coupling;

    cs_glob_luma_n_couplings++;
  }

  luma_coupling->dim = dim;
  luma_coupling->ref_axis = ref_axis;

  luma_coupling->luma_name = NULL;

  if (luma_name != NULL) {
    BFT_MALLOC(luma_coupling->luma_name, strlen(luma_name) + 1, char);
    strcpy(luma_coupling->luma_name, luma_name);
  }
  else {
    BFT_MALLOC(luma_coupling->luma_name, 1, char);
    luma_coupling->luma_name[0] = '\0';
  }
  
  if (vars_in != NULL) {
	  BFT_MALLOC(luma_coupling->vars_in, strlen(vars_in) + 1, char);
	  strcpy(luma_coupling->vars_in, vars_in);
  }
  else {
	  bft_printf(_("\nVariables to write to CS not detected."));
      bft_printf_flush();
	  
	  BFT_MALLOC(luma_coupling->vars_in, 1, char);
	  luma_coupling->vars_in[0] = '\0';
  }
  
  printf("CS: vars_out %s \n", vars_out);
  if (vars_out != NULL) {
	  BFT_MALLOC(luma_coupling->vars_out, strlen(vars_out) + 1, char);
	  strcpy(luma_coupling->vars_out, vars_out);
  }
  else {
	  bft_printf(_("\nVariables to read to CS not detected."));
      bft_printf_flush();
	  
	  BFT_MALLOC(luma_coupling->vars_out, 1, char);
	  luma_coupling->vars_out[0] = '\0';
  }

  /* Selection criteria  */

  luma_coupling->n_b_locations = 0;
  luma_coupling->n_v_locations = 0;
  luma_coupling->b_location_ids = NULL;
  luma_coupling->v_location_ids = NULL;

  luma_coupling->faces = NULL;
  luma_coupling->cells = NULL;

  luma_coupling->allow_nearest = allow_nonmatching;
  luma_coupling->tolerance = tolerance;
  luma_coupling->verbosity = verbosity;
  luma_coupling->visualization = visualization;

  /* Initialize communicators */

#if defined(HAVE_MPI)

  luma_coupling->comm = MPI_COMM_NULL;
  luma_coupling->n_luma_ranks = 0;
  luma_coupling->luma_root_rank = -1;

#endif

  return  luma_coupling;
}



/*----------------------------------------------------------------------------*/

END_C_DECLS
