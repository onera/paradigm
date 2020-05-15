#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include</stck/jcoulet/dev/dev-Test/json-c/install/include/json-c/json.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_dmesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/
  #define MAXBUFLEN 10000000 // Taille max (en char --> wc file.dat) du ficher json
/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
	         int           *method,
           char         **distri_dir)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-distri_dir") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *distri_dir = (argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/**
 *
 * \Write data in a file in order to rebuild a CGNS
 *
 */
static void _dumpJsonData
(
 const char  *filename,
 const int    iblock,
 const int    n_vtx,
 const int    n_face,
 const int    n_cell,
 const int    n_face_bound,
 const int    n_face_join,
 const int    n_total_part,
 double      *vtx_coord,
 int         *face_vtx_idx,
 PDM_g_num_t *face_vtx,
 PDM_g_num_t *face_cell,
 int         *cell_face_idx,
 int         *cell_face,
 int         *face_bound_idx,
 PDM_g_num_t *face_bound,
 int         *face_join_idx,
 PDM_g_num_t *face_join,
 int         *face_part_bound_proc_idx,
 int         *face_part_bound_part_idx,
 int         *face_part_bound

 )
{
  FILE * fp;
  fp = fopen(filename, "w");
  char endcomma = ' ';

  fprintf(fp, "{\n");
  fprintf(fp, "  \"nb_vtx\" : %d,\n", n_vtx);
  fprintf(fp, "  \"nb_face\" : %d,\n", n_face);
  fprintf(fp, "  \"nb_cell\" : %d,\n", n_cell);
  fprintf(fp, "  \"nb_face_bound\" : %d,\n", n_face_bound);
  fprintf(fp, "  \"nb_face_join\" : %d,\n", n_face_join);

  fprintf(fp, "  \"vtx_coord\" :\n  [");
  for (int k=0; k<n_vtx; k++)
  {
    fprintf(fp, "%c\n    [%12.5e,  %12.5e,  %12.5e]", endcomma, vtx_coord[3*k], vtx_coord[3*k+1], vtx_coord[3*k+2]);
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"NGonConnectivity\" :\n  [");
  for (int k=0; k<n_face; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=face_vtx_idx[k]; i<face_vtx_idx[k+1]-1; i++)
      fprintf(fp, "%d,  ", face_vtx[i]);
    fprintf(fp, "%d]", face_vtx[face_vtx_idx[k+1]-1]);
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"NGonParentElement\" :\n  [");
  for (int k=0; k<n_face; k++)
  {
    fprintf(fp, "%c\n    [%d, %d]", endcomma, face_cell[2*k], face_cell[2*k+1]);
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  if (cell_face != NULL)
  {
    endcomma = ' ';
    fprintf(fp, "  \"n_faceConnectivity\" :\n  [");
    for (int k=0; k<n_cell; k++)
    {
      fprintf(fp, "%c\n    [", endcomma);
      for (int i=cell_face_idx[k]; i<cell_face_idx[k+1]-1; i++)
        fprintf(fp, "%d,  ", cell_face[i]);
      fprintf(fp, "%d]", cell_face[cell_face_idx[k+1]-1]);
      endcomma = ',';
    }
    fprintf(fp, "\n  ],\n");
  }

  endcomma = ' ';
  fprintf(fp, "  \"face_bounds\" :\n  [");
  for (int k=0; k<n_face_bound; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=face_bound_idx[k]; i<face_bound_idx[k+1]-1; i++)
      fprintf(fp, "%d,  ", face_bound[i]);
    // Cas particulier ou le groupe n'a aucun élément
    if (face_bound_idx[k+1] > face_bound_idx[k])
      fprintf(fp, "%d]", face_bound[face_bound_idx[k+1]-1]);
    else
      fprintf(fp, "]");
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"face_joins\" :\n  [");
  for (int k=0; k<n_face_join; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=face_join_idx[k]; i<face_join_idx[k+1]-1; i++)
      fprintf(fp, "%d,  ", face_join[i]);
    // Cas particulier ou le groupe n'a aucun élément
    if (face_join_idx[k+1] > face_join_idx[k])
      fprintf(fp, "%d]", face_join[face_join_idx[k+1]-1]);
    else
      fprintf(fp, "]");
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"PartBoundGroups\" :\n  {");
  for (int iMatch=0; iMatch < n_total_part; iMatch++) //Local number of parts
  {
    if (face_part_bound_part_idx[iMatch+1] != face_part_bound_part_idx[iMatch])
    {
      int lastindex = -1;
      fprintf(fp, "%c\n    \"Match%i\" :\n    {", endcomma, iMatch);

      fprintf(fp, "\n      \"PointList\" : [");
      for (int k=face_part_bound_part_idx[iMatch]; k < face_part_bound_part_idx[iMatch+1]-1; k++)
        fprintf(fp, "%i, ", face_part_bound[4*k] );
      lastindex = 4*(face_part_bound_part_idx[iMatch+1]-1);
      fprintf(fp, "%i]", face_part_bound[lastindex]);

      fprintf(fp, ",\n      \"PointListDonor\" : [");
          for (int k=face_part_bound_part_idx[iMatch]; k < face_part_bound_part_idx[iMatch+1]-1; k++)
        fprintf(fp, "%i, ", face_part_bound[4*k+3] );
      lastindex = 4*(face_part_bound_part_idx[iMatch+1]-1)+3;
      fprintf(fp, "%i]", face_part_bound[lastindex]);

      char targetname[25];
      int procOpp = face_part_bound[4*face_part_bound_part_idx[iMatch] + 1];
      int partOpp = face_part_bound[4*face_part_bound_part_idx[iMatch] + 2];
      snprintf(targetname,sizeof(targetname),"P%dB%dN%d",procOpp, iblock, partOpp-1);
      fprintf(fp, ",\n      \"targetPart\" : \"%s\"", targetname);
      fprintf(fp, "\n    }");
    endcomma = ',';
    }
  }
  fprintf(fp, "\n  }\n");

  fprintf(fp, "}");
  fclose(fp);
}

static int _readJsonNumberOfBlocks
(
 const char *filename
)
{
  FILE *fp;
  char buffer[MAXBUFLEN];
  fp = fopen(filename, "r");
  fread(buffer, sizeof(char), MAXBUFLEN, fp);
  fclose(fp);

  struct json_object *parsed_json;
  struct json_object *intdata;
  int nbZones = -1;

  parsed_json = json_tokener_parse(buffer);
  json_object_object_get_ex(parsed_json, "nbZones", &intdata);
  nbZones = json_object_get_int(intdata);
  return nbZones;
}

static void _readJsonBlock
(
  const char  *filename,
  const int    blockId,
  int         *meshIds,
  int         *zoneIds
)
{
  PDM_printf("Parsing json data from %s\n", filename);
  FILE *fp;
  char buffer[MAXBUFLEN];
  fp = fopen(filename, "r");
  fread(buffer, sizeof(char), MAXBUFLEN, fp);
  fclose(fp);

  struct json_object *parsed_json;
  struct json_object *zonedata;
  struct json_object *intdata;
  struct json_object *array;
  struct json_object *arraydata;

  parsed_json = json_tokener_parse(buffer);
  // json_object_object_get_ex(parsed_json, "nbZones", &intdata);
  // PDM_printf("Nb of zones in json is %d\n\n", json_object_get_int(intdata));
  json_object_object_get_ex(parsed_json, "Zones", &zonedata);
  zonedata = json_object_array_get_idx(zonedata, blockId);

  int nb_vtx, nb_face, nb_cell, nb_bound, nb_join;
  int zone_gid = -1;
  json_object_object_get_ex(zonedata, "nb_vtx", &intdata);
  nb_vtx = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nb_face", &intdata);
  nb_face = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nb_cell", &intdata);
  nb_cell = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nb_boundary", &intdata);
  nb_bound = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbInterface", &intdata);
  nb_join = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "zone_gid", &intdata);
  zone_gid = json_object_get_int(intdata);

  // PDM_printf("Found in json : nb_vtx = %d, nb_face = %d, nb_cell = %d, nb_bound = %d, \
  nb_join = %d for zone_gid %d\n", nb_vtx, nb_face, nb_cell, nb_bound, nb_join, zone_gid);

  json_object_object_get_ex(zonedata, "vtx_coord", &array);
  assert(json_object_array_length(array) == nb_vtx);
  double * vtx_coord = (double *) malloc(3*nb_vtx * sizeof(double));
  for (int i = 0; i < nb_vtx; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    vtx_coord[3*i] = json_object_get_double(json_object_array_get_idx(arraydata, 0));
    vtx_coord[3*i+1] = json_object_get_double(json_object_array_get_idx(arraydata, 1));
    vtx_coord[3*i+2] = json_object_get_double(json_object_array_get_idx(arraydata, 2));
  }

  json_object_object_get_ex(zonedata, "NGonConnectivity", &array);
  assert(json_object_array_length(array) == nb_face);
  int * face_vtx_idx = (int *) malloc((nb_face+1) * sizeof(int));
  face_vtx_idx[0] = 0;
  //First gets size of face_vtx and fills face_vtx_idx
  for (int i = 0; i < nb_face; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    face_vtx_idx[i+1] = json_object_array_length(arraydata) + face_vtx_idx[i];
  }
  PDM_g_num_t * face_vtx = (PDM_g_num_t *) malloc(face_vtx_idx[nb_face] * sizeof(PDM_g_num_t));
  //Now fill face_vtx
  int iface = 0;
  for (int i = 0; i < nb_face; i++) {
    arraydata = json_object_array_get_idx(array, i);
    for (int k = 0; k < json_object_array_length(arraydata); k++) {
      face_vtx[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  json_object_object_get_ex(zonedata, "NGonParentElement", &array);
  assert(json_object_array_length(array) == nb_face);
  PDM_g_num_t * face_cell = (PDM_g_num_t *) malloc(2*nb_face * sizeof(PDM_g_num_t));
  for (int i = 0; i < nb_face; i++) {
    arraydata = json_object_array_get_idx(array, i);
    face_cell[2*i] = json_object_get_int(json_object_array_get_idx(arraydata, 0));
    face_cell[2*i+1] = json_object_get_int(json_object_array_get_idx(arraydata, 1));
  }

  //First get size of face_bound and fills face_bound_idx
  int * face_bound_idx = (int *) malloc((nb_bound + 1) * sizeof(int));
  face_bound_idx[0] = 0;
  json_object_object_get_ex(zonedata, "BoundGroups", &array);
  assert(json_object_array_length(array) == nb_bound);
  for (int i = 0; i < nb_bound; i++) {
    arraydata = json_object_array_get_idx(array, i);
    face_bound_idx[i+1] = json_object_array_length(arraydata) + face_bound_idx[i];
  }
  //Now fill face_bound
  PDM_g_num_t * face_bound = (PDM_g_num_t *) malloc(face_bound_idx[nb_bound] * sizeof(PDM_g_num_t));
  iface = 0;
  for (int i = 0; i < nb_bound; i++) {
    arraydata = json_object_array_get_idx(array, i);
    for (int k = 0; k < json_object_array_length(arraydata); k++) {
      face_bound[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  //First get size of face_join and fills face_join_idx
  int * face_join_idx = (int *) malloc((nb_join + 1) * sizeof(int));
  face_join_idx[0] = 0;
  json_object_object_get_ex(zonedata, "InterfaceGroups", &array);
  assert(json_object_array_length(array) == nb_join);
  for (int i = 0; i < nb_join; i++) {
    arraydata = json_object_array_get_idx(array, i);
    json_object_object_get_ex(arraydata, "PointList", &arraydata);
    face_join_idx[i+1] = json_object_array_length(arraydata) + face_join_idx[i];
  }
  //Now fill face_join
  PDM_g_num_t * face_join = (PDM_g_num_t *) malloc(face_join_idx[nb_join] * sizeof(PDM_g_num_t));
  int * joinGIds = (int *) malloc(2*nb_join * sizeof(int));
  iface = 0;
  for (int i = 0; i < nb_join; i++) {
    arraydata = json_object_array_get_idx(array, i);
    json_object_object_get_ex(arraydata, "sourceGId", &intdata);
    joinGIds[2*i] = json_object_get_int(intdata);
    json_object_object_get_ex(arraydata, "targetGId", &intdata);
    joinGIds[2*i+1] = json_object_get_int(intdata);
    json_object_object_get_ex(arraydata, "PointList", &arraydata);
    for (int k = 0; k < json_object_array_length(arraydata); k++) {
      face_join[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  int dmesh_id = PDM_dmesh_create(nb_cell, nb_face, nb_vtx, nb_bound, nb_join);
  PDM_dmesh_set(dmesh_id, vtx_coord, face_vtx_idx, face_vtx, face_cell,
                face_bound_idx, face_bound, joinGIds, face_join_idx, face_join);
  meshIds[blockId] = dmesh_id;
  zoneIds[blockId] = zone_gid;
  //Fuite mémoire -> les données allouées ici sont perdues
}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length  = 1.;
  int                n_part   = 1;
  char               *distri_dir = strdup("\0");

#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             (int *) &method,
             &distri_dir);

  /*
   *  Init
   */

  int i_rank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int              cubeid;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  int nbzone = 1;
  int mpart_id = -1;
  int multipart_size = 1;
  int *dmesh_ids;             // For a given block, store the identifier of the corresponding blockdata (size = nzone)
  int *dblock_ids;            // For a given block, store the globalId of the parent zone

  if (strcmp(distri_dir, "\0") == 0) // No distributed data provided -> generate cube (only 1 zone)
  {
    //Alloue les pointeurs de la structure interne dface_cell, dface_vtx_idx, dface_vtx, dvtx_coord, dface_group_idx, dface_group
    PDM_dcube_gen_init(&cubeid,
                        comm,
                        n_vtx_seg,
                        length,
              		      0.,
  		                  0.,
  		                  0.);

    PDM_dcube_gen_dim_get(cubeid,
                           &n_face_group,
                           &dn_cell,
                           &dn_face,
                           &dn_vtx,
                           &dface_vtxL,
                           &dFaceGroupL);

    PDM_dcube_gen_data_get(cubeid,
                            &dface_cell,
                            &dface_vtx_idx,
                            &dface_vtx,
                            &dvtx_coord,
                            &dface_group_idx,
                            &dface_group);

    dmesh_ids  = (int *) malloc(nbzone * sizeof(int));
    dblock_ids = (int *) malloc(nbzone * sizeof(int));
    int dmesh_id = -1;
    dmesh_id = PDM_dmesh_create(dn_cell, dn_face, dn_vtx, n_face_group, 0);
    PDM_dmesh_set(dmesh_id, dvtx_coord, dface_vtx_idx, dface_vtx, dface_cell, dface_group_idx, dface_group, NULL, NULL, NULL);
    dmesh_ids[0] = dmesh_id;
    dblock_ids[0] = 1;

    // Dump arrays
    char filename[25];
    snprintf(filename,sizeof(filename),"distributed_P%d.json",i_rank);
    _dumpJsonData(filename, 0, dn_vtx, dn_face, dn_cell, n_face_group, 0, 0,
                  dvtx_coord, dface_vtx_idx, dface_vtx, dface_cell, NULL, NULL, dface_group_idx, dface_group,
                  NULL, NULL, NULL, NULL, NULL);
  }
  else  // -> parse Json data
  {
    char filename[25];
    snprintf(filename,sizeof(filename),"distributed_P%d.json",i_rank);
    char path[50];
    strcpy(path, distri_dir);
    strcat(path, filename);
    // 1. Get number of zones known by each zone
    nbzone = _readJsonNumberOfBlocks(path);
    PDM_printf("[%i] Found %d blocks in distributed data\n", i_rank, nbzone);
    dmesh_ids  = (int *) malloc(nbzone * sizeof(int));
    dblock_ids = (int *) malloc(nbzone * sizeof(int));

    // 2. Read the data for each zone and store it in arrays
    for (int iblock=0; iblock < nbzone; iblock++)
    {
      _readJsonBlock(path, iblock, dmesh_ids, dblock_ids);
      // PDM_dmesh_dims_get(dmesh_id, &dn_cell, &dn_face, &dn_vtx, &n_face_group);
      // PDM_dmesh_data_get(dmesh_id, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell, &dface_group_idx, &dface_group);
    }
    // 3. Compute the multizone size (total number of different zones)
    int maxzone_gid = -1;
    for (int k = 0; k < nbzone; k++)
      maxzone_gid = (dblock_ids[k] > maxzone_gid) ? dblock_ids[k] : maxzone_gid;
    PDM_MPI_Allreduce(&maxzone_gid, &multipart_size, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  }

  int * n_part_array = (int *) malloc((multipart_size) * sizeof(int));
  for (int k=0; k<multipart_size; k++)
    n_part_array[k] = n_part;
  n_part_array[0] = (i_rank == 0) ? 0 : 1;
  n_part_array[1] = (i_rank == 0) ? 1 : 1;
   mpart_id = PDM_multipart_create(multipart_size, n_part_array, PDM_FALSE, method, comm);
   PDM_printf("From exe : created a multipart object, id is %i \n", mpart_id);
   for (int iblock=0; iblock < nbzone; iblock++)
   {
    // PDM_printf("[%i] check -- iblock %i : gid %i meshid %i\n", i_rank, iblock, dblock_ids[iblock], dmesh_ids[iblock]);
    PDM_multipart_register_block(mpart_id, dblock_ids[iblock]-1, dmesh_ids[iblock]);
   }
   PDM_MPI_Barrier(comm);



  PDM_multipart_run_ppart(mpart_id);

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_multipart_time_get(mpart_id,
                         multipart_size-1,
                         &elapsed,
                         &cpu,
                         &cpu_user,
                         &cpu_sys);

  if (i_rank == 0)
  {
    PDM_printf("[%i]   - elapsed total                    : %12.5e\n", i_rank, elapsed[0]);
    PDM_printf("[%i]   - elapsed building graph           : %12.5e\n", i_rank, elapsed[1]);
    PDM_printf("[%i]   - elapsed splitting graph          : %12.5e\n", i_rank, elapsed[2]);
    PDM_printf("[%i]   - elapsed building mesh partitions : %12.5e\n", i_rank, elapsed[3]);

    PDM_printf("[%i]   - cpu total                        : %12.5e\n", i_rank, cpu[0]);
    PDM_printf("[%i]   - cpu building graph               : %12.5e\n", i_rank, cpu[1]);
    PDM_printf("[%i]   - cpu splitting graph              : %12.5e\n", i_rank, cpu[2]);
    PDM_printf("[%i]   - cpu building mesh partitions     : %12.5e\n", i_rank, cpu[3]);

    PDM_printf("[%i]   - cpu_user total                   : %12.5e\n", i_rank, cpu_user[0]);
    PDM_printf("[%i]   - cpu_user building graph          : %12.5e\n", i_rank, cpu_user[1]);
    PDM_printf("[%i]   - cpu_user splitting graph         : %12.5e\n", i_rank, cpu_user[2]);
    PDM_printf("[%i]   - cpu_user building mesh partitions: %12.5e\n", i_rank, cpu_user[3]);

    PDM_printf("[%i]   - cpu_sys total                    : %12.5e\n", i_rank, cpu_sys[0]);
    PDM_printf("[%i]   - cpu_sys building graph           : %12.5e\n", i_rank, cpu_sys[1]);
    PDM_printf("[%i]   - cpu_sys splitting graph          : %12.5e\n", i_rank, cpu_sys[2]);
    PDM_printf("[%i]   - cpu_sys building mesh partitions : %12.5e\n", i_rank, cpu_sys[3]);
  }

  for (int iblock = 0; iblock < multipart_size; iblock++) {
    for (int i_part = 0; i_part < n_part_array[iblock]; i_part++) {

      int n_cell;
      int n_face;
      int n_face_part_bound;
      int n_vtx;
      int n_proc;
      int n_total_part;
      int scell_face;
      int sface_vtx;
      int sface_bound;
      int n_face_bound;
      int sface_join;
      int n_face_join;

      PDM_multipart_part_dim_get(mpart_id,
                                 iblock,
                                 i_part,
                                 &n_cell,
                                 &n_face,
                                 &n_face_part_bound,
                                 &n_vtx,
                                 &n_proc,
                                 &n_total_part,
                                 &scell_face,
                                 &sface_vtx,
                                 &sface_bound,
                                 &n_face_bound,
                                 &sface_join,
                                 &n_face_join);

      int          *cell_tag;
      int          *cell_face_idx;
      int          *cell_face;
      PDM_g_num_t  *cell_ln_to_gn;
      int          *face_tag;
      int          *face_cell;
      int          *face_vtx_idx;
      int          *face_vtx;
      PDM_g_num_t  *face_ln_to_gn;
      int          *face_part_bound_proc_idx;
      int          *face_part_bound_part_idx;
      int          *face_part_bound;
      int          *vtx_tag;
      double       *vtx;
      PDM_g_num_t  *vtx_ln_to_gn;
      int          *face_bound_idx;
      int          *face_bound;
      PDM_g_num_t  *face_bound_ln_to_gn;
      int          *face_join_idx;
      int          *face_join;
      PDM_g_num_t  *face_join_ln_to_gn;


      PDM_multipart_part_val_get(mpart_id,
                                 iblock,
                                 i_part,
                                 &cell_tag,
                                 &cell_face_idx,
                                 &cell_face,
                                 &cell_ln_to_gn,
                                 &face_tag,
                                 &face_cell,
                                 &face_vtx_idx,
                                 &face_vtx,
                                 &face_ln_to_gn,
                                 &face_part_bound_proc_idx,
                                 &face_part_bound_part_idx,
                                 &face_part_bound,
                                 &vtx_tag,
                                 &vtx,
                                 &vtx_ln_to_gn,
                                 &face_bound_idx,
                                 &face_bound,
                                 &face_bound_ln_to_gn,
                                 &face_join_idx,
                                 &face_join,
                                 &face_join_ln_to_gn);

      char filename[25];
      snprintf(filename,sizeof(filename),"partitionned_P%dB%dN%d.json",i_rank, iblock, i_part);
      _dumpJsonData(filename, iblock, n_vtx, n_face, n_cell, n_face_bound, n_face_join, n_total_part, vtx,
                    face_vtx_idx, face_vtx, face_cell, cell_face_idx, cell_face,face_bound_idx, face_bound,
                    face_join_idx, face_join, face_part_bound_proc_idx, face_part_bound_part_idx, face_part_bound);

    }
  }

  PDM_multipart_free(mpart_id);

  for (int iblock = 0; iblock<nbzone; iblock++)
    PDM_dmesh_free(dmesh_ids[iblock]);
  free(dmesh_ids);
  free(dblock_ids);
  free(n_part_array);

  if (strcmp(distri_dir, "\0") == 0) // No distributed data provided
  {
  //Desalloue les pointeurs de la structure interne dface_cell, dface_vtx_idx, dface_vtx, dvtx_coord, dface_group_idx, dface_group
    PDM_dcube_gen_free(cubeid);
  }

  PDM_MPI_Finalize();

  return 0;
}
