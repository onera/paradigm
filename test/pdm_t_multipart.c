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
 const int    n_faceBound,
 const int    n_faceJoin,
 const int    n_total_part,
 double      *VtxCoord,
 int         *face_vtx_idx,
 PDM_g_num_t *face_vtx,
 PDM_g_num_t *face_cell,
 int         *cell_face_idx,
 int         *cell_face,
 int         *FaceBoundIdx,
 PDM_g_num_t *FaceBound,
 int         *FaceJoinIdx,
 PDM_g_num_t *FaceJoin,
 int         *face_part_bound_proc_idx,
 int         *face_part_bound_part_idx,
 int         *face_part_bound

 )
{
  FILE * fp;
  fp = fopen(filename, "w");
  char endcomma = ' ';

  fprintf(fp, "{\n");
  fprintf(fp, "  \"nbVtx\" : %d,\n", n_vtx);
  fprintf(fp, "  \"nbFace\" : %d,\n", n_face);
  fprintf(fp, "  \"nbCell\" : %d,\n", n_cell);
  fprintf(fp, "  \"nbFaceBound\" : %d,\n", n_faceBound);
  fprintf(fp, "  \"nbFaceJoin\" : %d,\n", n_faceJoin);

  fprintf(fp, "  \"VtxCoord\" :\n  [");
  for (int k=0; k<n_vtx; k++)
  {
    fprintf(fp, "%c\n    [%12.5e,  %12.5e,  %12.5e]", endcomma, VtxCoord[3*k], VtxCoord[3*k+1], VtxCoord[3*k+2]);
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
  fprintf(fp, "  \"FaceBounds\" :\n  [");
  for (int k=0; k<n_faceBound; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=FaceBoundIdx[k]; i<FaceBoundIdx[k+1]-1; i++)
      fprintf(fp, "%d,  ", FaceBound[i]);
    // Cas particulier ou le groupe n'a aucun élément
    if (FaceBoundIdx[k+1] > FaceBoundIdx[k])
      fprintf(fp, "%d]", FaceBound[FaceBoundIdx[k+1]-1]);
    else
      fprintf(fp, "]");
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"FaceJoins\" :\n  [");
  for (int k=0; k<n_faceJoin; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=FaceJoinIdx[k]; i<FaceJoinIdx[k+1]-1; i++)
      fprintf(fp, "%d,  ", FaceJoin[i]);
    // Cas particulier ou le groupe n'a aucun élément
    if (FaceJoinIdx[k+1] > FaceJoinIdx[k])
      fprintf(fp, "%d]", FaceJoin[FaceJoinIdx[k+1]-1]);
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

  int nbVtx, nbFace, nbCell, nbBound, nbJoin;
  int zoneGId = -1;
  json_object_object_get_ex(zonedata, "nbVtx", &intdata);
  nbVtx = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbFace", &intdata);
  nbFace = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbCell", &intdata);
  nbCell = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbBoundary", &intdata);
  nbBound = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbInterface", &intdata);
  nbJoin = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "ZoneGId", &intdata);
  zoneGId = json_object_get_int(intdata);

  // PDM_printf("Found in json : nbVtx = %d, nbFace = %d, nbCell = %d, nbBound = %d, \
  nbJoin = %d for zoneGId %d\n", nbVtx, nbFace, nbCell, nbBound, nbJoin, zoneGId);

  json_object_object_get_ex(zonedata, "VtxCoord", &array);
  assert(json_object_array_length(array) == nbVtx);
  double * vtxCoord = (double *) malloc(3*nbVtx * sizeof(double));
  for (int i = 0; i < nbVtx; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    vtxCoord[3*i] = json_object_get_double(json_object_array_get_idx(arraydata, 0));
    vtxCoord[3*i+1] = json_object_get_double(json_object_array_get_idx(arraydata, 1));
    vtxCoord[3*i+2] = json_object_get_double(json_object_array_get_idx(arraydata, 2));
  }

  json_object_object_get_ex(zonedata, "NGonConnectivity", &array);
  assert(json_object_array_length(array) == nbFace);
  int * face_vtx_idx = (int *) malloc((nbFace+1) * sizeof(int));
  face_vtx_idx[0] = 0;
  //First gets size of face_vtx and fills face_vtx_idx
  for (int i = 0; i < nbFace; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    face_vtx_idx[i+1] = json_object_array_length(arraydata) + face_vtx_idx[i];
  }
  PDM_g_num_t * face_vtx = (PDM_g_num_t *) malloc(face_vtx_idx[nbFace] * sizeof(PDM_g_num_t));
  //Now fill face_vtx
  int iface = 0;
  for (int i = 0; i < nbFace; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    for (int k = 0; k < json_object_array_length(arraydata); k++)
    {
      face_vtx[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  json_object_object_get_ex(zonedata, "NGonParentElement", &array);
  assert(json_object_array_length(array) == nbFace);
  PDM_g_num_t * face_cell = (PDM_g_num_t *) malloc(2*nbFace * sizeof(PDM_g_num_t));
  for (int i = 0; i < nbFace; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    face_cell[2*i] = json_object_get_int(json_object_array_get_idx(arraydata, 0));
    face_cell[2*i+1] = json_object_get_int(json_object_array_get_idx(arraydata, 1));
  }

  //First get size of faceBound and fills faceBoundIdx
  int * faceBoundIdx = (int *) malloc((nbBound + 1) * sizeof(int));
  faceBoundIdx[0] = 0;
  json_object_object_get_ex(zonedata, "BoundGroups", &array);
  assert(json_object_array_length(array) == nbBound);
  for (int i = 0; i < nbBound; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    faceBoundIdx[i+1] = json_object_array_length(arraydata) + faceBoundIdx[i];
  }
  //Now fill faceBound
  PDM_g_num_t * faceBound = (PDM_g_num_t *) malloc(faceBoundIdx[nbBound] * sizeof(PDM_g_num_t));
  iface = 0;
  for (int i = 0; i < nbBound; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    for (int k = 0; k < json_object_array_length(arraydata); k++)
    {
      faceBound[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  //First get size of faceJoin and fills faceJoinIdx
  int * faceJoinIdx = (int *) malloc((nbJoin + 1) * sizeof(int));
  faceJoinIdx[0] = 0;
  json_object_object_get_ex(zonedata, "InterfaceGroups", &array);
  assert(json_object_array_length(array) == nbJoin);
  for (int i = 0; i < nbJoin; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    json_object_object_get_ex(arraydata, "PointList", &arraydata);
    faceJoinIdx[i+1] = json_object_array_length(arraydata) + faceJoinIdx[i];
  }
  //Now fill faceJoin
  PDM_g_num_t * faceJoin = (PDM_g_num_t *) malloc(faceJoinIdx[nbJoin] * sizeof(PDM_g_num_t));
  int * joinGIds = (int *) malloc(2*nbJoin * sizeof(int));
  iface = 0;
  for (int i = 0; i < nbJoin; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    json_object_object_get_ex(arraydata, "sourceGId", &intdata);
    joinGIds[2*i] = json_object_get_int(intdata);
    json_object_object_get_ex(arraydata, "targetGId", &intdata);
    joinGIds[2*i+1] = json_object_get_int(intdata);
    json_object_object_get_ex(arraydata, "PointList", &arraydata);
    for (int k = 0; k < json_object_array_length(arraydata); k++)
    {
      faceJoin[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  int dmeshId = PDM_dmesh_create(nbCell, nbFace, nbVtx, nbBound, nbJoin);
  PDM_dmesh_set(dmeshId, vtxCoord, face_vtx_idx, face_vtx, face_cell,
                faceBoundIdx, faceBound, joinGIds, faceJoinIdx, faceJoin);
  meshIds[blockId] = dmeshId;
  zoneIds[blockId] = zoneGId;
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
  PDM_g_num_t *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int              cubeid;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  int nbzone = 1;
  int mpartId = -1;
  int multipartSize = 1;
  int *dmeshIds;             // For a given block, store the identifier of the corresponding blockdata (size = nzone)
  int *dblockIds;            // For a given block, store the globalId of the parent zone

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

    dmeshIds  = (int *) malloc(nbzone * sizeof(int));
    dblockIds = (int *) malloc(nbzone * sizeof(int));
    int dmeshId = -1;
    dmeshId = PDM_dmesh_create(dn_cell, dn_face, dn_vtx, n_face_group, 0);
    PDM_dmesh_set(dmeshId, dvtx_coord, dface_vtx_idx, dface_vtx, dface_cell, dface_group_idx, dface_group, NULL, NULL, NULL);
    dmeshIds[0] = dmeshId;
    dblockIds[0] = 1;

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
    dmeshIds  = (int *) malloc(nbzone * sizeof(int));
    dblockIds = (int *) malloc(nbzone * sizeof(int));

    // 2. Read the data for each zone and store it in arrays
    for (int iblock=0; iblock < nbzone; iblock++)
    {
      _readJsonBlock(path, iblock, dmeshIds, dblockIds);
      // PDM_dmesh_dims_get(dmeshId, &dn_cell, &dn_face, &dn_vtx, &n_face_group);
      // PDM_dmesh_data_get(dmeshId, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell, &dface_group_idx, &dface_group);
    }
    // 3. Compute the multizone size (total number of different zones)
    int maxZoneGId = -1;
    for (int k = 0; k < nbzone; k++)
      maxZoneGId = (dblockIds[k] > maxZoneGId) ? dblockIds[k] : maxZoneGId;
    PDM_MPI_Allreduce(&maxZoneGId, &multipartSize, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  }

  int * n_partArray = (int *) malloc((multipartSize) * sizeof(int));
  for (int k=0; k<multipartSize; k++)
    n_partArray[k] = n_part;
  n_partArray[0] = (i_rank == 0) ? 0 : 1;
  n_partArray[1] = (i_rank == 0) ? 1 : 1;
   mpartId = PDM_multipart_create(multipartSize, n_partArray, PDM_FALSE, method, comm);
   PDM_printf("From exe : created a multipart object, id is %i \n", mpartId);
   for (int iblock=0; iblock < nbzone; iblock++)
   {
    // PDM_printf("[%i] check -- iblock %i : gid %i meshid %i\n", i_rank, iblock, dblockIds[iblock], dmeshIds[iblock]);
    PDM_multipart_register_block(mpartId, dblockIds[iblock]-1, dmeshIds[iblock]);
   }
   PDM_MPI_Barrier(comm);



  PDM_multipart_run_ppart(mpartId);

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_multipart_time_get(mpartId,
                         multipartSize-1,
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

  for (int iblock = 0; iblock < multipartSize; iblock++) {
    for (int i_part = 0; i_part < n_partArray[iblock]; i_part++) {

      int n_cell;
      int n_face;
      int n_face_part_bound;
      int n_vtx;
      int n_proc;
      int n_total_part;
      int scell_face;
      int sface_vtx;
      int sFaceBound;
      int n_faceBound;
      int sFaceJoin;
      int n_faceJoin;

      PDM_multipart_part_dim_get(mpartId,
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
                                 &sFaceBound,
                                 &n_faceBound,
                                 &sFaceJoin,
                                 &n_faceJoin);

      int          *cell_tag;
      int          *cell_face_idx;
      int          *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int          *face_tag;
      int          *face_cell;
      int          *face_vtx_idx;
      int          *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int          *face_part_bound_proc_idx;
      int          *face_part_bound_part_idx;
      int          *face_part_bound;
      int          *vtx_tag;
      double       *vtx;
      PDM_g_num_t *vtx_ln_to_gn;
      int          *faceBoundIdx;
      int          *faceBound;
      PDM_g_num_t *faceBoundLNToGN;
      int          *faceJoinIdx;
      int          *faceJoin;
      PDM_g_num_t *faceJoinLNToGN;


      PDM_multipart_part_val_get(mpartId,
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
                                 &faceBoundIdx,
                                 &faceBound,
                                 &faceBoundLNToGN,
                                 &faceJoinIdx,
                                 &faceJoin,
                                 &faceJoinLNToGN);

      char filename[25];
      snprintf(filename,sizeof(filename),"partitionned_P%dB%dN%d.json",i_rank, iblock, i_part);
      _dumpJsonData(filename, iblock, n_vtx, n_face, n_cell, n_faceBound, n_faceJoin, n_total_part, vtx,
                    face_vtx_idx, face_vtx, face_cell, cell_face_idx, cell_face,faceBoundIdx, faceBound,
                    faceJoinIdx, faceJoin, face_part_bound_proc_idx, face_part_bound_part_idx, face_part_bound);

    }
  }

  PDM_multipart_free(mpartId);

  for (int iblock = 0; iblock<nbzone; iblock++)
    PDM_dmesh_free(dmeshIds[iblock]);
  free(dmeshIds);
  free(dblockIds);
  free(n_partArray);

  if (strcmp(distri_dir, "\0") == 0) // No distributed data provided
  {
  //Desalloue les pointeurs de la structure interne dface_cell, dface_vtx_idx, dface_vtx, dvtx_coord, dface_group_idx, dface_group
    PDM_dcube_gen_free(cubeid);
  }

  PDM_MPI_Finalize();

  return 0;
}
