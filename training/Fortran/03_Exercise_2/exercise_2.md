---
jupytext:
  text_representation:
    extension: '.md'
    format_name: myst
    format_version: '0.7'
    jupytext_version: 1.4.0+dev
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Exercise 2 : Localization of a point cloud inside a mesh

+++

*(Load custom magics)*

```{code-cell} ipython3
import os, sys
module_path = os.path.abspath(os.path.join('../../utils'))
if module_path not in sys.path:
    sys.path.append(module_path)
```

```{code-cell}
%reload_ext visu_magics
%reload_ext code_magics
```

Your job is to fill the code cells left blank using the API referenced [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_doc_pretty/user_manual/prepro_algo/index.html#python-api).

+++

## Use modules

```{code-cell}
%%code_block -p exercise_2 -i 1
#include "pdm_configf.h"

program exercise_2

  use pdm
  use pdm_pointer_array
  use pdm_generate_mesh
  use pdm_mesh_location
  use pdm_part_to_part
  use iso_c_binding
  use pdm_writer_wrapper

  implicit none

  include "mpif.h"
```

+++

## Declare variables

```{code-cell}
%%code_block -p exercise_2 -i 2
  !---------------------------------------------------------------
  integer,              parameter    :: comm = MPI_COMM_WORLD

  integer(pdm_g_num_s), parameter    :: src_n_vtx_seg = 10
  integer(c_int),       parameter    :: src_n_part = 1

  integer(pdm_g_num_s), parameter    :: tgt_n_vtx_seg = 10
  integer(c_int),       parameter    :: tgt_n_part = 1
  double precision,     parameter    :: tgt_xmin = 0.3d0
  double precision,     parameter    :: tgt_ymin = 0.3d0

  logical,              parameter    :: nodal = .false.

  integer(pdm_l_num_s),      pointer :: src_n_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: src_n_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: src_n_face(:)     => null()
  type(PDM_pointer_array_t), pointer :: src_vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: src_edge_vtx      => null()
  type(PDM_pointer_array_t), pointer :: src_face_edge_idx => null()
  type(PDM_pointer_array_t), pointer :: src_face_edge     => null()
  type(PDM_pointer_array_t), pointer :: src_face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: src_vtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: src_edge_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: src_face_ln_to_gn => null()

  integer(pdm_l_num_s),      pointer :: tgt_n_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: tgt_n_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: tgt_n_face(:)     => null()
  type(PDM_pointer_array_t), pointer :: tgt_vtx_coord     => null()
  type(PDM_pointer_array_t), pointer :: tgt_edge_vtx      => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_edge_idx => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_edge     => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_vtx      => null()
  type(PDM_pointer_array_t), pointer :: tgt_vtx_ln_to_gn  => null()
  type(PDM_pointer_array_t), pointer :: tgt_edge_ln_to_gn => null()
  type(PDM_pointer_array_t), pointer :: tgt_face_ln_to_gn => null()

  type(c_ptr)                        :: mesh_loc = C_NULL_PTR
  integer(pdm_g_num_s),      pointer :: vtx_ln_to_gn(:)  => null()
  integer(pdm_g_num_s),      pointer :: face_ln_to_gn(:) => null()
  double precision,          pointer :: vtx_coord(:,:)   => null()
  integer(pdm_l_num_s),      pointer :: face_edge_idx(:) => null()
  integer(pdm_l_num_s),      pointer :: face_edge(:)     => null()
  integer(pdm_l_num_s),      pointer :: face_vtx(:)      => null()
  integer(pdm_l_num_s),      pointer :: edge_vtx(:)      => null()

  type(PDM_pointer_array_t), pointer :: src_send_field1  => null()
  type(PDM_pointer_array_t), pointer :: tgt_recv_field1  => null()
  type(PDM_pointer_array_t), pointer :: src_send_field2  => null()
  type(PDM_pointer_array_t), pointer :: tgt_recv_field2  => null()
  double precision,          pointer :: field1(:)        => null()
  double precision,          pointer :: field2(:)        => null()

  integer                            :: n_pts, elt_n_vtx
  integer(pdm_l_num_s),      pointer :: src_to_tgt_idx(:)            => null()
  integer(pdm_g_num_s),      pointer :: points_gnum(:)               => null()
  double precision,          pointer :: points_coords(:,:)           => null()
  double precision,          pointer :: points_uvw(:,:)              => null()
  integer(pdm_l_num_s),      pointer :: points_weights_idx(:)        => null()
  double precision,          pointer :: points_weights(:)            => null()
  double precision,          pointer :: points_dist2(:)              => null()
  double precision,          pointer :: points_projected_coords(:,:) => null()
  integer(pdm_l_num_s),      pointer :: cell_vtx_idx(:)              => null()
  integer(pdm_l_num_s),      pointer :: cell_vtx(:)                  => null()

  type(c_ptr)                        :: ptp = C_NULL_PTR
  integer(c_int)                     :: request1, request2

  integer                            :: n_located
  integer(pdm_l_num_s),      pointer :: located(:)   => null()
  integer                            :: n_unlocated
  integer(pdm_l_num_s),      pointer :: unlocated(:) => null()

  integer(pdm_g_num_s),      pointer :: location(:)           => null()
  double precision,          pointer :: dist2(:)              => null()
  double precision,          pointer :: projected_coords(:,:) => null()
  double precision                   :: error(2)
  type(my_field_t)                   :: tgt_visu_fields(3)
  double precision,          pointer :: visu_field1(:)        => null()
  double precision,          pointer :: visu_field2(:)        => null()
  double precision,          pointer :: visu_field3(:)        => null()


  integer :: i_rank, ierr
  integer :: i_part, i_elt, i_pt, i_vtx, vtx_id, i
  !---------------------------------------------------------------
```

## Initialize MPI

```{code-cell}
%%code_block -p exercise_2 -i 3
  ! Initialize MPI
  call mpi_init(ierr)
  call mpi_comm_rank(comm, i_rank, ierr)

```

## Generate a partitioned "source" mesh

```{code-cell}
%%code_block -p exercise_2 -i 4
  ! Generate source mesh
  call pdm_generate_mesh_rectangle_ngon(comm,                         &
                                        PDM_MESH_NODAL_POLY_2D,       &
                                        0.d0,                         &
                                        0.d0,                         &
                                        0.d0,                         &
                                        1.d0,                         &
                                        1.d0,                         &
                                        src_n_vtx_seg,                &
                                        src_n_vtx_seg,                &
                                        src_n_part,                   &
                                        PDM_SPLIT_DUAL_WITH_PARMETIS, &
                                        src_n_vtx,                    &
                                        src_n_edge,                   &
                                        src_n_face,                   &
                                        src_vtx_coord,                &
                                        src_edge_vtx,                 &
                                        src_face_edge_idx,            &
                                        src_face_edge,                &
                                        src_face_vtx,                 &
                                        src_vtx_ln_to_gn,             &
                                        src_edge_ln_to_gn,            &
                                        src_face_ln_to_gn,            &
                                        0.4d0)
```

## Generate a partitioned "target" mesh
We will use its vertices as a point cloud.

```{code-cell}
%%code_block -p exercise_2 -i 5
  ! Generate target mesh
  call pdm_generate_mesh_rectangle_ngon(comm,                        &
                                        PDM_MESH_NODAL_QUAD4,        &
                                        tgt_xmin,                    &
                                        tgt_ymin,                    &
                                        0.d0,                        &
                                        1.d0,                        &
                                        1.d0,                        &
                                        tgt_n_vtx_seg,               &
                                        tgt_n_vtx_seg,               &
                                        tgt_n_part,                  &
                                        PDM_SPLIT_DUAL_WITH_HILBERT, &
                                        tgt_n_vtx,                   &
                                        tgt_n_edge,                  &
                                        tgt_n_face,                  &
                                        tgt_vtx_coord,               &
                                        tgt_edge_vtx,                &
                                        tgt_face_edge_idx,           &
                                        tgt_face_edge,               &
                                        tgt_face_vtx,                &
                                        tgt_vtx_ln_to_gn,            &
                                        tgt_edge_ln_to_gn,           &
                                        tgt_face_ln_to_gn)
```

## Create the `Mesh Location` object

```{code-cell}
%%code_block -p exercise_2 -i 6
  ! Create Mesh Location object
  call pdm_mesh_location_create(mesh_loc, &
                                1,        &
                                comm,     &
                                PDM_OWNERSHIP_KEEP)
```

## Set the target point cloud

```{code-cell}
%%code_block -p exercise_2 -i 7
  ! Set target point cloud
  call pdm_mesh_location_n_part_cloud_set(mesh_loc, &
                                          0,        &
                                          tgt_n_part)

  do i_part = 1, tgt_n_part
    call pdm_pointer_array_part_get(tgt_vtx_ln_to_gn, &
                                    i_part-1,         &
                                    vtx_ln_to_gn)

    call pdm_pointer_array_part_get(tgt_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)

    call pdm_mesh_location_cloud_set(mesh_loc,          &
                                     0,                 &
                                     i_part-1,          &
                                     tgt_n_vtx(i_part), &
                                     vtx_coord,         &
                                     vtx_ln_to_gn)
  enddo
```

## Set the source mesh
Here you have essentially two options :
- you can either define the mesh with "nodal" connectivity (i.e. Finite-Element style)
- or with "descending" connectivity (i.e. Finite-Volume style)

```{code-cell}
%%code_block -p exercise_2 -i 8
  ! Set source mesh
  call pdm_mesh_location_mesh_n_part_set(mesh_loc, src_n_part)

  do i_part = 1, src_n_part

    call pdm_pointer_array_part_get(src_vtx_ln_to_gn, &
                                    i_part-1,         &
                                    vtx_ln_to_gn)

    call pdm_pointer_array_part_get(src_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)


    call pdm_pointer_array_part_get(src_face_ln_to_gn, &
                                    i_part-1,          &
                                    face_ln_to_gn)

    call pdm_pointer_array_part_get(src_face_edge_idx, &
                                    i_part-1,          &
                                    face_edge_idx)


    if (nodal) then
      call pdm_pointer_array_part_get(src_face_vtx, &
                                      i_part-1,     &
                                      face_vtx)

      call pdm_mesh_location_nodal_part_set_2d(mesh_loc,           &
                                               i_part-1,           &
                                               src_n_face(i_part), &
                                               face_edge_idx,      &
                                               face_vtx,           &
                                               face_ln_to_gn,      &
                                               src_n_vtx(i_part),  &
                                               vtx_coord,          &
                                               vtx_ln_to_gn)
    else
      call pdm_pointer_array_part_get(src_face_edge, &
                                      i_part-1,      &
                                      face_edge)

      call pdm_pointer_array_part_get(src_edge_vtx, &
                                      i_part-1,     &
                                      edge_vtx)

      call pdm_mesh_location_part_set_2d(mesh_loc,           &
                                         i_part-1,           &
                                         src_n_face(i_part), &
                                         face_edge_idx,      &
                                         face_edge,           &
                                         face_ln_to_gn,      &
                                         src_n_edge(i_part), &
                                         edge_vtx,           &
                                         src_n_vtx(i_part),  &
                                         vtx_coord,          &
                                         vtx_ln_to_gn)

    endif
  enddo
```

## Set some optional parameters

```{code-cell}
%%code_block -p exercise_2 -i 9
  ! Set the geometric tolerance (optional)
  call pdm_mesh_location_tolerance_set(mesh_loc, 1.d-3)

  ! Set the location preconditioning method (optional)
  call pdm_mesh_location_method_set(mesh_loc, PDM_MESH_LOCATION_OCTREE)
```

## Compute the localization
```{code-cell}
%%code_block -p exercise_2 -i 10
  ! Compute location
  call pdm_mesh_location_compute(mesh_loc)

  ! Dump elapsed and CPU times
  call pdm_mesh_location_dump_times(mesh_loc)
```

## Results

Now that we have located the target points in the source mesh, we can exchange data between the two.
To complete this exercise, we will interpolate two fields from the source mesh to the target cloud:
  1. a cell-based field: we can simply use the face global ids for such a field, and check it matches the location data.
  2. a node-based field: we can use the node coordinates.

First, compute the spatially interpolated fields on the source side.
For the first field, the interpolation is straightforward : the target value is simply the same as the host source.
The second field interpolation is trickier as you will need the cell->vertex connectivity built during the location computation to link the interpolation weights to the appropriate source nodes.

### Interpolate the first field (cell-based)
```{code-cell}
%%code_block -p exercise_2 -i 11
  ! Interpolate first field (cell-based)
  call pdm_pointer_array_create(src_send_field1, &
                                src_n_part,      &
                                PDM_TYPE_DOUBLE)
  do i_part = 1, src_n_part

    call pdm_mesh_location_points_in_elt_get(mesh_loc,           &
                                             0,                  &
                                             i_part-1,           &
                                             src_to_tgt_idx,     &
                                             points_gnum,        &
                                             points_coords,      &
                                             points_uvw,         &
                                             points_weights_idx, &
                                             points_weights,     &
                                             points_dist2,       &
                                             points_projected_coords)

    call pdm_pointer_array_part_get(src_face_ln_to_gn, &
                                    i_part-1,          &
                                    face_ln_to_gn)

    n_pts = src_to_tgt_idx(src_n_face(i_part)+1)

    allocate(field1(n_pts))
    do i_elt = 1, src_n_face(i_part)
      do i_pt = src_to_tgt_idx(i_elt)+1, src_to_tgt_idx(i_elt+1)
        field1(i_pt) = face_ln_to_gn(i_elt)
      enddo
    enddo

    call pdm_pointer_array_part_set(src_send_field1, &
                                    i_part-1,        &
                                    field1)

  enddo
```

### Interpolate the second field (node-based)

```{code-cell}
%%code_block -p exercise_2 -i 12
  ! Interpolate second field (node-based)
  call pdm_pointer_array_create(src_send_field2, &
                                src_n_part,      &
                                PDM_TYPE_DOUBLE)
  do i_part = 1, src_n_part

    call pdm_mesh_location_points_in_elt_get(mesh_loc,           &
                                             0,                  &
                                             i_part-1,           &
                                             src_to_tgt_idx,     &
                                             points_gnum,        &
                                             points_coords,      &
                                             points_uvw,         &
                                             points_weights_idx, &
                                             points_weights,     &
                                             points_dist2,       &
                                             points_projected_coords)

    call pdm_pointer_array_part_get(src_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)

    call pdm_mesh_location_cell_vertex_get(mesh_loc,     &
                                           i_part-1,     &
                                           cell_vtx_idx, &
                                           cell_vtx)

    n_pts = src_to_tgt_idx(src_n_face(i_part)+1)

    allocate(field2(n_pts))
    do i_elt = 1, src_n_face(i_part)
      do i_pt = src_to_tgt_idx(i_elt)+1, src_to_tgt_idx(i_elt+1)
        field2(i_pt) = 0.d0

        elt_n_vtx = cell_vtx_idx(i_elt+1) - cell_vtx_idx(i_elt)

        if (points_weights_idx(i_pt+1) - points_weights_idx(i_pt) /= elt_n_vtx) then
          print *, "Error elt_n_vtx"
          stop
        endif

        do i_vtx = 1, elt_n_vtx
          vtx_id = cell_vtx(cell_vtx_idx(i_elt) + i_vtx)
          field2(i_pt) = field2(i_pt) + vtx_coord(1,vtx_id) * points_weights(points_weights_idx(i_pt) + i_vtx)
        enddo

      enddo
    enddo

    call pdm_pointer_array_part_set(src_send_field2, &
                                    i_part-1,        &
                                    field2)

  enddo
```

### Exchange the interpolated fields from source to target

Now, use the [`Part to Part`](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_doc_pretty/user_manual/comm_graph/ptp.html) object to exchange the interpolated fields from the source mesh to the target cloud.
This `Part to Part` object was built when computing the location and can be accessed from the `Mesh Location` object.

```{code-cell}
%%code_block -p exercise_2 -i 13
  ! Get Part to Part object (it is now owned by the user)
  call pdm_mesh_location_part_to_part_get(mesh_loc, &
                                          0,        &
                                          ptp,      &
                                          PDM_OWNERSHIP_USER)

  ! Initiate exchange of first field
  request1 = -1
  call pdm_part_to_part_iexch(ptp,                                            &
                              PDM_MPI_COMM_KIND_P2P,                          &
                              PDM_STRIDE_CST_INTERLACED,                      &
                              PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2, &
                              1,                                              &
                              null(),                                         &
                              src_send_field1,                                &
                              null(),                                         &
                              tgt_recv_field1,                                &
                              request1)

  ! Initiate exchange of second field
  request2 = -1
  call pdm_part_to_part_iexch(ptp,                                            &
                              PDM_MPI_COMM_KIND_P2P,                          &
                              PDM_STRIDE_CST_INTERLACED,                      &
                              PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2, &
                              1,                                              &
                              null(),                                         &
                              src_send_field2,                                &
                              null(),                                         &
                              tgt_recv_field2,                                &
                              request2)

  ! Finalize both exchanges
  call pdm_part_to_part_iexch_wait(ptp, request1)
  call pdm_part_to_part_iexch_wait(ptp, request2)
```

## Check the interpolated received on the target side

Finally, visualize the interpolated target fields.
(Beware of unlocated points!)

```{code-cell}
%%code_block -p exercise_2 -i 14
  ! Check and visualize the interpolated target fields
  do i = 1, 3
    call pdm_pointer_array_create(tgt_visu_fields(i)%pa, &
                                  tgt_n_part,            &
                                  PDM_TYPE_DOUBLE)
  enddo
  tgt_visu_fields(1)%name = "field1"
  tgt_visu_fields(2)%name = "field1"
  tgt_visu_fields(3)%name = "is_located"


  do i_part = 1, tgt_n_part
    n_located = pdm_mesh_location_n_located_get(mesh_loc, &
                                                0,        &
                                                i_part-1)

    call pdm_mesh_location_located_get(mesh_loc, &
                                       0,        &
                                       i_part-1, &
                                       located)

    n_unlocated = pdm_mesh_location_n_unlocated_get(mesh_loc, &
                                                    0,        &
                                                    i_part-1)

    call pdm_mesh_location_unlocated_get(mesh_loc, &
                                         0,        &
                                         i_part-1, &
                                         unlocated)


    call pdm_mesh_location_point_location_get(mesh_loc, &
                                              0,        &
                                              i_part-1, &
                                              location, &
                                              dist2,    &
                                              projected_coords)

    call pdm_pointer_array_part_get(tgt_vtx_coord,             &
                                    i_part-1,                  &
                                    PDM_STRIDE_CST_INTERLACED, &
                                    3,                         &
                                    vtx_coord)

    call pdm_pointer_array_part_get(tgt_vtx_ln_to_gn, &
                                    i_part-1,         &
                                    vtx_ln_to_gn)

    call pdm_pointer_array_part_get(tgt_recv_field1, &
                                    i_part-1,        &
                                    field1)

    call pdm_pointer_array_part_get(tgt_recv_field2, &
                                    i_part-1,        &
                                    field2)

    allocate(visu_field1(tgt_n_vtx(i_part)), &
             visu_field2(tgt_n_vtx(i_part)), &
             visu_field3(tgt_n_vtx(i_part)))

    do i = 1, n_unlocated
      vtx_id = unlocated(i)
      visu_field1(vtx_id) = -1.d0
      visu_field2(vtx_id) = -1.d0
      visu_field3(vtx_id) =  0.d0
    enddo

    do i = 1, n_located
      vtx_id = located(i)
      error(1) = abs(field1(i) - location(i))
      error(2) = abs(field2(i) - vtx_coord(1,vtx_id))
      if (minval(error) > 1.e-9) then
        print *, "!! error vtx", vtx_ln_to_gn(vtx_id), " :", error
      endif

      visu_field1(vtx_id) = field1(i)
      visu_field2(vtx_id) = field2(i)
      visu_field3(vtx_id) = 1.d0

      call pdm_pointer_array_part_set(tgt_visu_fields(1)%pa, &
                                      i_part-1,              &
                                      visu_field1)

      call pdm_pointer_array_part_set(tgt_visu_fields(2)%pa, &
                                      i_part-1,              &
                                      visu_field2)

      call pdm_pointer_array_part_set(tgt_visu_fields(3)%pa, &
                                      i_part-1,              &
                                      visu_field3)

    enddo
  enddo


  call writer_wrapper(comm,              &
                      "visu",            &
                      "src_mesh",        &
                      src_n_part,        &
                      src_n_vtx,         &
                      src_vtx_coord,     &
                      src_vtx_ln_to_gn,  &
                      src_n_face,        &
                      src_face_edge_idx, &
                      src_face_vtx,      &
                      src_face_ln_to_gn)


  call writer_wrapper(comm,              &
                      "visu",            &
                      "tgt_mesh",        &
                      tgt_n_part,        &
                      tgt_n_vtx,         &
                      tgt_vtx_coord,     &
                      tgt_vtx_ln_to_gn,  &
                      tgt_n_face,        &
                      tgt_face_edge_idx, &
                      tgt_face_vtx,      &
                      tgt_face_ln_to_gn, &
                      vtx_field=tgt_visu_fields)
```


## Free memory
```{code-cell}
%%code_block -p exercise_2 -i 15
  ! Free memory
  call pdm_mesh_location_free(mesh_loc)

  call pdm_part_to_part_free(ptp)

  ! TODO: free everything...
```

## Finalize
```{code-cell}
%%code_block -p exercise_2 -i 16
  call mpi_finalize(ierr)

end program exercise_2
```


## Run the code
Moment of truth!

```{code-cell}
%merge_code_blocks -l fortran -p exercise_2 -n 2 -c
```


## Visualize the results

```{code-cell}
%%visualize -nl -sv
visu/SRC_MESH.case
visu/TGT_MESH.case : is_located
```
