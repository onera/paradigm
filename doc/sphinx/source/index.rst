**********************
ParaDiGM documentation
**********************

**ParaDiGM** (*Parallel Distributed General Mesh*) is a parallel computational geometry library developed at ONERA.
It provides a progressive framework, which consists of a set of low-, mid- and high-level services usable by developers of scientific computing software.

.. The library is written in C but also has Fortran and Python/Numpy APIs.


###############
Getting started
###############


.. container:: toc-cards

  .. card:: Installation
    :link: installation
    :link-type: ref
    :img-top: ../../images/index_installation.png

    A guide to install and configure ParaDiGM



  .. card:: General presentation
    :link: general
    :link-type: ref

    ParaDiGM's key concepts, terminology, conventions and philosophy


.. toctree::
   :caption: Getting Started
   :maxdepth: 1
   :hidden:

   getting_started/installation
   getting_started/general



|

###########
User manual
###########

.. container:: toc-cards

  .. card:: Simple mesh generation
    :link: simple_mesh_gen
    :link-type: ref

    Automatic generation of meshes with simple parametrable shapes



  .. card:: Partitioning
    :link: partitioning
    :link-type: ref

    Graph partitioning, connectivity reconstruction, partition extension and local renumbering



  .. card:: Mesh Structures
    :link: struct
    :link-type: ref

    Mesh data structures



  .. card:: Communication graphs
    :link: comm_graph
    :link-type: ref

    High-level capabilities for building and operating generic communication graphs



  .. card:: Global numbering
    :link: gnum
    :link-type: ref

    Global ID generation and manipulation



  .. card:: Pre-/co-/post-processing
    :link: prepro_algo
    :link-type: ref

    Parrallel geometric and topological algorithms for pre-/co-/post-processing



  .. card:: Parallel I/O
    :link: io
    :link-type: ref

    Parallel read and write



  .. card:: Mesh Transformation
    :link: mesh_transform
    :link-type: ref

    Algorithms performing geometrical or topological transformations on meshes



  .. card:: Miscellaneous utils
    :link: util
    :link-type: ref

    Utilities for basic operations



.. toctree::
   :caption: User manual
   :maxdepth: 1
   :hidden:

   user_manual/simple_mesh_gen/index
   user_manual/partitioning/index
   user_manual/comm_graph/index
   user_manual/gnum/index
   user_manual/prepro_algo/index
   user_manual/mesh_transform/index
   user_manual/io/index
   user_manual/struct/index
   user_manual/util/index
   user_manual/faq

|

################
Developer manual
################


.. container:: toc-cards

  .. card:: Coding rules & guidelines
    :link: coding_rules
    :link-type: ref


.. toctree::
   :caption: Developer manual
   :maxdepth: 1
   :hidden:

   developer_manual/api/index
   developer_manual/coding_rules
   developer_manual/mesh_adaptation/index





.. toctree::
  :maxdepth: 1
  :caption: Appendix
  :hidden:

  changelog
  license


|
