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

# Introduction

Welcome to the **ParaDiGM** library introductory day !

The aim of the day is to give you : 
 - An overview of **ParaDiGM**'s features
 - An understanding of **ParaDiGM**'s advanced parallelism concepts
 - An understanding of how **ParaDiGM** can be called in your software written in C/C++, Fortran or Python.

This is an interactive course, so don't hesitate to interrupt us to ask your questions.

The training will take place in three stages:
- General presentation of **ParaDiGM**
- Description of the abstract concept of inter-partition communication graphs
- Focus on two functionalities through interactive exercises :
    - Exercise 1 : Mesh partitioning
    - Exercise 2 : Localization of a point cloud inside a mesh

# **ParaDiGM** highlights

## Origins of **ParaDiGM**

From 2009 to 2015, various HPC libraries dedicated to different themes were written. **ppart** for parallel graph partitioning, **pario** for parallel I/O as a wrapping to MPI-IO, **CWIPI** for HPC coupling... These libraries shared many principles and structures and were complementary. The idea of merging these libraries into a single one emerged in 2016.  **ParaDiGM** was born!

As the **CWIPI** library was already well used in the academic and industrial communities, it was retained, but its entire algorithmic core was replaced by **ParaDiGM**.

## Objectives

An efficient parallel algorithm takes much longer to write and validate than a sequential one. A simple sequential operation can become very difficult if you want to maintain good *load* and *memory* balancing during all algorithm steps.

Mettre image de la planche des fonctionnalites illustrees ....

**ParaDiGM** aims to offer a set of efficient services to simplify the writing of massively parallel distributed numerical simulation software, from reading data files to writing results. Numerical codes are generally based on a discretization of the study domain which can take the form of an unstructured mesh.
**ParaDiGM** only offers some services for unstructured meshes.<span style="color:red">???</span>

## API

The API has been designed to be as user-friendly and intuitive as possible. All functionalities are in object form. To use a feature, the user **creates** an object and then **provides** the necessary data in the form of arrays and scalars. Once all the data has been provided, the user executes the feature **compute**. As with data, results are **retrieved** in the form of arrays or **ParaDiGM** objects. All concrete results required by the user are obtained directly by dedicated functions, without the need for abstract intermediate structures. User manipulation of abstract **ParaDiGM** objects is reduced to a strict minimum.

The native API is in C and can be used in C/C++ software.

Two other APIs are available. The first one is in Python/Numpy, automatically generated from Cython, and the second one is in Fortran using Fortran's *iso-c-binding* interfaces. These APIs are not simply direct interfaces to C functions. They have been designed to be user-friendly and intuitive in each language.
<span style="color:red">*(OK mais préciser que c'est pas encore au point, et que le retour des utilisateurs est le bienvenu pour converger vers une API optimale.)*</span>
The Python API reinforces the notion of objects, and results are provided in the form of dictionaries. The Fortran API takes Fortran pointers as input/output and not `c_ptr` types with which Fortran developers are unfamiliar. The C API contains pointer arrays. This notion is not defined in Fortran. When giving data or retrieving results in this form, the user must use the Fortran **pdm_pointer_array** intermediate structure.


## License

**ParaDiGM** is licensed under LGPL V3.0.

## Diffusion

diffusion interne + github

## Releases

The latest stable version is 2.4.0, released on November ??, 2023.

The first stable version 1.0.0 was released on March 22, 2017.

Minor versions are released every 3 to 4 months and main versions are released every 2 to 4 years.

Backward compatibility is not guaranteed between two major versions. 

API compatibility is guaranteed between two minor versions, except for new beta functionalities. 

## Man power

Eco-système :  : développeur de code scientifique -> parler du man power (DAAA/DMPE)

## Organisation

GitLab git externe
Licence
GitHub
Documentation Sphinx
Différence ParaDiGM et ParaDiGMa -> mentionner extension

## Application examples

Nicolas avec MoDeTheC : modérniser existant
Sonics : créer un nouveau code
Maia: outil pre-post plus haut niveau
Cedre
CWIPI

## Installation Instructions

### Basic Installation

<span style="color:red">*(vraiment utile de détailler, autant pointer vers doc/README, non?)*</span>

>**cmake .**

>**make**

>**make install**

### CMake general options

> **cmake . -D\<option1_name\>=\<option1_value\> ... -D\<optionn_name\>=\<optionn_value\>** with options :  

 - **CMAKE\_INSTALL\_PREFIX=\<prefix\>** : Installation directory path

 - **PDM_ENABLE_Fortran=<ON | OFF> (default : OFF)** : Enable Fortran interface

 - **PDM_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)** : Enable python interface
      
      If a simple autodetection fails, you can use these options to find Python :

        - Python_ROOT_DIR=<path> 
        - Python_LIBRARY=<path>
        - Python_INCLUDE_DIR=<path>
        - Python_EXECUTABLE=<path>
        
      Refer to FindPython in the CMake documentation for more information.
      shared libraries are necessary for python interface (CWP_ENABLE_SHARED=ON)

 - **PDM_ENABLE_SHARED=<ON | OFF> (default : ON)** : Enable shared libraries

 - **PDM_ENABLE_STATIC=<ON | OFF> (default : ON)** : Enable static libraries

 - **PDM_ENABLE_PARMETIS=<ON | OFF> (default : ON)** : Enable [ParMETIS](https://github.com/KarypisLab/ParMETIS) library (parallel graph partitioning)

      If a simple autodetection fails, you can use PARMETIS_DIR=\<path\> and METIS_DIR=\<path\> options       

      CMake looks for :

        - parmetis.h and metis.h includes
        - parmetis and metis libraries
  
      To link shared libraries, ParMETIS has to be compiled with "-fPIC" option. **ParaDiGM** is compatible with a 32-bit or 64-bit installation.

 - **PDM_ENABLE_PTSCOTCH=<ON | OFF> (default : ON)** : Enable [PTSCOTCH](https://gitlab.inria.fr/scotch/scotch) library (parallel graph partitioning) :
      If a simple autodetection fails, you can use these options to find PTSCOTCH :
        PTSCOTCH_DIR=\<path\>

     CMake looks for :

        - ptscotch.h include file
        - scotch, scotcherr, ptscotch, ptscotcherr libraries
  
     To link shared libraries, PTSCOTCH has to be compiled with "-fPIC" and SCOTCH_PTHREAD_MPI=OFF. *ParaDiGM* is compatible with a 32-bit or 64-bit installation.


 - **PDM_ENABLE_LONG_G_NUM= <ON | OFF> (default : ON)** : Enable long global numbering

        - ON : PDM_g_num_t type is "long int"
        - OFF : PDM_g_num_t type is "int"

 - **PDM_ENABLE_DOC= <ON | OFF> (default : OFF)** : Enable Documentation

     prerequis (sphinx, sphinx fortran, ...)

### CMake compiler options

> **CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...**

or 

> use the following CMake options:

 - **CMAKE_C_COMPILER=\<C compiler\>**
 - **CMAKE_CXX_COMPILER=\<CXX compiler\>**
 - **CMAKE_Fortran_COMPILER=\<Fortran compiler\>**

### CMake MPI options:

 - **MPI_C_COMPILER=\<C mpi wrapper\>**
 - **MPI_CXX_COMPILER=\<CXX mpi wrapper\>**
 - **MPI_Fortran_COMPILER=\<Fortran mpi wrapper\>**

   If a simple autodetection fails, you can use these options to find MPI:

    - **MPI_\<lang\>_LIBRARIES**
    - **MPI_\<lang\>_INCLUDE_PATH**
  
   Refer to FindMPI in the CMake documentation for more informations.

## Concepts and definition

### Mesh

Most computational methods rely on a mesh for the spatial discretization of partial differential equations.
A mesh is composed of entities of different dimensions. The following terminology is used in **ParaDiGM** :
- **cells**: 3D entities such as tetrahedra, pyramids, prisms, hexahedra or arbitrary polyhedra ;
- **faces**: 2D entities such as triangles, quadrangles or arbitrary (simply connected) polygons ;
- **edges**: 1D entities (segment between two points) ;
- **vertices** (*shortened as "**vtx**"*): points defined by their Cartesian coordinates $(x, y, z)$.

A mesh can either be *structured* or *unstructured*.

Structured meshes are typically made of blocks, each one arranged in a regular grid.
Adjacency relations between the mesh entities are therefore implicit : cell $C_{i,j,k}$ is adjacent to cells $C_{i-1,j,k}$, $C_{i+1,j,k}$, $C_{i,j-1,k}$, and so on...

Unstructured meshes, however, require an explicit description of the connectivity between mesh entities.

The entities and connectivities of interest depend on the numerical method.
For example, Finite Element methods typically only require the cell$\to$vertex connectivity (and face$\to$vtx for boundary faces).
On the other hand, cell-centered Finite Volume methods generally require the cell$\to$face and face$\to$vtx connectivities.
Other method, such as node-centered Finite Volume methods may also require the connectivities relative to the edges.

In **ParaDiGM** all connectivities are stored as *1-based*, possibly signed, *flat* arrays.
Because each entity $A$ may be connected to a variable number of entities $B$, an **index** is necessary to access the array $\texttt{connect}$ representing the connectivity $A \to B$.
This index is an array $\texttt{connect\_idx}$ of length $n_A + 1$ which contains the ranges, i.e. the entities $B$ connected to $A_i$ are given by $ \texttt{connect}[j]$, for $j \in \left[ \texttt{connect\_idx}[i], \texttt{connect\_idx}[i+1] \right)$.
The first element in the index array is always zero, and the last element is the length of the connectivity array.
<span style="color:red">*(à adapter en Fortran)*</span>

Let's take a look at a simple example to illustrate this notion:

<img align="left" width="120" height="120" style="margin: 0px 30px 0px 30px;" src="house.svg">

<br/>
Here we have a simple mesh composed of 3 faces and 9 vertices.

The face$\to$vtx connectivity and its index are
\begin{flalign}
  \texttt{face\_vtx\_idx} & =  [0, 4, 12, 15]&&\\\nonumber
  \texttt{face\_vtx}      & =  [{\color{red}2, 3, 6, 5}, {\color{green}1, 2, 5, 6, 3, 4, 8, 7}, {\color{blue}7, 8, 9}]&&
\end{flalign}

<!-- <code>
  face_vtx_idx = [0, 4, 12, 15]
  face_vtx     = [<span style="color:red">2, 3, 6, 5</span>, <span style="color:green;">1, 2, 5, 6, 3, 4, 8, 7</span>, <span style="color:blue">7, 8, 9</span>]
</code> -->

<!-- $$face_vtx = [\underbrace{2, 3, 6, 5}_{face 1}, \quad \underbrace{1, 2, 5, 6, 3, 4, 8, 7}_{face 2}, \quad \underbrace{7, 8, 9}_{face 3} ]$$
 -->

<br/><br/>


*Note that edges can only have two endpoints, so the index for the edge$\to$vtx is useless.*


Vertices are described by the $3 \times n_\mathrm{vtx}$ *flat* array of the Cartesian coordinates. The coordinates are stored in an *interlaced* fashion:
$\left(x_0, y_0, z_0, x_1, y_1, z_1, \ldots \right)$.
<span style="color:red">*(à adapter en Fortran)*</span>

*Note that in **ParaDiGM** coordinates are always assumed to be three-dimensional, even for 2D, planar meshes.*


#### Groups?

<span style="color:red">*TODO : expliquer notion de groupe?*</span>


### Global/absolute numbering

### Parallel distribué MPI

ex : génération de gnum (pas un exercice mais montrer du code) -> exposer graph de communication

### Distributed blocks of an array in absolute numbering 

### 

Part_to_part essentiel pour transmettre les résultats avec les données
-> faire un point intéractif pour bien comprendre le part_to_part : deux nuages de points : dessiner un graphe (nuage 8 et 10 points)
 -- analogie (facteur et grand-parents)
 -- les gens dans la salle sont des partitions qui s'échangent un ballon
 -- exercice Python un peu manuel avec des gnums

## Features overview

Manière dont c'est regroupé dans la documentation
-> dire les points sur lesquels ont va se focaliser pendant le tp : multipart, (bonus : part_extension), localisation, (bonus : part_to_part)

# Exercise 0

cf autre .md
