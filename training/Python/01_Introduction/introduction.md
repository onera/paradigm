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

Welcome to the ParaDiGM library introductory day !

The aim of the day is to give you : 
 - An overview of functionality
 - An understanding of ParaDiGM's advanced parallelism concepts
 - An understanding of how ParaDiGM can be called in your software written in C/C++, Fortran or Python.

This is an interactive course, so don't hesitate to interrupt us to ask your questions.

The training will take place in three stages:
- General presentation of ParaDiGM
- Description of the abstract concept of inter-partition communication graphs
- Focus on two functionalities through interactive exercises :
    - Exercice 1 : Mesh partitioning
    - Exercice 2 : Localization of a point cloud inside a mesh

# ParaDiGM highlights

## Origins of ParaDiGM

From 2009 to 2015, various HPC libraries dedicated to different themes were written. ppart for parallel graph partitioning, pario for parallel I/O as an overlay to MPI-IO, CWIPI for HPC coupling... These libraries shared many principles and structures and were complementary. The idea of merging these libraries into a single one emerged in 2016.  ParaDiGM was born!

As the CWIPI library was already well used in the academic and industrial communities, it was retained, but its entire algorithmic core was replaced by ParaDiGM.

## Goal

An efficient parallel algorithm takes much longer to write and validate than a sequential one. A simple sequential operation can become very difficult if you want to maintain good load balancing and memory balancing during all algorithm steps.

Mettre image de la planche des fonctionnalites illustrees ....

ParaDiGM aims to offer a set of efficient services to simplify the writing of massively parallel distributed simulation software, from reading data files to writing results. Numerical codes are generally based on a discretization of the study domain which can take the form of an unstructured mesh. ParaDiGM offers some services for unstructured meshes. Most ParaDiGM services are based on unstructured mesh data.       

## Licence

ParaDiGM is licensed under LGPL V3.0

## Diffusion

diffusion interne + github

## Releases

The latest stable version is 2.4.0, released on November ??, 2023.

The first stable version 1.0.0 was released on March 22, 2017.

Minor versions are released every 3 to 4 months and main versions are released every 2 to 4 years

Backward compatibility is not guaranteed between two major versions. 

API compatibility is guaranteed between two minor versions, except for new beta functionalities. 

## Man power

Eco-système : Cible de la bibliothèque : développeur de code scientifique -> parler du man power (DAAA/DMPE)
Nicolas avec ModeTech : modérniser existant
Sonics : créer un nouveau code
Maia: outil pre-post plus haut niveau

## Organisation

GitLab
Licence
GitHub
Documentation Sphinx
Différence ParaDiGM et ParaDiGMa -> mentionner extension

## Installation Instructions

### Basic Installation

>**cmake .**

>**make**

>**make install**

### CMake general options

> **cmake . -D\<option1_name\>=\<option1_value\> ... -D\<optionn_name\>=\<optionn_value\>** with options :  

 - **CMAKE\_INSTALL\_PREFIX=\<prefix\>** : Installation directory path

 - **PDM_ENABLE_Fortran=<ON | OFF> (default : OFF)** : Enable fortran interface

 - **PDM_ENABLE_PYTHON_BINDINGS=<ON | OFF> (default : OFF)** : Enable python interface
      
      If a simple autodetection fails, you can use these options to find Python :

        - Python_ROOT_DIR=<path> 
        - Python_LIBRARY=<path>
        - Python_INCLUDE_DIR=<path>
        - Python_EXECUTABLE=<path>
        
      Refere to FindPython in the CMake documentation for more informations.
      shared libraries are necessary for python interface (CWP_ENABLE_SHARED=ON)

 - **PDM_ENABLE_SHARED=<ON | OFF> (default : ON)** : Enable shared libraries

 - **PDM_ENABLE_STATIC=<ON | OFF> (default : ON)** : Enable static libraries

 - **PDM_ENABLE_PARMETIS=<ON | OFF> (default : ON)** : Enable [ParMETIS](https://github.com/KarypisLab/ParMETIS) library (parallel graph partition)

      If a simple autodetection fails, you can use PARMETIS_DIR=\<path\> and METIS_DIR=\<path\> options       

      To link shared libraries, ParMETIS has to be compiled with "-fPIC" option. ParaDiGM is compatible with a 32bit or 64bit version.

      CMake looks for :

        - parmetis.h and metis.h includes
        - parmetis and metis libraries

 - **PDM_ENABLE_PTSCOTCH=<ON | OFF> (default : ON)** : Enable [PTSCOTCH](https://gitlab.inria.fr/scotch/scotch) library (parallel graph partition) :
      If a simple autodetection fails, you can use these options to find PTSCOTCH :
        PTSCOTCH_DIR=<path>

     To link shared libraries, PTSCOTCH has to be compiled with "-fPIC" and SCOTCH_PTHREAD_MPI=OFF. ParaDiGM is compatible with a 32bit or 64bit version.

     CMake looks for :

        - ptscotch.h include file
        - scotch, scotcherr, ptscotch, ptscotcherr libraries

 - **PDM_ENABLE_LONG_G_NUM= <ON | OFF> (default : ON)** : Enable long global number

        - ON : PDM_g_num_t type is "long int"
        - OFF : PDM_g_num_t type is "int"

 - **PDM_ENABLE_DOC= <ON | OFF> (default : OFF)** : Enable Docuementation

     prerequis (sphinx, sphinx fortran, ...)

### CMake compiler options

> **CC=<C compiler> CXX=<CXX compiler> FC=<Fortran compiler> cmake ...**

or use the following cmake options
    CMAKE_C_COMPILER=<C compiler>
    CMAKE_CXX_COMPILER=<CXX compiler>
    CMAKE_Fortran_COMPILER=<Fortran compiler>

### CMake MPI options

    - MPI_C_COMPILER=<C mpi wrapper>
    - MPI_CXX_COMPILER=<CXX mpi wrapper>
    - MPI_Fortran_COMPILER=<Fortran mpi wrapper>

If a simple autodetection fails, you can use these options to find MPI :

    - MPI_<lang>_LIBRARIES
    - MPI_<lang>_INCLUDE_PATH

Refere to FindMPI in the CMake documentation for more informations.

## Technique

Reprendre les slides de Julien Coulet

### Structures de données

Qu'est-ce qu'un maillage non structuré? Comment on le représente? Présenter la logique des tableaux avec index.

### Parallel distribué MPI

### Numérotation absolu

ex : génération de gnum (pas un exercice mais montrer du code) -> exposer graph de communication

### Bloc/Partition

Part_to_part essentiel pour transmettre les résultats avec les données
-> faire un point intéractif pour bien comprendre le part_to_part : deux nuages de points : dessiner un graphe (nuage 8 et 10 points)
 -- analogie (facteur et grand-parents)
 -- les gens dans la salle sont des partitions qui s'échangent un ballon
 -- exercice Python un peu manuel avec des gnums

## Eventail des fonctionalités

Manière dont c'est regroupé dans la documentation
-> dire les points sur lesquels ont va se focaliser pendant le tp : multipart, (bonus : part_extension), localisation, (bonus : part_to_part)

# Exercice 0

cf autre .md
