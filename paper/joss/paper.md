---
title: 'ParaDiGM: Parallel Distributed General Mesh'
tags:
  - Parallel computing
  - Computational geometry
  - Mesh
  - MPI

authors:
  - name: Eric Quémerais
    orcid: 0009-0007-2018-6358
    affiliation: "1" 

  - name: Bastien Andrieu
    orcid: 0009-0000-9937-0244
    affiliation: "1" 

  - name: Bruno Maugars
#    orcid: 
    affiliation: "2" 

  - name: Karmijn Hoogveld
    orcid: 0009-0007-8251-1152
    affiliation: "1" 

  - name: Clement Benazet
    orcid: 0009-0007-8251-1152
    affiliation: "2" 

  - name: Berenger Berthoul
#    orcid: 
    affiliation: "2" 

  - name: Nicolas Dellinger
    orcid: 0000-0002-0843-6509
    affiliation: "3" 


affiliations:
 - index: 1
   name: DMPE, ONERA, Université Paris-Saclay, 92320 Châtillon, France
 - index: 2
   name: DAAA, ONERA, Université Paris-Saclay, 92320 Châtillon, France
 - index: 3
   name: DMPE, ONERA, Université de Toulouse, 31000 Toulouse, France
   
date: 27 january, 2026

bibliography: paper.bib

---

# ParaDiGM: Parallel Distributed General Mesh

## Summary 

**ParaDiGM** (**Para**llel **Di**stributed **G**eneral **M**esh) is an open-source C library (LGPL) designed to overcome the technological bottlenecks associated with geometric data handling and mesh manipulation in massively parallel numerical simulations. As High-Performance Computing (HPC) enters the exascale era, managing meshes exceeding 100 billion elements has become a primary hurdle. ParaDiGM surmounts these limitations through a fully distributed architecture capable of scaling beyond 10,000 cores.

Originally branched from the **CWIPI** [@Quemerais2026] coupling library, ParaDiGM provides high-performance geometric services—such as parallel wall distance computation, point cloud location, and dynamic repartitioning—to **ONERA**’s production suite, including **CEDRE** [@Refloch2011], **elsA** [@Cambier2011], **SoNiCS**[@lienhardt2025], and **MoDeTheC**[@Dellinger2024]. By generalizing the **Partitioned** and **Distributed View** concepts, it ensures that every stage of the computation remains memory-balanced. With native APIs for **C/C++**, **Fortran**, and **Python/NumPy**, ParaDiGM acts as a versatile middleware for aerospace research and large-scale computational physics.

## Statement of Need

In the era of exascale computing, numerical simulation faces a critical shift. In **Computational Fluid Dynamics (CFD)**, mesh management can no longer rely on centralized or semi-distributed approaches. Integrating distributed mesh capabilities into existing production solvers often encounters major architectural and performance barriers.

### From Static to Dynamic and Repetitive Geometry
Geometric algorithms, once confined to pre-processing, are now required repeatedly within the solver's main execution loop. Modern studies involve evolving meshes, driven either by **moving bodies** or **Dynamic Mesh Adaptation (AMR)**. Whether computing wall distances for turbulence modeling or performing mesh repartitioning, these operations must deliver turnaround times of the same **order of magnitude as a solver iteration**. Failing to meet this constraint critically penalizes the overall simulation wall-clock time.

### The Challenge of Dynamic Load Balancing
A further bottleneck lies in data distribution: the initial partition provided to geometric algorithms is typically the one optimized for physical computation, which is often **sub-optimal for geometric operations**. Without internal dynamic load-balancing, these operations create computation imbalances and memory bottlenecks, compromising the **execution of simulations** at a very large scale.

### Limitations of Current Solutions
Existing geometric frameworks frequently impose constraints that hinder their adoption:
* **Rigid Data Structures:** Unlike heavy frameworks (Trilinos [@Heroux2005], Arcane [@Grospellier2017]), **ParaDiGM** does not enforce complex object hierarchies. This **non-intrusive** "Progressive Framework" approach allows adoption without a deep refactoring of the host solver.
* **Monolithic Data Models:** Unlike platforms such as *Salome* (MED format) [@Ribes2007], **ParaDiGM** relies on simple **CSR arrays**, avoiding software overhead and memory peaks critical for optimized solvers.

### The ParaDiGM Approach: A Progressive Middleware
**ParaDiGM** is a **Software Development Kit (SDK)** rather than a standalone application. Unlike tools like **Gmsh**[@Geuzaine2009], it is not intended for mesh generation from **CAD** models. Its role begins **as soon as a first mesh is obtained**: it provides solvers with the low-level functions needed to distribute, partition, and manipulate discretized meshes in a parallel fashion.

!
! Ajouter l'image des fonctionnalités
!

Its philosophy is built on three pillars:
1.  **Interoperability and Agnosticism:** Ensuring full compatibility with legacy languages (Fortran, C) and modern environments (Python/NumPy) for a wide variety of numerical methods (Finite Volumes, Finite Elements, SPH).
2.  **Hybrid Partitioning Strategy:** ParaDiGM unifies third-party solutions (e.g., **PT-Scotch**[@Chevalier2008], **ParMetis**[@Karypsis1997]) while supplementing them with native high-performance algorithms like **Space Filling Curves (SFC)** for frequent repartitioning.
3.  **Total Distribution and Geometric Performance:** By ensuring homogeneous load distribution, the framework has enabled the development of highly efficient algorithms, such as **distributed point cloud location**[@Andrieu2026], achieving performance levels compatible with solver iteration frequencies.

## Ongoing Work and Future Challenges

A major challenge for modern numerical simulation is the transition from CPU-based algorithms to implementations optimized for **GPU architectures**. To this end, recent work has been conducted to port ParaDiGM’s most critical components to these accelerators [@Cazalbou2024].

These developments, intended to be industrialized and integrated into a future version, focus on:
* **Hybrid Parallel Octrees:** The construction and traversal of search trees (octrees) are generally the most memory-intensive and time-consuming steps in geometric algorithms.
* **CPU/GPU Optimization:** By offloading these massive data structures to GPUs, ParaDiGM aims to drastically reduce turnaround times for location operations while maintaining the load balance essential for exascale performance.

This evolution will ensure that ParaDiGM remains a state-of-the-art infrastructure, capable of taking full advantage of the next generation of heterogeneous supercomputers.

# Acknowledgements 

The authors would like to thank the following people for their contributions to the development, testing, and dissemination of the ParaDiGM library within the community:

Sébastien Bourasseau, Nicolas Lantos, Niels Guilbert, Alain Hervault, Julien Magnenet, Lucas Manueco, Lionel Matuszewski, Yacine Mezemate, Bertrand Michel, Christina Paulin, Christophe Peyret, and Julien Vanharen (ONERA, France); Jérôme Esbrat, Victor Pacotte, and Bruno Peres (EOLEN); Mickael Philit (Safran, France).

The authors also wish to thank BPI France, the Directorate General for Civil Aviation (DGAC), and the General Scientific Directorate of ONERA for their financial support.

# Author contributions 

The contributions to this software are listed according to the CRediT taxonomy:

- **E. Quémerais**: Conceptualization, Methodology, Software, Validation, Writing – Review & Editing, Project Administration, Funding Acquisition, Supervision  
- **B. Andrieu**: Software, Validation, Writing – Original Draft  
- **C. Benazet**: Software, Validation
- **B. Berthoul**: Software
- **N. Dellinger**: Software
- **K. Hoogveld**: Software, Validation, Writing – Original Draft
- **B. Maugars**: Software, Validation  

# References



