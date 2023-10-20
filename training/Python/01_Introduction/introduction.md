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

# Inntroduction

=> formation intéractive, posezzzzzz vos questions !!

## Historique

Développements similaires dans plusieurs codes (CEDRE, SPACE, CWIPI) : partitionnement, IO, localisation (géométrie)
Parallèle complexe intérêt de multualiser => détails pypart, pario ...(agrégation pour centraliser)

## Sucess story

Objectif paralléliser CEDRE et les IO (date?) : pas seulement solveur aussi pre-post processing (gains qualitatifs?)
Distance aux parois (gains quantitatifs?)
Besoin bibliothèque géométrique efficace

Eco-système : Cible de la bibliothèque : développeur de code scientifique
Nicolas avec ModeTech : modérniser existant
Sonics : créer un nouveau code
Maia: outil pre-post plus haut niveau

## Organisation

GitLab
Licence
GitHub
Documentation Sphinx
Différence ParaDiGM et ParaDiGMa -> mentionner extension

## Installation

-> pas leur faire faire mais penser à les ajouter sur le GitLab : les faire aller sur GitLab où est-ce que c'est décrit comment faire
-> démo en direct ?
-> build.sh le montrer voir faire utiliser ?
-> mentionner les dépendances : MPI, Scotch, Metis (32 ou 64 bit peu importe)

Installation depuis GitLab
Comment générer la documentation ?

=> pas s'éterniser, c'est standard !

## Technique

Reprendre les slides de Julien Coulet

### Parallel distribué MPI

### Numérotation absolu

ex : génération de gnum (pas un exercice mais montrer du code)

### Bloc/Partition

## Eventail des fonctionalités

Manière dont c'est regroupé dans la documentation
-> dire les points sur lesquels ont va se focaliser pendant le tp : multipart, (bonus : part_extension), localisation, (bonus : part_to_part)

# Exercice 0

cf autre .md
