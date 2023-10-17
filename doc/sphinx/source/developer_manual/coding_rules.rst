.. _coding_rules:

#########################
Coding rules & guidelines
#########################

Structure
---------

The API of each algorithm is in 5 steps :
- create : initialize the structure associated to the feature
- set : provide the data necessary for the algorithm to the structure
- compute : run the algorithm
- get : retrieve the output data
- free : free the memory of the feature structure

When data needs to be exchanged at the end of the algorithm, a getter has to be provided on the ``part_to_part``
associated to the feature. Data related to the exchange can be getted directly in the ``part_to_part``.

In any case, a feature should take into account a mesh with several partition (``n_part > 1``). Adk J. Coulet and B. Maugars, if
you want to know whether you should handle several zones.

Indentation
-----------

C, Fortran and Python sources use 2 **spaces** for indentation.


Naming
------

Function, variable and C struct names should follow the ``snake_case`` convention.
Python class names follow the ``CamelCase`` convention.

.. Function names should respect the following templates:
.. ...

Signature
---------

In a functions signature, always put the *something*_idx argument before the *something*.

Documentation
-------------

Add with an array is 0 or 1-based in the function signature documentation. If it is 0-based, refer to it as ID rather than number.

