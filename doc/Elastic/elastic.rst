
.. module:: sgt


 
Elastic Constants
=================
This tool allows for calculating the elastic constants. It is taking into account the symmetry.

Example
-------- 
This snipped sets ub the distortions and performs the calculations.

 .. literalinclude:: elastic.py
 
The results get saved to a file that can be loaded at a later point for further analyses:

 .. literalinclude:: elasticanalyse.py 
 
Interface
---------



.. autoclass::  ElasticDistortion
   :members:


 