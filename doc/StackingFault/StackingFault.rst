.. module:: StackingFault
==============
Stacking Fault
==============

The program has been designed to make proper super-cells for generalized-stacking-fault (GSF) studies. 
At the moment, it can be used for studying the GSF in {111} plane of fcc structure. 
It reads an initial unit-cell, and produces a super-cell containing an arbitrary number of slabs,
 every one shifted along [-211] direction with respect to its neighbors.
 
 .. literalinclude:: StackingFault.py
 
 The resulting structure is:
 
 .. image::  sc.png 

