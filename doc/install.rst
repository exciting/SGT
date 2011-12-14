
.. module:: sgt

=============
What is SGT
=============

Collection of Structure generation tools

.. _installation:
 
==========
Installation
==========

The :mod:`sgt` 
package is a python module. Please folow the installation steps thoroghly to satisfy all dependencies.

Requirements
-------------
:mod:`sgt`  builds on the Atomic Simulation Environment `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`_ 
which in turn needs `numpy <http://numpy.scipy.org/>`_. 
Additionaly some parts need  
`matplotlib <http://matplotlib.sourceforge.net/>`_

Installation Steps
------------------

Numpy and Matplotlib are available as convenient packages which are easy to install. 
For windows and Mac OS X the Enthought Python Distribution 
(`EPD <http://www.enthought.com/products/epd.php>`_) installs python with numpy and matplotlib conviniently. 
Under Linux your package manager will probably have the appropriate packages.

The other packages have to be downloaded manually and installed to the python path. 
In order to do so chose one location to unpack all the packages. e.g. :file:`~/lib/python`
 
To install ASE please download the latest stable release from the
`ASE homepage <https://wiki.fysik.dtu.dk/ase/index.html>`_.
The installation requires to unpack the sources an add the location of the ASE module 
to the python path. see :ref:`pythonpath`
 

.. _pythonpath:

Set Python Path
-----------------
The python path is the environment variable that lists the directories in which the python interpreter searches 
when it is asked to load a module.

For Basch users: add the following to your :file:`.bashrc` (change if unpack locations differ)::

 export PYTHONPATH=${HOME}/lib/python/ase/trunk:${PYTHONPATH}
 export PYTHONPATH=${HOME}/lib/python/structure:${PYTHONPATH}

For csh or tsh users add the following to the :file:`.cshrc` file::

 setenv PYTHONPATH ${HOME}/lib/python/ase/trunk:${PYTHONPATH}
 setenv PYTHONPATH ${HOME}/lib/python/structure:${PYTHONPATH}

Getting started
----------------

The documentation includes examples that should now work as is. see :ref:`toc`.

