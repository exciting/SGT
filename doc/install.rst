
.. module:: sgt

=============
What is SGT
=============

Collection of Structure generation tools

.. _installation:
 
==========
Installation
==========

<<<<<<< HEAD
The :mod:`sgt` 
package is a Python module. Please follow the installation steps thoroughly to satisfy all dependencies.
=======
The SGT 
package is a python module. Please follow the installation steps thoroughly to satisfy all dependencies.
>>>>>>> 4ee1f73e5616062df8938597e3c4685a072425da

Requirements
-------------
SGT  builds on the Atomic Simulation Environment `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`_ 
which in turn needs `numpy <http://numpy.scipy.org/>`_. 
Additionally some parts need  
`matplotlib <http://matplotlib.sourceforge.net/>`_

Installation Steps
------------------

Numpy and Matplotlib are available as convenient packages which are easy to install. 
For windows and Mac OS X the Enough Python Distribution 
<<<<<<< HEAD
(`EPD <http://www.enthought.com/products/epd.php>`_) installs Python with Numpy and Matplotlib conveniently. 
=======
(`EPD <http://www.enthought.com/products/epd.php>`_) installs python with numpy
 and matplotlib conveniently. 
>>>>>>> 4ee1f73e5616062df8938597e3c4685a072425da
Under Linux your package manager will probably have the appropriate packages.

The other packages have to be downloaded manually and installed to the Python path. 
In order to do so chose one location to unpack all the packages. e.g. :file:`~/lib/python`
 
To install ASE please download the latest stable release from the
`ASE home page <https://wiki.fysik.dtu.dk/ase/index.html>`_.
The installation requires to unpack the sources an add the location of the ASE module 
<<<<<<< HEAD
to the Python path. see :ref:`pythonpath`
 
=======
to the python path. see :ref:`pythonpath`

Finally, SGT can be installed. `Download the Package here, <https://github.com/exciting/SGT/zipball/master>`_ 
unpack it, and set the python path :ref:`pythonpath`
>>>>>>> 4ee1f73e5616062df8938597e3c4685a072425da

.. _pythonpath:

Set Python Path
-----------------
The python path is the environment variable that lists the directories in which the python interpreter searches 
when it is asked to load a module.

For Bash users: add the following to your :file:`.bashrc` (change if unpack locations differ)::

 export PYTHONPATH=${HOME}/lib/python/ase/trunk:${PYTHONPATH}
 export PYTHONPATH=${HOME}/lib/python/structure:${PYTHONPATH}

For csh or tsh users add the following to the :file:`.cshrc` file::

 setenv PYTHONPATH ${HOME}/lib/python/ase/trunk:${PYTHONPATH}
 setenv PYTHONPATH ${HOME}/lib/python/structure:${PYTHONPATH}

Getting started
----------------

The documentation includes examples that should now work as is. see :ref:`toc`.

