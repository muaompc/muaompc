************
Installation
************

Dependencies
============

The following packages are required:

* `Python <http://www.python.org/download/releases/>`_ interpreter.
  This code has been fully tested with Python versions 3.2.3, 3.4.3 and
  Python 2.7.3, 2.7.6.

* `NumPy <http://www.numpy.org>`_, tested with versions 1.6.2, 1.10.4.
  It basically manages the linear
  algebra operations, and some extra features are used.

* `SciPy <http://www.scipy.org>`_, tested with versions 0.10.
  It is used for some tiny features are used.

* `PyParsing <http://pyparsing.wikispaces.com>`_, tested with versions 2.1.4.
  It is used for parsing the problem.

* A C89/C90 compiler.
  To compile the generated code, a C/C++ compiler that supports C89/C90
  or later standards is required.

Optional packages are:

* `Cython <http://cython.org/>`_ to compile the Python interface to the
  generated C-code.

Building and installing
=======================

``muaompc`` installation is made directly from `source code <http://ifatwww.et.uni-magdeburg.de/syst/mpctool/#download>`_.

Install from source in Linux and Mac OS X
-----------------------------------------

Linux and OS X users typically have all required (and most optional) packages already installed.
To install ``muaompc``, switch to the directory where you unpacked ``muaompc``
(you should find a file called ``setup.py`` in that directory) and in a terminal type::

   python setup.py install --user --force

The ``--user`` option indicates that no administrator privileges are required.
The ``--force`` option will overwritte old files from previous installations (if any).
Alternatively, for a global installation, type::

   sudo python setup.py install --force

And that is all about installing the package.
The rest of this document will show you how to use it.


Install from source in Windows Systems
--------------------------------------

For Windows users, we
recommend installing the `Anaconda <https://www.continuum.io/downloads#_windows>`_ platform,
as it contains all of the necessary python packages.

For a full installation of ``muaompc`` do the following:

* Install Anaconda. 

* Open an ``Ananconda Prompt``, switch to the directory where you unpacked ``muaompc`` (the one containing the file ``setup.py``), and type::

   python setup.py install --force

