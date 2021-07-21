#!/usr/bin/python
"""muAO-MPC: microcontroller applications online model predictive control.

Automatic code generation for model predictive control
for embedded systems,
with a particular focus on microcontroller applications.
"""
# this file is based on numpy's setup.py, no copyright infringement intended. 
DOCLINES = __doc__.split("\n")

import os
from setuptools import setup, find_packages

CLASSIFIERS = [
"Development Status :: 4 - Beta",
"Intended Audience :: Science/Research",
"Intended Audience :: Developers",
"License :: OSI Approved :: BSD License",
"Programming Language :: C",
"Programming Language :: Python",
"Programming Language :: Python :: 3",
"Topic :: Scientific/Engineering",
"Operating System :: OS Independent",
]

NAME                = 'muaompc'
MAINTAINER          = "Pablo Zometa"
MAINTAINER_EMAIL    = "muao-mpc@ovgu.de"
DESCRIPTION         = DOCLINES[0]
LONG_DESCRIPTION    = "\n".join(DOCLINES[2:])
URL                 = "http://ifatwww.et.uni-magdeburg.de/syst/mpctool/"
LICENSE             = 'BSD'
AUTHOR              = "Pablo Zometa, et.al."
AUTHOR_EMAIL        = "muao-mpc@ovgu.de"
PLATFORMS           = ["Windows", "Linux", "Mac OS-X"]
MAJOR               = 1
MINOR               = 0
MICRO               = 0
STATUS              = 'beta'
ISRELEASED          = True
VERSION             = '%d.%d.%d-%s' % (MAJOR, MINOR, MICRO, STATUS)

def write_version_py(filename='muaompc/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM MUAO-MPC SETUP.PY
version = '%(version)s'
"""
    with open(filename, 'w') as f:
        f.write(cnt % {'version': VERSION})

def get_ldt_package_data():
    # this is easier done using setuptools or the like, but for compatibility
    # we use only the standard library's distutils
    ldt_base_dirs = ['muaompc/_ldt/codegen/common',
            'muaompc/_ldt/codegen/mpc',
            'muaompc/_ldt/codegen/solver/fgm',
            'muaompc/_ldt/codegen/solver/alm',
            'muaompc/_ldt/codegen/solver/pbm',
            'muaompc/_ldt/codegen/former/cvp',]
    ldt_data_dir = 'tmpl'
    p_dir = dict()
    p_data = dict()

    for ldt_base in ldt_base_dirs:
        ldt_data = []
        for dirname, dirnames, filenames in os.walk(ldt_base + 
                                                    '/' + ldt_data_dir):
            for filename in filenames:
                rel_dirname = dirname.replace(ldt_base + '/', '')
                ldt_data.append(rel_dirname + '/' + filename)

        ldt_base_pkg = ldt_base.replace('/', '.')
        p_dir[ldt_base_pkg] = ldt_base
        p_data[ldt_base_pkg] = ldt_data
    return p_dir, p_data

def get_data_files():
    basedir = 'examples'
    f_data = []
    for dirname, dirnames, filenames in os.walk(basedir):
        for filename in filenames:
            f_data.append(dirname + '/' + filename)
    return basedir, f_data

def setup_package():
    p_dir = dict()
    p_data = dict()

    write_version_py()

    p_ldt_dir, p_ldt_data = get_ldt_package_data()
    f_dir, f_data = get_data_files()
    ps_ldt = ['muaompc',
              'muaompc._ldt', 'muaompc._ldt.parse', 'muaompc._ldt.mpc',
              'muaompc._ldt.codegen', 'muaompc._ldt.codegen.common',
              'muaompc._ldt.codegen.mpc',
              'muaompc._ldt.codegen.former','muaompc._ldt.codegen.solver',
              'muaompc._ldt.codegen.solver.fgm',
              'muaompc._ldt.codegen.solver.alm',
              'muaompc._ldt.codegen.solver.pbm',
              'muaompc._ldt.codegen.former.cvp']
    p_dir.update(p_ldt_dir)
    p_data.update(p_ldt_data)
    packages = ps_ldt

    setup(name=NAME,
          version=VERSION,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          url=URL,
          license=LICENSE,
          classifiers=CLASSIFIERS,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          platforms=PLATFORMS,
          install_requires=['cython', 'numpy', 'scipy', 'pyparse'],
          # test_suite='tests',
          packages = packages,
          package_dir = p_dir,
          package_data = p_data,
          data_files = [(f_dir, f_data)],
          )

if __name__ ==  '__main__':
    setup_package()
