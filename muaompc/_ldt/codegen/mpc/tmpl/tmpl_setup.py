from distutils.core import setup
from distutils.extension import Extension
from sys import path

from Cython.Build import cythonize

path.append('./src/cython/')
from {prefix}{former}setup import ext_modules as former_ext_modules
from {prefix}{solver}setup import ext_modules as solver_ext_modules
ext_modules = former_ext_modules + solver_ext_modules

setup(name='{prefix}', packages=['{prefix}'],
        ext_package='{prefix}',
        ext_modules=cythonize(ext_modules),
        package_dir={{'{prefix}': 'src'}},
        package_data={{'{prefix}': ['mpc.pickle']}},
        )

