from distutils.core import setup
from distutils.extension import Extension

ext_modules = [Extension("stocmpccvp",
    ["src/cython/stocmpccvp.pyx", "src/cython/stocmpcCcvp.pyx",
        "src/stocmpccvp.c",
    "src/stocmpcdynmem.c",
    "src/stocmpccvpdynmem.c",
    "src/cjson.c",  "src/stocmpcmtxops.c"],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"]),
    Extension("stocmpcCcvp",
    ["src/cython/stocmpcCcvp.pyx",
    ],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"])]

