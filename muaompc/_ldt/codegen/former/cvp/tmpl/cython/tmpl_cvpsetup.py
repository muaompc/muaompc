from setuptools import Extension

ext_modules = [Extension("{prefix}cvp",
    ["src/cython/{prefix}cvp.pyx", 
        "src/{prefix}cvp.c",
        "src/{prefix}dynmem.c",
        "src/{prefix}cvpdynmem.c",
        "src/cjson.c",  "src/{prefix}mtxops.c"],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"]),
    Extension("{prefix}Ccvp",
    ["src/cython/{prefix}Ccvp.pyx",
    ],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"])]

