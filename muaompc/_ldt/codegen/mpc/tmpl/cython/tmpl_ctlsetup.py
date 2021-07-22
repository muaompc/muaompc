from setuptools import Extension

ext_modules = [Extension("{prefix}ctl",
    ["src/cython/{prefix}ctl.pyx",
        "src/{prefix}ctldynmem.c","src/{prefix}ctl.c",
        "src/{prefix}{solver}dynmem.c","src/{prefix}{solver}.c",
        "src/{prefix}{former}dynmem.c","src/{prefix}{former}.c",
        "src/cjson.c", "src/{prefix}dynmem.c",
        "src/{prefix}mtxops.c"],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"])]
