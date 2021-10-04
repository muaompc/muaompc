from setuptools import Extension

ext_modules = [Extension("{prefix}pbm",
    ["src/cython/{prefix}pbm.pyx",
        "src/{prefix}pbmdynmem.c","src/{prefix}pbm.c",
        "src/{prefix}cvpdynmem.c","src/{prefix}cvp.c",
        "src/cjson.c", "src/{prefix}dynmem.c",
        "src/{prefix}pbmsolve.c",
        "src/{prefix}mtxops.c"],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"])]
