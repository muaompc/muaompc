from distutils.extension import Extension

ext_modules = [Extension("{prefix}alm",
    ["src/cython/{prefix}alm.pyx",
        "src/{prefix}almdynmem.c","src/{prefix}alm.c",
        "src/{prefix}fgmdynmem.c","src/{prefix}fgm.c",
        "src/{prefix}cvpdynmem.c","src/{prefix}cvp.c",
        "src/cjson.c", "src/{prefix}dynmem.c",
        "src/{prefix}mtxops.c"],
    include_dirs = ["src/include"], library_dirs = ["src"],
    libraries=["m"])]
