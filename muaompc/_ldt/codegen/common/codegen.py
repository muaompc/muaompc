"""Generate C code common to code generation classes."""

from muaompc._ldt.codegen.codegen import BaseCodeGenerator


class CCodeGenerator(BaseCodeGenerator, object):

    """Generate C code common to code generation classes.

    This class is intended to be used as is, to generate C code,
    or be subclassed to generate interfaces to the C code.
    """

    def __init__(self, prbname, prefix, destdir):
        """Assign prefix to self and set the C code paths.

        :param str prefix: the text to be prepended to the names of the generated
        files and to the internal code global identifiers,
        e.g. names of functions, structures,
        defines, enums, etc. (like a namespace)
        """
        self.prbname = prbname
        self.prefix = prefix
        self.destdir = destdir
        self.path = self._get_paths('codegen.common')

    def generate_code(self, real_t):
        self._make_destination_dir_tree()
        self._copy_static_src(['cjson.c'])
        self._generate_c_code(real_t)

    def _generate_c_code(self, real_t):
        self._generate_arithmetic_header(real_t)
        self._replace_prefix('mtxops.c')
        self._replace_prefix('mtxops.h')
        self._replace_prefix('dynmem.c')
        self._replace_prefix('dynmem.h')
        self._replace_prefix('math.c')
        self._replace_prefix('math.h')

        self._replace_prefix('Makefile.mk')

    def _generate_arithmetic_header(self, real_t):
        fname = 'arithmetic.h'
        tmpl = self._get_tmpl(fname)
        fmt = tmpl.format(real_t=real_t)
        self._write_file(fmt, fname)
