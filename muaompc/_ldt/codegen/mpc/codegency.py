"""
"""
# 19.02.2015 This is just a simple template, maybe of use in the future
import shutil
import os

from muaompc._ldt.codegen.codegen import BaseCodeGenerator as BCG


class CythonCodeGenerator(BCG, object):

    def __init__(self, prefix, path):
        self.prefix = prefix
        self.path = path

    def generate_code(self):
        self._make_destination_dir_tree()
        self._generate_cython_code()

    def _generate_cython_code(self):
        self._replace_prefix('setup.py', subdir='cython')
