"""Generate C code to form a condensed problem with vector parameters."""

import numpy as np

from muaompc._ldt.codegen.codegen import BaseCodeGenerator as BCG
from muaompc._ldt.codegen.codegen import BaseDataGenerator as BDG


PMETRIC_TERMS = ['val', 'fac0', 'aux']


class CCodeGenerator(BCG, object):

    """Generate C code to form a condensed problem with vector parameters.

    This class is intended to be used as is, to generate C code,
    or be subclassed to generate interfaces to the C code.
    """

    def __init__(self, cvp, ccg):
        """Assign the arguments to self, and set the C code paths.

        This depends on common C-code already being generated.

        :param cvp: an optimization problem representation in the form
        of a condensed qp.
        :type cvp: muaompc._ldt.codegen.former.cvp.cvp.CondensedVectorialParameters
        :param ccg: contains the prefix and the destination directory of the
        generated common C-code.
        :type ccg: muaompc._ldt.codegen.common.codegen.CCodeGenerator
        """
        self.cvp = cvp
        self.prefix = ccg.prefix
        self.destdir = ccg.destdir
        self.prbname = ccg.prbname
        self.path = self._get_paths('codegen.former.cvp')

    def generate_code(self):
        self._make_destination_dir_tree()
        self._generate_c_code()

    def _generate_c_code(self):
        self._generate_c_header('cvp')
        self._generate_c_body('cvp')
        self._generate_c_dynamic_allocation('cvpdynmem')

    def _generate_c_header(self, fname):
        fname += '.h'
        def_flags = self._get_def_flags()
        par_inds = self._get_enum_param_c_tmpl(self.cvp.smb.par)
        prb_terms = self._get_prb_terms_c_tmpl()
        prb_terms += self._get_prb_socc_terms_c_tmpl()
        parameters_terms = self._get_parameters_terms_c_tmpl(self.cvp.par)
        tmp = dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   def_flags=def_flags,
                   par_inds=par_inds,
                   prb_terms=prb_terms,
                   parameters_terms=parameters_terms)
        self._replace_dict(tmp, self.prefix, fname)

    def _generate_c_body(self, fname):
        fname += '.c'
        copy_parameters = self._get_copy_parameters_c_tmpl(self.cvp.par)
        tmp = dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   copy_parameters=copy_parameters)
        self._replace_dict(tmp, self.prefix, fname)

    def _generate_c_dynamic_allocation(self, fname):
        fname_c = fname+'.c'
        init_data = self._get_init_data_c_tmpl()
        init_struct_prb = self._get_init_prb_c_tmpl()
        init_struct_prb += self._get_init_prb_socc_c_tmpl()
        alloc_parameters = self._get_alloc_parameters_c_tmpl(self.cvp.par)
        tmp = dict(prefix=self.prefix,
                   PREFIX=self.prefix.upper(),
                   alloc_parameters=alloc_parameters,
                   init_struct_prb=init_struct_prb,
                   init_data=init_data)
        self._replace_dict(tmp, self.prefix, fname_c)
        # Header file
        fname_h = fname+'.h'
        self._replace_prefix(fname_h)

    def _eval_func_over_cvp(self, func):
        fmt = ''
        for name in self.cvp.par:
            fmt += func('par', name)
        for name, par_fac in self.cvp.pmetric.items():
            fmt += func('pmetric', name, par_fac)
        return fmt

    def _get_def_flags(self):
        flags = self.cvp._get_flags()
        fmt = ''
        for name, flag in flags.items():
            fmt += '#define %s_CVP_PRB_%s %s\n' % (self.prefix.upper(), name.upper(), flag)
        return fmt

    def _get_enum_param_c_tmpl(self, smb_par):
        ntabs = 1
        fmt = ''
        for k, par in enumerate(smb_par):
            fmt += self._get_term_enum_c_tmpl(par.name, k, ntabs)
        return fmt

    def _get_init_data_c_tmpl(self):
        def _init_data_term(key, name, par_fac=None):
            ntab = 1
            tab = '    '
            tabs = ntab * tab
            idx = (self.prefix+'_'+name).upper()
            term = 'cvp->%s[%s]' % (key, idx)
            check_ret = tabs+'if (%s_DYNMEM_OK != ret) {return ret;}\n' % (
                self.prefix.upper())
            fmt = ''
            if key == 'par':
                fmt += tabs+'ret = %s_get_json_term(%s, data, "%s", "%s");\n' % (
                    self.prefix, term, key, name)
                fmt += check_ret

            return fmt
        return self._eval_func_over_cvp(_init_data_term)

    def _get_prb_terms_c_tmpl(self):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        fmt = ''
        for term_name in self.cvp._get_prb_terms():
            term = 'struct %s_term *%s' % (self.prefix, term_name)
            fmt += tabs+'%s;\n' % (term)
        return fmt

    def _get_prb_socc_terms_c_tmpl(self):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        fmt = ''
        if len(self.cvp.soccs.struct['uppers']) != 0:
            fmt += tabs+'uint32_t *socc_num;\n'
            fmt += tabs+'struct %s_cvp_prb_socc **socc;\n' % (self.prefix)
        return fmt

    def _get_parameters_terms_c_tmpl(self, terms):
        fmt = ''
        tab = '    '
        for name in terms:
            fmt += 1*tab + 'real_t *' + name + ';\n'
        return fmt

    def _get_copy_parameters_c_tmpl(self, terms):
        fmt = ''
        tab = '    '
        for name in terms:
            idx = (self.prefix+'_'+name).upper()
            assign = 'ctl_par[%s].data = parameters->%s' % (idx, name)
            fmt += 1*tab + assign + ';\n'
        return fmt

    def _get_init_prb_c_tmpl(self):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        fmt = ''
        for term_name in self.cvp._get_prb_terms():
            idx = (self.prefix+'_'+term_name).upper()
            term_lhs = 'cvp->prb->%s' % (term_name)
            term_rhs = 'cvp->pmetric[%s]' % (idx)
            fmt += tabs+'%s = %s->val;\n' % (term_lhs, term_rhs)

        return fmt

    def _get_prb_socc_terms_c_tmpl(self):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        fmt = ''
        if len(self.cvp.soccs.struct['uppers']) != 0:
            fmt += tabs+'uint32_t *socc_num;\n'
            fmt += tabs+'struct %s_cvp_prb_socc **socc;\n' % (self.prefix)
        return fmt

    def _get_init_prb_socc_c_tmpl(self):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        fmt = ''
        if len(self.cvp.soccs.struct['uppers']) != 0:
            fmt += tabs+'cvp->prb->socc_num = cvp->socc_num;'
        return fmt

    def _get_alloc_parameters_c_tmpl(self, terms):
        fmt = ''
        tab = '    '
        for name in terms:
            idx = (self.prefix+'_'+name).upper()
            assign = 'parameters->%s = cvp->par[%s]->data' % (name, idx)
            fmt += 1*tab + assign + ';\n'
        return fmt


class CVPDataGenerator(BDG, object):
    """Generate json data for a condensed vectorial parameters formulation.
    """

    def __init__(self, cvpcg):
        """Assign the arguments to self.

        :param cvp: an instance of condensed vectorial parameters formulation
        :type cvp: muaompc._ldt.codegen.former.cvp.cvp.CondensedVectorialParameters

        """
        self.cvp = cvpcg.cvp
        self.smb = cvpcg.cvp.smb
        self.path = cvpcg.path
        self.prefix = cvpcg.prefix

    def generate_data(self, num, dataname):
        """Write a json file with the CVP values for the given data.

        Given numeric values for the symbols found in the general
        problem description, create and write into json format
        the matrices for the condensed vectorial parameters formulation.

        :param dict num: numeric values for the symbols (dict keys)
        in the optimization problem
        :param str prefix: the text to prepend to the name of the generated
        data file. The name of the file is then 'prefix_cvp.json'
        :return: the data written to the json file
        :rtype: dict
        """

        data = self._get_data(num)
        pfname = self._get_full_fname_prefix(dataname) + '.json'
        self._write_json_data(pfname, data)

    def _get_data(self, num):
        cprb = self.cvp.form_condensed_prb(num)
        par = self._get_par(num)
        pmetric = self._get_pmetric(cprb['pmetric'])
        optvar_len = self._get_optvar_len(num)
        lagmul_len = cprb['lagmul_len']
        socc = self._get_socc(cprb['socc'])
        return dict(par=par, pmetric=pmetric,
                    socc=socc,
                    optvar=optvar_len, lagmul=lagmul_len)

    def _get_pmetric(self, pmetric_):
        pmetric_t = dict()
        for name, pmetric in pmetric_.items():
            pmetric_t[name] = self._get_pmetric_terms(pmetric)
        #FIXME: does this work also for the list of pmetric in iec?
        return pmetric_t

    def _get_pmetric_terms(self, pmetric):
        # terms affine wrt the parameters
        fac_terms = []
        par_names = []
        par_id = []
        for name, fac in pmetric.items():
            if name == 'fac0':
                fac0_term = self._get_constant_mtx_term(fac)
                continue

            is_not_zero_mtx = not np.allclose(fac, np.zeros(fac.shape))
            if is_not_zero_mtx:
                fac_terms.append(self._get_constant_mtx_term(fac))
                par_names.append(name)
                for id_, par in enumerate(self.cvp.smb.par):
                    if name == par.name:
                        par_id.append(id_)
        # column vector val and fac0 have the same number of rows as any factor
        fac_num = len(fac_terms)
        val_term = fac0_term
        aux_term = self._get_constant_mtx_term(np.zeros(
            (fac0_term['rows'], fac0_term['cols'])))

        return dict(fac0=fac0_term, fac=fac_terms, aux=aux_term, par=par_names,
                    par_id=par_id, val=val_term, fac_num=fac_num)

    def _get_socc(self, socc):
        socc_pmetric = [self._get_pmetric(pmetric)
                for pmetric in socc['pmetric']]
        socc_num = len(socc_pmetric)
        dsocc = dict(socc_num=socc_num,
                pmetric=socc_pmetric)
        return dsocc

class CDataGenerator(CCodeGenerator, CVPDataGenerator, object):

    def __init__(self, cvpcg):
        """Assign the arguments to self.

        :param ccg: an instance of C code generator of a
        condensed vectorial parameters formulation
        :type cvp: muaompc._ldt.codegen.former.cvp.codegen.CCodeGenerator

        """
        self.cvp = cvpcg.cvp
        self.smb = cvpcg.cvp.smb
        self.path = cvpcg.path
        self.prefix = cvpcg.prefix
        self.dict_ = None  # Assignment made during data generation

    def generate_data(self, num, dataname):
        self.dict_ = dict(prefix=self.prefix,
                          PREFIX=self.prefix.upper(),
                          dataname=dataname,
                          DATANAME=dataname.upper())
        data = self._get_data(num)
        self._write_static_data(data, dataname)

    def _write_static_data(self, data, dataname):
        self._generate_h_struct(data, dataname)
        self._generate_c_struct(data, dataname)

    def _generate_h_struct(self, data, dataname):
        fname = 'cvpdata'
        fname += '.h'
        prefix = self._get_full_fname_prefix(dataname)
        tmp = self.dict_
        self._replace_dict(tmp, prefix, fname)

    def _generate_c_struct(self, data, dataname):
        fname = 'cvpdata'
        fname += '.c'
        alloc_terms = self._get_alloc_terms_c_tmpl(data, dataname)
        alloc_socc_terms = self._get_alloc_socc_terms_c_tmpl(data, dataname)
        init_terms = self._get_init_terms_c_tmpl(data, dataname)
        init_socc_terms = self._get_init_socc_terms_c_tmpl(data, dataname)
        init_prb = self._get_init_prb_c_tmpl(dataname)
        init_prb += self._get_init_prb_socc_c_tmpl(data, dataname)
        init_parameters = self._get_init_parameters_c_tmpl(self.cvp.par)
        prefix = self._get_full_fname_prefix(dataname)
        pname = self.prefix+dataname
        tmp = dict(alloc_terms=alloc_terms,
                   alloc_socc_terms=alloc_socc_terms,
                   init_terms=init_terms,
                   init_socc_terms=init_socc_terms,
                   init_prb=init_prb,
                   init_parameters=init_parameters,
                   pname=pname)
        tmp.update(self.dict_)
        self._replace_dict(tmp, prefix, fname)

    def _assign_term_data(self, pname, key, name, data, tabs, par_aux=None,
                          is_par_name=False):
        if par_aux is None:
            term = data[key][name]
        else:
            if is_par_name:
                term = data[key][name]['fac'][par_aux]
                par_name = data[key][name]['par'][par_aux]
                pname += '_'+par_name
            else:
                term = data[key][name][par_aux]
                pname += '_'+par_aux

        fmt = self._assign_term_Cdata(pname, term, tabs)
        return fmt

    def _assign_term_Cdata(self, pname, term, tabs):
        fmt = ''
        fmt += tabs+'real_t %s_data[] = %s;\n' % (
                pname, self._np2Carray(term['data']))
        fmt += tabs+'struct %s_term %s = {' % (self.prefix, pname)
        fmt += '%s, ' % (term['rows'])
        fmt += '%s, ' % (term['cols'])
        fmt += '%s_data};\n' % (pname)
        return fmt

    def _get_alloc_terms_c_tmpl(self, data, dataname):
        def _alloc_terms_term(key, name, par_fac=None, data=data,
                              dataname=dataname):
            ntab = 0
            tab = '    '
            tabs = ntab * tab
            pname = self.prefix+dataname+'_'+name
            fmt = ''  # alloc_terms
            if par_fac is None:
                fmt += self._assign_term_data(pname, key, name, data, tabs)
            elif data[key][name]['fac_num'] == 0:
                # constant term, only val is and fac_num are important
                if name in self.cvp._get_prb_terms():
                    fmt += self._assign_term_data(pname, key, name, data, tabs,
                                              'val')
                else:
                    fmt += tabs+'struct %s_term %s_val = {0, 0, NULL};\n' % (
                        self.prefix, pname)

                fmt += tabs+'struct %s_term %s_aux = {0, 0, NULL};\n' % (
                        self.prefix, pname)
                fmt += tabs+'struct %s_term %s_fac0 = {0, 0, NULL};\n' % (
                        self.prefix, pname)
                fmt += tabs+'struct %s_term **%s_fac = NULL;\n' % (
                        self.prefix, pname)
                fmt += tabs+'struct %s_term **%s_par = NULL;\n' % (
                        self.prefix, pname)
                fmt += tabs+'struct %s_pmetric %s;\n' % (self.prefix, pname)
                fmt += tabs+'uint32_t %s_fac_num = %d;\n' % (
                    pname, data[key][name]['fac_num'])
                fmt += tabs+'uint32_t *%s_par_id = NULL;\n' % (pname)
            else:
                n = data[key][name]['fac_num']

                for term_name in PMETRIC_TERMS:
                    fmt += self._assign_term_data(pname, key, name, data, tabs,
                                                  term_name)
                fmt += tabs+'struct %s_term *%s_fac[%d];\n' % (
                    self.prefix, pname, n)
                fmt += tabs+'struct %s_term *%s_par[%d];\n' % (
                    self.prefix, pname, n)
                fmt += tabs+'struct %s_pmetric %s;\n' % (self.prefix, pname)
                for j, par_name in enumerate(data[key][name]['par']):
                    fmt += self._assign_term_data(pname, key, name, data,
                                                  tabs, j, True)

                fmt += tabs+'uint32_t %s_fac_num = %d;\n' % (
                    pname, data[key][name]['fac_num'])
                fmt += tabs+'uint32_t %s_par_id[] = %s;\n' % (
                    pname, self._np2Carray(data[key][name]['par_id']))
            return fmt

        fmt = self._eval_func_over_cvp(_alloc_terms_term)
        fmt += '    '+'struct %s_cvp_prb %s_prb;\n' % (self.prefix,
                                                       self.prefix+dataname)
        return fmt

    def _get_alloc_socc_terms_c_tmpl(self, data, dataname):
        #FIXME: this function is very similar to  _get_alloc_terms_c_tmpl
        ntab = 0
        tab = '    '
        tabs = ntab * tab
        pname = self.prefix+dataname
        fmt = ''  # alloc_terms
        fmt += tabs+'uint32_t %s_socc_num = %d;\n' % (
                    pname, data['socc']['socc_num'])
        for k, pmetric in enumerate(data['socc']['pmetric']):
            for name in pmetric:
                tname = pname+'_'+name+str(k)
                if pmetric[name]['fac_num'] == 0:
                    if name in ['Wm', 'wvT']:  # terms assumed constant
                        pmtname = tname+'_val'
                        fmt += self._assign_term_Cdata(pmtname, pmetric[name]['val'], tabs)
                    else:
                        fmt += tabs+'struct %s_term %s_val = {0, 0, NULL};\n' % (
                            self.prefix, tname)

                    fmt += tabs+'struct %s_term %s_aux = {0, 0, NULL};\n' % (
                            self.prefix, tname)
                    fmt += tabs+'struct %s_term %s_fac0 = {0, 0, NULL};\n' % (
                            self.prefix, tname)
                    fmt += tabs+'struct %s_term **%s_fac = NULL;\n' % (
                            self.prefix, tname)
                    fmt += tabs+'struct %s_term **%s_par = NULL;\n' % (
                            self.prefix, tname)
                    fmt += tabs+'struct %s_pmetric %s;\n' % (self.prefix, tname)
                    fmt += tabs+'uint32_t %s_fac_num = %d;\n' % (
                        tname, pmetric[name]['fac_num'])
                    fmt += tabs+'uint32_t *%s_par_id = NULL;\n' % (tname)

                else:
                    par_num = len(pmetric[name]['par'])
                    for term_name in PMETRIC_TERMS:
                        pmtname = tname+'_'+term_name
                        fmt += self._assign_term_Cdata(pmtname, pmetric[name][term_name], tabs)
                    fmt += tabs+'struct %s_term *%s_fac[%d];\n' % (
                        self.prefix, tname, par_num)
                    fmt += tabs+'struct %s_term *%s_par[%d];\n' % (
                        self.prefix, tname, par_num)
                    fmt += tabs+'struct %s_pmetric %s;\n' % (self.prefix, tname)
                    for j, par_name in enumerate(pmetric[name]['par']):
                        term = pmetric[name]['fac'][j]
                        parname = tname+'_'+par_name
                        fmt += self._assign_term_Cdata(parname, term, tabs)

                    fmt += tabs+'uint32_t %s_fac_num = %d;\n' % (
                        tname, pmetric[name]['fac_num'])
                    fmt += tabs+'uint32_t %s_par_id[] = %s;\n' % (
                        tname, self._np2Carray(pmetric[name]['par_id']))

        for k in range(data['socc']['socc_num']):
            fmt += tabs+'struct %s_cvp_prb_socc %s_prb_socc%d;\n' % (self.prefix, pname, k)
            fmt += tabs+'struct %s_cvp_socc %s_socc%d;\n' % (self.prefix, pname, k)

        if data['socc']['socc_num'] == 0:
            fmt += tabs+'struct %s_cvp_socc **%s_socc = NULL;\n' % (self.prefix, pname)
            fmt += tabs+'struct %s_cvp_prb_socc *%s_prb_socc = NULL;\n' % (self.prefix, pname)
        else:
            fmt += tabs+'struct %s_cvp_socc *%s_socc[%d];\n' % (self.prefix, pname, data['socc']['socc_num'])
            fmt += tabs+'struct %s_cvp_prb_socc *%s_prb_socc[%d];\n' % (self.prefix, pname, data['socc']['socc_num'])

        return fmt

    def _get_init_terms_c_tmpl(self, data, dataname):
        def _init_terms_term(key, name, par_fac=None, data=data,
                             dataname=dataname):
            ntab = 1
            tab = '    '
            tabs = ntab * tab
            pname = self.prefix+dataname+'_'+name
            idx = (self.prefix+'_'+name).upper()
            term = 'cvp->%s[%s]' % (key, idx)
            fmt = ''  # init_terms
            if par_fac is None:
                fmt += tabs+'%s = &%s;\n' % (term, pname)
            else:
                fmt += tabs+'%s = &%s;\n' % (term, pname)
                fmt += tabs+'%s->fac_num = &%s_fac_num;\n' % (term, pname)
                fmt += tabs+'%s->par_id = %s_par_id;\n' % (term, pname)
                fmt += tabs+'%s->val = &%s_val;\n' % (term, pname)
                fmt += tabs+'%s->aux = &%s_aux;\n' % (term, pname)
                fmt += tabs+'%s->fac0 = &%s_fac0;\n' % (term, pname)
                fmt += tabs+'%s->par = %s_par;\n' % (term, pname)
                fmt += tabs+'%s->fac = %s_fac;\n' % (term, pname)

                for j, par_name in enumerate(data[key][name]['par']):
                    pname_j = pname+'_'+par_name
                    fmt += tabs+'%s->fac[%d] = &%s;\n' % (term, j, pname_j)
                    par_term = self.prefix+dataname+'_'+par_name
                    fmt += tabs+'%s->par[%d] = &%s;\n' % (term, j, par_term)
            return fmt

        return self._eval_func_over_cvp(_init_terms_term)

    def _get_init_socc_terms_c_tmpl(self, data, dataname):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        pname = self.prefix+dataname
        fmt = ''  # init_terms
        fmt += tabs+'cvp->socc_num = &%s_socc_num;\n' % (pname)
        for k in range(data['socc']['socc_num']):
            fmt += tabs+'cvp->socc[%d] = &%s_socc%d;\n' % (k, pname, k)

        for k, pmetric in enumerate(data['socc']['pmetric']):
            for name in pmetric:
                idx = (self.prefix+'_SOCC_'+name).upper()
                term = 'cvp->socc[%d]->pmetric[%s]' % (k, idx)
                tname = pname+'_'+name+str(k)
                fmt += tabs+'%s = &%s;\n' % (term, tname)
                fmt += tabs+'%s->fac_num = &%s_fac_num;\n' % (term, tname)
                fmt += tabs+'%s->par_id = %s_par_id;\n' % (term, tname)
                for term_name in PMETRIC_TERMS:
                    pmtname = tname+'_'+term_name
                    fmt += tabs+'%s->%s = &%s;\n' % (term, term_name, pmtname)
                fmt += tabs+'%s->par = %s_par;\n' % (term, tname)
                fmt += tabs+'%s->fac = %s_fac;\n' % (term, tname)

                for j, par_name in enumerate(pmetric[name]['par']):
                    parname = tname+'_'+par_name
                    fmt += tabs+'%s->fac[%d] = &%s;\n' % (term, j, parname)
                    par_term = pname+'_'+par_name
                    fmt += tabs+'%s->par[%d] = &%s;\n' % (term, j, par_term)

        return fmt

    def _get_init_prb_c_tmpl(self, dataname):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        fmt = ''
        for term_name in self.cvp._get_prb_terms():
            pname = self.prefix+dataname+'_'+term_name
            term = 'cvp->prb->%s' % (term_name)
            fmt += tabs+'%s = &%s_val;\n' % (term, pname)
        return fmt

    def _get_init_prb_socc_c_tmpl(self, data, dataname):
        ntab = 1
        tab = '    '
        tabs = ntab * tab
        pname = self.prefix+dataname
        fmt = ''
        if data['socc']['socc_num'] != 0:
            fmt += 'cvp->prb->socc_num = &%s_socc_num;' % (pname)
            fmt += 'cvp->prb->socc = %s_prb_socc;' % (pname)
        for k in range(data['socc']['socc_num']):
            fmt += tabs+'cvp->prb->socc[%d] = &%s_prb_socc%d;\n' % (k, pname, k)
        for k, pmetric in enumerate(data['socc']['pmetric']):
            for name in pmetric:
                tname = pname+'_'+name+str(k)
                fmt += tabs+'cvp->prb->socc[%d]->%s = &%s_val;\n' % (k, name, tname)
        return fmt

    def _get_init_parameters_c_tmpl(self, terms):
        fmt = ''
        tab = '    '
        for name in terms:
            idx = (self.prefix+'_'+name).upper()
            assign = 'parameters->%s = cvp->par[%s]->data' % (name, idx)
            fmt += 1*tab + assign + ';\n'
        return fmt
