classdef {prefix}ctl < handle
    properties(SetAccess=private)
        data
    end

    properties(Dependent)
        parameters
        conf
        prb
        u_opt
    end

    methods(Static, Access=private)
        ctl = {prefix}_matlab_create_ctl(fname)
        ctl = {prefix}_matlab_form_problem(ctl)
        ctl = {prefix}_matlab_solve_problem(ctl)
    end

    methods
        function self = {prefix}ctl(fname)
            self.data = self.{prefix}_matlab_create_ctl(fname);
        end

        function form_problem(self)
            self.data = self.{prefix}_matlab_form_problem(self.data);
        end

        function solve_problem(self)
            self.data = self.{prefix}_matlab_solve_problem(self.data);
        end

        function u_opt = get.u_opt(self)
            u_opt = self.data.u_opt;
        end

        function set.u_opt(self, u_opt)
            self.data.u_opt = u_opt;
        end

        function parameters = get.parameters(self)
            parameters = self.data.parameters;
        end

        function set.parameters(self, parameters)
            self.data.parameters = parameters;
        end

        function conf = get.conf(self)
            conf = self.data.conf;
        end

        function set.conf(self, conf)
            self.data.conf = conf;
        end

        function prb = get.prb(self)
            prb = self.data.prb;
        end

        function set.prb(self, prb)
            self.data.prb = prb;
        end

    end
end
