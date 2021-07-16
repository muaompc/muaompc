class MPC(object):
    def __init__(self, prbname, ccg, ddg, sdg, solver, former):
        self.prbname = prbname
        self.prefix = ccg.prefix
        self.destdir = ccg.destdir
        self.path = ccg.path
        self.ddg = ddg  # dynamic data gen
        self.sdg = sdg  # static data gen
        self.solver = solver
        self.former = former
