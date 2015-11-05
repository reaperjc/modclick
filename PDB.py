import os





class PDB:
    def __init__(self, pdb):
        self.pdb     = os.path.abspath(pdb)
        self.abspath = os.path.abspath(pdb)
        self.partes  = self.abspath.split("/")
        self.name    = self.partes[-1][:-4]
        self.path    = "/".join(self.partes[:-1] + ["", ])

