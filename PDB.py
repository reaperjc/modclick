import os
import functions


class PDB:
    def __init__(self, pdb, inputs):
        self.pdb = os.path.abspath(pdb)
        self.abspath = os.path.abspath(pdb)

        partes = self.abspath.split("/")

        self.name = partes[-1][:-4]
        self.path = "/".join(partes[:-1] + ["", ])

        self.original_lines = []


        #check if pdb file exist!!!
        functions.check(self.pdb)


        # self.modifyPDBAtoms(inputs)


    def modifyPDBAtoms(self, inputs):
        '''writes new PDB file with atomtypes modifications. The new file has the "Mod" word before .pdb extension.
        If clickOutput flag is True, the pdb file will be modified with normal
        IUPAC names for a standard PDB file'''

        d_keywordRecord = {"ATOM": None, "HETATM": None, "ANISOU": None, "SIGUIJ": None}

        l_newPDB = []
        for line in open(self.pdb):
            if line[:6].strip() in d_keywordRecord:
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                chain = line[21]
                resnum = line[22:26].strip()

                # save line in original_lines list
                self.original_lines.append(line)
                # first we can modify a line just if the table has the right resname and atomname
                if inputs.d_table.has_key(resname):
                    if inputs.d_table[resname].has_key(atomname):
                        if inputs.f:
                            # if -f flag is given, we use only the line with residue specified in d_pdbResnumTable
                            if self.name in inputs.d_pdbResnumTable and \
                                chain in inputs.d_pdbResnumTable[self.name] and \
                                resnum in inputs.d_pdbResnumTable[self.name][chain]:

                                l_newPDB.append(line[:12] + functions.formatAtomName(inputs.d_table[resname][atomname]) + line[16:])
                            else:
                                # ~ l_newPDB.append(line[:12] + formatAtomName(d_table[resname][atomname]) + line[16:])
                                l_newPDB.append(line)
                        else:
                            l_newPDB.append(line[:12] + functions.formatAtomName(inputs.d_table[resname][atomname]) + line[16:])
                    else:
                        l_newPDB.append(line)
                else:
                    l_newPDB.append(line)
            else:
                l_newPDB.append(line)

        #write modified pdb
        f = open("%sMod.pdb" % self.pdb[:-4], "w")

        for i in l_newPDB:
            f.write(i)
        f.close()