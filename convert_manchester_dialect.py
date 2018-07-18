"""
This file converts files from the manchester SBML dialect

"""

import os, re, shutil
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import cbmpy

def convert_ooe():
    m1 = ['ooeni_BAA_331_PSU1.xml', re.compile('OEOE_.*?_i')]
    m2 = ['BMID000000141913.xml', re.compile('LEUM_.*?_i')]
    m3 = ['BMID000000141757.xml', re.compile('LSL_.*?_i')]
    m4 = ['BMID000000140363.xml', re.compile('LSA.*?_i')]
    m5 = ['BMID000000140939.xml', re.compile('LJ.*?_i')]
    m6 = ['BMID000000141178.xml', re.compile('LBA.*?_i')]
    libSet = 'ooe'
    models = [m1, m2, m3, m4, m5, m6]
    return libSet, models, True


def convertManchesterDialect(models, modelDir):
    for m, gm in models:
        print(m)
        mFin = open(os.path.join(modelDir, m), 'r')
        sbmlS = mFin.read()
        mFin.close()
        genes = genes0 = re.findall(gm, sbmlS)
        genes = [g.strip()[:-2] for g in genes if g.strip().endswith('_i')]
        done = []
        for g_ in range(len(genes0)):
            if genes0[g_] not in done:
                sbmlS = sbmlS.replace(genes0[g_], genes[g_])
                done.append(genes0[g_])
        m1 = os.path.join(modelDir, m.replace('.xml','.1.xml'))
        mFout = open(m1, 'w')
        mFout.write(sbmlS)
        mFout.flush()
        mFout.close()
        badmod = cbmpy.readCOBRASBML(m1, delete_intermediate=True)
        m2 = os.path.join(modelDir, m.replace('.xml','.fbc.xml'))
        badmod.parameters = []
        for fb in badmod.flux_bounds:
            fb.__sbo_term__ = ''
            fb.setName('')
        glbl = badmod.getGeneLabels()
        delSpecies = []
        print(len(badmod.species))
        for s in badmod.species:
            if not cbmpy.CBCommon.checkChemFormula(s.getChemFormula(), quiet=True):
                s.setChemFormula('')
            if s.getId() in glbl:
                if len(s.isReagentOf()) == 0:
                    delSpecies.append(s.getId())
        for s in delSpecies:
            badmod.deleteSpecies(s)
        print(len(badmod.species))
        for r in badmod.reactions:
            r._modifiers_ = []
        os.remove(m1)
        #cbmpy.writeSBML3FBC(badmod, m2)
        cbmpy.writeSBML3FBCV2(badmod, m2)
        bettermod = cbmpy.readSBML3FBC(m2)






if __name__ == '__main__':
    # setup data directories
    dataDirBase = os.path.join(cDir, 'data')
    modelDir = os.path.join(cDir, 'lib_model')

    convert_all = [convert_ooe()]

    # build the models
    for libSet, models, useV2 in convert_all:
        dataDir = os.path.join(dataDirBase, libSet)
        print('\n*****\nProcessing...\n*****\n{}\n'.format(models))
        convertManchesterDialect(models, dataDir)





