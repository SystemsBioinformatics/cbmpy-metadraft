"""
This file generates the complete default set of MetaDraft template models
"""


import os
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
from biotools import createSeqplusModel

def build_bigg1():
    #ecoli
    eco = ('eco', "100908EcoliMG1655.gbff", 'Ecoli_iAF1260.glc.l3.xml')
    # #H Pylori
    hpyl = ('hpyl', 'Hpylori_iIT341-NC_000915.gbk', 'Hpylori_iIT341.xml')
    # #M Barkeri
    mbar = ('mbar', ["Mbarkeri_iAF692-NC_007355.gbk", "Mbarkeri_iAF692-p-NC_007349.gbk"], 'Mbarkeri_iAF692.xml')
    # #M tuberculosis
    mtub = ('mtub', 'Mtuberculosis_iNJ661-NC_000962.gbk', 'Mtuberculosis_iNJ661.xml')
    # #S aureus
    saur = ('saur', ['Saureus_iSB619-NC_002745.gbk', 'Saureus_iSB619-p-NC_003140.gbk'], 'Saureus_iSB619.xml')
    libSet = 'bigg1'
    models = [eco, hpyl, mtub, saur]
    return libSet, models, False

def build_bigg2():
    #BiGG2
    bar = ('bar', 'iAF692-(CP000099.1).gbk', 'iAF692.xml')
    clos = ('clos', 'iHN637-(CP000099.1).gbk ', 'iHN637.xml')
    eco = ('eco', 'iJO1366-(NC_000913.3).gbk', 'iJO1366.xml')
    bsub = ('bsub', 'iYO844-(AL009126.3).gbk', 'iYO844.xml')
    sty = ('sty', 'STM_v1_0-(AE006468.1).gbk', 'STM_v1_0.xml')
    #sdy = ('sdy', 'iSDY_1059-(CP000034.1).gbk', 'iSDY_1059.xml')
    sce = ('sce', ["iMM904-(GCF_000146045_2).gbff"], 'iMM904.xml')
    psd = ('psd', ["iJN746-(NC_002947_4).gbk"], 'iJN746.xml')
    geo = ('geo', 'iAF987-(CP000148.1).gbk', 'iAF987.xml')
    kle = ('kle', 'iYL1228-(CP000647.1).gbk', 'iYL1228.xml')
    hpy = ('hpy', 'iIT341-(NC_000915.1).gbk', 'iIT341.xml')
    saur = ('saur', 'iSB619-(NC_002745.2).gbk', 'iSB619.xml')
    shig = ('shig', 'iS_1188-(AE014073.1).gbk', 'iS_1188.xml')
    #mus = ('mus', 'iMM1415-(GCF_000001635.25).gbff', 'iMM1415.xml')
    libSet = 'bigg2'
    models = [bar, clos, eco, bsub, sty, sce, psd, geo, kle, hpy, saur, shig]
    ## There is something funny with this sequence file, model genes are marked as unknown
    ## need to check template generator
    #models += [mus]
    return libSet, models, False

def build_lab():
    #LAB models
    llac1363 = ('llac1363', 'Lactococcus_lactis_MG1363.gbk', 'Lactococcus_lactis_MG1363.xml')
    lplawcfs1 = ('lplawcfs1', 'Lactobacillus_plantarum_WCFS1.gbk', 'Lactobacillus_plantarum_WCFS1.xml')
    libSet = 'lab'
    models = [llac1363, lplawcfs1]
    return libSet, models, False

def build_yeast_B1_v2():
    # needs FBCV2 to work around libSBML bug
    scer = ('scer', ["Saccharomyces_cerevisiae_S288c-(BK006935.2).gbff", "Saccharomyces_cerevisiae_S288c-(BK006935.2)-mito.gbff"], 'Scerevisiae_iND750.xml')
    libSet = 'bigg1'
    models = [scer]
    return libSet, models, True

def build_synecocystis_v2():
    # needs FBCV2 to work around libSBML bug
    syn = ('syn', 'iJN678-(BA000022.2).gbk', 'iJN678.xml')
    libSet = 'bigg2'
    models = [syn]
    return libSet, models, True

def build_ooe_v2():
    psu1 = ('psu1', 'ooeni_BAA_331_PSU1-(GCF_000014385.1_ASM1438v1).gbff', 'ooeni_BAA_331_PSU1.fbc.xml')
    leum = ('leum', 'BMID000000141913(GCF_000014445.1_ASM1444v1_genomic).gbff', 'BMID000000141913.fbc.xml')
    lsl = ('lsl', 'BMID000000141757(GCF_000008925.1_ASM892v1).gbff', 'BMID000000141757.fbc.xml')
    lsa = ('lsa', 'BMID000000140363(GCF_000026065.1_ASM2606v1).gbff', 'BMID000000140363.fbc.xml')
    lj = ('lj', 'BMID000000140939(GCF_000008065.1_ASM806v1).gbff', 'BMID000000140939.fbc.xml')
    lba = ('lba', 'BMID000000141178(GCF_000011985.1_ASM1198v1).gbff', 'BMID000000141178.fbc.xml')
    libSet = 'ooe'
    models = [psu1, leum, lsl, lsa, lj, lba]

    return libSet, models, True



if __name__ == '__main__':
    # setup data directories
    dataDirBase = os.path.join(cDir, 'data')
    modelDir = os.path.join(cDir, 'lib_model')
    genedb = os.path.join(cDir, 'dbx', '_metadraft_genedb.sql')

    # create a list of template models to create
    # still not good; build_synecocystis_v2()
    #build_all = [build_bigg1(), build_bigg2(), build_lab(), build_yeast_B1_v2()]
    build_all = [build_ooe_v2()]

    # build the models
    for libSet, models, useV2 in build_all:
        dataDir = os.path.join(dataDirBase, libSet)
        print('\n*****\nProcessing...\n*****\n{}\n'.format(models))
        print(createSeqplusModel(models, dataDir, modelDir, libSet, gene_db=genedb, add_cobra_annot=False, useV2=useV2, compress_output=False))



