"""
This file generates the complete default set of MetaDraft template models
"""


import os
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
from biotools import createSeqplusModel

def build_bigg1():
    "Build BiGG1"
    eco = ('eco', "100908EcoliMG1655.gbff", 'Ecoli_iAF1260.glc.l3.xml')
    hpyl = ('hpyl', 'Hpylori_iIT341-NC_000915.gbk', 'Hpylori_iIT341.xml')
    mbar = ('mbar', ["Mbarkeri_iAF692-NC_007355.gbk", "Mbarkeri_iAF692-p-NC_007349.gbk"], 'Mbarkeri_iAF692.xml')
    mtub = ('mtub', 'Mtuberculosis_iNJ661-NC_000962.gbk', 'Mtuberculosis_iNJ661.xml')
    saur = ('saur', ['Saureus_iSB619-NC_002745.gbk', 'Saureus_iSB619-p-NC_003140.gbk'], 'Saureus_iSB619.xml')
    libSet = 'bigg1'
    models = [eco, hpyl, mtub, saur]
    return libSet, models, True

def build_lab():
    "LAB models"
    # older versions of the models
    ##llac1363 = ('llac1363', 'Lactococcus_lactis_MG1363.gbk', 'Lactococcus_lactis_MG1363.xml')
    ##lplawcfs1 = ('lplawcfs1', 'Lactobacillus_plantarum_WCFS1.gbk', 'Lactobacillus_plantarum_WCFS1.xml')
    ##lplawcfs1 = ('lplawcfs1', 'Lactobacillus_plantarum_WCFS1.gbk', 'Lactobacillus_plantarum_WCFS1_includingChargeChemForm.xml')
    # newer version supplied by Sebastian``
    lla = ('lla', 'iLLA516.gb', 'iLLA516.xml')
    lpl = ('lpl', 'iLPL728.gb', 'iLPL728.xml')
    libSet = 'lab'
    models = [lla, lpl]
    return libSet, models, True

def build_bigg2_original():
    "Build BiGG2 original set"
    mba = ('mba', 'iAF692-(CP000099.1).gbk', 'iAF692.xml')
    clj = ('clj', 'iHN637-(CP000099.1).gbk ', 'iHN637.xml')
    eco = ('eco', 'iJO1366-(NC_000913.3).gbk', 'iJO1366.xml')
    bsu = ('bsu', 'iYO844-(AL009126.3).gbk', 'iYO844.xml')
    sty = ('sty', 'STM_v1_0-(AE006468.1).gbk', 'STM_v1_0.xml')
    #sdy = ('sdy', 'iSDY_1059-(CP000034.1).gbk', 'iSDY_1059.xml') # suspicious
    sce = ('sce', ["iMM904-(GCF_000146045_2).gbff"], 'iMM904.xml')
    psd = ('psd', ["iJN746-(NC_002947_4).gbk"], 'iJN746.xml')
    gme = ('gme', 'iAF987-(CP000148.1).gbk', 'iAF987.xml')
    kpn = ('kpn', 'iYL1228-(CP000647.1).gbk', 'iYL1228.xml')
    hpy = ('hpy', 'iIT341-(NC_000915.1).gbk', 'iIT341.xml')
    saur = ('saur', 'iSB619-(NC_002745.2).gbk', 'iSB619.xml')
    shig = ('shig', 'iS_1188-(AE014073.1).gbk', 'iS_1188.xml')
    scer = ('scer', ["Saccharomyces_cerevisiae_S288c-(BK006935.2).gbff", "Saccharomyces_cerevisiae_S288c-(BK006935.2)-mito.gbff"], 'iND750.xml')
    libSet = 'bigg2'
    models = [mba, clj, eco, bsu, sty, sce, psd, gme, kpn, hpy, saur, shig, scer] # sdy
    return libSet, models, True


def build_bigg2_new():
    eco9 = ('eco9', 'iECABU_c1320.gb', 'iECABU_c1320.xml')
    ssp = ('ssp', 'iJN678.gb', 'iJN678.xml')
    eco5 = ('eco5', 'iBWG_1329.gb', 'iBWG_1329.xml')
    eco28 = ('eco28', 'iECSP_1301.gb', 'iECSP_1301.xml')
    eco51 = ('eco51', 'ic_1306.gb', 'ic_1306.xml')
    eco46 = ('eco46', 'iUMNK88_1353.gb', 'iUMNK88_1353.xml')
    eco12 = ('eco12', 'iECDH10B_1368.gb', 'iECDH10B_1368.xml')
    eco_core = ('eco_core', 'e_coli_core.gb', 'e_coli_core.xml')
    eco15 = ('eco15', 'iECED1_1282.gb', 'iECED1_1282.xml')
    eco29 = ('eco29', 'iECUMN_1333.gb', 'iECUMN_1333.xml')
    eco3 = ('eco3', 'iAPECO1_1312.gb', 'iAPECO1_1312.xml')
    eco33 = ('eco33', 'iETEC_1333.gb', 'iETEC_1333.xml')
    eco24 = ('eco24', 'iECP_1309.gb', 'iECP_1309.xml')
    eco17 = ('eco17', 'iECIAI1_1343.gb', 'iECIAI1_1343.xml')
    sau = ('sau', 'iSB619.gb', 'iSB619.xml')
    eco20 = ('eco20', 'iECO103_1326.gb', 'iECO103_1326.xml')
    eco14 = ('eco14', 'iECD_1391.gb', 'iECD_1391.xml')
    eco22 = ('eco22', 'iECO26_1355.gb', 'iECO26_1355.xml')
    eco23 = ('eco23', 'iECOK1_1307.gb', 'iECOK1_1307.xml')
    sdy = ('sdy', 'iSDY_1059.gb', 'iSDY_1059.xml')
    sbo1 = ('sbo1', 'iSbBS512_1146.gb', 'iSbBS512_1146.xml')
    eco30 = ('eco30', 'iECW_1372.gb', 'iECW_1372.xml')
    sce1 = ('sce1', 'iND750.gb', 'iND750.xml')
    eco8 = ('eco8', 'iEC55989_1330.gb', 'iEC55989_1330.xml')
    cgr = ('cgr', 'iCHOv1.gb', 'iCHOv1.xml')
    sbo = ('sbo', 'iSBO_1134.gb', 'iSBO_1134.xml')
    eco41 = ('eco41', 'iJR904.gb', 'iJR904.xml')
    eco45 = ('eco45', 'iUMN146_1321.gb', 'iUMN146_1321.xml')
    eco38 = ('eco38', 'iEcolC_1368.gb', 'iEcolC_1368.xml')
    eco35 = ('eco35', 'iEcE24377_1341.gb', 'iEcE24377_1341.xml')
    sfl2 = ('sfl2', 'iSFxv_1172.gb', 'iSFxv_1172.xml')
    eco50 = ('eco50', 'iZ_1308.gb', 'iZ_1308.xml')
    eco32 = ('eco32', 'iEKO11_1354.gb', 'iEKO11_1354.xml')
    sce = ('sce', 'iMM904.gb', 'iMM904.xml')
    eco37 = ('eco37', 'iEcSMS35_1347.gb', 'iEcSMS35_1347.xml')
    eco39 = ('eco39', 'iG2583_1286.gb', 'iG2583_1286.xml')
    eco44 = ('eco44', 'iNRG857_1313.gb', 'iNRG857_1313.xml')
    eco16 = ('eco16', 'iECH74115_1262.gb', 'iECH74115_1262.xml')
    eco7 = ('eco7', 'iEC042_1314.gb', 'iEC042_1314.xml')
    eco43 = ('eco43', 'iML1515.gb', 'iML1515.xml')
    eco26 = ('eco26', 'iECSE_1348.gb', 'iECSE_1348.xml')
    mmu = ('mmu', 'iMM1415.gb', 'iMM1415.xml')
    eco49 = ('eco49', 'iY75_1357.gb', 'iY75_1357.xml')
    ppu = ('ppu', 'iJN746.gb', 'iJN746.xml')
    eco31 = ('eco31', 'iECs_1301.gb', 'iECs_1301.xml')
    lla = ('lla', 'iNF517.gb', 'iNF517.xml')
    eco47 = ('eco47', 'iUTI89_1310.gb', 'iUTI89_1310.xml')
    ype = ('ype', 'iPC815.gb', 'iPC815.xml')
    sfl1 = ('sfl1', 'iSF_1195.gb', 'iSF_1195.xml')
    eco27 = ('eco27', 'iECSF_1327.gb', 'iECSF_1327.xml')
    eco40 = ('eco40', 'iJO1366.gb', 'iJO1366.xml')
    eco21 = ('eco21', 'iECO111_1330.gb', 'iECO111_1330.xml')
    sfl = ('sfl', 'iSFV_1184.gb', 'iSFV_1184.xml')
    eco18 = ('eco18', 'iECIAI39_1322.gb', 'iECIAI39_1322.xml')
    eco19 = ('eco19', 'iECNA114_1301.gb', 'iECNA114_1301.xml')
    sfl3 = ('sfl3', 'iS_1188.gb', 'iS_1188.xml')
    eco1 = ('eco1', 'iAF1260.gb', 'iAF1260.xml')
    eco13 = ('eco13', 'iECDH1ME8569_1439.gb', 'iECDH1ME8569_1439.xml')
    tma = ('tma', 'iLJ478.gb', 'iLJ478.xml')
    eco11 = ('eco11', 'iECB_1328.gb', 'iECB_1328.xml')
    eco4 = ('eco4', 'iB21_1397.gb', 'iB21_1397.xml')
    eco42 = ('eco42', 'iLF82_1304.gb', 'iLF82_1304.xml')
    eco2 = ('eco2', 'iAF1260b.gb', 'iAF1260b.xml')
    mtu = ('mtu', 'iNJ661.gb', 'iNJ661.xml')
    eco6 = ('eco6', 'iE2348C_1286.gb', 'iE2348C_1286.xml')
    eco48 = ('eco48', 'iWFL_1372.gb', 'iWFL_1372.xml')
    hpy = ('hpy', 'iIT341.gb', 'iIT341.xml')
    sel = ('sel', 'iJB785.gb', 'iJB785.xml')
    eco25 = ('eco25', 'iECS88_1305.gb', 'iECS88_1305.xml')
    sso = ('sso', 'iSSON_1240.gb', 'iSSON_1240.xml')
    eco10 = ('eco10', 'iECBD_1354.gb', 'iECBD_1354.xml')
    eco36 = ('eco36', 'iEcHS_1320.gb', 'iEcHS_1320.xml')
    eco34 = ('eco34', 'iEcDH1_1363.gb', 'iEcDH1_1363.xml')
    libSet = 'bigg2'
    models = [eco9, eco5, eco51, eco46, eco12, eco_core, eco15, eco29, eco3, eco33, eco24,\
              eco17, sau, eco20, eco14, eco22, eco23, sdy, sbo1, eco30, sce1, eco8, sbo, eco41, eco45,\
              eco38, eco35, sfl2, eco50, eco32, sce, eco37, eco39, eco44, eco16, eco7, eco43, eco26, \
              ppu, lla, eco47, ype, sfl1, eco27, eco40, eco21, sfl, eco18, eco19, sfl3, eco1, \
              tma, eco11, eco4, eco42, eco2, mtu, eco6, eco48, hpy, sel, eco25, sso, eco10, eco36, eco34,
              eco28, eco13, eco49]
    #models = [cgr, ssp, eco31] # incomplete genbank file, incoherent gene identifiers, usable identifers that should be fixed in source model``
    #models = [sbo, eco33] # ,  test
    return libSet, models, True

"""
Note on the gene vs. locus_tag issue, it is now a GenBank best practice to assign locus_tags to genes,
as these are unique and registered in BioProject https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
So models without locus_tags are probably old, what I will do is if metadraft finds a locus tag it will use it otherwise it will
use the gene name (these are not unique except in the model namespace!).
"""


if __name__ == '__main__':
    # setup data directories
    dataDirBase = os.path.join(cDir, 'data')
    modelDir = os.path.join(cDir, 'lib_model')
    genedb = os.path.join(cDir, 'dbx', '_metadraft_genedb.sql')

    # create a list of template models to create
    # still not good; build_synecocystis_v2()

    build_all = [build_lab(), build_bigg1(), build_bigg2_original(), build_bigg2_new()]
    #build_all = [build_bigg2_new()]

    # build the models
    for libSet, models, useV2 in build_all:
        dataDir = os.path.join(dataDirBase, libSet)
        print('\n*****\nProcessing...\n*****\n{}\n'.format(models))
        print(createSeqplusModel(models, dataDir, modelDir, libSet, gene_db=genedb, add_cobra_annot=False, useV2=useV2, compress_output=False))



"""
##def build_ooe_v2():
    ##"OOE models - not part of standard metadraft"
    ##psu1 = ('psu1', 'ooeni_BAA_331_PSU1-(GCF_000014385.1_ASM1438v1).gbff', 'ooeni_BAA_331_PSU1.fbc.xml')
    ##leum = ('leum', 'BMID000000141913(GCF_000014445.1_ASM1444v1_genomic).gbff', 'BMID000000141913.fbc.xml')
    ##lsl = ('lsl', 'BMID000000141757(GCF_000008925.1_ASM892v1).gbff', 'BMID000000141757.fbc.xml')
    ##lsa = ('lsa', 'BMID000000140363(GCF_000026065.1_ASM2606v1).gbff', 'BMID000000140363.fbc.xml')
    ##lj = ('lj', 'BMID000000140939(GCF_000008065.1_ASM806v1).gbff', 'BMID000000140939.fbc.xml')
    ##lba = ('lba', 'BMID000000141178(GCF_000011985.1_ASM1198v1).gbff', 'BMID000000141178.fbc.xml')
    ##libSet = 'ooe'
    ##models = [psu1, leum, lsl, lsa, lj, lba]

    ##return libSet, models, True
"""
