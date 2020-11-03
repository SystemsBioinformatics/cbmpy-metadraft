"""
MetaToolkit: MetaDraft
======================

MetaDraft: for the reconstruction of Genome Scale Models

MetaToolkit: MetaDraft (https://github.com/SystemsBioinformatics/cbmpy-metadraft)
Copyright (C) 2016-2019 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier PhD
Contact email: b.g.olivier@vu.nl

"""

from __future__ import division, print_function
from __future__ import absolute_import

import os, sys, random, copy, time, math, json, xlwt
import sqlite3

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
dataDir = os.path.join(cDir, 'data')

try:
    import Bio
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
except ImportError:
    print(
        '\n####\nERROR: please install Biopython.\n- conda install biopython\n or\n- pip install biopython\n'
    )
    os.sys.exit(1)

# Compatible with both pre- and post Biopython 1.78:
oldbiopython = True
try:
    from Bio.Alphabet import ProteinAlphabet

    print(
        '\n####\nINFO: please update your BioPython ({}) to 1.78 or newer.\n- conda update biopython\n or\n- pip update biopython\n'.format(
            Bio.__version__
        )
    )
except ImportError:
    oldbiopython = False
    print('You are using BioPython {}'.format(Bio.__version__))

import cbmpy

# backwards compatability
cbm = cbmpy


def web_NCBIblastq(seq, options=None):
    """
    Uses NCBI web API to perform a blast on the input sequences, returns BLAST XML

     - *seq* a string formatted as a FASTA file
     - *options* a dictionary of options, None uses default

      - *method* [default='blastp'] the blast method to use
      - *dbase* [default='nr'] the NCBI database
      - *hit_num* [default=10] the number of hits to return
      - *matrix* [default='BLOSUM62'] scoring matrix to use

    """
    optionsx = {'method': 'blastp', 'dbase': 'nr', 'hit_num': 10, 'matrix': 'BLOSUM62'}

    if options is not None:
        for o in options:
            if o in optionsx:
                optionsx[o] = options[o]

    print('\nContacting NCBI ...')
    bres = NCBIWWW.qblast(
        optionsx['method'],
        optionsx['dbase'],
        seq,
        hitlist_size=optionsx['hit_num'],
        matrix_name=optionsx['matrix'],
    )
    print('\ndone.')
    return bres


def readGenBankFullFile(fname):
    """
    Read in a full genbank file and extract all annotation. Returns three dictionaries containing CDS, Gene and Other information.

     - *fname* the full filename of the input filename in GenBank full format (typically *.gbk or *.gb)

    """
    print("\nINFO: This method only works with single sequence files\n")

    SEQ = SeqIO.read(fname, 'genbank')

    CDS_features = {}
    GENE_features = {}
    MISC_features = {}
    for fe_ in SEQ.features:
        if fe_.type == 'CDS':
            CDS_features[fe_.qualifiers['locus_tag'][0]] = {
                'ref_db': fe_.ref_db,
                'id': fe_.id,
                'qualifiers': fe_.qualifiers,
                'strand': fe_.strand,
                'ref': fe_.ref,
                'sub_features': fe_.sub_features,
            }
        elif fe_.type == 'gene':
            GENE_features[fe_.qualifiers['locus_tag'][0]] = {
                'ref_db': fe_.ref_db,
                'id': fe_.id,
                'qualifiers': fe_.qualifiers,
                'strand': fe_.strand,
                'ref': fe_.ref,
                'sub_features': fe_.sub_features,
                'cds': False,
            }
        else:
            MISC_features[fe_.type] = fe_

    for g_ in GENE_features:
        if g_ in CDS_features:
            GENE_features[g_]['cds'] = True

    return CDS_features, GENE_features, MISC_features


def readFASTAFile(fname):
    """
    Read in a full genbank file and extract all annotation. Returns three dictionaries containing CDS, Gene and Other information.

     - *fname* the full filename of the input filename in GenBank full format (typically *.gbk or *.gb)

    """

    fasta_sequences = SeqIO.parse(fname, 'fasta')

    CDS_features = {}
    GENE_features = {}
    MISC_features = {}

    for fasta in fasta_sequences:
        CDS_features[fasta.id] = {
            'ref_db': None,
            'id': fasta.id,
            'qualifiers': None,
            'strand': None,
            'ref': None,
            'sub_features': None,
        }

    ##SEQ = SeqIO.read(fname, 'genbank')
    ##for fe_ in SEQ.features:
    ##if fe_.type == 'CDS':
    ##CDS_features[fe_.qualifiers['locus_tag'][0]] = {'ref_db' : fe_.ref_db,
    ##'id' : fe_.id,
    ##'qualifiers' : fe_.qualifiers,
    ##'strand' : fe_.strand,
    ##'ref' : fe_.ref,
    ##'sub_features' : fe_.sub_features
    ##}
    ##elif fe_.type == 'gene':
    ##GENE_features[fe_.qualifiers['locus_tag'][0]] = {'ref_db' : fe_.ref_db,
    ##'id' : fe_.id,
    ##'qualifiers' : fe_.qualifiers,
    ##'strand' : fe_.strand,
    ##'ref' : fe_.ref,
    ##'sub_features' : fe_.sub_features,
    ##'cds' : False
    ##}
    ##else:
    ##MISC_features[fe_.type] = fe_

    ##for g_ in GENE_features:
    ##if g_ in CDS_features:
    ##GENE_features[g_]['cds'] = True

    return CDS_features, GENE_features, MISC_features


def checkModelLocusTags(sbml, genbank, allow_gene_names=False):
    """
    Checks the gene identifiers (assuming they are locus tags) against a genbank file of the same organism

    - *sbml* the model SBML (*.xml) file
    - *genbank* the associated GenBank (*.gbk) full file_s) (including CDS annotations and sequences)
    - *allow_gene_names* allow gene names, non-unique as gene identifiers if locus tags are not present. USE WITH CAUTION!!!

    """

    cmod = cbmpy.readSBML3FBC(sbml)

    cntr = 0
    if type(genbank) == str:
        genbank = [genbank]

    gprMap = {}
    gprMapAnnot = {}
    gnoprMap = {}
    fileNum = 0
    no_locus_tag = []
    for G in genbank:
        print('CheckModelLocusTags is processing: {}'.format(G))
        for seq_record in SeqIO.parse(G, "genbank"):
            # print(seq_record.id)
            cntr += 1
        if cntr > 1:
            print("INFO: Multiple sequences encountered in file: {}".format(G))
            # raise RuntimeError, "\nMutltiple sequences encountered in file: {}".format(G)
        # print(repr(seq_record.seq))
        # print(len(seq_record))

        # add all the annotations from the genbank record(s) into the model? Use first record
        if fileNum == 0:

            gprMapAnnot[seq_record.id] = seq_record.annotations.copy()
            try:
                for r_ in gprMapAnnot[seq_record.id]['references']:
                    if r_.pubmed_id != '':
                        cmod.addMIRIAMannotation(
                            'isDescribedBy', 'PubMed', r_.pubmed_id
                        )
                    gprMapAnnot[seq_record.id].pop('references')
                    cmod.setAnnotation('genbank_id', seq_record.id)
                    cmod.setAnnotation('genbank_name', seq_record.name)
            except KeyError:
                print('checkModelLocusTags: no references')

            # global features
            features = [f.qualifiers for f in seq_record.features if f.type == 'source']
            if len(features) > 0:
                features = features[0]
            for f_ in features:
                if f_ == 'db_xref':
                    for ff_ in features[f_]:
                        if ff_.startswith('taxon:'):
                            cmod.setAnnotation(
                                'genbank_taxon_id', ff_.replace('taxon:', '')
                            )
                            break
                cmod.setAnnotation('genbank_{}'.format(f_), features[f_])
            for r_ in gprMapAnnot[seq_record.id]:
                cmod.setAnnotation(
                    'genbank_{}'.format(r_), gprMapAnnot[seq_record.id][r_]
                )
            fileNum += 1

        GBFile = open(G, 'r')

        USING_GENE_NAME = False
        for cds in Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile):
            if cds.seq != None:
                if cds.name != '<unknown name>':
                    gprMap[cds.name] = cds
                # this is dangerous as there is no defined 1:1 mapping
                elif 'gene' in cds.annotations and allow_gene_names:
                    gprMap[cds.annotations['gene']] = cds
                    USING_GENE_NAME = True
                    no_locus_tag.append(cds.annotations['gene'])
                else:
                    gnoprMap[cds.name] = cds
            else:
                gnoprMap[cds.name] = cds
        GBFile.close()

    oldLoTags = []
    # for g_ in cmod.getGeneIds():
    for g_ in cmod.getGeneLabels():
        if g_ not in gprMap:
            # print(g_)
            if g_ != 'None':
                oldLoTags.append(g_)

    # cbmpy.CBTools.pprint.pprint(gprMap)

    oldtags = [
        gprMap[a].annotations['old_locus_tag'].split(' ') + [a]
        for a in gprMap
        if 'old_locus_tag' in gprMap[a].annotations
    ]
    oldtags2 = [
        gnoprMap[a].annotations['old_locus_tag'].split(' ') + [a]
        for a in gnoprMap
        if 'old_locus_tag' in gnoprMap[a].annotations
    ]

    old2new = {}
    old2new2 = {}

    for x in oldtags:
        for y in x:
            if y != x[-1]:
                old2new[y] = x[-1]

    for x in oldtags2:
        for y in x:
            if y != x[-1]:
                old2new2[y] = x[-1]

    print('\n\nChecking locus tags\n===================')
    updated = {}
    unknown = []
    noseq = []
    good = []
    F = open(sbml.replace('.xml', '.seqcheck.csv'), 'w')
    # if cmod.__FBC_VERSION__ < 2:
    # geneIDs = cmod.getGeneIds()
    # else:
    # geneIDs = cmod.getGeneLabels()
    geneIDs = cmod.getGeneLabels()
    for g_ in geneIDs:
        if g_ not in gprMap:
            if g_ in old2new:
                print('old {} --> new {}'.format(g_, old2new[g_]))
                updated[g_] = old2new[g_]
                F.write('UPDATED,{},{}\n'.format(g_, old2new[g_]))
            elif g_ in old2new2:
                print('NoSEQ: old {} --> new {}'.format(g_, old2new2[g_]))
                noseq.append((g_, old2new2[g_]))
                F.write('NOSEQ,{},{}\n'.format(g_, old2new2[g_]))
            else:
                print('UNKNOWN gene: {}'.format(g_))
                unknown.append(g_)
                F.write('UNKNOWN,{}\n'.format(g_))
        else:
            good.append(g_)
    F.close()
    if USING_GENE_NAME:
        print(
            '\nINFO: model contained {} genes without locus tags. For these genes \
the /gene name was used instead. This gene name may not be unique, please check model!!'.format(
                len(no_locus_tag)
            )
        )
        # print(no_locus_tag)

    # cbmpy.CBTools.storeObj(gprMap, sbml.replace('.xml', '.seqplus.dat'))

    return good, updated, noseq, unknown, cmod, gprMap, gprMapAnnot, no_locus_tag


def addSeqAnnotation(emod, gpr, good, updated, update_tags=True, no_locus_tag=None):
    # now we annotate the reactions with the sequences and optionally update all the gene
    # names to match the current GenBank locus tags
    geneAnnot = {}

    if no_locus_tag is None:
        no_locus_tag = []

    for gp in emod.gpr:
        # if emod.__FBC_VERSION__ < 2:
        # geneIDs = gp.getGeneIds()
        # else:
        geneIDs = gp.getGeneLabels()
        for g_ in geneIDs:
            if g_ in good:
                if g_ in no_locus_tag:
                    print(
                        '-> {}, {}, {}\n{}'.format(
                            g_, gpr[g_].id, gpr[g_].name, gpr[g_]
                        )
                    )

                emod.getReaction(gp.getProtein()).setAnnotation(
                    'gbank_seq_{}'.format(g_), '{}'.format(str(gpr[g_].seq))
                )
                gp.setAnnotation('gbank_seq', '{}'.format(str(gpr[g_].seq)))
                gp.setAnnotation('gbank_id', '{}'.format(str(gpr[g_].id)))
                gp.setAnnotation('gbank_name', '{}'.format(str(gpr[g_].name)))
                gp.setAnnotation(
                    'gbank_description', '{}'.format(str(gpr[g_].description))
                )
                if g_ not in geneAnnot:
                    geneAnnot[g_] = gp.getAnnotations()
            elif g_ in updated and update_tags:
                print('Updating gene identifier: {} --> {}'.format(g_, updated[g_]))
                emod.getReaction(gp.getProtein()).setAnnotation(
                    'gbank_seq_{}'.format(updated[g_]),
                    '{}'.format(str(gpr[updated[g_]].seq)),
                )
                gp.setAnnotation('gbank_seq', '{}'.format(str(gpr[updated[g_]].seq)))
                gp.setAnnotation('gbank_id', '{}'.format(str(gpr[updated[g_]].id)))
                gp.setAnnotation('gbank_name', '{}'.format(str(gpr[updated[g_]].name)))
                gp.setAnnotation(
                    'gbank_description', '{}'.format(str(gpr[updated[g_]].description))
                )
                if updated[g_] not in geneAnnot:
                    geneAnnot[updated[g_]] = gp.getAnnotations()

    # print(geneAnnot)
    if update_tags and len(updated) > 0:
        if emod.__FBC_VERSION__ < 2:
            # updateGeneIdsFBCv1(emod, updated, geneAnnot, annotation_key='GENE ASSOCIATION', replace_existing=True)
            updateGeneIdsFBCv2(
                emod,
                updated,
                geneAnnot,
                replace_existing=True,
                no_locus_tag=no_locus_tag,
            )
        else:
            updateGeneIdsFBCv2(
                emod,
                updated,
                geneAnnot,
                replace_existing=True,
                no_locus_tag=no_locus_tag,
            )


def updateGeneIdsFBCv2(
    emod, updated, gene_annotation, replace_existing=True, no_locus_tag=None
):
    """
    Update gene locus tags from updated dictionary. If this fails it tries some standard annotation
    keys: GENE ASSOCIATION, GENE_ASSOCIATION, gene_association, gene association.

     - *emod* the model instance
     - *updated* a dictionary with updated gene lables
     - *gene_annotation* update the gene annotation (if necessary)
     - *replace_existing* [default=True] replace existing annotations, otherwise only new ones are added

    """
    if no_locus_tag is None:
        no_locus_tag = []
    for of_ in updated:
        G = emod.getGeneByLabel(of_)
        G.setLabel(updated[of_])


def updateGeneIdsFBCv1(
    emod,
    updated,
    gene_annotation,
    annotation_key='GENE ASSOCIATION',
    replace_existing=True,
):
    """
    Update gene locus tags from updated dictionary. If this fails it tries some standard annotation
    keys: GENE ASSOCIATION, GENE_ASSOCIATION, gene_association, gene association.

     - *emod* the model instance
     - *updated* a dictionary with updated gene lables
     - *gene_annotation* update the gene annotation (if necessary)
     - *annotation_key* the annotation dictionary key that holds the gene association for the protein/enzyme
     - *replace_existing* [default=True] replace existing annotations, otherwise only new ones are added

    """

    if replace_existing:
        for g in emod.getGeneIds():
            emod.deleteGene(g)
        for gpr in [g.getid() for g in emod.gpr]:
            emod.deleteGPRAssociation(gpr)
        print(len(emod.genes), len(emod.gpr))
        emod.__genes_idx__ = []
    gid = name = None
    g0 = len(emod.genes)
    gpr0 = len(emod.gpr)
    ga_keys = []
    for r_ in emod.getReactionIds():
        GA = None
        R = emod.getReaction(r_)
        ##  print r.annotation
        if annotation_key in R.annotation:
            GA = annotation_key
        elif 'GENE ASSOCIATION' in R.annotation:
            GA = 'GENE ASSOCIATION'
        elif 'GENE_ASSOCIATION' in R.annotation:
            GA = 'GENE_ASSOCIATION'
        elif 'gene_association' in R.annotation:
            GA = 'gene_association'
        elif 'gene association' in R.annotation:
            GA = 'gene association'
        if GA != None:
            GPR = R.getAnnotation(GA)
            if GPR != None:
                for l_ in updated:
                    if l_ in GPR:
                        GPR = GPR.replace(l_, updated[l_])
            R.setAnnotation(GA, GPR)
            emod.createGeneProteinAssociation(r_, GPR, gid, name, update_idx=False)
            if GA not in ga_keys:
                ga_keys.append(GA)
    print('INFO: used key(s) \'{}\''.format(ga_keys))
    emod.__updateGeneIdx__()
    for g_ in emod.genes:
        if g_.getLabel() in gene_annotation:
            # print(g_)
            g_.annotation = gene_annotation[g_.getLabel()]
    print(
        'INFO: Added {} new genes and {} associations to model'.format(
            len(emod.genes) - g0, len(emod.gpr) - gpr0
        )
    )


def createSequence(modrefseq, filtered_ids=None, description=None):
    emod = cbmpy.readSBML3FBC(modrefseq)
    if len(emod.genes) == 0:
        emod.createGeneAssociationsFromAnnotations()
    geneseq = {}
    if description == None:
        secDescr = os.path.split(modrefseq)[-1]
    else:
        secDescr = description
    if filtered_ids == None:
        filtered_ids = emod.getReactionIds()
    for r_ in emod.reactions:
        if r_.getId() in filtered_ids:
            for a in r_.annotation:
                if a.startswith('gbank_seq_'):
                    newName = a.replace('gbank_seq_', '')
                    if newName not in geneseq:
                        if oldbiopython:
                            seq = Bio.Seq.Seq(
                                r_.getAnnotation(a), Bio.Alphabet.ProteinAlphabet()
                            )
                        else:
                            seq = Bio.Seq.Seq(r_.getAnnotation(a))

                        geneseq[newName] = Bio.SeqRecord.SeqRecord(
                            seq,
                            id=newName,
                            name=newName,
                            description=secDescr,
                            annotations={"molecule_type": "protein"},
                        )
                        # print('Adding {}'.format(a))
    print('{} genseqs added'.format(len(geneseq)))
    return geneseq


def writeFASTA(fname, sequences, paranoid_style=True):

    F = open(fname, 'w')
    cntr = 0
    # print(len(sequences))
    for seq in sequences:
        if paranoid_style:
            sequences[seq].description = ''

        F.write(sequences[seq].format('fasta'))
        cntr += 1
        if cntr >= 500:
            F.flush()
            cntr = 0
    F.flush()
    F.close()
    print('FASTA sequence file created: {}'.format(fname))


def createMegaGenomeBasic(fname, fdir, input_files=None):
    """
    Searches fdir for .fasta files and concatenates them to creates a megagenome

    - *fname* the output filename
    - *fdir* the directory to scan for .fasta files
    - *input_files* custom order of input files

    """
    catlist = []
    if input_files == None:
        for f_ in os.listdir(fdir):
            if f_.endswith('.fasta'):
                catlist.append(os.path.join(fdir, f_))
    else:
        for f_ in input_files:
            if os.path.exists(os.path.join(fdir, f_)):
                catlist.append(os.path.join(fdir, f_))
            else:
                print('File {} does not exist'.format(os.path.join(fdir, f_)))

    Fout = open(fname, 'w')
    for f_ in catlist:
        F = open(f_, 'r')
        Fout.write(F.read())
        Fout.flush()
        F.close()
    Fout.close()
    print('MegaGenome file created as: {}'.format(fname))


def createMetaProteome(fname, link, optimized=True, paranoid_style=True, K=0.13):
    """
    Creates a MetaProteome, optionally uses linkDict organisms filtered_ids::

     - *fname* the output filename
     - *link* a linkDict
     - *oidList* list of ordered organism keys
     - *optimized* use phylogenetically filtered id's (if present)
     - *paranoid_style* output FASTA files for use with inparanoid
     - *K* [default=0.13] Altschul coefficient used to calculate sequence search space size: K*n*m

    """
    # debug
    # global sequences

    metaprotdat = {}

    F = open(fname, 'w')
    outseq = []
    total_prot_length = 0
    # protein_lengths = {}
    for o_ in link:
        # print(o_)
        if not o_.startswith('__'):
            if 'filtered_id' in link[o_] and optimized:
                print('createMetaProteome [{}] is using filtered ids'.format(o_))
                sequences = createSequence(
                    link[o_]['sbml_out'],
                    filtered_ids=link[o_]['filtered_id'],
                    description=None,
                )
            else:
                print('createMetaProteome [{}] is using raw ids'.format(o_))
                sequences = createSequence(
                    link[o_]['sbml_out'], filtered_ids=None, description=None
                )

            for seq in sequences:
                total_prot_length += len(sequences[seq])
                # protein_lengths[seq] = [len(sequences[seq])]
                if paranoid_style:
                    sequences[seq].description = ''
                F.write(sequences[seq].format('fasta'))
            outseq.append(sequences)
            F.flush()
    F.close()

    # gone for now ... may be useful for debug purposes
    # cbmpy.CBTools.storeObj(outseq, fname.replace('.fasta',''))

    # """
    # search space = n(seq_len)*m(db_len)*K(Altschul)
    # """
    # Km = K*total_prot_length
    # for p_ in protein_lengths:
    # search_space = protein_lengths[p_][0]*Km
    # protein_lengths[p_].append(search_space)
    # protein_lengths[p_].append(math.log(search_space, 2))

    metaprotdat = {'file_name': fname, 'total_length': total_prot_length}
    #'protein_lengths' : protein_lengths,
    #'K' : K}
    if '__metaproteome__' in link:
        link['__metaproteome__'].update(metaprotdat)
    else:
        link['__metaproteome__'] = metaprotdat
    # create report structure
    link['__metaproteome__']['reports'] = {}
    link['__metaproteome__']['reports']['genes'] = {}
    link['__metaproteome__']['reports']['genes']['unmatched'] = []
    link['__metaproteome__']['reports']['genes']['selected'] = []
    link['__metaproteome__']['reports']['genes']['unselected'] = []
    link['__metaproteome__']['reports']['reactions'] = {}
    link['__metaproteome__']['reports']['reactions']['selected'] = []
    link['__metaproteome__']['reports']['reactions']['unselected'] = []
    link['__metaproteome__']['reports']['metabolites'] = {}
    link['__metaproteome__']['reports']['metabolites']['selected'] = []
    link['__metaproteome__']['reports']['metabolites']['unselected'] = []
    link['__metaproteome__']['notes'] = {}

    print('MegaGenome file created as: {}'.format(fname))


# IDMAP0 = None
# IDMAP1 = None
def idFilter(link, oidList):
    """
    Updates the linkDict with ID's filtered in oidList order or "phylogenetic ordering":

     - *link* a linkDict
     - *oidList* list of ordered organism keys

    """
    idmap = []
    for o_ in oidList:
        keys = list(link[o_]['reaction2gene'].keys())
        idmap.append(set(keys))

    idx = copy.deepcopy(idmap[0])
    print(len(idx), ','.join([str(len(l_)) for l_ in idmap]))

    for x_ in range(1, len(idmap)):
        idmap[x_].difference_update(idx)
        idx = idx.union(idmap[x_])
        print(len(idx), ','.join([str(len(l_)) for l_ in idmap]))

    for o_ in range(len(oidList)):
        link[oidList[o_]]['filtered_id'] = list(idmap[o_])


def connectSQLDB(fname):
    """
    Opens a SQLite DB as a CBMPy-DB instance. Returns the instance or None

     - *fname* the sql file name

    """
    if not os.path.exists(fname):
        return None
    db = cbmpy.CBNetDB.DBTools()
    print('INITDB: connecting geneDB')
    db.connectSQLiteDB(fname)
    return db


def addGeneInformationToDB(gbkf, gendb, table, idxc):
    """
    Extracts CDS information from a GenBank file (*.gb, *.gbk *.gbff) and adds it to a gene database

    - *gbkf* a GenBank file (*.gb, *.gbk, *. gbff) containing CDS annotations
    - *gendb* a CBMPy DBTools (SQLite) database that has the suggested format:

     - 'id TEXT PRIMARY KEY', 'pid TEXT', 'type TEXT', 'annotation TEXT', 'db_xref TEXT', 'notes TEXT', 'rdf TEXT', 'sbo TEXT', 'seq TEXT'

    - *table* the table name e.g. 'GENES'
    - *idxc* the index column name e.g. 'id'

    Example call for MetaDraft GeneDB: addGeneInformationToDB(gb_file, gene_db, 'GENES', 'id')

    """
    # print(gbkf)
    assert os.path.exists(gbkf)
    GBFile = open(gbkf, 'r')
    GBcds = Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
    cntr = 0
    # ltags = []
    for cds in GBcds:
        if cds.seq != None and not gendb.checkEntryInColumn(table, idxc, cds.name):
            # ltags.append(cds.name)
            ADD2DB = True
            data = {
                'id': cds.name,
                'pid': cds.id,
                'annotation': json.dumps(cds.annotations),
                #'seq' : str(cds.seq),
                'seq': '',
                'type': 'cds',
                'db_xref': gendb.URLEncode(str(cds.dbxrefs)),
            }

            if data['id'] == '<unknown name>':
                ADD2DB = False
                # adding dodgy non-LT genes is not feasible as no 1:1 mapping is possible
                # if 'gene' in cds.annotations:
                # data['id'] = cds.annotations['gene']
                # else:
                # ADD2DB = False
            else:
                ADD2DB = True
                try:
                    gendb.insertData(table, data, commit=False)
                    cntr += 1
                except sqlite3.IntegrityError:
                    print('ERROR: gene id {} already exists')

    gendb.db_cursor.connection.commit()
    GBFile.close()
    ## turns out this is not needed as db_xrefs can be obtained from cds annotations, left here just in case
    # try:
    # genome = Bio.SeqIO.read(gbkf, 'genbank')
    # for f in genome.features:
    # if 'locus_tag' in f.qualifiers and f.qualifiers['locus_tag'][0] in ltags and f.type == 'gene' and not \
    # gendb.checkEntryInColumn(table, 'db_xref', f.qualifiers['locus_tag'][0]):
    # if 'db_xref' in f.qualifiers:
    # db_xref = ','.join(f.qualifiers['db_xref'])
    # gendb.updateData(table, idxc, f.qualifiers['locus_tag'][0], {'db_xref' : '\"{}\"'.format(db_xref)}, commit=False)
    # gendb.db_cursor.connection.commit()
    # except ValueError:
    # print('INFO: could not add extra gene annotations for file: {}'.format(gbkf))

    # gendb.dumpTableToTxt(table, 'sqldump.txt')

    print('INFO: addGeneInformationToDB committed {} new records'.format(cntr))


def createBasicFASTAfromFile(
    gbkf, ext_replace='.in.fasta', gene_prefix=None, remove_pseudo=True
):
    """
    Extracts sequence information from a FASTA (*.fasta) or GenBank file (*.gb, *.gbk *.gbff) and writes a
    simplified FASTA output.

    - *gbkf* FASTA (*.fasta) or GenBank file (*.gb, *.gbk, *. gbff) containing CDS annotations
    - *ext_replace* [default='.in.fasta'] replace the extension (as above) with this
    - *gene_prefix* [default=None] prefix gene names with this in the fasta file
    - *remove_pseudo* [default=True] Remove pseudo genes from output sequence if they are encoded in GenBank format [pseudo=true]
    any pseudo genes removed will be written to <gbkf>.pseudo.txt
    """

    if type(gbkf) == str:
        gbkf = [gbkf]

    proteins = {}
    pseudo_genes = {}
    cntr = 0
    for fasta in gbkf:
        psugenes = []
        if fasta.endswith('.fasta') or fasta.endswith('.faa'):
            if cntr == 0:
                if fasta.endswith('.fasta'):
                    outF = fasta.replace('.fasta', '.in.fasta')
                elif fasta.endswith('.faa'):
                    outF = fasta.replace('.faa', '.in.fasta')
                cntr += 1
            for seq_record in SeqIO.parse(fasta, "fasta"):
                # remove pseudo genes from output file. This is based on a GenBank FASTA file format.
                if remove_pseudo and '[pseudo=true]' in seq_record.description:
                    psugenes.append(seq_record.description)
                    continue
                seq_record.description = ''
                if gene_prefix is not None:
                    seq_record.id = gene_prefix + seq_record.id
                proteins[seq_record.id] = seq_record
            # print(psugenes)
            print('Pseudo genes:', len(psugenes))

        elif fasta.endswith('.gbk') or fasta.endswith('.gb') or fasta.endswith('.gbff'):
            if cntr == 0:
                outF = (
                    fasta.replace('.gbk', '.in.fasta')
                    .replace('.gbff', '.in.fasta')
                    .replace('.gb', '.in.fasta')
                )
                cntr += 1
            GBFile = open(fasta, 'r')
            GBcds = Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
            for cds in GBcds:
                if cds.seq is not None:
                    if cds.name != '<unknown name>':
                        cds.id = cds.name
                    else:
                        cds.id = (
                            cds.id.replace('.', '_').replace(':', '_').replace('-', '_')
                        )
                    cds.description = ''
                    if gene_prefix is not None:
                        cds.id = gene_prefix + cds.id
                    if cds.name != '<unknown name>':
                        proteins[cds.name] = cds
                    else:
                        proteins[cds.id] = cds
            GBFile.close()
        else:
            raise RuntimeError('ERROR: Unknown file: {}'.format(fasta))
        if remove_pseudo and len(psugenes) > 0:
            pseudo_genes[fasta] = psugenes

    writeFASTA(outF, proteins)

    pseudo_cntr = 0
    # write pseudo_genes to file
    if remove_pseudo and len(pseudo_genes) > 0:
        # print(pseudo_genes)
        with open(outF + '.pseudo.txt', 'w') as F:
            for s in pseudo_genes:
                F.write('#{}\n'.format(s))
                for g in pseudo_genes[s]:
                    F.write('{}\n'.format(g))
                    pseudo_cntr += 1

    print('\nPseudogenes: {}'.format(pseudo_cntr))
    print('Proteins: {}\n'.format(len(proteins)))

    return outF


def createSeqplusModel(
    modlist,
    data_dir,
    model_dir,
    lib_set,
    gene_db=None,
    add_cobra_annot=False,
    useV2=True,
    compress_output=False,
):
    """
    This function takes pairs of models and genbank files, updates the locus tags in the model (if necessary),
    merges GenBank annotation, including AA sequences, into the SBML model::

     - *modlist* a list of SBML/Genbank (id, GenBank, SBML) tuples for example: `('eco', "100908EcoliMG1655.gbk", 'Ecoli_iAF1260.glc.l3.xml')`

      - *id* short model identifier
      - *genbank* one (as a string) or more than one (as a list) of full format GenBank files with CDS annotation
      - *sbml* the SBML L3 FBC V1 model that contains GeneAssociation information

     - *data_dir* the input data directory
     - *model_dir* output seqplus models
     - *lib_set* the library set (subdirectory of data)
     - *gene_db* [default=None] the local _metadraft_genedb database
     - *add_cobra_annot* [default=False] add COBRAML notes annotation to new SBML file
     - *useV2* use FBC version 2
     - *compress_output* [default=False] create ZIP compressed SBML files

    Output is as follows

     - the annotated models
     - a JSON formatted link dictionary for reassociating genes, proteins and models.

    """

    if model_dir == None:
        model_dir = data_dir
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    oidList = []
    taxIdList = []

    for oid, fgb, fmod in modlist:
        print('Processing:', oid, fgb, fmod)
        if '-' in oid and '-' in lib_set:
            print(
                '\nERROR: \"-\" not allowed in library name ({}) or organism shortcut ({})!'.format(
                    lib_set, oid
                )
            )
            continue
        lib_set = lib_set.replace('-', '')
        oid = oid.replace('-', '')
        oid = '{}-{}'.format(lib_set, oid)
        oidList.append(oid)
        linkDict = {}
        linkDict['__idx__'] = {}
        linkDict[oid] = {}
        LD = linkDict[oid]
        LD['genbank_in'] = fgb
        LD['sbml_in'] = fmod
        LD['data_path'] = data_dir

        fmod = os.path.join(data_dir, fmod)
        if type(fgb) == str:
            fgb = [os.path.join(data_dir, fgb)]
        else:
            fgb = [os.path.join(data_dir, f) for f in fgb]

        # add CDS gene information to database
        if gene_db is not None and os.path.exists(gene_db):
            DB = connectSQLDB(gene_db)
            for f_ in fgb:
                addGeneInformationToDB(f_, DB, 'GENES', 'id')
            DB.closeDB()
            del DB, f_

        # First we see if we can map the gene association "genes" to the genbank locus tags
        (
            good,
            updated,
            noseq,
            unknown,
            emod,
            gpr,
            gprannot,
            no_locus_tag,
        ) = checkModelLocusTags(fmod, fgb)
        print(
            '\ngood: {}\nupdated: {}\nnoseq: {}\nunknown: {}\n'.format(
                len(good), updated, noseq, unknown
            )
        )
        # Now we annotate the model with the genbank sequences
        addSeqAnnotation(
            emod, gpr, good, updated, update_tags=True, no_locus_tag=no_locus_tag
        )
        # if emod.__FBC_VERSION__ < 2:
        # LD['gene2reaction'] = emod.getAllProteinGeneAssociations()
        # LD['reaction2gene'] = emod.getAllGeneProteinAssociations()
        # else:
        LD['gene2reaction'] = emod.getAllProteinGeneAssociations(use_labels=True)
        LD['reaction2gene'] = emod.getAllGeneProteinAssociations(use_labels=True)

        for g_ in LD['gene2reaction']:
            linkDict['__idx__'][g_] = oid

        print(len(emod.gpr))
        # return LD, emod

        # ... and taxon data.
        # TODO update to add MIRIAM hasTaxon identifier
        taxId = ''
        if emod.hasAnnotation('genbank_taxon_id'):
            taxId = emod.getAnnotation('genbank_taxon_id')
            taxIdList.append(taxId)
        LD['taxon_id'] = taxId

        # this is almost "certain" never to happen but just in case .... check for reactions without genes
        for k_ in list(LD['reaction2gene'].keys()):
            if len(LD['reaction2gene'][k_]) == 0:
                print(
                    'INFO: Reaction has no genes associated:',
                    k_,
                    LD['reaction2gene'][k_],
                )
                LD['reaction2gene'].pop(k_)
                # time.sleep(1)

        if not compress_output:
            sbmlfname = os.path.split(fmod)[-1].replace('.xml', '.seqplus.xml')
        else:
            sbmlfname = os.path.split(fmod)[-1].replace('.xml', '.seqplus.gz.xml')
        sbmloutalt = os.path.join(
            model_dir,
            '({})-({}'.format(oid, sbmlfname.replace('.seqplus.', ').seqplus.')),
        )
        # sbmloutalt = '({})-({})-({}'.format(lib_set, oid, sbmlfname.replace('.seqplus.',').seqplus.'))
        if not useV2:
            cbmpy.writeSBML3FBC(
                emod,
                sbmloutalt,
                directory=model_dir,
                add_cobra_annot=add_cobra_annot,
                xoptions={'zip_model': compress_output},
            )
        else:
            cbmpy.writeSBML3FBCV2(
                emod,
                sbmloutalt,
                directory=model_dir,
                add_cobra_annot=add_cobra_annot,
                zip_model=compress_output,
            )
        LD['sbml_out'] = sbmloutalt
        LD['sbml_out_generic'] = sbmloutalt

        Fj = open(
            os.path.join(
                model_dir, sbmloutalt.replace('.seqplus.xml', '.seqplus.json')
            ),
            'w',
        )
        json.dump(linkDict, Fj, indent=1, separators=(',', ': '))
        Fj.close()
        # print(sbmloutalt)
        # print(os.path.join(model_dir, sbmloutalt.replace('.seqplus.xml', '.seqplus.json')))

    return oidList


def createSeqplusModelMetaIdx(fmod, oid, oclass, metadraft_lib_model):
    """
    Create a MetaDraft template model index file from a seqplus model file::

     - *fmod* the seqplus model file e.g. 'EstherDB.xml'
     - *oid* the model unique shortname e.g. 'edb'
     - *class* the model category e.g. 'vu'
     - *metadraft_lib_model* the target MetaDraft lib_model directory

     'edb', 'vu', "<path>"

    """

    import cbmpy

    dmod = cbmpy.readSBML3FBC(fmod)
    fgb = os.path.abspath(fmod)
    input_path, fmod = os.path.split(fgb)
    # dmod.createGeneAssociationsFromAnnotations()
    oclass = oclass.replace('-', '')
    oid = oid.replace('-', '')
    oid = '{}-{}'.format(oclass, oid)
    if fmod.endswith('.xml'):
        fmod = fmod[:-4]
    fmod = '({})-({}).seqplus'.format(oid, fmod)
    print(fmod)

    linkDict = {}
    linkDict[oid] = {}
    linkDict["__idx__"] = {}
    LD = linkDict[oid]
    LD['genbank_in'] = fgb
    LD['sbml_in'] = fgb
    LD['data_path'] = input_path
    LD['gene2reaction'] = dmod.getAllProteinGeneAssociations()
    for g_ in LD['gene2reaction']:
        linkDict['__idx__'][g_] = oid
    LD['reaction2gene'] = dmod.getAllGeneProteinAssociations()
    LD['taxon_id'] = "unknown"
    LD['sbml_out'] = os.path.join(metadraft_lib_model, "{}.xml".format(fmod))
    LD['sbml_out_generic'] = os.path.join(metadraft_lib_model, "{}.xml".format(fmod))
    LD['fasta_out'] = None

    Fj = open(os.path.join(input_path, '{}-link.json'.format(fmod)), 'w')
    json.dump(linkDict, Fj, indent=1, separators=(',', ': '))
    Fj.close()

    cbmpy.writeSBML3FBC(
        dmod,
        os.path.join(input_path, fmod + '.xml'),
        add_cbmpy_annot=True,
        add_cobra_annot=False,
        add_groups=False,
    )


def createCDSdb(list_of_files):
    """
    Creates an Excel workbook with CDS information from a set of genbank files. This is formatted (aka) a "CDS database"
    that is used for proteomics/metagenomics

     - *list_of_files* a real list of genbank files that contain CDS information

    """
    output = {}
    for gbkf in list_of_files:
        output[gbkf] = []
        GBFile = open(gbkf, 'r')
        GBcds = Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
        for cds in GBcds:
            if cds.seq != None:
                data = {
                    'pid': cds.id,
                    'product': cds.annotations['product'],
                    'seq': str(cds.seq),
                }
            output[gbkf].append(data)
        GBFile.close()

    wb = xlwt.Workbook(encoding='utf-8')
    for f in output:
        ws = wb.add_sheet(f)
        for p in range(len(output[f])):
            ws.write(p, 0, output[f][p]['pid'])
            ws.write(p, 1, output[f][p]['product'])
            ws.write(p, 2, output[f][p]['seq'])

    wb.save('{}.xls'.format('cds_db'))
    print('\nCDS database file saved as: cds_db.xls\n')
