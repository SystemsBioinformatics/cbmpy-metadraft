import time

def gene_report_template(ver, ngene, ifile, db, nmg, selg, nselg, annotG, cp):
    return """
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>MetaDraft Gene Report</title>
</head>
<body>
 <h2 id="page-top">MetaDraft Gene Report</h2>
 <h3>Analysis</h3>
 <p>
 <strong>Input fasta:</strong> {}<br/>
 <strong>Metaproteome used:</strong> {}<br/>
 <strong>Report date:</strong> {}<br/>
 <strong>MetaDraft version:</strong> {}
 </p>
 <!--<h3>Notes</h3>
 <p></p>-->
 <h3>Genes ({})</h3>
 <h4>Input genes not matched</h4>
 <p>{}</p>
 <h4>Genes selected</h4>
 <table width="80%" border="0" cellspacing="1" cellpadding="1">
   <tbody>{}</tbody>
 </table>
 <h4>Genes not selected</h4>
 <table width="80%" border="0" cellspacing="1" cellpadding="1">
   <tbody>{}</tbody>
</table>
 <h3>Annotation</h3>
 {}
 <p>&nbsp;</p>
 {}
</body>
</html>""".format(ifile, db, time.strftime('%y-%m-%d'), ver, ngene, nmg, selg, nselg, annotG, cp)

def reaction_report_template(ver, nreact, ifile, db, selr, nselr, annot, cp):
    return """
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>MetaDraft Reaction Report</title>
</head>
<body>
 <h2 id="page-top">MetaDraft Reaction Report</h2>
 <h3>Analysis</h3>
 <p>
 <strong>Input fasta:</strong> {}<br/>
 <strong>Metaproteome used:</strong> {}<br/>
 <strong>Report date:</strong> {}<br/>
 <strong>MetaDraft version:</strong> {}
 </p>
 <!--<h3>Notes</h3>
 <p></p>-->
 <h3>Reactions ({})</h3>
 <h4>Reactions selected</h4>
 <table width="80%" border="0" cellspacing="1" cellpadding="1">
   <tbody>{}</tbody>
 </table>
 <h4>Reactions not selected</h4>
 <table width="80%" border="0" cellspacing="1" cellpadding="1">
   <tbody>{}</tbody>
</table>
 <h3>Annotation</h3>
 {}
 <p>&nbsp;</p>
 {}
</body>
</html>""".format(ifile, db, time.strftime('%y-%m-%d'), ver, nreact, selr, nselr, annot, cp)

def metabolite_report_template(ver, nmetab, ifile, db, selm, nselm, annot, cp):
    return """
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>MetaDraft Metabolite Report</title>
</head>
<body>
 <h2 id="page-top">MetaDraft Metabolite Report</h2>
 <h3>Analysis</h3>
 <p>
 <strong>Input fasta:</strong> {}<br/>
 <strong>Metaproteome used:</strong> {}<br/>
 <strong>Report date:</strong> {}<br/>
 <strong>MetaDraft version:</strong> {}
 </p>
 <!--<h3>Notes</h3>
 <p></p>-->
 <h3>Metabolites ({})</h3>
 <h4>Metabolites selected</h4>
 <table width="80%" border="0" cellspacing="1" cellpadding="1">
   <tbody>{}</tbody>
 </table>
 <h3>Annotation</h3>
 {}
 <p>&nbsp;</p>
 {}
</body>
</html>""".format(ifile, db, time.strftime('%y-%m-%d'), ver, nmetab, selm, annot, cp)

def combine_index_template(sbmlf, xlf):
    return """
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>MetaDraft Report</title>
</head>
 <body>
  <h3 id="section-reports">Model Reports</h3>
  <ol>
   <li><a href="./1_summary_report.html">Summary report</a></li>
   <li><a href="./2_gene_report.html">Gene report</a></li>
   <li><a href="./3_reaction_report.html">Reaction report</a></li>
   <li><a href="./4_metabolite_report.html">Metabolite report</a></li>
  </ol>
  <h3 id="section-files">Included Files</h3>
  <ul>
   <li><a href="./{}">SBML 3 FBC</a> version of the reconstruction</li>
   <li><a href="./{}">Excel</a> format summary of the reconstruction</li>
  </ul>
 </body>
</html>
""".format(sbmlf, xlf)

def combine_manifest_file(sbmlf, xlf):
    MFstr = ''
    MFstr += '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">\n'
    MFstr += ' <content location="." format="http://identifiers.org/combine.specifications/omex"/>\n'
    MFstr += ' <content location="./metadata.rdf" format="http://identifiers.org/combine.specifications/omex-metadata"/>\n'
    MFstr += ' <content location="./{}" format="http://identifiers.org/combine.specifications/sbml.level-3.version-1"/>\n'.format(sbmlf)
    MFstr += ' <content location="./{}" format="http://mediatypes.appspot.com/application/vnd.ms-excel"/>\n'.format(xlf)
    MFstr += ' <content location="./index.html" format="text/html"/>\n'
    MFstr += ' <content location="./1_summary_report.html" format="text/html"/>\n'
    MFstr += ' <content location="./2_gene_report.html" format="text/html"/>\n'
    MFstr += ' <content location="./3_reaction_report.html" format="text/html"/>\n'
    MFstr += ' <content location="./4_metabolite_report.html" format="text/html"/>\n'
    MFstr = '<?xml version="1.0" encoding="utf-8"?>\n{}\n</omexManifest>\n'.format(MFstr)
    return MFstr

def combine_metadata_file(vc_given='MetaDraft', vc_family='Software', vc_email='', vc_org='', scTime=None):
    if scTime is None:
        scTime = time.strftime('%y-%m-%d %H:%M')
    MDstr = '<?xml version="1.0" encoding="UTF-8"?>\n'
    MDstr += '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n'
    MDstr += '    xmlns:dcterms="http://purl.org/dc/terms/"\n'
    MDstr += '    xmlns:vCard="http://www.w3.org/2006/vcard/ns#"\n'
    MDstr += '    xmlns:bqmodel="http://biomodels.net/models-qualifiers">\n'
    MDstr += ' <rdf:Description rdf:about=".">\n'
    MDstr += ' <dcterms:creator>\n'
    MDstr += ' <rdf:Bag>\n'
    MDstr += '  <rdf:li rdf:parseType="Resource">\n'
    MDstr += '   <vCard:hasName rdf:parseType="Resource">\n'
    MDstr += '    <vCard:family-name>{}</vCard:family-name>\n'.format(vc_family)
    MDstr += '    <vCard:given-name>{}</vCard:given-name>\n'.format(vc_given)
    MDstr += '   </vCard:hasName>\n'
    MDstr += '   <vCard:hasEmail rdf:resource="{}" />\n'.format(vc_email)
    MDstr += '   <vCard:organization-name>\n'
    MDstr += '      {}\n'.format(vc_org)
    MDstr += '   </vCard:organization-name>\n'
    MDstr += '  </rdf:li>\n'
    MDstr += ' </rdf:Bag>\n'
    MDstr += ' </dcterms:creator>\n'
    MDstr += '   <dcterms:created rdf:parseType="Resource">\n'
    MDstr += '    <dcterms:W3CDTF>{}</dcterms:W3CDTF>\n'.format(scTime)
    MDstr += '   </dcterms:created>\n'
    MDstr += '   <dcterms:modified rdf:parseType="Resource">\n'
    MDstr += '    <dcterms:W3CDTF>{}</dcterms:W3CDTF>\n'.format(scTime)
    MDstr += '   </dcterms:modified>\n'
    MDstr += ' </rdf:Description>\n'
    MDstr += '</rdf:RDF> \n'
    return MDstr
