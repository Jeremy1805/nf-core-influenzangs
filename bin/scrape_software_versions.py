#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'nf-core/influenzangs': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'TrimGalore': ['v_trim_galore.txt', r"version (\S+)"],
    'PicardMarkDuplicates': ['v_picard_mark_dup.txt', r"([\d\.]+)-SNAPSHOT"],
    'PicardCreateSequenceDictionary': ['v_picard_dict.txt', r"([\d\.]+)-SNAPSHOT"],
    'GATK4': ['v_gatk.txt', r"The Genome Analysis Toolkit \(GATK\) v(\S+)"],
    'python': ['v_python.txt', r"Python (\S+)"],
    'VarScan': ['v_varscan.txt', r"VarScan v(\S+)"],
    'bcftools': ['v_bcftools.txt', r"bcftools (\S+)"],
    'bedtools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'ClustalW': ['v_clustal.txt', r"CLUSTAL (\S+) Multiple Sequence Alignments"],
    'R': ['v_R.txt', r"R version (\S+) "],
    'pyfaidx': ['v_pyfaidx.txt', r"(\S+)"],
    'BBMap': ['v_bbmap.txt', r"BBMap version (\S+)"],
    'bamutils': ['v_bamutils.txt', r"Version: (\S+);"],
    'bc': ['v_bc.txt', r"bc (\S+)"]
}
results = OrderedDict()
results['nf-core/influenzangs'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/influenzangs Software Versions'
section_href: 'https://github.com/nf-core/influenzangs'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
