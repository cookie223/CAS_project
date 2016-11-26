import urllib2
import re

"""
This file parses and downloads sequence and domain markers of Cas9
"""


def parse_data():
    # list_page = urllib2.urlopen('http://www.ebi.ac.uk/interpro/protein/G3ECR1/similar-proteins')
    list_page = urllib2.urlopen('http://www.ebi.ac.uk/interpro/protein/Q07766/similar-proteins') # Cas 5
    protein_pattern = re.compile('href="/interpro/protein/([A-Z0-9]+)"')
    proteins = protein_pattern.findall(list_page.read())

    sequence_base_url = 'http://www.ebi.ac.uk/interpro/protein/{}?export=fasta'
    markers_base_url = 'http://www.ebi.ac.uk/interpro/protein/{}'
    marker_patter = re.compile('''entryAcs=(.*?)&start=([0-9]+)&end=([0-9]+)''')

    with open('data/full_sequence_Cas5/protein_list', mode='w') as f:
        for protein in proteins:
            f.write(protein)
            f.write('\n')

    for protein in proteins:
        with open('data/full_sequence_Cas5/{}.fasta'.format(protein), mode='w') as f:
            f.write(urllib2.urlopen(sequence_base_url.format(protein)).read())
        pattern_web = urllib2.urlopen(markers_base_url.format(protein)).read()
        markers = marker_patter.findall(pattern_web)
        markers.sort(key=lambda x: int(x[1]))
        with open('data/full_sequence_Cas5/{}.domains'.format(protein), mode='w') as f:
            for marker in markers:
                f.write(marker[0])
                f.write(' ')
                f.write(marker[1])
                f.write(' ')
                f.write(marker[2])
                f.write('\n')

if __name__ == '__main__':
    parse_data()
