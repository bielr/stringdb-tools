from xml.etree import ElementTree
from urllib.request import urlopen
from gzip import GzipFile
import stringdb
import sys
import pdb




def create_column(cursor):
    cursor.execute("""
        alter table mapping.uniprot
        add column swiss_prot boolean default false;
        """);

def parse_xml(cursor, f):
    xml = iter(ElementTree.iterparse(f, events=('start', 'end')))
    _, root = next(xml)

    for event, elem in xml:
        if event != 'end' and elem.tag == '{http://uniprot.org/uniprot}entry':
            uniprot_ac = elem.find('{http://uniprot.org/uniprot}accession')

            if uniprot_ac is not None:
                print(f'set {uniprot_ac.text}')

                cursor.execute("""
                    update mapping.uniprot set swiss_prot = true where uniprot_ac = %(uniprot_ac)s;
                    """,
                    {'uniprot_ac': uniprot_ac.text})

            elem.clear()
            root.clear()

if __name__ == '__main__':
    with stringdb.connect_to_docker() as conn:
        with conn.cursor() as cursor:
            print('connected')

            create_column(cursor)

            if len(sys.argv) == 1:
                with urlopen('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz') as ftp_conn:
                    with GzipFile(fileobj=ftp_conn) as decompressed_f:
                        parse_xml(cursor, decompressed_f)

            elif sys.argv[1] == '-':
                parse_xml(cursor, sys.stdin)

            else:
                parse_xml(cursor, sys.argv[1])

            conn.commit()
