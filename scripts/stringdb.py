import csv
import itertools
import networkx as nx
#import obonet
import pronto
import sys


def connect_to_localhost():
    import psycopg2
    return psycopg2.connect(host='localhost', port=5432, user='stringdb', password='stringdb', dbname='stringdb')

def connect_to_docker():
    import psycopg2
    import docker

    docker_client = docker.from_env()

    containers = docker_client.containers.list(
        filters = {'label': [
            'com.docker.compose.project=stringdb',
            'com.docker.compose.service=stringdb'
        ]})

    assert len(containers) == 1

    container = containers[0]
    network_settings = container.attrs['NetworkSettings']

    host = network_settings['Networks']['stringdb-net']['IPAddress']
    env = dict(e.split('=', 1) for e in container.attrs['Config']['Env'])

    return psycopg2.connect(host=host, port=5432, user=env['POSTGRES_USER'], password=env['POSTGRES_PASSWORD'], dbname='stringdb')


def get_explicit_annotations(cursor, go_id, evidence_codes):
    cursor.execute("""
        select
          string_id
        from
          mapping.gene_ontology
        where
          go_id = %(go_id)s
          and
          evidence_code in %(evidence_codes)s;
        """,
        {'go': go_id, 'evidence_codes': evidence_codes})

    return [string_id for string_id, in cursor.fetchall()]

def count_annotations(cursor, go_is_a_g, go_id, evidence_codes):
    gos = tuple(nx.descendants(go_is_a_g, go_id)) + (go_id,)

    cursor.execute("""
        select
          count(distinct string_id)
        from
          mapping.gene_ontology
        where
          go_id in %(gos)s
          and
          evidence_code in %(evidence_codes)s;
        """,
        {'gos': gos, 'evidence_codes': evidence_codes})

    [(cnt,)] = cursor.fetchall()
    return cnt

def get_explicit_prot_annotations(cursor, prot_string_id, evidence_codes):
    cursor.execute("""
        select
          go_id
        from
          mapping.gene_ontology
        where
          string_id = %(string_id)s
          and
          evidence_code in %(evidence_codes)s;
        """,
        {'string_id': prot_string_id, 'evidence_codes': evidence_codes})

    return [go_id for go_id, in cursor.fetchall()]

def get_prot_annotations(cursor, go_is_a_g, prot_string_id, evidence_codes):
    explicit_anns = get_explicit_prot_annotations(cursor, prot_string_id, evidence_codes)
    anns = set(explicit_anns)

    for go_id in explicit_anns:
        anns.update(nx.descendants(go_is_a_g, go_id))

    return anns

def get_species_names(cursor):
    cursor.execute("""
        select
          species_id,
          official_name
        from
          items.species
        """)

    return cursor.fetchall()

def get_species_prots(cursor, species_id):
    cursor.execute("""
        select
          protein_id
        from
          items.proteins
        where
          species_id = %(species_id)s
        """,
        {'species_id': species_id})

    return [string_id for string_id, in cursor.fetchall()]

def get_species_prots_uniprot(cursor, species_id, min_identity):
    cursor.execute("""
        select
          protein_id string_id,
          array(
            select
              uniprot_ac
            from
              mapping.uniprot
            where
              species_id = %(species_id)s
              and
              string_id = protein_id
              and
              identity >= %(min_identity)s)
            uniprot_ac
        from
          items.proteins
        where
          species_id = %(species_id)s;
        """,
        {'species_id': species_id, 'min_identity': min_identity})

    return cursor.fetchall()

def get_uniprot_species_id(cursor, species_id, uniprot_acs):
    cursor.execute("""
        select
          uniprot_ac,
          string_id
        from
          mapping.uniprot
        where
          species_id = %(species_id)s
          and
          uniprot_ac in %(uniprot_acs)s;
        """,
        {'species_id': species_id, 'uniprot_acs': tuple(uniprot_acs)})

    return cursor.fetchall()
