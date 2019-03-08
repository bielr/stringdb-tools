import csv
import pandas as pd

import stringdb


def connect_to_localhost(*, dbname='stringdb_virus', user='helen'):
    return stringdb.connect_to_localhost(dbname=dbname, user=user)

def connect_through_docker_network(*, dbname='stringdb_virus', user='helen'):
    return stringdb.connect_through_docker_network(dbname=dbname, user=user)

def connect_to_docker(*, dbname='stringdb_virus', project='stringdb-virus', service='stringdb-virus'):
    return stringdb.connect_to_docker(dbname=dbname, project=project, service=service)


# species-species interactions only
def get_species_interactions(cursor, species_ids, score_types):
    cursor.execute("""
        select
          node_type_a,
          node_type_b,
          node_id_a,
          node_id_b,
          evidence_scores[i][1] score_type,
          evidence_scores[i][2] evidence_score
        from (
          select
            node_type_a,
            node_type_b,
            node_id_a,
            node_id_b,
            evidence_scores,
            generate_subscripts(evidence_scores, 1) as i
          from
            network.node_node_links
          where
            node_type_a in %(species_ids)s
            and
            node_type_a = node_type_b
        ) as indexed_scores
        where
          evidence_scores[i][1] in %(score_types)s;
        """,
        {'species_ids': tuple(species_ids),
         'score_types': tuple(score_types)})

    return pd.DataFrame(
            cursor.fetchall(),
            columns=['node_type_a', 'node_type_b', 'node_id_a', 'node_id_b', 'score_type', 'evidence_score'])


# species-species interactions only
def get_protein_interactions(cursor, species_id, protein_ids, score_types):
    if protein_ids:
        cursor.execute("""
            select
              node_type_a,
              node_type_b,
              node_id_a,
              node_id_b,
              evidence_scores[i][1] score_type,
              evidence_scores[i][2] evidence_score
            from (
              select
                node_type_a,
                node_type_b,
                node_id_a,
                node_id_b,
                evidence_scores,
                generate_subscripts(evidence_scores, 1) as i
              from
                network.node_node_links
              where
                node_type_b = %(species_id)s
                and
                node_type_a = %(species_id)s
                and
                node_id_a in %(protein_ids)s
                and
                node_id_b in %(protein_ids)s
            ) as indexed_scores
            where
              evidence_scores[i][1] in %(score_types)s;
            """,
            {'species_id': species_id,
             'protein_ids': tuple(protein_ids),
             'score_types': tuple(score_types)})

        records = cursor.fetchall()

    else:
        records = []

    return pd.DataFrame(
            records,
            columns=['node_type_a', 'node_type_b', 'node_id_a', 'node_id_b', 'score_type', 'evidence_score'])


# host-host not included
# virus-virus not included
# virus-host included
# host-virus included
def get_all_virus_host_interactions(cursor, host_species_id, score_types):
    cursor.execute("""
        select
          node_type_a,
          node_type_b,
          node_id_a,
          node_id_b,
          evidence_scores[i][1] score_type,
          evidence_scores[i][2] evidence_score
        from (
          select
            node_type_a,
            node_type_b,
            node_id_a,
            node_id_b,
            evidence_scores,
            generate_subscripts(evidence_scores, 1) as i
          from
            network.node_node_links
          where
            ( node_type_b = %(host_species_id)s
              or
              node_type_a = %(host_species_id)s )
            and
            node_type_a != node_type_b
        ) as indexed_scores
        where
          evidence_scores[i][1] in %(score_types)s;
        """,
        {'host_species_id': host_species_id,
         'score_types': tuple(score_types)})

    return pd.DataFrame(
            cursor.fetchall(),
            columns=['node_type_a', 'node_type_b', 'node_id_a', 'node_id_b', 'score_type', 'evidence_score'])


def get_viruses_for_host(cursor, host_species_id, score_types):
    cursor.execute("""
        select distinct(node_type_a) node_type from
          network.node_node_links
        where
          node_type_b = %(host_species_id)s
          and
          node_type_a != node_type_b

        union

        select distinct(node_type_b) node_type from
          network.node_node_links
        where
          node_type_a = %(host_species_id)s
          and
          node_type_a != node_type_b
        """,
        {'host_species_id': host_species_id})

    return [vid for vid, in cursor.fetchall()]


# host-host not included
# virus-virus not included
# virus-host included
# host-virus included
def get_virus_host_interactions(cursor, host_species_id, virus_species_id, score_types):
    cursor.execute("""
        select
          node_type_a,
          node_type_b,
          node_id_a,
          node_id_b,
          evidence_scores[i][1] score_type,
          evidence_scores[i][2] evidence_score
        from (
          select
            node_type_a,
            node_type_b,
            node_id_a,
            node_id_b,
            evidence_scores,
            generate_subscripts(evidence_scores, 1) as i
          from
            network.node_node_links
          where
            ( node_type_a = %(host_species_id)s and node_type_b = %(virus_species_id)s )
            or
            ( node_type_b = %(host_species_id)s and node_type_a = %(virus_species_id)s )
        ) as indexed_scores
        where
          evidence_scores[i][1] in %(score_types)s;
        """,
        {'host_species_id': host_species_id,
         'virus_species_id': virus_species_id,
         'score_types': tuple(score_types)})

    return pd.DataFrame(
            cursor.fetchall(),
            columns=['node_type_a', 'node_type_b', 'node_id_a', 'node_id_b', 'score_type', 'evidence_score'])


# virus-host included
# host-virus included
# host-host included (virus-host neighbors only)
# virus-virus included (virus-host neighbors only)
def get_all_virus_host_networks(cursor, host_species_id, score_types):
    vh = get_all_virus_host_interactions(cursor, host_species_id, score_types)
    virus_species_ids = frozenset(vh.node_type_a) | frozenset(vh.node_type_b)
    virus_species_ids = tuple(virus_species_ids - {host_species_id})
    vv = get_species_interactions(cursor, virus_species_ids, score_types)

    return pd.concat([vh, vv])


# virus-host included
# host-virus included
# host-host included (virus-host neighbors only)
# virus-virus included (virus-host neighbors only)
def get_virus_host_network(cursor, host_species_id, virus_species_id, score_types):
    vh = get_virus_host_interactions(cursor, host_species_id, virus_species_id, score_types)

    node_ids = pd.concat([
        vh.loc[:, ['node_type_a','node_id_a']].rename(columns={'node_type_a':'node_type', 'node_id_a':'node_id'}),
        vh.loc[:, ['node_type_b','node_id_b']].rename(columns={'node_type_b':'node_type', 'node_id_b':'node_id'})
    ])

    vv = get_protein_interactions(
            cursor,
            virus_species_id,
            node_ids.loc[node_ids.node_id == virus_species_id].node_id.to_numpy(),
            score_types)

    hh = get_protein_interactions(
            cursor,
            host_species_id,
            node_ids.loc[node_ids.node_id == host_species_id].node_id.to_numpy(),
            score_types)

    return pd.concat([vh,vv,hh])


# maps node_id_a from string id's to string external id's
#      node_id_b from string id's to string external id's
def network_with_external_ids(cursor, network):
    orig_columns = network.columns

    external_ids = stringdb.get_prots_external_ids(cursor, frozenset(network.node_id_a) | frozenset(network.node_id_b))
    external_ids = external_ids.set_index('string_id')

    network = pd.concat([
            network.set_index('node_id_a'),
            external_ids.rename(columns={'external_id': 'node_id_a'})
        ],
        axis=1, join='inner')

    network = pd.concat([
            network.set_index('node_id_b'),
            external_ids.rename(columns={'external_id': 'node_id_b'})
        ],
        axis=1, join='inner')

    return network.loc[:, orig_columns]


def get_species_official_names(cursor, species_ids=None):
    if species_id is None:
        cursor.execute("""
            select
              species_id, official_name
            from
              items.species;
            """)
    else:
        cursor.execute("""
            select
              species_id, official_name
            from
              items.species
            where
              species_id in %(species_ids)s;
            """,
            {'species_ids': tuple(species_ids)})

    return pd.DataFrame(
        cursor.fetchall(),
        columns=['species_id', 'official_name'])
