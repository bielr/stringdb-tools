#!/bin/bash

set -e

script_dir="$(dirname "$0")"
stringdb_dir="$script_dir/.."

source "$script_dir/python3.6.5.venv/bin/activate"

# Create tmp_gene_ontology from all_go_knowledge_explicit

docker-compose exec -T stringdb \
    psql -1 stringdb stringdb -c "
        create schema if not exists mapping;

        create table mapping.tmp_gene_ontology (
            species_id       integer  not null,
            string_name      text     not null,
            common_name      text     not null,
            go_id            text     not null,
            go_name          text     not null,
            evidence_source  text     not null,
            evidence_code    text     not null,
            confidence_score smallint not null
        );"

wget "https://string-db.org/mapping_files/gene_ontology_mappings/all_go_knowledge_explicit.tsv.gz" -O- --tries=0 | \
    gzip -d |                                                                                                      \
    docker-compose exec -T stringdb                                                                                \
        psql -1 stringdb stringdb -c                                                                               \
        "\\copy mapping.tmp_gene_ontology from STDIN delimiter E'\\t' csv;"


# Create gene_ontology_update from go_tools.py (alternatives) of distinct go_id's in tmp_gene_ontology

docker-compose exec -T stringdb \
    psql stringdb stringdb -c "
        create table mapping.gene_ontology_update (
            go_old_id text,
            go_new_id text
        );"

docker-compose exec -T stringdb                                                                         \
    psql stringdb stringdb -c                                                                           \
    "\\copy (
        select
          distinct(go_id)
        from
          mapping.tmp_gene_ontology
    ) to STDOUT with csv delimiter E'\\t';
    " |                                                                                                 \
    PYTHONPATH="$stringdb_dir/../semantic_similarity:$stringdb_dir/../geneontology/scripts:$PYTHONPATH" \
        python "$script_dir/go_tools.py" alternatives - - |                                             \
    docker-compose exec -T stringdb                                                                     \
        psql stringdb stringdb -c                                                                       \
        "\\copy mapping.gene_ontology_update from STDIN delimiter E'\\t' csv;"


# Merge tmp_gene_ontology and gene_ontology_update, then create indices

docker-compose exec -T stringdb \
    psql -1 stringdb stringdb -c "
        create table mapping.gene_ontology as
          select
            tmp.species_id       species_id,
            prots.protein_id     string_id,
            tmp.string_name      string_name,
            tmp.common_name      common_name,
            tmp.go_id            go_old_id,
            tmp.go_name          go_old_name,
            upd.go_new_id        go_id,
            tmp.evidence_source  evidence_source,
            tmp.evidence_code    evidence_code,
            tmp.confidence_score confidence_score
          from
            mapping.tmp_gene_ontology tmp
            inner join
              mapping.gene_ontology_update upd
            on
              upd.go_old_id = tmp.go_id
            left join
              items.proteins prots
            on
              prots.species_id = tmp.species_id
              and
              prots.protein_external_id = tmp.species_id || '.' || tmp.string_name
          order by
            tmp.species_id, prots.protein_id, upd.go_new_id;

        drop table mapping.tmp_gene_ontology;
        drop table mapping.gene_ontology_update;

        create index
          si_gene_ontology_string_id
        on
          mapping.gene_ontology
        using
          brin (string_id)
        with
          (pages_per_range = 4);

        create index
          si_gene_ontology_species_id
        on
          mapping.gene_ontology
        using
          brin (species_id)
        with
          (pages_per_range = 128);"

deactivate
