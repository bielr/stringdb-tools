#!/bin/bash

set -e

docker-compose exec -T stringdb \
    psql -1 stringdb stringdb -c "
        create schema if not exists mapping;

        create table mapping.tmp_gene_ontology_full (
            species_id       integer  not null,
            string_name      text     not null,
            common_name      text     not null,
            go_id            text     not null,
            go_name          text     not null,
            evidence_source  text     not null,
            evidence_code    text     not null,
            confidence_score smallint not null
        );"

wget "https://version-10-5.string-db.org/mapping_files/gene_ontology_mappings/all_go_knowledge_full.tsv.gz" -O- --tries=0 | \
    gzip -d |                                                                                                      \
    docker-compose exec -T stringdb                                                                                \
        psql -1 stringdb stringdb -c                                                                               \
        "\\copy mapping.tmp_gene_ontology_full from STDIN delimiter E'\\t' csv;"

docker-compose exec -T stringdb \
    psql -1 stringdb stringdb -c "
        create table mapping.gene_ontology_full as
          select
            tmp.species_id       species_id,
            prots.protein_id     string_id,
            tmp.string_name      string_name,
            tmp.common_name      common_name,
            tmp.go_id            go_id,
            tmp.go_name          go_name,
            tmp.evidence_source  evidence_source,
            tmp.evidence_code    evidence_code,
            tmp.confidence_score confidence_score
          from
            mapping.tmp_gene_ontology_full tmp
            left join
              items.proteins prots
            on
              prots.species_id = tmp.species_id
              and
              prots.protein_external_id = tmp.species_id || '.' || tmp.string_name
          order by
            tmp.species_id, prots.protein_id, tmp.go_id;

        drop table mapping.tmp_gene_ontology_full;

        create index
          si_gene_ontology_full_string_id
        on
          mapping.gene_ontology_full
        using
          brin (string_id)
        with
          (pages_per_range = 4);

        create index
          si_gene_ontology_full_species_id
        on
          mapping.gene_ontology_full
        using
          brin (species_id)
        with
          (pages_per_range = 128);"
