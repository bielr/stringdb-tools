#!/bin/bash

set -e

docker-compose exec -T stringdb \
    psql -1 stringdb stringdb -c "
        create schema if not exists mapping;

        create table mapping.tmp_uniprot (
            species_id  integer  not null,
            uniprot_ac  text     not null,
            uniprot_id  text     not null,
            string_name text     not null,
            identity    real     not null,
            bit_score   real     not null
        );"

wget 'https://version-10-5.string-db.org/mapping_files/uniprot_mappings/full_uniprot_2_string.04_2015.tsv.gz' -O- --tries=0 | \
    gzip -d |                                                                             \
    tail -n +2 |                                                                          \
    sed -e 's/|/\t/' |                                                                    \
                                                                                          \
    docker-compose exec -T stringdb                                                       \
        psql -1 stringdb stringdb -c                                                      \
        "\\copy mapping.tmp_uniprot from STDIN delimiter E'\\t' csv;"

docker-compose exec -T stringdb \
    psql -1 stringdb stringdb -c "
        create table mapping.uniprot as
          select
            tmp.species_id       species_id,
            prots.protein_id     string_id,
            tmp.string_name      string_name,
            tmp.uniprot_ac       uniprot_ac,
            tmp.uniprot_id       uniprot_id,
            tmp.identity         \"identity\",
            tmp.bit_score        bit_score
          from
            mapping.tmp_uniprot tmp
            left join
              items.proteins prots
            on
              prots.species_id = tmp.species_id
              and
              prots.protein_external_id = tmp.species_id || '.' || tmp.string_name
          order by
            tmp.species_id, prots.protein_id, tmp.uniprot_ac;

        drop table mapping.tmp_uniprot;

        create index
          si_uniprot_string_id
        on
          mapping.uniprot
        using
          brin (string_id)
        with
          (pages_per_range = 4);

        create index
          si_uniprot_species_id
        on
          mapping.uniprot
        using
          brin (species_id)
        with
          (pages_per_range = 128);

        create index
          si_uniprot_uniprot_ac
        on
          mapping.uniprot
        using
          btree (uniprot_ac)
        with
          (fillfactor = 100);"
