#!/bin/zsh

script_dir="$(dirname "$0")"

wget https://stringdb-static.org/download/homology_schema.v10.5.sql.gz -O- --tries=0 | \
    gzip -d | \
    "$script_dir/import_stringdb.sh" -

docker-compose exec -T stringdb \
    psql stringdb stringdb -c "
    drop index homology.pi_blast_data;

    alter table homology.blast_data
      drop start_a,
      drop end_a,
      drop start_b,
      drop end_b,
      drop size_b;

    vacuum homology.blast_data;

    alter table homology.blast_data
      add column species_id_a integer;

    update homology.blast_data
    set
      species_id_a = species_id
    from
      items.proteins
    where
      protein_id = protein_id_a;

    vacuum homology.blast_data;

    create index si_blast_data_species
    on
      homology.blast_data
    using
      brin (species_id_a, species_id_b)
    with (pages_per_range=128);

    create index si_blast_data_protein_id_a
    on
      homology.blast_data
    using
      brin (protein_id_a)
    with (pages_per_range=128);

    create index si_blast_data_protein_id_b
    on
      homology.blast_data
    using
      brin (protein_id_b)
    with (pages_per_range=128);
    "
