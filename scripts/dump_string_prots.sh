#!/usr/bin/zsh

taxo_name="$1"
taxo_id="$2"

out_file_prefix="/opt/shared/dump"

docker-compose exec stringdb sh -c "mkdir -p '$out_file_prefix'"
# mkdir -p "$out_file_prefix"

if [ $# = 2 ]
then
    out_file="$out_file_prefix/$taxo_name-prots.tab"

    docker-compose exec stringdb psql stringdb stringdb -c \
    "\\copy (
        select
            protein_id,
            species_id,
            regexp_replace(protein_external_id, '^[0-9]*\.', '')
                protein_name
        from
            items.proteins
        where
            species_id = $taxo_id
    ) to '$out_file' with csv header delimiter E'\\t';
    "
else
    echo invalid number of arguments
fi

if [ -f "$out_file" ]
then
    docker-compose exec stringdb chmod a+rw "$out_file"
fi
