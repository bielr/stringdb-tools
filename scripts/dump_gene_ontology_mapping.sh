#!/usr/bin/zsh

taxo_name="$1"
taxo_id="$2"

out_file_prefix="/opt/shared/dump"

docker-compose exec stringdb sh -c "mkdir -p '$out_file_prefix'"
# mkdir -p "$out_file_prefix"

if [ $# = 2 ]
then
    out_file="$out_file_prefix/$taxo_name-go-mapping.tab"

    docker-compose exec stringdb psql stringdb stringdb -c \
    "\\copy (
        select
          string_id,
          string_name,
          go_id,
          evidence_source,
          evidence_code,
          confidence_score
        from
          mapping.gene_ontology
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
