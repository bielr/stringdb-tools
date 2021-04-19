#!/usr/bin/zsh

id_column="$1"

out_file_prefix="/opt/shared/dump"

docker-compose exec stringdb sh -c "mkdir -p '$out_file_prefix'"

if [ $# = 1 ]
then
    out_file="$out_file_prefix/$id_column-swissprot-list.tab"

    docker-compose exec stringdb psql stringdb stringdb -c \
    "\\copy (
        select
            $id_column
        from
            mapping.uniprot
        group by
            $id_column
        having
            bool_and(swiss_prot)
    ) to '$out_file' with csv header delimiter E'\\t';
    "
else
    echo invalid number of arguments
fi

if [ -f "$out_file" ]
then
    docker-compose exec stringdb chmod a+rw "$out_file"
fi
