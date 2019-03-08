#!/bin/bash

script_dir="$(dirname "$0")"

wget https://version-10-5.string-db.org/download/items_schema.v10.5.sql.gz -O- --tries=0 | \
    gzip -d | \
    "$script_dir/import_stringdb.sh" -
