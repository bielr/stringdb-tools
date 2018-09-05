#!/bin/bash

script_dir="$(dirname "$0")"

wget https://stringdb-static.org/download/network_schema.v10.5.sql.gz -O- --tries=0 | \
    gzip -d |                                                                         \
    "$script_dir/import_stringdb.sh" - &&                                             \
                                                                                      \
    docker-compose exec -T stringdb                                                   \
        psql -1 stringdb stringdb -c "
        create index
          si_node_node_links_node_type_b
        on
          network.node_node_links
        using
          brin (node_type_b)
        with
          (pages_per_range = 128);"
