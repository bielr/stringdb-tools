#!/bin/zsh

die() {
    echo "$@"
    exit 1
}

for dump_name in "$@"
do
    if [ "$dump_name" = - ]
    then
        dump_path=/dev/stdin
    else
        script_dir="$(dirname "$(readlink -f "$0")")"
        dump_path="$script_dir/../$dump_name"
        test -f "$dump_path" || die "error: $dump_path does not exist"
    fi

    echo "importing..."

    docker-compose exec -T stringdb \
        psql -1 stringdb stringdb < "$dump_path"
done
