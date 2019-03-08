#!/bin/bash

set -e

die() {
    echo "$@"
    exit 1
}

(python --version | grep -q 'Python 3') || die "Python 3 environment required"

script_dir="$(dirname "$0")"

"$script_dir/import_stringdb_items.sh"
"$script_dir/import_stringdb_network.sh"
"$script_dir/import_stringdb_homology.sh"

"$script_dir/import_all_go_knowledge_explicit.sh"
"$script_dir/import_full_uniprot_2_string.sh"
python "$script_dir/import_swissprot.py"
