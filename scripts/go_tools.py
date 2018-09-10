from os import path
import csv
import itertools
import networkx as nx
import numpy as np
import sys

import stringdb
import geneontology as godb
import semantic_similarity


def get_curated_frequencies_path():
    script_dir = path.dirname(path.realpath(__file__))
    stringdb_dir = path.dirname(script_dir)

    # '/data/pinaweb/stringdb/extra/go_curated_frequencies.tab'
    return path.join(stringdb_dir, 'extra', 'go_curated_frequencies.tab')

def load_go_list(fd):
    return [go_id.strip() for go_id in fd.readlines() if not go_id.isspace()]



def open_arg_file(path, mode):
    if path == '-':
        if mode == 'r':
            return sys.stdin
        else:
            return sys.stdout
    else:
        return open(path, mode)

if __name__ == '__main__':
    cmd = sys.argv[1]

    go_onto = godb.load_go_obo()
    go_is_a_g = godb.onto_rel_graph(go_onto)
    evidence_codes = godb.get_curated_evidence_codes()

    assert len(list(nx.weakly_connected_components(go_is_a_g))) == 3

    if cmd == 'alternatives':
        go_list_path = sys.argv[2]
        output_mapping_path = sys.argv[3]

        go_list = load_go_list(open_arg_file(go_list_path, 'r'))
        go_missing = [go for go in go_list if go not in go_is_a_g]
        assert all(godb.is_obsolete(go, go_onto) for go in go_missing)
        assert all(godb.is_disconnected(go, go_onto) for go in go_missing)

        go_alt_id_g = godb.onto_alt_id_graph(go_onto, go_is_a_g)

        go_no_alt_id       = [go for go in go_missing if list(godb.find_alternatives(go, alt_id_g=go_alt_id_g)) == [go]]
        go_no_valid_alt_id = [go for go in go_missing if not godb.find_valid_alternatives(go, alt_id_g=go_alt_id_g, rel_g=go_is_a_g)]

        writer = csv.writer(open_arg_file(output_mapping_path, 'w+'), delimiter='\t')

        writer.writerows((go, alt_id) for go in go_list
                                      for alt_id in godb.find_valid_alternatives(go, alt_id_g=go_alt_id_g, rel_g=go_is_a_g)[:1])

        print(f'skipped {len(go_no_valid_alt_id)} terms with no valid alternative id\'s', file=sys.stderr)

    elif cmd == 'curated-frequencies':
        counts = []
        source = sys.argv[2]
        output_counts_path = sys.argv[3]

        if source == 'stringdb':
            conn = stringdb.connect()
            cursor = conn.cursor()

            count_annotations = lambda go: \
                stringdb.count_annotations(
                    cursor = cursor,
                    go_is_a_g = go_is_a_g,
                    go_id = go,
                    evidence_codes = evidence_codes)

        elif source == 'geneontology':
            conn = godb.connect_to_docker()
            cursor = conn.cursor()
            count_annotations = lambda go: \
                godb.count_protein_annotations(cursor, go, evidence_codes=evidence_codes)

        else:
            print(f'unknown source db of GO annotations: {source}', file=sys.stderr)
            sys.exit(1)

        for i, go in enumerate(go_is_a_g.nodes()):
            go_cnt = count_annotations(go)

            print(f'{i}/{nx.number_of_nodes(go_is_a_g)} count {go} = {go_cnt}', file=sys.stderr)
            counts.append((go, go_cnt))

        writer = csv.writer(open_arg_file(output_counts_path, 'w+'), delimiter='\t')
        writer.writerows(counts)

    elif cmd == 'dissim-matrix':
        source = sys.argv[2]
        species_id = int(sys.argv[3])
        output_path = sys.argv[4]

        if source == 'stringdb':
            with stringdb.connect() as string_conn:
                with string_conn.cursor() as string_cursor:
                    annotations = [(string_id, []) for string_id in stringdb.get_species_prots(string_cursor, species_id)]
                    unannotated_cnt = 0

                    for i, (string_id, gos) in enumerate(annotations):
                        if i % 100 == 0:
                            print(f'\rretrieving annotations ({i}/{len(annotations)})... ', end='', file=sys.stderr)

                        if not gos:
                            unannotated_cnt += 1
                        else:
                            gos.extend(stringdb.get_prot_annotations(string_cursor, go_is_a_g, string_id, evidence_codes))

                    print(f'found {unannotated_cnt} unannotated proteins', file=sys.stderr)

        ic = semantic_similarity.init_ic(get_curated_frequencies_path())
        pair_dissim = lambda go1,  go2:  semantic_similarity.get_mica_dissim(go_is_a_g, ic, go1, go2)
        agg_dissim  = lambda gos1, gos2: semantic_similarity.agg_bma_min(dissim, gos1, gos2)
        dissim      = lambda gos1, gos2: min(ns_wise_comparisons(go_onto, agg_dissim, gos1, gos2))

        dissim_mat = np.zeros((len(annotations), len(annotations)))

        for i, (string_id1, gos1) in enumerate(annotations):
            for j, (string_id2, gos2) in enumerate(annotations):
                if i <= j:
                    dissim_mat[i,j] = agg_dissim(gos1, gos2)
                else:
                    dissim_mat[i,j] = dissim_mat[j,i]

        with open(output_path, 'w+') as out_f:
            writer = csv.writer(out_f, delimiter='\t')
            for i in range(len(annotations)):
                for j in range(len(annotations)):
                    writer.writerow(annotations[i][0], annotations[j][0], dissim_mat[i,j])

