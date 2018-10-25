from os import path
import csv
import itertools
import math
import networkx as nx
import numpy as np
import sys
import time

import stringdb
import geneontology as godb
import semantic_similarity as semsim


def get_curated_frequencies_path():
    script_dir = path.dirname(path.realpath(__file__))
    stringdb_dir = path.dirname(script_dir)

    # '/data/pinaweb/stringdb/extra/go_curated_frequencies.tab'
    return path.join(stringdb_dir, 'extra', 'go_curated_frequencies.tab')


def load_go_list(fd):
    return [go_id.strip() for go_id in fd.readlines() if not go_id.isspace()]


def init_default_hrss(agg = semsim.agg_bma_max):
    go_onto = godb.load_go_obo()
    go_is_a_g = godb.onto_rel_graph(go_onto)
    ic = semsim.init_ic(go_onto, get_curated_frequencies_path())

    return semsim.HRSS(agg=agg, onto=go_onto, rel_g=go_is_a_g, ic=ic)

def init_default_mica_dissim(agg = semsim.agg_bma_min):
    onto = godb.load_go_obo()
    go_is_a_g = godb.onto_rel_graph(go_onto)
    ic = semsim.init_ic(go_onto, get_curated_frequencies_path())

    return semsim.MICADissim(agg=agg, onto=onto, rel_g=go_is_a_g, ic=ic)


def open_arg_file(path, mode):
    if path == '-':
        if mode == 'r':
            return sys.stdin
        else:
            return sys.stdout
    else:
        return open(path, mode)


def find_valid_go_alternatives(go_list, go_onto, go_is_a_g, evidence_codes):
    go_missing = [go for go in go_list if go not in go_is_a_g]
    assert all(godb.is_obsolete(go, go_onto) for go in go_missing)
    assert all(godb.is_disconnected(go, go_onto) for go in go_missing)

    go_alt_id_g = godb.onto_alt_id_graph(go_onto, go_is_a_g)

    go_no_alt_id       = [go for go in go_missing if list(godb.find_alternatives(go, alt_id_g=go_alt_id_g)) == [go]]
    go_no_valid_alt_id = [go for go in go_missing if not godb.find_valid_alternatives(go, alt_id_g=go_alt_id_g, rel_g=go_is_a_g)]

    for go in go_list:
        for alt_id in godb.find_valid_alternatives(go, alt_id_g=go_alt_id_g, rel_g=go_is_a_g)[:1]:
            yield go, alt_id

    print(f'skipped {len(go_no_valid_alt_id)} terms with no valid alternative id\'s', file=sys.stderr)


def make_annotation_counter(source, go_is_a_g, evidence_codes):
    if source == 'stringdb':
        conn = stringdb.connect_to_docker()
        cursor = conn.cursor()

        return lambda go: \
            stringdb.count_annotations(cursor=cursor, go_is_a_g=go_is_a_g, go_id=go, evidence_codes=evidence_codes)

    elif source == 'geneontology':
        conn = godb.connect_to_docker()
        cursor = conn.cursor()

        return lambda go: \
            godb.count_protein_annotations(cursor, go, evidence_codes=evidence_codes)

    else:
        raise ValueError(f'unknown source db of GO annotations: {source}')


def get_curated_frequencies(source, go_is_a_g, evidence_codes):
    count_annotations = make_annotation_counter(source, go_is_a_g, evidence_codes)

    for i, go in enumerate(go_is_a_g.nodes()):
        go_cnt = count_annotations(go)

        print(f'{i}/{nx.number_of_nodes(go_is_a_g)} count {go} = {go_cnt}', file=sys.stderr)
        yield go, go_cnt


def get_all_annotations_for_species(source, species_id, evidence_codes):
    if source == 'stringdb':
        with stringdb.connect_to_docker() as string_conn:
            with string_conn.cursor() as string_cursor:
                prots = stringdb.get_species_prots(string_cursor, species_id)

                annotations = []

                for i, string_id in enumerate(prots):
                    if i % 100 == 0:
                        print(f'\rretrieving annotations ({i}/{len(prots)})... ', end='', file=sys.stderr)

                    prot_gos = stringdb.get_explicit_prot_annotations(string_cursor, string_id, evidence_codes)

                    if prot_gos:
                        annotations.append((string_id, prot_gos))

                print('done', file=sys.stderr)

    else:
        raise ValueError(f'unknown source db of GO annotations: {source}')

    print(' > found',  len(prots) - len(annotations), 'unannotated proteins', file=sys.stderr)
    ann_lens = np.array(sorted([len(gos) for prot, gos in annotations]))
    print(' > annotation length percentiles (excluding unannotated):', [np.percentile(ann_lens, p) for p in range(0, 101, 10)], file=sys.stderr)

    return annotations


def classify_annotations_by_namespace(annotations, onto):
    for prot, gos in annotations:
        yield (prot,
                [go for go in gos if 'biological_process' in onto[go].other['namespace']],
                [go for go in gos if 'cellular_component' in onto[go].other['namespace']],
                [go for go in gos if 'molecular_function' in onto[go].other['namespace']])

def count_annotations_by_namespace(annotations, onto):
    for prot, bgos, cgos, mgos in classify_annotations_by_namespace(annotations, onto):
        yield prot, len(bgos), len(cgos), len(mgos)

def compute_comparison_matrix(cmpobj, annotations, namespace=None):
    comparison_mat = np.zeros((len(annotations), len(annotations)))

    if namespace is None:
        compare = cmpobj.compare
    else:
        compare = lambda gos1, gos2: cmpobj.compare_for_namespace(namespace, gos1, gos2)

    for i, (prot_id1, gos1) in enumerate(annotations):
        if i % 10 == 0:
            print(f'\rcomputing {namespace or ""} matrix entries ({i}/{len(annotations)} rows)... ', end='', file=sys.stderr)

        for j, (prot_id2, gos2) in enumerate(annotations):
            if i <= j:
                comparison_mat[i,j] = compare(gos1, gos2)
            else:
                comparison_mat[i,j] = comparison_mat[j,i]

    print(f'done', file=sys.stderr)
    return comparison_mat


def compute_comparison_matrices(cmpobj, annotations):
    comparison_mats = np.zeros((3, len(annotations), len(annotations)))
    comparison_mats[0] = compute_comparison_matrix(cmpobj, annotations, 'biological_process')
    comparison_mats[1] = compute_comparison_matrix(cmpobj, annotations, 'cellular_component')
    comparison_mats[2] = compute_comparison_matrix(cmpobj, annotations, 'molecular_function')

    return comparison_mats


def write_annotation_comparison_matrix(annotations, comparison_mat, out_f):
    writer = csv.writer(out_f, delimiter='\t')
    for i in range(len(annotations)):
        for j in range(len(annotations)):
            if comparison_mat[i,j] != 0:
                writer.writerow((annotations[i][0], annotations[j][0], comparison_mat[i,j]))

def write_annotation_comparison_matrices(annotations, comparison_mat, out_f):
    writer = csv.writer(out_f, delimiter='\t')
    writer.writerow(('protein1', 'protein2', 'biological_process', 'cellular_component', 'molecular_function'))

    for i in range(len(annotations)):
        for j in range(len(annotations)):
            if np.any(comparison_mat[:,i,j] != 0):
                writer.writerow((annotations[i][0], annotations[j][0],
                    comparison_mat[0,i,j], comparison_mat[1,i,j], comparison_mat[2,i,j], comparison_mat[3,i,j]))


if __name__ == '__main__':
    cmd = sys.argv[1]

    go_onto = godb.load_go_obo()
    go_is_a_g = godb.onto_rel_graph(go_onto)
    evidence_codes = godb.get_curated_evidence_codes()

    assert len(list(nx.weakly_connected_components(go_is_a_g))) == 3

    if cmd == 'alternatives':
        go_list_path = sys.argv[2]
        output_path = sys.argv[3]

        go_list = load_go_list(open_arg_file(go_list_path, 'r'))
        alternatives = find_valid_go_alternatives(go_list, go_onto, go_is_a_g, evidence_codes)

        writer = csv.writer(open_arg_file(output_path, 'w+'), delimiter='\t')
        writer.writerows(alternatives)

    elif cmd == 'curated-frequencies':
        source = sys.argv[2]
        output_path = sys.argv[3]

        freqs = get_curated_frequencies(source, go_onto, go_is_a_g, evidence_codes)

        writer = csv.writer(open_arg_file(output_path, 'w+'), delimiter='\t')
        writer.writerows(freqs)

    elif cmd == 'namespace-ann-counts':
        source = sys.argv[2]
        species_id = int(sys.argv[3])
        output_path = sys.argv[4]

        annotations = get_all_annotations_for_species(source, species_id, evidence_codes)
        annotations = count_annotations_by_namespace(annotations, go_onto)

        writer = csv.writer(open_arg_file(output_path, 'w+'), delimiter='\t')
        writer.writerow(('protein', 'biological_process', 'cellular_component', 'molecular_function'))
        writer.writerows(annotations)

    elif cmd == 'semsim-matrix':
        measure = sys.argv[2]
        source = sys.argv[3]
        species_id = int(sys.argv[4])
        output_path = sys.argv[5]

        annotations = get_all_annotations_for_species(source, species_id, evidence_codes)

        if measure == 'mica-dissim':
            cmpobj = init_default_mica_dissim()
        elif measure == 'hrss':
            cmpobj = init_default_hrss()
        else:
            raise ValueError(f'unknown (dis)similarity measure: {measure}')

        mats = compute_comparison_matrices(cmpobj, annotations)

        print('writing results...', file=sys.stderr)
        write_annotation_comparison_matrices(annotations, mats, open_arg_file(output_path, 'w+'))
