library(rlang)
library(igraph)
library(tidyverse)

load_g = function(filename_noext) {
    g_path = sprintf('shared/dump/%s.tab', filename_noext)
    edges = read_tsv(g_path, col_names=TRUE)

    cat("score quantiles:\n")
    quantiles = edges$evidence_score %>% quantile(probs = seq(0, 100, 5)/100)
    print(quantiles)

    edges %>%
        filter(evidence_score >= mean(quantiles[[4]])) %>%
        select(-evidence_score) %>%
        graph_from_data_frame(directed=FALSE) %>%
        igraph::simplify()
}

analyze_g = function(g) {
    res = list(whole=list(), big_comp=list())
    res$whole$vcount = vcount(g)
    res$whole$ecount = ecount(g)

    cat("vcount", res$whole$vcount, "ecount", res$whole$ecount, "\n")

    comps_g = components(g)
    res$whole$component_sizes = comps_g$csize

    cat("connected components:", comps_g$no, "\n")
    print(sort(comps_g$csize))

    dec_g = decompose(g)
    max_i = which.max(lapply(dec_g, vcount))
    big_g = dec_g[[max_i]]

    res$big_comp$vcount = vcount(big_g)
    res$big_comp$ecount = ecount(big_g)
    cat("big component: vcount", res$big_comp$vcount, "ecount", res$big_comp$ecount, "\n")

    res$whole$degree_table = table(degree(g, V(g)))
    res$big_comp$degree_table = table(degree(big_g, V(big_g)))
    plot(res$big_comp$degree_table)

    res$big_comp$diameter = diameter(big_g)
    res$big_comp$radius   = radius(big_g)
    cat("diameter", res$big_comp$diameter, "radius", res$big_comp$radius, "\n")

    art = articulation_points(big_g)
    cat("# articulation points:", length(art), "\n")

    ecc_art = eccentricity(big_g, art)

    cat("eccentricity table", "\n")
    print(table(ecc_art))

    min_ecc = min(ecc_art)

    for (min_ecc_node in names(ecc_art)[ecc_art == min_ecc]) {
        cat("components removing", min_ecc_node, "\n")
        print(components(delete_vertices(big_g, min_ecc_node))$csize)
    }

    list(big_comp_g=big_g, stats=res)
}

reindex_mapping = function(mapping) {
    match(mapping, unique(mapping))
}

subgraph_quot = function(g, subg_vs) {
    mapping = 1:vcount(g)
    found = FALSE
    g_vs_name = V(g)$name
    subg_vs_name = subg_vs$name

    for (i in 1:vcount(g)) {
        if (g_vs_name[[i]] %in% subg_vs_name) {
            if (!found) {
                collapse_idx = i
                found = TRUE
            }
            mapping[[i]] = collapse_idx
        } else {
            mapping[[i]] = i
        }
    }

    qg = contract.vertices(g, mapping=reindex_mapping(mapping), vertex.attr.comb='first')
    simplify(qg)
}

process_mcl_clusters = function(big_comp, mcl) {
    membership = reindex_mapping(mcl$Cluster)
    names(membership) = V(big_comp)$name

    cluster_names = sorted(unique(unname(membership)))
    clusters = lapply(cluster_names, function(cluster_id) V(big_comp)[mcl$Cluster == cluster_id])
    names(clusters) = cluster_names

    list(
        mcl = mcl,
        clusters = clusters,
        membership = mcl$Cluster+1)
}

get_mcl_clusters = function(big_comp) {
    mcl = MCL::mcl(big_comp, addLoops=FALSE)

    process_mcl_clusters(big_comp, mcl)
}

process_fastgreedy_clusters = function(big_comp, fg) {
    cluster_names = sort(unique(membership(fg)))
    clusters = lapply(cluster_names, function(cluster_id) V(big_comp)[membership(fg) == cluster_id])
    names(clusters) = cluster_names

    list(fg = fg,
         clusters = clusters,
         membership = membership(fg))
}

get_fastgreedy_clusters = function(big_comp) {
    fg = cluster_fast_greedy(big_comp)

    process_fastgreedy_clusters(big_comp, fg)
}

index_to_name = function(lst) {
    names(lst) = 1:length(lst)
    lst
}

clusters_to_df = function(clusters) {
    data_frame(cluster_name = names(clusters),
               cluster = clusters)
}

unlist_df = function(dfs, idx_name) {
    is = if (length(names(dfs)) == 0) 1:length(dfs) else names(dfs)

    idx_name = quo_name(enquo(idx_name))
    bind_rows(lapply(is, function(i) mutate(dfs[[i]], !! idx_name := i)))
}

tidy_clustering = function(clusters, col) {
    col = quo_name(enquo(col))

    lapply(clusters, function(cluster) data_frame(!!col := cluster$name)) %>%
        unlist_df(cluster_id) %>%
        select(cluster_id, !! col)
}

intersect_clusters = function(clusters1, clusters2) {
    tidyr::crossing(
        tibble(cluster_name1 = names(clusters1)),
        tibble(cluster_name2 = names(clusters2))) %>%
      rowwise() %>%
      mutate(csize1      = length(clusters1[[cluster_name1]]),
             csize2      = length(clusters2[[cluster_name2]]),
             inters_size = length(intersect(clusters1[[cluster_name1]],
                                            clusters2[[cluster_name2]]))) %>%
      ungroup() %>%
      filter(inters_size > 0 & inters_size < pmin(csize1, csize2)) %>%
      arrange(desc(inters_size))
      # mutate(prop1 = inters_size/csize1,
      #        prop2 = inters_size/csize2,
      #        prop = inters_size / pmin(csize1, csize2)) %>%
      # arrange(prop)
}

cluster_adjacency = function(g, clust) {
    adjs = lapply(clust$clusters, function(vs) adjacent_vertices(g, vs))

    adj_df = list()

    for (clust_name in names(adjs)) {
        adj = adjs[[clust_name]]

        adj_df[[clust_name]] =
            data_frame(src_clust=clust_name, src=names(adj), tgt=adj) %>%
                unnest(tgt = lapply(tgt, names)) %>%
                mutate(tgt_clust = as.character(clust$membership[tgt])) %>%
                select(src_clust, tgt_clust, src, tgt) %>%
                filter(src_clust != tgt_clust)
    }

    bind_rows(adj_df)
}


load_uniprot_mapping = function(species_name) {
    read_tsv(sprintf('shared/dump/%s-uniprot-mapping.tab', species_name), col_names=TRUE)
}

uniprot_top_matches = function(uniprot_mapping) {
    uniprot_mapping %>%
        group_by(string_id) %>%
        arrange(desc(identity), desc(bit_score)) %>%
        do(head(., 1)) %>%
        ungroup() %>%
        select(string_id, string_name, uniprot_ac) %>%
        arrange(string_id, uniprot_ac)
}

load_go_mapping = function(species_name) {
    read_tsv(sprintf('shared/dump/%s-go-mapping.tab', species_name), col_names=TRUE)
}

map_string_id = function(data, mapping, col=string_id, uniprot_top_match=TRUE) {
    col = enquo(col)
    col_name = quo_name(col)
    col_str = as.character(col_name)

    if (is.vector(data)) {
        data = data_frame(!!col_name := data)
    }

    if (is.character(data[[col_str]])) {
        mapping = mapping %>%
           mutate(string_id = as.character(string_id))
    } else if (is.integer(data[[col_str]])) {
        mapping = mapping %>%
           mutate(string_id = as.integer(string_id))
    }

    join_by = 'string_id'
    names(join_by) = col_str

    if (col_str != 'string_id') {
        suffix = paste0('.', col_str)
        new_names = list()

        for (name in names(mapping)) {
            if (name == 'string_id') {
                new_names[[name]] = name
            } else {
                new_names[[paste0(name, suffix)]] = name
            }
        }

        mapping = rename(mapping, !!!new_names)
    }

    left_join(data, mapping, by=join_by)
}
