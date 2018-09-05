#!/usr/bin/zsh

taxo_name="$1"
taxo_id="$2"
evtype="$3"
scorethreshold="$4"

out_file_prefix="/opt/shared/dump"

docker-compose exec stringdb sh -c "mkdir -p '$out_file_prefix'"
# mkdir -p "$out_file_prefix"


if [ $# = 2 ]
then
    out_file="$out_file_prefix/$taxo_name.tab"

    docker-compose exec stringdb psql stringdb stringdb -c \
    "\\copy (
        select
          node_id_a,
          node_id_b,
          combined_score,
          array(select unnest(evidence_scores[:][1])) evidence_score_types
        from
          network.node_node_links
        where
          node_type_b = $taxo_id
    ) to '$out_file' with csv header delimiter E'\\t';
    "
elif [ $# = 3 ]
then
    out_file="$out_file_prefix/$taxo_name$evtype.tab"

    docker-compose exec stringdb psql stringdb stringdb -c \
    "\\copy (
        select
          node_id_a,
          node_id_b,
          evidence_scores[i][2] evidence_score
        from (
          select
            node_id_a,
            node_id_b,
            evidence_scores,
            generate_subscripts(evidence_scores, 1) as i
          from
            network.node_node_links
          where
            node_type_b = $taxo_id
        ) as indexed_scores
        where
          evidence_scores[i][1] = $evtype
    ) to '$out_file' with csv header delimiter E'\\t';
    "
elif [ $# = 4 ]
then
    out_file="$out_file_prefix/$taxo_name$evtype-threshold$scorethreshold.tab"

    docker-compose exec stringdb psql stringdb stringdb -c \
    "\\copy (
        select
          node_id_a,
          node_id_b,
          evidence_scores[i][2] evidence_score
        from (
          select
            node_id_a,
            node_id_b,
            evidence_scores,
            generate_subscripts(evidence_scores, 1) as i
          from
            network.node_node_links
          where
            node_type_b = $taxo_id
        ) as indexed_scores
        where
          evidence_scores[i][1] = $evtype
          and
          evidence_scores[i][2] >= $scorethreshold
    ) to '$out_file' with csv header delimiter E'\\t';
    "
fi

if [ -f "$out_file" ]
then
    docker-compose exec stringdb chmod a+rw "$out_file"
fi
