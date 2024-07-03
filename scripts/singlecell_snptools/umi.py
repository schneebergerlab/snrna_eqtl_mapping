from copy import copy
from collections import defaultdict
import itertools as it


def edit_dist(umi1, umi2):
    return sum(i != j for i, j in zip(umi1, umi2))


def umi_dedup(umi_counts):
    edges = defaultdict(set)
    nodes = sorted(umi_counts, key=umi_counts.__getitem__, reverse=True)
    for umi_i, umi_j in it.combinations(nodes, r=2):
        if edit_dist(umi_i, umi_j) <= 1:
            if umi_counts[umi_i] >= (2 * umi_counts[umi_j] + 1):
                edges[umi_i].add(umi_j)
    deduped_umi_counts = copy(umi_counts)
    for parent in nodes:
        for child in edges[parent]:
            try:
                deduped_umi_counts[parent] += deduped_umi_counts.pop(child)
            except KeyError:
                # umi already merged to a different parent
                continue
    return deduped_umi_counts


def umi_dedup_hap(umi_counts):
    edges = defaultdict(set)
    nodes = sorted(umi_counts, key=lambda k: sum(umi_counts[k].values()), reverse=True)
    for umi_i, umi_j in it.combinations(nodes, r=2):
        if edit_dist(umi_i, umi_j) <= 1:
            umi_i_count = sum(umi_counts[umi_i].values())
            umi_j_count = sum(umi_counts[umi_j].values())
            if umi_i_count >= (2 * umi_j_count + 1):
                edges[umi_i].add(umi_j)
    deduped_umi_counts = copy(umi_counts)
    for parent in nodes:
        for child in edges[parent]:
            try:
                deduped_umi_counts[parent].update(deduped_umi_counts.pop(child))
            except KeyError:
                # umi already merged to a different parent
                continue
    return deduped_umi_counts