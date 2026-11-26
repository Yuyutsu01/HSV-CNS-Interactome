"""
Phase 4 — CNS-Weighted Interactome Completion
Gaps filled:
  1. Rebuild interactome with expanded Phase 2 data (627 high-confidence hits)
  2. Network perturbation analysis (heat diffusion, simplified HotNet2)
  3. Community/module detection (label propagation + greedy modularity)
  4. Hub & bottleneck analysis
  5. Network robustness (targeted vs random node removal)
  6. Re-run hypothesis test H1 with expanded data
  7. Updated publication plots
"""

import os, sys, json, csv, math
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict, Counter
from scipy import stats
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome"
PHASE1_DIR = os.path.join(BASE_DIR, "phase1")
PHASE2_DIR = os.path.join(BASE_DIR, "phase2")
PHASE3_DIR = os.path.join(BASE_DIR, "phase3")
PHASE4_DIR = os.path.join(BASE_DIR, "phase4")
RESULTS_DIR = os.path.join(PHASE4_DIR, "results")
PLOT_DIR = os.path.join(PHASE4_DIR, "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

sys.path.insert(0, PHASE3_DIR)
from interactome_builder import (
    DOMAIN_MOTIF_MAP, GO_ANNOTATIONS, CNS_GO_TERMS,
    REACTOME_PATHWAYS, CNS_EXPRESSION,
    HYDROPHOBICITY, CHARGE, HELIX_PROPENSITY,
)


def build_expanded_interactome(hits_df, disorder_data):
    """Build bipartite interactome from expanded Phase 2 high-confidence hits."""
    G = nx.Graph()

    for _, row in hits_df.iterrows():
        viral_protein = row["protein"]
        category = row["category"]
        targets_str = str(row.get("targets", ""))
        targets = [t.strip() for t in targets_str.split(",") if t.strip()]
        motif_class = row["motif_class"]
        ml_score = float(row.get("ml_score", 0))
        pathway = str(row.get("pathway", ""))

        if not G.has_node(viral_protein):
            pdata = disorder_data.get(viral_protein, {})
            G.add_node(viral_protein, bipartite=0, node_type="viral",
                       category=category,
                       idr_fraction=pdata.get("idr_fraction", 0),
                       length=pdata.get("length", 0))

        for host_gene in targets:
            if not host_gene:
                continue
            brain_tpm = CNS_EXPRESSION.get(host_gene, 0)
            cns_specificity = np.log1p(brain_tpm) / np.log1p(350)

            if not G.has_node(host_gene):
                go_terms = GO_ANNOTATIONS.get(host_gene, [])
                cns_go_score = 0
                for term_str in go_terms:
                    go_id = term_str.split("|")[0]
                    if go_id in CNS_GO_TERMS:
                        cns_go_score = max(cns_go_score, CNS_GO_TERMS[go_id][1])

                G.add_node(host_gene, bipartite=1, node_type="host",
                           brain_tpm=brain_tpm,
                           cns_specificity=round(cns_specificity, 3),
                           cns_go_score=cns_go_score,
                           go_terms="; ".join(go_terms))

            edge_key = (viral_protein, host_gene)
            if G.has_edge(*edge_key):
                edata = G.edges[edge_key]
                edata["n_slims"] += 1
                edata["max_ml_score"] = max(edata["max_ml_score"], ml_score)
                edata["motif_classes"].append(motif_class)
                edata["pathways"].add(pathway)
            else:
                cns_edge_weight = ml_score * cns_specificity
                domain_info = DOMAIN_MOTIF_MAP.get(motif_class, {})
                G.add_edge(viral_protein, host_gene,
                           n_slims=1,
                           max_ml_score=ml_score,
                           cns_edge_weight=round(cns_edge_weight, 4),
                           motif_classes=[motif_class],
                           pathways={pathway},
                           domain=domain_info.get("domain", "unknown"),
                           binding_mode=domain_info.get("binding_mode", "unknown"))
    return G


def compute_network_metrics(G):
    """Full topology metrics including hub and bottleneck scores."""
    degree_cent = nx.degree_centrality(G)
    between_cent = nx.betweenness_centrality(G)
    closeness_cent = nx.closeness_centrality(G)

    try:
        eigenvector_cent = nx.eigenvector_centrality(G, max_iter=500)
    except nx.PowerIterationFailedConvergence:
        eigenvector_cent = {n: 0 for n in G.nodes()}

    metrics = {}
    for node, ndata in G.nodes(data=True):
        m = {
            "node_type": ndata.get("node_type", "unknown"),
            "category": ndata.get("category", ""),
            "degree": G.degree(node),
            "degree_centrality": round(degree_cent[node], 4),
            "betweenness_centrality": round(between_cent[node], 4),
            "closeness_centrality": round(closeness_cent[node], 4),
            "eigenvector_centrality": round(eigenvector_cent.get(node, 0), 4),
        }

        if ndata.get("node_type") == "viral":
            neighbors = list(G.neighbors(node))
            host_neighbors = [n for n in neighbors if G.nodes[n].get("node_type") == "host"]
            cns_weighted_degree = sum(
                G.edges[node, h].get("cns_edge_weight", 0) for h in host_neighbors
            )
            mean_host_cns = np.mean([
                G.nodes[h].get("cns_specificity", 0) for h in host_neighbors
            ]) if host_neighbors else 0

            pathways_targeted = set()
            for h in host_neighbors:
                pathways_targeted.update(G.edges[node, h].get("pathways", set()))

            m.update({
                "cns_weighted_degree": round(cns_weighted_degree, 3),
                "mean_host_cns_specificity": round(mean_host_cns, 3),
                "n_pathways": len(pathways_targeted),
                "pathways": sorted(pathways_targeted),
                "n_host_targets": len(host_neighbors),
            })
        else:
            m.update({
                "brain_tpm": ndata.get("brain_tpm", 0),
                "cns_specificity": ndata.get("cns_specificity", 0),
                "cns_go_score": ndata.get("cns_go_score", 0),
            })

        metrics[node] = m

    return metrics


def hub_bottleneck_analysis(G, metrics):
    """Identify hubs (high degree) and bottlenecks (high betweenness)."""
    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]
    host_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "host"]

    deg_threshold_viral = np.percentile([metrics[n]["degree"] for n in viral_nodes], 75) if viral_nodes else 0
    betw_threshold_viral = np.percentile([metrics[n]["betweenness_centrality"] for n in viral_nodes], 75) if viral_nodes else 0

    deg_threshold_host = np.percentile([metrics[n]["degree"] for n in host_nodes], 75) if host_nodes else 0
    betw_threshold_host = np.percentile([metrics[n]["betweenness_centrality"] for n in host_nodes], 75) if host_nodes else 0

    results = {"viral_hubs": [], "viral_bottlenecks": [], "viral_hub_bottlenecks": [],
               "host_hubs": [], "host_bottlenecks": [], "host_hub_bottlenecks": []}

    for n in viral_nodes:
        is_hub = metrics[n]["degree"] >= deg_threshold_viral
        is_bottleneck = metrics[n]["betweenness_centrality"] >= betw_threshold_viral
        entry = {"protein": n, "degree": metrics[n]["degree"],
                 "betweenness": metrics[n]["betweenness_centrality"],
                 "category": metrics[n].get("category", "")}
        if is_hub and is_bottleneck:
            results["viral_hub_bottlenecks"].append(entry)
        elif is_hub:
            results["viral_hubs"].append(entry)
        elif is_bottleneck:
            results["viral_bottlenecks"].append(entry)

    for n in host_nodes:
        is_hub = metrics[n]["degree"] >= deg_threshold_host
        is_bottleneck = metrics[n]["betweenness_centrality"] >= betw_threshold_host
        entry = {"protein": n, "degree": metrics[n]["degree"],
                 "betweenness": metrics[n]["betweenness_centrality"],
                 "brain_tpm": metrics[n].get("brain_tpm", 0)}
        if is_hub and is_bottleneck:
            results["host_hub_bottlenecks"].append(entry)
        elif is_hub:
            results["host_hubs"].append(entry)
        elif is_bottleneck:
            results["host_bottlenecks"].append(entry)

    for key in results:
        results[key].sort(key=lambda x: -x["degree"])

    return results


def heat_diffusion_perturbation(G, metrics, n_permutations=1000):
    """
    Simplified HotNet2-style heat diffusion.
    Seed heat at each viral protein; diffuse through weighted edges.
    Compare observed heat distribution to random permutations.
    """
    alpha = 0.5
    n_steps = 10

    def diffuse_heat(G, seed_node, alpha, steps):
        heat = {n: 0.0 for n in G.nodes()}
        heat[seed_node] = 1.0
        for _ in range(steps):
            new_heat = {}
            for n in G.nodes():
                retained = alpha * heat[n]
                neighbors = list(G.neighbors(n))
                if not neighbors:
                    new_heat[n] = retained
                    continue
                total_weight = sum(G.edges[n, nb].get("cns_edge_weight", 0.01) for nb in neighbors)
                incoming = 0
                for nb in neighbors:
                    w = G.edges[n, nb].get("cns_edge_weight", 0.01)
                    nb_neighbors = list(G.neighbors(nb))
                    nb_total_w = sum(G.edges[nb, nn].get("cns_edge_weight", 0.01) for nn in nb_neighbors)
                    if nb_total_w > 0:
                        incoming += (1 - alpha) * heat[nb] * (w / nb_total_w)
                new_heat[n] = retained + incoming
            total = sum(new_heat.values())
            if total > 0:
                heat = {n: v / total for n, v in new_heat.items()}
            else:
                heat = new_heat
        return heat

    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]
    host_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "host"]

    perturbation_results = {}
    for viral in viral_nodes:
        heat = diffuse_heat(G, viral, alpha, n_steps)
        host_heat = {n: heat[n] for n in host_nodes}
        cns_host = [n for n in host_nodes if G.nodes[n].get("cns_go_score", 0) >= 2.0]
        non_cns_host = [n for n in host_nodes if G.nodes[n].get("cns_go_score", 0) < 2.0]

        cns_heat_sum = sum(host_heat.get(n, 0) for n in cns_host)
        non_cns_heat_sum = sum(host_heat.get(n, 0) for n in non_cns_host)

        random_cns_sums = []
        rng = np.random.RandomState(42)
        for _ in range(n_permutations):
            perm_cns = rng.choice(host_nodes, size=len(cns_host), replace=False)
            random_cns_sums.append(sum(host_heat.get(n, 0) for n in perm_cns))

        z_score = (cns_heat_sum - np.mean(random_cns_sums)) / max(np.std(random_cns_sums), 1e-10)
        p_value = np.mean([1 if rs >= cns_heat_sum else 0 for rs in random_cns_sums])

        top_heated = sorted(host_heat.items(), key=lambda x: -x[1])[:10]

        perturbation_results[viral] = {
            "category": G.nodes[viral].get("category", ""),
            "cns_heat_sum": round(cns_heat_sum, 6),
            "non_cns_heat_sum": round(non_cns_heat_sum, 6),
            "cns_heat_fraction": round(cns_heat_sum / max(cns_heat_sum + non_cns_heat_sum, 1e-10), 4),
            "z_score": round(z_score, 3),
            "permutation_p": round(p_value, 4),
            "top_heated_hosts": [{"host": h, "heat": round(v, 6)} for h, v in top_heated],
        }

    return perturbation_results


def community_detection(G):
    """Detect communities using greedy modularity and label propagation."""
    communities_greedy = list(nx.community.greedy_modularity_communities(G))
    communities_lpa = list(nx.community.label_propagation_communities(G))

    modularity_greedy = nx.community.modularity(G, communities_greedy)

    community_map_greedy = {}
    for idx, comm in enumerate(communities_greedy):
        for node in comm:
            community_map_greedy[node] = idx

    community_summary = []
    for idx, comm in enumerate(communities_greedy):
        comm_nodes = list(comm)
        viral_in_comm = [n for n in comm_nodes if G.nodes[n].get("node_type") == "viral"]
        host_in_comm = [n for n in comm_nodes if G.nodes[n].get("node_type") == "host"]

        pathways = set()
        for v in viral_in_comm:
            for h in host_in_comm:
                if G.has_edge(v, h):
                    pathways.update(G.edges[v, h].get("pathways", set()))

        mean_cns = np.mean([G.nodes[h].get("cns_specificity", 0) for h in host_in_comm]) if host_in_comm else 0

        community_summary.append({
            "community_id": idx,
            "size": len(comm_nodes),
            "viral_proteins": sorted(viral_in_comm),
            "host_proteins": sorted(host_in_comm),
            "n_viral": len(viral_in_comm),
            "n_host": len(host_in_comm),
            "pathways": sorted(pathways),
            "mean_cns_specificity": round(mean_cns, 3),
        })

    community_summary.sort(key=lambda x: -x["size"])

    return {
        "method": "greedy_modularity + label_propagation",
        "n_communities_greedy": len(communities_greedy),
        "n_communities_lpa": len(communities_lpa),
        "modularity_greedy": round(modularity_greedy, 4),
        "community_map": community_map_greedy,
        "communities": community_summary,
    }


def network_robustness(G):
    """Test network robustness: targeted vs random removal."""
    def largest_component_fraction(graph):
        if graph.number_of_nodes() == 0:
            return 0
        return max(len(c) for c in nx.connected_components(graph)) / graph.number_of_nodes()

    n_original = G.number_of_nodes()
    original_lcc = largest_component_fraction(G)

    degree_sorted = sorted(G.nodes(), key=lambda n: G.degree(n), reverse=True)
    betw = nx.betweenness_centrality(G)
    betw_sorted = sorted(G.nodes(), key=lambda n: betw[n], reverse=True)

    fractions = np.arange(0, 0.55, 0.05)
    targeted_degree = []
    targeted_betw = []
    random_avg = []

    for frac in fractions:
        n_remove = int(frac * n_original)

        G_td = G.copy()
        to_remove = degree_sorted[:n_remove]
        G_td.remove_nodes_from(to_remove)
        targeted_degree.append(largest_component_fraction(G_td))

        G_tb = G.copy()
        to_remove_b = betw_sorted[:n_remove]
        G_tb.remove_nodes_from(to_remove_b)
        targeted_betw.append(largest_component_fraction(G_tb))

        rand_vals = []
        rng = np.random.RandomState(42)
        for _ in range(20):
            G_r = G.copy()
            rand_remove = rng.choice(list(G_r.nodes()), size=min(n_remove, G_r.number_of_nodes()), replace=False)
            G_r.remove_nodes_from(rand_remove)
            rand_vals.append(largest_component_fraction(G_r))
        random_avg.append(np.mean(rand_vals))

    return {
        "fractions_removed": [round(f, 2) for f in fractions],
        "lcc_targeted_degree": [round(v, 4) for v in targeted_degree],
        "lcc_targeted_betweenness": [round(v, 4) for v in targeted_betw],
        "lcc_random_avg": [round(v, 4) for v in random_avg],
        "original_lcc_fraction": round(original_lcc, 4),
    }


def hypothesis_test_expanded(metrics, disorder_data):
    """H1 with expanded interactome data."""
    target_scores = []
    control_scores = []

    for protein, m in metrics.items():
        if m["node_type"] != "viral":
            continue
        pdata = disorder_data.get(protein, {})
        boundaries = pdata.get("idr_boundaries", [])
        idr_aa = sum(e - s + 1 for s, e in [tuple(b) for b in boundaries])
        if idr_aa == 0:
            idr_aa = max(1, int(pdata.get("idr_fraction", 0.1) * pdata.get("length", 100)))

        slim_density = m["degree"] / max(idr_aa, 1) * 100
        cns_frac = m.get("mean_host_cns_specificity", 0)
        cns_weighted_degree = m.get("cns_weighted_degree", 0)
        n_pathways = m.get("n_pathways", 0)

        composite_h1 = (
            0.4 * min(cns_weighted_degree / 15, 1.0) +
            0.3 * min((slim_density * cns_frac) / 25, 1.0) +
            0.3 * min(n_pathways / 15, 1.0)
        )

        entry = {
            "protein": protein,
            "category": m.get("category", ""),
            "degree": m["degree"],
            "idr_aa": idr_aa,
            "slim_density": round(slim_density, 2),
            "cns_interaction_frac": round(cns_frac, 3),
            "cns_targeting_score": round(slim_density * cns_frac, 3),
            "cns_weighted_degree": round(cns_weighted_degree, 3),
            "pathway_diversity": n_pathways,
            "composite_h1_score": round(composite_h1, 4),
            "eigenvector_centrality": m.get("eigenvector_centrality", 0),
            "betweenness_centrality": m.get("betweenness_centrality", 0),
        }

        if m.get("category") == "target":
            target_scores.append(entry)
        else:
            control_scores.append(entry)

    t_vals = [e["composite_h1_score"] for e in target_scores]
    c_vals = [e["composite_h1_score"] for e in control_scores]

    if len(t_vals) >= 2 and len(c_vals) >= 2:
        u_stat, u_pval = stats.mannwhitneyu(t_vals, c_vals, alternative="greater")
    else:
        u_stat, u_pval = float("nan"), float("nan")

    pooled_std = np.std(t_vals + c_vals)
    effect_size = (np.mean(t_vals) - np.mean(c_vals)) / pooled_std if pooled_std > 0 else 0

    if len(t_vals) >= 2 and len(c_vals) >= 2:
        t_stat, t_pval = stats.ttest_ind(t_vals, c_vals, alternative="greater")
    else:
        t_stat, t_pval = float("nan"), float("nan")

    n_perm = 10000
    rng = np.random.RandomState(42)
    all_vals = t_vals + c_vals
    obs_diff = np.mean(t_vals) - np.mean(c_vals)
    perm_diffs = []
    for _ in range(n_perm):
        rng.shuffle(all_vals)
        perm_diffs.append(np.mean(all_vals[:len(t_vals)]) - np.mean(all_vals[len(t_vals):]))
    perm_p = np.mean([1 if d >= obs_diff else 0 for d in perm_diffs])

    # Post-hoc power analysis
    sys.path.insert(0, BASE_DIR)
    from utils.power_analysis import full_power_report
    power = full_power_report(len(t_vals), len(c_vals), effect_size)
    print(f"  H1 power analysis: achieved power={power['achieved_power']:.1%}, "
          f"MDE={power['min_detectable_effect_d_at_80pct']}")
    print(f"  {power['interpretation']}")

    return {
        "targets": sorted(target_scores, key=lambda x: -x["composite_h1_score"]),
        "controls": sorted(control_scores, key=lambda x: -x["composite_h1_score"]),
        "mann_whitney_U": round(float(u_stat), 2) if not np.isnan(u_stat) else None,
        "mann_whitney_p": round(float(u_pval), 4) if not np.isnan(u_pval) else None,
        "welch_t_stat": round(float(t_stat), 3) if not np.isnan(t_stat) else None,
        "welch_t_p": round(float(t_pval), 4) if not np.isnan(t_pval) else None,
        "permutation_p": round(perm_p, 4),
        "effect_size_cohens_d": round(effect_size, 3),
        "target_mean": round(np.mean(t_vals), 4),
        "control_mean": round(np.mean(c_vals), 4),
        "target_std": round(np.std(t_vals), 4),
        "control_std": round(np.std(c_vals), 4),
        "observed_diff": round(obs_diff, 4),
        "power_analysis": power,
    }


def pathway_enrichment(G):
    """GO + Reactome enrichment per viral protein."""
    results = {}
    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]
    all_host = set(n for n, d in G.nodes(data=True) if d.get("node_type") == "host")

    for viral in viral_nodes:
        host_targets = set(n for n in G.neighbors(viral) if G.nodes[n].get("node_type") == "host")

        go_counts = defaultdict(int)
        for gene in host_targets:
            for term_str in GO_ANNOTATIONS.get(gene, []):
                go_id, go_name = term_str.split("|")
                go_counts[(go_id, go_name)] += 1

        bg_go = defaultdict(int)
        for gene in all_host:
            for term_str in GO_ANNOTATIONS.get(gene, []):
                go_id, go_name = term_str.split("|")
                bg_go[(go_id, go_name)] += 1

        go_enrichment = []
        for (go_id, go_name), count in go_counts.items():
            bg_count = bg_go.get((go_id, go_name), 0)
            if bg_count < 2:
                continue
            table = [[count, bg_count - count],
                     [len(host_targets) - count, len(all_host) - len(host_targets) - (bg_count - count)]]
            table = [[max(0, c) for c in row] for row in table]
            try:
                odds, pval = stats.fisher_exact(table, alternative="greater")
            except ValueError:
                continue
            go_enrichment.append({
                "go_id": go_id, "go_name": go_name, "count": count,
                "total_targets": len(host_targets), "bg_count": bg_count,
                "odds_ratio": round(odds, 2), "p_value": round(pval, 4),
            })
        go_enrichment.sort(key=lambda x: x["p_value"])

        reactome = []
        for pw_name, pw_data in REACTOME_PATHWAYS.items():
            pw_genes = set(pw_data["genes"])
            overlap = host_targets & pw_genes
            if not overlap:
                continue
            bg_overlap = all_host & pw_genes
            table = [[len(overlap), len(bg_overlap) - len(overlap)],
                     [len(host_targets) - len(overlap), len(all_host) - len(host_targets) - (len(bg_overlap) - len(overlap))]]
            table = [[max(0, c) for c in row] for row in table]
            try:
                odds, pval = stats.fisher_exact(table, alternative="greater")
            except ValueError:
                continue
            reactome.append({
                "pathway": pw_name, "reactome_id": pw_data["id"],
                "overlap_genes": sorted(overlap), "count": len(overlap),
                "pathway_size": len(pw_genes),
                "odds_ratio": round(odds, 2), "p_value": round(pval, 4),
            })
        reactome.sort(key=lambda x: x["p_value"])

        results[viral] = {
            "n_host_targets": len(host_targets),
            "go_enrichment": go_enrichment[:15],
            "reactome_enrichment": reactome[:10],
        }

    return results


# ═══════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════

def plot_expanded_network(G, output_path):
    """Network visualization with community coloring."""
    fig, ax = plt.subplots(1, 1, figsize=(20, 16))

    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]
    host_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "host"]

    pos = nx.spring_layout(G, k=2.0, iterations=150, seed=42)

    target_viral = [n for n in viral_nodes if G.nodes[n].get("category") == "target"]
    control_viral = [n for n in viral_nodes if G.nodes[n].get("category") == "control"]

    host_sizes = [max(80, G.nodes[n].get("brain_tpm", 10) * 2.5) for n in host_nodes]
    host_colors = [G.nodes[n].get("cns_specificity", 0) for n in host_nodes]

    nx.draw_networkx_nodes(G, pos, nodelist=host_nodes, node_size=host_sizes,
                           node_color=host_colors, cmap=plt.cm.YlOrRd,
                           vmin=0, vmax=1, alpha=0.6, ax=ax)

    nx.draw_networkx_nodes(G, pos, nodelist=target_viral, node_size=900,
                           node_color="#e74c3c", node_shape="s", alpha=0.9, ax=ax,
                           edgecolors="black", linewidths=2)
    nx.draw_networkx_nodes(G, pos, nodelist=control_viral, node_size=500,
                           node_color="#3498db", node_shape="s", alpha=0.7, ax=ax,
                           edgecolors="black", linewidths=1)

    edge_weights = [G.edges[e].get("cns_edge_weight", 0.01) for e in G.edges()]
    max_w = max(edge_weights) if edge_weights else 1
    target_edges = [(u, v) for u, v in G.edges()
                    if G.nodes[u].get("category") == "target" or G.nodes[v].get("category") == "target"]
    control_edges = [(u, v) for u, v in G.edges() if (u, v) not in target_edges]

    nx.draw_networkx_edges(G, pos, edgelist=target_edges,
                           width=[max(0.3, G.edges[e].get("cns_edge_weight", 0.01) / max_w * 3) for e in target_edges],
                           alpha=0.5, edge_color="#e74c3c", ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=control_edges,
                           width=[max(0.2, G.edges[e].get("cns_edge_weight", 0.01) / max_w * 2) for e in control_edges],
                           alpha=0.2, edge_color="#95a5a6", ax=ax)

    labels = {}
    for n in viral_nodes:
        labels[n] = n
    for n in host_nodes:
        if G.degree(n) >= 4 or G.nodes[n].get("brain_tpm", 0) > 80:
            labels[n] = n
    nx.draw_networkx_labels(G, pos, labels, font_size=7, font_weight="bold", ax=ax)

    ax.set_title("HSV-1 SLiM-Mediated CNS Interactome (Expanded, 627 high-confidence SLiMs)",
                 fontsize=14, fontweight="bold")
    legend_elements = [
        mpatches.Patch(color="#e74c3c", label="Neurovirulence targets (ICP34.5, US3, pUL37)"),
        mpatches.Patch(color="#3498db", label="Control proteins"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#ff6b35",
                   markersize=10, label="High CNS expression"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_robustness(robustness, output_path):
    """Plot network robustness curves."""
    fig, ax = plt.subplots(figsize=(8, 6))
    fracs = robustness["fractions_removed"]
    ax.plot(fracs, robustness["lcc_targeted_degree"], "r-o", label="Targeted (degree)", linewidth=2, markersize=5)
    ax.plot(fracs, robustness["lcc_targeted_betweenness"], "b-s", label="Targeted (betweenness)", linewidth=2, markersize=5)
    ax.plot(fracs, robustness["lcc_random_avg"], "g--^", label="Random (avg 20 trials)", linewidth=2, markersize=5)
    ax.set_xlabel("Fraction of nodes removed", fontsize=12)
    ax.set_ylabel("Largest connected component fraction", fontsize=12)
    ax.set_title("Network Robustness: Targeted vs Random Attack", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_heat_diffusion(perturbation, output_path):
    """Plot heat diffusion results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    proteins = sorted(perturbation.keys(), key=lambda x: -perturbation[x]["cns_heat_fraction"])
    categories = [perturbation[p]["category"] for p in proteins]
    cns_fracs = [perturbation[p]["cns_heat_fraction"] for p in proteins]
    z_scores = [perturbation[p]["z_score"] for p in proteins]
    colors = ["#e74c3c" if c == "target" else "#3498db" for c in categories]

    ax = axes[0]
    bars = ax.barh(range(len(proteins)), cns_fracs, color=colors, alpha=0.7)
    ax.set_yticks(range(len(proteins)))
    ax.set_yticklabels(proteins, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("CNS heat fraction", fontsize=10)
    ax.set_title("A. Heat Diffusion: CNS Pathway Concentration", fontsize=11, fontweight="bold")

    ax = axes[1]
    bars = ax.barh(range(len(proteins)), z_scores, color=colors, alpha=0.7)
    ax.set_yticks(range(len(proteins)))
    ax.set_yticklabels(proteins, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(1.96, color="red", linestyle="--", alpha=0.5, label="z=1.96 (p<0.05)")
    ax.set_xlabel("z-score (vs permuted CNS labels)", fontsize=10)
    ax.set_title("B. CNS Heat Enrichment Significance", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_hypothesis_expanded(h1, output_path):
    """Updated H1 hypothesis test figure."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    t_vals = [e["composite_h1_score"] for e in h1["targets"]]
    c_vals = [e["composite_h1_score"] for e in h1["controls"]]

    ax = axes[0]
    bp = ax.boxplot([t_vals, c_vals], patch_artist=True, widths=0.5)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Neurovirulence\n(n={len(t_vals)})", f"Control\n(n={len(c_vals)})"])
    bp["boxes"][0].set_facecolor("#e74c3c"); bp["boxes"][0].set_alpha(0.5)
    bp["boxes"][1].set_facecolor("#3498db"); bp["boxes"][1].set_alpha(0.5)
    for e in h1["targets"]:
        ax.plot(1, e["composite_h1_score"], "o", color="#e74c3c", markersize=10, zorder=5)
        ax.annotate(e["protein"], (1.15, e["composite_h1_score"]), fontsize=8)
    for e in h1["controls"]:
        ax.plot(2, e["composite_h1_score"], "o", color="#3498db", markersize=8, zorder=5)
        ax.annotate(e["protein"], (2.15, e["composite_h1_score"]), fontsize=7)

    p_str = f"MW p={h1['mann_whitney_p']}, perm p={h1['permutation_p']}" if h1["mann_whitney_p"] else "N/A"
    ax.set_title(f"A. CNS-Targeting Score\n({p_str})", fontsize=10, fontweight="bold")
    ax.set_ylabel("Composite H1 Score")

    ax = axes[1]
    all_prots = h1["targets"] + h1["controls"]
    all_prots.sort(key=lambda x: -x["composite_h1_score"])
    names = [p["protein"] for p in all_prots]
    colors = ["#e74c3c" if p["category"] == "target" else "#3498db" for p in all_prots]
    scores = [p["composite_h1_score"] for p in all_prots]
    y = np.arange(len(names))
    ax.barh(y, scores, color=colors, alpha=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Composite H1 Score")
    ax.set_title("B. Protein Ranking", fontsize=11, fontweight="bold")

    ax = axes[2]
    for e in h1["targets"]:
        ax.scatter(e["cns_weighted_degree"], e["pathway_diversity"],
                   c="#e74c3c", s=150, zorder=5, edgecolors="black")
        ax.annotate(e["protein"], (e["cns_weighted_degree"] + 0.1, e["pathway_diversity"] + 0.2), fontsize=8)
    for e in h1["controls"]:
        ax.scatter(e["cns_weighted_degree"], e["pathway_diversity"],
                   c="#3498db", s=100, zorder=5, edgecolors="gray", alpha=0.7)
        ax.annotate(e["protein"], (e["cns_weighted_degree"] + 0.1, e["pathway_diversity"] + 0.2), fontsize=7)
    ax.set_xlabel("CNS-weighted degree")
    ax.set_ylabel("Pathway diversity")
    ax.set_title("C. CNS Interaction vs Pathway Breadth", fontsize=11, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    print("=" * 80)
    print("PHASE 4 — CNS-Weighted Interactome (Expanded Data)")
    print("=" * 80)

    with open(os.path.join(PHASE1_DIR, "results", "disorder_summary_v3.json")) as f:
        disorder_data = {d["protein"]: d for d in json.load(f)}

    hits_path = os.path.join(PHASE2_DIR, "results", "slim_hits_high_confidence.csv")
    hits_df = pd.read_csv(hits_path, encoding="utf-8")
    print(f"\nLoaded {len(hits_df)} high-confidence SLiM hits from expanded Phase 2")

    # 1. Build expanded interactome
    print("\n[1] Building expanded CNS-weighted interactome...")
    G = build_expanded_interactome(hits_df, disorder_data)
    viral_count = sum(1 for _, d in G.nodes(data=True) if d.get("node_type") == "viral")
    host_count = sum(1 for _, d in G.nodes(data=True) if d.get("node_type") == "host")
    print(f"  Network: {G.number_of_nodes()} nodes ({viral_count} viral, {host_count} host), "
          f"{G.number_of_edges()} edges")

    # 2. Network metrics
    print("\n[2] Computing network metrics...")
    metrics = compute_network_metrics(G)
    print(f"\n  {'Protein':<15} {'Cat':<8} {'Deg':>4} {'CNS-WD':>8} {'MeanCNS':>8} {'#PW':>4} {'Betw':>8} {'Eigen':>8}")
    print("  " + "-" * 72)
    for p in sorted(metrics.keys()):
        m = metrics[p]
        if m["node_type"] != "viral":
            continue
        print(f"  {p:<15} {m.get('category',''):<8} {m['degree']:>4} "
              f"{m.get('cns_weighted_degree',0):>8.3f} {m.get('mean_host_cns_specificity',0):>8.3f} "
              f"{m.get('n_pathways',0):>4} {m['betweenness_centrality']:>8.4f} "
              f"{m['eigenvector_centrality']:>8.4f}")

    # 3. Hub & bottleneck analysis
    print("\n[3] Hub & bottleneck analysis...")
    hub_results = hub_bottleneck_analysis(G, metrics)
    for key, entries in hub_results.items():
        if entries:
            print(f"  {key}: {', '.join(e['protein'] for e in entries)}")

    # 4. Heat diffusion perturbation
    print("\n[4] Heat diffusion perturbation analysis (1000 permutations)...")
    perturbation = heat_diffusion_perturbation(G, metrics, n_permutations=1000)
    print(f"\n  {'Protein':<15} {'Cat':<8} {'CNS_frac':>9} {'z-score':>8} {'perm_p':>8}")
    print("  " + "-" * 55)
    for p in sorted(perturbation.keys(), key=lambda x: -perturbation[x]["cns_heat_fraction"]):
        d = perturbation[p]
        print(f"  {p:<15} {d['category']:<8} {d['cns_heat_fraction']:>9.4f} "
              f"{d['z_score']:>8.3f} {d['permutation_p']:>8.4f}")

    # 5. Community detection
    print("\n[5] Community detection...")
    communities = community_detection(G)
    print(f"  Greedy modularity: {communities['n_communities_greedy']} communities, Q={communities['modularity_greedy']:.4f}")
    print(f"  Label propagation: {communities['n_communities_lpa']} communities")
    for comm in communities["communities"][:5]:
        viral_str = ", ".join(comm["viral_proteins"]) if comm["viral_proteins"] else "none"
        print(f"  Community {comm['community_id']}: {comm['size']} nodes "
              f"(viral: {viral_str}, CNS={comm['mean_cns_specificity']:.3f})")

    # 6. Network robustness
    print("\n[6] Network robustness analysis...")
    robustness = network_robustness(G)
    print(f"  At 20% removal - targeted(deg): {robustness['lcc_targeted_degree'][4]:.4f}, "
          f"targeted(betw): {robustness['lcc_targeted_betweenness'][4]:.4f}, "
          f"random: {robustness['lcc_random_avg'][4]:.4f}")

    # 7. Pathway enrichment
    print("\n[7] Pathway enrichment...")
    enrichment = pathway_enrichment(G)
    for protein in ["ICP34.5", "US3", "pUL37"]:
        if protein in enrichment:
            data = enrichment[protein]
            print(f"\n  {protein} ({data['n_host_targets']} host targets):")
            for r in data["reactome_enrichment"][:5]:
                sig = "*" if r["p_value"] < 0.05 else ""
                print(f"    {r['pathway']:<25} p={r['p_value']:.4f}{sig} "
                      f"({r['count']}/{r['pathway_size']}: {', '.join(r['overlap_genes'])})")

    # 8. Hypothesis test
    print("\n[8] Hypothesis test H1 (expanded data)...")
    h1 = hypothesis_test_expanded(metrics, disorder_data)
    print(f"\n  H1: Neurovirulence proteins have higher CNS-targeting scores")
    print(f"  Target mean:  {h1['target_mean']:.4f} +/- {h1['target_std']:.4f}")
    print(f"  Control mean: {h1['control_mean']:.4f} +/- {h1['control_std']:.4f}")
    print(f"  Effect size (Cohen's d): {h1['effect_size_cohens_d']:.3f}")
    print(f"  Mann-Whitney U: {h1['mann_whitney_U']}, p={h1['mann_whitney_p']}")
    print(f"  Welch t: {h1['welch_t_stat']}, p={h1['welch_t_p']}")
    print(f"  Permutation p: {h1['permutation_p']}")

    print(f"\n  {'Protein':<15} {'Cat':<8} {'Deg':>4} {'Density':>8} {'CNS-WD':>8} {'#PW':>4} {'H1':>8}")
    print("  " + "-" * 60)
    for e in h1["targets"] + sorted(h1["controls"], key=lambda x: -x["composite_h1_score"]):
        print(f"  {e['protein']:<15} {e['category']:<8} {e['degree']:>4} "
              f"{e['slim_density']:>8.2f} {e['cns_weighted_degree']:>8.3f} "
              f"{e['pathway_diversity']:>4} {e['composite_h1_score']:>8.4f}")

    # 9. Plots
    print("\n[9] Generating plots...")
    plot_expanded_network(G, os.path.join(PLOT_DIR, "cns_interactome_expanded.png"))
    print(f"  Interactome network")
    plot_robustness(robustness, os.path.join(PLOT_DIR, "network_robustness.png"))
    print(f"  Robustness")
    plot_heat_diffusion(perturbation, os.path.join(PLOT_DIR, "heat_diffusion.png"))
    print(f"  Heat diffusion")
    plot_hypothesis_expanded(h1, os.path.join(PLOT_DIR, "hypothesis_test_h1_expanded.png"))
    print(f"  Hypothesis test")

    # 10. Save all results
    print("\n[10] Saving results...")

    edges_data = []
    for u, v, d in G.edges(data=True):
        edges_data.append({
            "viral_protein": u if G.nodes[u].get("node_type") == "viral" else v,
            "host_protein": v if G.nodes[v].get("node_type") == "host" else u,
            "n_slims": d.get("n_slims", 1),
            "max_ml_score": d.get("max_ml_score", 0),
            "cns_edge_weight": d.get("cns_edge_weight", 0),
            "domain": d.get("domain", ""),
            "binding_mode": d.get("binding_mode", ""),
            "motif_classes": "; ".join(d.get("motif_classes", [])),
            "pathways": "; ".join(d.get("pathways", set())),
        })
    pd.DataFrame(edges_data).to_csv(os.path.join(RESULTS_DIR, "interactome_edges_expanded.csv"), index=False)

    metrics_list = []
    for node, m in metrics.items():
        entry = {"node": node}
        entry.update({k: v for k, v in m.items() if k not in ("pathways", "viral_partners")})
        if "pathways" in m:
            entry["pathways"] = "; ".join(m["pathways"])
        metrics_list.append(entry)
    pd.DataFrame(metrics_list).to_csv(os.path.join(RESULTS_DIR, "network_metrics_expanded.csv"), index=False)

    with open(os.path.join(RESULTS_DIR, "hypothesis_test_h1_expanded.json"), "w") as f:
        json.dump(h1, f, indent=2, default=str)

    with open(os.path.join(RESULTS_DIR, "heat_diffusion_perturbation.json"), "w") as f:
        json.dump(perturbation, f, indent=2, default=str)

    comm_save = {k: v for k, v in communities.items() if k != "community_map"}
    with open(os.path.join(RESULTS_DIR, "community_detection.json"), "w") as f:
        json.dump(comm_save, f, indent=2, default=str)

    with open(os.path.join(RESULTS_DIR, "hub_bottleneck_analysis.json"), "w") as f:
        json.dump(hub_results, f, indent=2, default=str)

    with open(os.path.join(RESULTS_DIR, "network_robustness.json"), "w") as f:
        json.dump(robustness, f, indent=2, default=str)

    with open(os.path.join(RESULTS_DIR, "pathway_enrichment_expanded.json"), "w") as f:
        json.dump(enrichment, f, indent=2, default=str)

    G_export = G.copy()
    for u, v, d in G_export.edges(data=True):
        d["pathways"] = "; ".join(d.get("pathways", set()))
        d["motif_classes"] = "; ".join(d.get("motif_classes", []))
    nx.write_graphml(G_export, os.path.join(RESULTS_DIR, "interactome_expanded.graphml"))

    summary = {
        "data_source": "slim_hits_high_confidence.csv (627 hits, expanded ELM + ML scoring)",
        "network": {"nodes": G.number_of_nodes(), "edges": G.number_of_edges(),
                     "viral": viral_count, "host": host_count},
        "hypothesis_test": {
            "mann_whitney_p": h1["mann_whitney_p"],
            "permutation_p": h1["permutation_p"],
            "effect_size": h1["effect_size_cohens_d"],
            "target_mean": h1["target_mean"],
            "control_mean": h1["control_mean"],
        },
        "communities": {"n_greedy": communities["n_communities_greedy"],
                         "modularity": communities["modularity_greedy"]},
        "perturbation_top_cns_viral": sorted(
            [(p, d["cns_heat_fraction"], d["z_score"]) for p, d in perturbation.items()],
            key=lambda x: -x[1])[:5],
    }
    with open(os.path.join(RESULTS_DIR, "phase4_complete_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\nPhase 4 complete! Results in {RESULTS_DIR}")


if __name__ == "__main__":
    main()
