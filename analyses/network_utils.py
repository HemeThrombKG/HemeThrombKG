# -*- coding: utf-8 -*-

"""Utils and visualization functions for HemeThrombKG and HemeKG 2.0."""

from collections import defaultdict
from itertools import permutations

import igraph
import networkx as nx
import pybel

COLORS_CATEGORY = {
    'literature': "thistle",  # purple
    'pathway': "lightblue",  # blue
    'common': "#F60000",  # red
    'other': "lightgray",
}
FUNCTIONS = ['Protein', 'Abundance', 'Gene', 'RNA']
BLACKLIST_FUNCTIONS = ['Complex', 'Composite', 'Reaction']


def get_bel_graph(file):
    """Load bel graphs from bel scripts."""
    return pybel.from_bel_script(
        f'../bel_files/{file}.bel',
        no_identifier_validation=True, allow_definition_failures=True,
    )


def get_cellular_node_paths(kg, cell_terms):
    """Get cell and cell-related terms from KG."""
    cell_terms_set = set()

    for node in kg:

        if node.function not in BLACKLIST_FUNCTIONS:
            if node.name in cell_terms:
                cell_terms_set.add(node)

    return cell_terms_set


def get_node_paths(kg, nodes):
    """Get paths between node pairs in KG."""
    connected_nodes = []

    for node1, node2 in permutations(nodes, 2):

        if not nx.has_path(kg, node1, node2):
            continue

        shortest_path = nx.shortest_path(kg, node1, node2)

        for node in shortest_path:
            connected_nodes.append(node)

    return connected_nodes


def get_nodes_by_region(kg, proteins_list):
    """Get nodes in subcellular or extracellular regions from KG."""
    node_set = set()

    for node in kg:

        if (node.function in FUNCTIONS) and (node.name in proteins_list):
            node_set.add(node)

        # Get proteins from complex nodes
        elif node.function == "Complex":

            if node.entity is not None:
                node_set.add(node)

            else:
                for member in node.members:
                    if (member.function in FUNCTIONS) and (member.name in proteins_list):
                        node_set.add(node)

        # Get proteins from reaction nodes
        elif node.function == "Reaction":
            for reactant in node.reactants:
                if (reactant.function in FUNCTIONS) and (reactant.name in proteins_list):
                    node_set.add(node)

            for product in node.products:
                if (product.function in FUNCTIONS) and (product.name in proteins_list):
                    node_set.add(node)

    return node_set


def normalize(values, bounds):
    """Normalize node sizes."""
    return [bounds['desired']['lower'] + (x - bounds['actual']['lower']) * (
        bounds['desired']['upper'] - bounds['desired']['lower']) / (
                bounds['actual']['upper'] - bounds['actual']['lower']) for x in values]


def get_node_counts(graph, intracellular=True):
    """Get frequencies of nodes in edges within a network."""
    node_dict = defaultdict(int)

    for source, target, data in graph.edges(data=True):

        if intracellular:

            if source.function not in BLACKLIST_FUNCTIONS:
                node_dict[source.name] += 1

            if target.function not in BLACKLIST_FUNCTIONS:
                node_dict[target.name] += 1

        else:

            if source.function in FUNCTIONS:
                node_dict[source.name] += 1

            elif source.function == 'Complex':

                for source_member in source.members:
                    node_dict[source_member.name] += 1

            if target.function in FUNCTIONS:
                node_dict[target.name] += 1

            elif target.function == 'Complex':

                for target_member in target.members:
                    node_dict[target_member.name] += 1

    return dict(sorted(node_dict.items(), key=lambda item: item[1], reverse=True))


def overlay_graphs(merge_graph, pathway_graph, literature_kg, db_proteins):
    node_to_category = {}

    for node in merge_graph:
        if node in literature_kg and node not in pathway_graph:
            node_to_category[node] = 'literature'

        elif node in pathway_graph and node not in literature_kg:
            node_to_category[node] = 'pathway'

        else:
            node_to_category[node] = 'common'

    cc = max(nx.algorithms.components.connected_components(merge_graph.to_undirected()), key=len)

    non_connected = {
        node
        for node in merge_graph
        if node not in cc
    }

    merge_graph.remove_nodes_from(non_connected)

    nodes = []
    names = []
    node_colors = []

    for node in merge_graph:

        nodes.append(node.as_bel())

        if node_to_category[node] == "common":

            if (node.function in FUNCTIONS) and (node.name in db_proteins):

                names.append(node.name)
                node_colors.append(COLORS_CATEGORY[node_to_category[node]])

            else:
                node_colors.append(COLORS_CATEGORY['other'])
                names.append('')

        else:
            node_colors.append(COLORS_CATEGORY[node_to_category[node]])
            names.append('')

    graph = igraph.Graph()
    graph.add_vertices(nodes)
    graph.vs["color"] = node_colors
    graph.vs["label"] = names

    for source, target in merge_graph.edges():
        graph.add_edge(
            source.as_bel(),
            target.as_bel(),
        )

    visual_style = {
        "vertex_size": 6,
        "vertex_frame_color": "#dddedf",
        "vertex_label_size": 10,
        "vertex_label_color": "#363636",
        "vertex_label_dist": 1,
        "edge_color": "#D3D3D3",
        "edge_width": 0.3,
    }

    return graph, visual_style


def render_graph(merge_graph, pathway_graph, literature_kg, db_proteins, normalize_vals_dict, intracellular=True):
    node_to_category = {}

    for node in merge_graph:
        if node in literature_kg and node not in pathway_graph:
            node_to_category[node] = 'literature'

        elif node in pathway_graph and node not in literature_kg:
            node_to_category[node] = 'pathway'

        else:
            node_to_category[node] = 'common'

    cc = max(nx.algorithms.components.connected_components(merge_graph.to_undirected()), key=len)

    non_connected = {
        node
        for node in merge_graph
        if node not in cc
    }

    merge_graph.remove_nodes_from(non_connected)

    if intracellular:
        node_dict = get_node_counts(merge_graph)
    else:
        node_dict = get_node_counts(merge_graph, intracellular=False)

    nodes = []
    names = []
    node_colors = []
    node_sizes = []

    for node in merge_graph:

        nodes.append(node.as_bel())

        if node_to_category[node] == "common":

            if (node.function in FUNCTIONS) and (node.name in db_proteins):

                names.append(node.name)
                node_colors.append(COLORS_CATEGORY[node_to_category[node]])
                node_sizes.append(node_dict[node.name])

            else:
                names.append('')
                node_colors.append(COLORS_CATEGORY['other'])
                node_sizes.append(1)

        else:
            names.append('')
            node_colors.append(COLORS_CATEGORY[node_to_category[node]])
            node_sizes.append(1)

    norm = normalize(
        node_sizes,
        normalize_vals_dict
    )

    sizes = []

    for name, size in zip(names, norm):

        if name == '':
            sizes.append(4)
        else:
            sizes.append(size)

    graph = igraph.Graph()
    graph.add_vertices(nodes)
    graph.vs["color"] = node_colors
    graph.vs["label"] = names

    for source, target in merge_graph.edges():
        graph.add_edge(
            source.as_bel(),
            target.as_bel(),
        )

    visual_style = {
        "vertex_size": sizes,
        "vertex_frame_color": "#dddedf",
        "vertex_label_size": 10,
        "vertex_label_color": "#363636",
        "vertex_label_dist": 1,
        "edge_color": "#D3D3D3",
        "edge_width": 0.3,
    }

    return graph, visual_style
