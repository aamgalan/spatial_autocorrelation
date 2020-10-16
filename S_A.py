#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""S_A.py: The main routine for computing the fast spatial autocorrelation."""

import scipy as sp, numpy as np
from scipy.spatial import Delaunay
from scipy.cluster.hierarchy import dendrogram, linkage
from disjoint_set import DisjointSet

__author__ = "Anar Amgalan"
__credits__ = ["Anar Amgalan", "Steven S. Skiena", "Lilianne R. Mujica-Parodi"]
__version__ = "1.0"
__maintainer__ = "Anar Amgalan"


def get_edges_from_linkage_matrix(linkage_Z, n_node):
    """turn the outputs of linkage into a merge order

    Parameters
    ----------
    linkage_Z : linkage_Z output of linkage

    n_node : int

    Returns
    -------
    edge :
        list of edges that constitute the ... ?

    cluster_to_node :
        something else ... I forget what it is

    """

    cluster_to_node = {}
    edges = []

    for i_cluster, cluster in enumerate(linkage_Z):

        if cluster[0] < n_node:
            cluster_to_node[i_cluster + n_node] = int(cluster[0])

        else:
            cluster_to_node[i_cluster + n_node] = cluster_to_node[int(cluster[0])]

        edges.append(
            [
                int(node) if node < n_node else cluster_to_node[node]
                for node in cluster[:2]
            ]
        )

    return edges, cluster_to_node


def get_merge_order(z, method="single"):
    """compute merge order from 2D points and return list of edges.

    Parameters
    ----------
    coordinates : [[], ... ]

    method : string
        methods from the list ['single', 'median', 'average', 'complete']
        the methods that scipy linkage provides
        the default is 'single' linkage

    Returns
    -------
    list of pair of node IDs

    """

    # TODO check coordinates

    linkage_Z = linkage(z, method=method)
    merge_order_method = get_edges_from_linkage_matrix(linkage_Z, len(z))[0]

    return merge_order_method


def get_mst_merge_order(coordinates):
    """given the coordinates in euclidean space,
    this returns the order in which edges added to join them

    """

    delaunay_edges = []

    for simplex in Delaunay(coordinates).simplices:

        delaunay_edges.append(sorted([simplex[0], simplex[1]]))
        delaunay_edges.append(sorted([simplex[0], simplex[2]]))
        delaunay_edges.append(sorted([simplex[1], simplex[2]]))
        
    # some diagnostic outputs, that aren't necessary anymore
    # print(len(delaunay_edges), ' edges in Delaunay triangulation')
    # delaunay_edges = np.array(sorted(delaunay_edges))[np.arange(0, len(delaunay_edges), 2)]
    # print(len(delaunay_edges), ' edges in Delaunay triangulation')

    delaunay_edges_weighted = []

    for delaunay_edge in delaunay_edges:

        delaunay_edges_weighted.append(
            [
                delaunay_edge[0],
                delaunay_edge[1],
                np.linalg.norm(
                    coordinates[delaunay_edge[0]] - coordinates[delaunay_edge[1]]
                ),
            ]
        )
        
    # some diagnostic outputs, that aren't necessary anymore
    # print(delaunay_edges_weighted[:10])
    # print(delaunay_edges_weighted[-10:])
    
    vertices = range(len(coordinates))
    g = graph(vertices, delaunay_edges_weighted)
    mst = g.kruskal()
    mst_order = [edge[:2] for edge in mst]

    return mst_order


class graph(object):
    """graph object, on which the MST is calculated"""

    def __init__(self, vertices, edges):

        self.vertices = vertices

        #
        self.edges = sorted(edges, key=lambda edge: edge[2])

    def kruskal(self):

        # TODO what are these counters?
        edge_i, edge_n = (0, 0)
        self.ds = DisjointSet()
        self.mst = []
        while edge_n < len(self.vertices) - 1:
            vertex_1, vertex_2, weight = self.edges[edge_i]
            edge_i += 1
            cluster_1 = self.ds.find(vertex_1)
            cluster_2 = self.ds.find(vertex_2)
            if cluster_1 != cluster_2:
                self.ds.union(cluster_1, cluster_2)
                self.mst.append([vertex_1, vertex_2, weight])
                edge_n += 1

        return self.mst


def SkienaA(z, merge_order, rescaled_scalar=True):
    """the main algorithm computing the S_A statistic
    from z (z_i in the publication) and merge_order

    Parameters
    ----------
    z : List of N real numbers

    merge_order : a legitimate merge order defined on the N vertices 
    
    rescaled_scalar : Boolean, whether to return a single real value in [-1, 1]


    """

    # TODO check that merge_order is legitimate and matches z in size

    n_V = len(z)
    variance_sum = 0.0
    variance_mean_running = [0.0]
    ds = DisjointSet()
    # vertex_to_cluster = dict(zip(range(n_V), range(n_V)))
    cluster_info = []
    for vertex_id, vertex_quantity in zip(range(n_V), z):
        cluster_info.append(
            {"size": 1, "mean": vertex_quantity, "variance": 0}  # 'id':vertex_id,
        )
    # print(cluster_info)
    for vertex_1, vertex_2 in merge_order:

        cluster_1 = ds.find(vertex_1)
        cluster_2 = ds.find(vertex_2)
        ds.union(cluster_1, cluster_2)

        cluster_combined = cluster_info[cluster_1]

        cluster_info_1 = cluster_info[cluster_1]
        cluster_info_2 = cluster_info[cluster_2]

        mean_12 = (
            cluster_info_1["size"] * cluster_info_1["mean"]
            + cluster_info_2["size"] * cluster_info_2["mean"]
        ) / (cluster_info_1["size"] + cluster_info_2["size"])

        mean_change_1 = mean_12 - cluster_info_1["mean"]
        mean_change_2 = mean_12 - cluster_info_2["mean"]

        cluster_variance_change_1 = (
            cluster_info_1["size"] * mean_change_1 * mean_change_1
        )
        cluster_variance_change_2 = (
            cluster_info_2["size"] * mean_change_2 * mean_change_2
        )
        cluster_info_1["variance"] += cluster_variance_change_1
        cluster_info_2["variance"] += cluster_variance_change_2
        variance_sum += cluster_variance_change_1 + cluster_variance_change_2

        cluster_size_new = cluster_info_1["size"] + cluster_info_2["size"]
        cluster_info_1["size"] = cluster_size_new
        cluster_info_2["size"] = cluster_size_new

        cluster_info_1["mean"] = mean_12
        cluster_info_2["mean"] = mean_12

        variance_mean_running.append(variance_sum / n_V)
        
        # some diagnostic outputs, that aren't necessary anymore
        # print(cluster_combined['size'], cluster_combined['variance'], np.sqrt(variance_sum/n_V))
        # print(cluster_info_2['size'], cluster_info_2['variance'])

    if rescaled_scalar == True:

        skienaA = np.mean(variance_mean_running) / variance_mean_running[-1]

        return 2.0 * (1.0 - skienaA) - 1.0

    return variance_mean_running
