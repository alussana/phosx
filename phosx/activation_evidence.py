#!/usr/bin/env python3

from os import path
import sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
import itertools

from phosx.utils import read_pssms, read_pssm_score_quantiles, hdf5_to_dict, sync_dataframe_with_names, decay_from_1
from phosx.pssms import quantile_scaling, max_pssm_scoring


class KinaseNetwork:
    def __init__(self):
        # Initialize an empty network with dictionaries for nodes and edges
        self.nodes = {}  # node_id -> attributes dictionary
        self.edges = {}  # (source_id, target_id) -> complementarity
        self.name_to_id = {}  # name -> node_id mapping for fast lookup
        
    def add_node(self, node_id, name, full_seq, kinase_domain_seq, aloop_seq, activity_score):
        """Add a node to the network with all required attributes"""
        self.nodes[node_id] = {
            'name': name,
            'full_seq': full_seq,
            'kinase_domain_seq': kinase_domain_seq,
            'aloop_seq': aloop_seq,
            'activity_score': activity_score
        }
        # Update the name-to-id mapping
        self.name_to_id[name] = node_id
        
    def add_edge(self, source_id, target_id, complementarity):
        """Add a directed edge between two nodes with complementarity score"""
        if source_id in self.nodes and target_id in self.nodes:
            self.edges[(source_id, target_id)] = complementarity
        else:
            raise ValueError("Both source and target nodes must exist in the network")
            
    def get_node(self, node_id):
        """Return attributes of a specific node by ID"""
        return self.nodes.get(node_id, None)
    
    def get_node_by_name(self, name):
        """Return node ID and attributes for a node with the given name"""
        node_id = self.name_to_id.get(name)
        if node_id is not None:
            return node_id, self.nodes[node_id]
        return None, None  # Return None, None if no node is found
    
    def get_edge(self, source_id, target_id):
        """Return complementarity of a specific edge"""
        return self.edges.get((source_id, target_id), None)
    
    def get_neighbors(self, node_id):
        """Return list of nodes that are targets of edges from this node"""
        return [target_id for (source, target_id) in self.edges.keys() 
                if source == node_id]
    
    def remove_node(self, node_id):
        """Remove a node and all its associated edges"""
        if node_id in self.nodes:
            # Remove from name_to_id mapping first
            node_name = self.nodes[node_id]['name']
            del self.name_to_id[node_name]
            # Remove the node
            del self.nodes[node_id]
            # Remove all edges involving this node
            self.edges = {k: v for k, v in self.edges.items() 
                         if k[0] != node_id and k[1] != node_id}
    
    def remove_edge(self, source_id, target_id):
        """Remove a specific edge"""
        if (source_id, target_id) in self.edges:
            del self.edges[(source_id, target_id)]


def compute_complementarity(k1:str, k2:str, assigned_substrates: pd.DataFrame):
    """
    Compute 1 - Jaccard index between two columns of a binary DataFrame.
    
    Parameters:
    k1 (str): Name of first column
    k2 (str): Name of second column
    assigned_substrates (pd.DataFrame): DataFrame with binary values (0s and 1s)
    
    Returns:
    float: 1 - Jaccard index value, or 0 if columns not found
    """
    # Check if the columns exist in the DataFrame
    if k1 not in assigned_substrates.columns or k2 not in assigned_substrates.columns:
        return 0.0

    # Get the sets of indices where each column has 1s
    set1 = set(assigned_substrates[assigned_substrates[k1] == 1].index)
    set2 = set(assigned_substrates[assigned_substrates[k2] == 1].index)
    
    # Calculate intersection and union
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    # Handle division by zero case
    if union == 0:
        return 0.0
    
    # Return 1 - Jaccard index
    return 1 - (intersection / union)


def compute_redundancy(k1:str, k2:str, assigned_substrates: pd.DataFrame):
    """
    Compute Jaccard index between two columns of a binary DataFrame.
    If the columns are the same or not found, return 0.0
    
    Parameters:
    k1 (str): Name of first column
    k2 (str): Name of second column
    assigned_substrates (pd.DataFrame): DataFrame with binary values (0s and 1s)
    
    Returns:
    float: Jaccard index value, or 0 if columns are the same or not found
    """
    # Check if columns are the same
    if k1 == k2:
        return 0.0
    
    # Check if the columns exist in the DataFrame
    if k1 not in assigned_substrates.columns or k2 not in assigned_substrates.columns:
        return 0.0

    # Get the sets of indices where each column has 1s
    set1 = set(assigned_substrates[assigned_substrates[k1] == 1].index)
    set2 = set(assigned_substrates[assigned_substrates[k2] == 1].index)
    
    # Calculate intersection and union
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    # Handle division by zero case
    if union == 0:
        return 0.0
    
    # Return Jaccard index
    return intersection / union


def build_kin_net(adj_df, assigned_substrates_df, activity_df, kinase_metadata_dict):
    """
    Build a KinaseNetwork from a pandas DataFrame adjacency matrix.
    
    Parameters:
    adj_df (pandas.DataFrame): Adjacency matrix where columns are source nodes,
                              rows are target nodes, and values are 0 or 1
    
    Returns:
    KinaseNetwork: A populated KinaseNetwork object
    """
    # Create an empty KinaseNetwork
    network = KinaseNetwork()
    
    # Get source nodes from columns and target nodes from index
    source_nodes = adj_df.columns.tolist()
    target_nodes = adj_df.index.tolist()
    
    # Create a set of all unique node names
    all_nodes = set(source_nodes).union(set(target_nodes))
    
    # Add all nodes to the network with default attributes
    for i, node_name in enumerate(all_nodes):
        activity = activity_df.loc[activity_df.index==node_name, "Activity Score"].values
        if len(activity) == 0:
            activity = np.nan
        else:
            activity = activity[0]
        network.add_node(
            node_id=i,  # Using integer IDs starting from 0
            name=node_name,
            full_seq=kinase_metadata_dict["full_seq"][node_name],
            kinase_domain_seq=kinase_metadata_dict["kinase_domain_seq"][node_name],
            aloop_seq=kinase_metadata_dict["aloop_seq"][node_name],
            activity_score=activity,
        )
    
    # Add edges based on the adjacency matrix
    for source_name in source_nodes:
        for target_name in target_nodes:
            if adj_df.loc[target_name, source_name] == 1:
                # Use get_node_by_name to get the IDs efficiently
                source_id, _ = network.get_node_by_name(source_name)
                target_id, _ = network.get_node_by_name(target_name)
                complementarity = compute_complementarity(source_name, target_name, assigned_substrates_df)
                network.add_edge(source_id, target_id, complementarity=complementarity)
    
    return network


def compute_activation_evidence(
    activity_df: pd.DataFrame,
    s_t_assigned_substrates_df: pd.DataFrame,
    y_assigned_substrates_df: pd.DataFrame,
    s_t_pssm_h5_file: str = str(path.join(path.dirname(__file__), "../phosx/data/S_T_PSSMs.h5")),
    s_t_pssm_score_quantiles_h5_file: str = str(path.join(path.dirname(__file__), "../phosx/data/S_T_PSSM_score_quantiles.h5")),
    y_pssm_h5_file: str = str(path.join(path.dirname(__file__), "../phosx/data/Y_PSSMs.h5")),
    y_pssm_score_quantiles_h5_file: str = str(path.join(path.dirname(__file__), "../phosx/data/Y_PSSM_score_quantiles.h5")),
    metadata_h5_file: str = str(path.join(path.dirname(__file__), "../phosx/data/kinase_metadata.h5")),
    s_t_quantile_threshold: int = 0.95,
    y_quantile_threshold: int = 0.90,
    redundacy_threshold: float = 0.75,
    decay_factor: float = 64,
    upregulation: bool = True,
    n_proc: int = 1,
    plot_figures: bool = False,
    out_plot_dir: str = "phosx_output",
    out_path=None,
):
    print("     Loading input objects    : ", file=sys.stderr, end="")
    s_t_pssm_df_dict = read_pssms(s_t_pssm_h5_file)
    s_t_pssm_bg_scores_df = read_pssm_score_quantiles(s_t_pssm_score_quantiles_h5_file)
    y_pssm_df_dict = read_pssms(y_pssm_h5_file)
    y_pssm_bg_scores_df = read_pssm_score_quantiles(y_pssm_score_quantiles_h5_file)
    kinase_metadata_dict = hdf5_to_dict(metadata_h5_file)
    aloop_df = pd.DataFrame.from_dict(kinase_metadata_dict["aloop_seq"], orient="index", columns=["Sequence"])
    print("DONE", file=sys.stderr)

    if upregulation:
        activity_series = activity_df.loc[activity_df["Activity Score"] >= 0, "Activity Score"]
    else:
        activity_series = -activity_df.loc[activity_df["Activity Score"] <= 0, "Activity Score"]
    
    kinase_list = list(activity_series.index)

    aloop_df = aloop_df.loc[[i in kinase_list for i in aloop_df.index]]
    seq_series = aloop_df["Sequence"].dropna()

    s_t_pssm_df_dict = {k: v for k, v in s_t_pssm_df_dict.items() if k in kinase_list}
    s_t_pssm_bg_scores_df = s_t_pssm_bg_scores_df.loc[:,[k in kinase_list for k in s_t_pssm_bg_scores_df.columns]]
    y_pssm_df_dict = {k: v for k, v in y_pssm_df_dict.items() if k in kinase_list}
    y_pssm_bg_scores_df = y_pssm_bg_scores_df.loc[:,[k in kinase_list for k in y_pssm_bg_scores_df.columns]]
    s_t_assigned_substrates_df = s_t_assigned_substrates_df.loc[:,[k in kinase_list for k in s_t_assigned_substrates_df.columns]]
    y_assigned_substrates_df = y_assigned_substrates_df.loc[:,[k in kinase_list for k in y_assigned_substrates_df.columns]]   

    # Score A-loops sequences with each Ser/Thr PSSM
    arg1 = [seq_series[i] for i in range(len(seq_series))]
    arg2 = [s_t_pssm_df_dict for i in range(len(seq_series))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(
            max_pssm_scoring,
            tqdm(
                zip(arg1, arg2),
                ncols=80,
                total=len(seq_series),
                desc=" Scoring A-loops w/ S/T PSSMs ",
            ),
        )
    s_t_pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    s_t_pssm_scoring_df.index = seq_series.index
    # quantile-scale the PSSM scores for each kinase
    s_t_sorted_bg_scores_dict = {kinase: np.sort(s_t_pssm_bg_scores_df[kinase].values) for kinase in s_t_pssm_bg_scores_df.columns}
    s_t_pssm_scoring_scaled01_df = s_t_pssm_scoring_df.apply(
        quantile_scaling,
        args=[s_t_sorted_bg_scores_dict],
        axis=0,
    )
    s_t_binarized_pssm_scores = s_t_pssm_scoring_scaled01_df.apply(lambda x: x.map(lambda y: 1 if y > s_t_quantile_threshold else 0 ))

    # Score A-loops sequences with each Tyr PSSM
    arg1 = [seq_series[i] for i in range(len(seq_series))]
    arg2 = [y_pssm_df_dict for i in range(len(seq_series))]
    with Pool(processes=n_proc) as pool:
        dfs_list = pool.starmap(
            max_pssm_scoring,
            tqdm(
                zip(arg1, arg2),
                ncols=80,
                total=len(seq_series),
                desc="  Scoring A-loops w/ Y PSSMs  ",
            ),
        )
    y_pssm_scoring_df = pd.concat(dfs_list, axis=1).T
    y_pssm_scoring_df.index = seq_series.index
    # quantile-scale the PSSM scores for each kinase
    y_sorted_bg_scores_dict = {kinase: np.sort(y_pssm_bg_scores_df[kinase].values) for kinase in y_pssm_bg_scores_df.columns}
    y_pssm_scoring_scaled01_df = y_pssm_scoring_df.apply(
        quantile_scaling,
        args=[y_sorted_bg_scores_dict],
        axis=0,
    )
    y_binarized_pssm_scores = y_pssm_scoring_scaled01_df.apply(lambda x: x.map(lambda y: 1 if y > y_quantile_threshold else 0 ))

    print("      Analyzing network       : ", file=sys.stderr, end="")
    # build kinase network (do not consider self-loops)
    adj_df = pd.concat([s_t_binarized_pssm_scores, y_binarized_pssm_scores], axis=1)
    adj_df = sync_dataframe_with_names(adj_df, kinase_list, fill_value=0)
    assigned_substrates_df = pd.concat([s_t_assigned_substrates_df, y_assigned_substrates_df]).replace(np.nan, 0)
    network = build_kin_net(adj_df, assigned_substrates_df, activity_df, kinase_metadata_dict)

    # get a, the activity vector (activity for each kinase based on phosphosites)
    a = activity_series.values
    A = np.diag(a) # make a diagonal matrix from the activity vector

    # get Dt, the transposed adjacency matrix (kinase --> A-loop-targeting kinase)
    Dt = adj_df.transpose().values
    np.fill_diagonal(Dt, 0) # remove self-loops
    
    # get C, the redundancy/correlation matrix (indicating which kinases have similar substrate sets)
    pairs = list(itertools.product(kinase_list, repeat=2)) # create all possible combinations of names (including same-name pairs)
    scores = [compute_redundancy(pair[0], pair[1], assigned_substrates_df) for pair in pairs] # compute redundancy scores for all pairs of kinases
    n = len(kinase_list)
    score_matrix = [scores[i:i+n] for i in range(0, len(scores), n)] # Reshape the scores into a square matrix
    df = pd.DataFrame(score_matrix, 
                     index=kinase_list, 
                     columns=kinase_list)
    df = df.round(3)
    df[df < redundacy_threshold] = 0
    df[df != 0] = 1
    C = df.values

    # get e, the upstream activation evidence vector (extent of upstream activation evidence for each kinase)
    E = Dt @ A
    e = np.max(E, axis=1) # max upstream activation evidence
    #E = np.diag(e) # make a diagonal matrix from the evidence vector

    # get B, the outer product of the activity vector and its reciprocal
    B = np.outer(a, 1/a).round(4)

    # get F, the outer product of the upstream activation evidence vector and its reciprocal
    with np.errstate(divide='ignore', invalid='ignore'):
        F = np.outer(e, 1/e).round(4)
    F[np.isinf(F)] = 1 # if competing kinase has no upstream activation evidence, we don't consider it
    F[np.isnan(F)] = 1 # if both the kinase and competing kinase have no upstream activation evidence, we don't consider it
    F[F > 1] = 1 # if the kinase has stronger upstream activation evidence than the competing kinase, we don't consider it

    # get decayed_B, [...]
    v_decay_from_1 = np.vectorize(decay_from_1)
    X = v_decay_from_1(B, decay_factor=decay_factor).round(4)

    # get M, the modifier matrix
    M = (1 - ((C * X * (1 - F)))).round(4)
    m =  np.min(M, axis=1) # choose modifiers for each kinase considering the strongest competing upstream activation evidence

    # get z, the modified activity vector
    z = a * m

    z_series = pd.Series(z, index=activity_series.index)
    if upregulation == False:
        z_series = -z_series

    print("DONE", file=sys.stderr)
    
    return z_series
