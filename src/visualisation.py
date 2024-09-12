"""
Module for creating and visualizing a graph representation of the protein in the lattice.
"""

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from typing import Dict


def create_lattice_graph(lattice) -> nx.Graph:
    """
    Create a graph representation of the lattice using NetworkX.

    Args:
        lattice: The lattice object.

    Returns:
        A NetworkX graph where nodes represent amino acids and edges represent
        connections between consecutive amino acids in the sequence.
        Hydrophobic amino acids are colored green and polar amino acids are colored red.
    """
    G = nx.Graph()  # Create a NetworkX graph
    
    # Add amino acid nodes with the 'type' attribute
    for aa in lattice.sequence.aa_coord:
        G.add_node(aa['chain_index'], type=aa['type'])

    # Add edges only between consecutive amino acids in the sequence
    for i in range(len(lattice.sequence.aa_coord) - 1):
        aa_current = lattice.sequence.aa_coord[i]
        aa_next = lattice.sequence.aa_coord[i + 1]
        
        # Add an edge between the current and the next amino acid in the sequence
        G.add_edge(aa_current['chain_index'], aa_next['chain_index'], connected=True)
    
    return G

def visualize_lattice_graph(lattice) -> None:
    """
    Visualize the lattice graph using Matplotlib.

    Args:
        lattice: The lattice object containing sequence and grid information.

    Returns:
        None
    """
    G = create_lattice_graph(lattice)

    # Define node positions for visualization
    pos = {node: (lattice.sequence.aa_coord[node - 1]['x'], lattice.sequence.aa_coord[node - 1]['y']) for node in G.nodes}

    # Define node colors based on hydrophobicity
    node_colors = ['green' if G.nodes[node]['type'] == 'H' else 'red' for node in G.nodes]

    # Define edge colors based on connection
    edge_colors = ['green' if G.edges[edge]['connected'] else 'black' for edge in G.edges]

    # Draw the graph
    nx.draw(G, pos, with_labels=True, node_color=node_colors, edge_color=edge_colors, node_size=1000, font_size=8)

    # Display the graph
    plt.show()
