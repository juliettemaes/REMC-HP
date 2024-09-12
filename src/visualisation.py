import networkx as nx
import matplotlib.pyplot as plt

class LatticeHPGraph:
    def __init__(self, lat):
        """
        Initialize the LatticeHPGraph class.
        
        :param lat: A dictionary containing two keys:
                    - 'lattice': A 2D array where positions are either 0 (empty) or a number representing an amino acid.
                    - 'sequence': A string of 'H' and 'P' representing the corresponding amino acids in the lattice.
        """
        self.lattice = lat['lattice']
        self.sequence = lat['sequence']
        self.graph = nx.Graph()

    def _get_neighbors(self, row, col):
        """
        Get the valid neighboring positions (up, down, left, right) of a given cell in the lattice.
        
        :param row: Row index of the current cell.
        :param col: Column index of the current cell.
        :return: A list of valid neighboring positions as tuples (row, col).
        """
        neighbors = []
        if row > 0:  # Up
            neighbors.append((row - 1, col))
        if row < len(self.lattice) - 1:  # Down
            neighbors.append((row + 1, col))
        if col > 0:  # Left
            neighbors.append((row, col - 1))
        if col < len(self.lattice[0]) - 1:  # Right
            neighbors.append((row, col + 1))
        return neighbors

    def create_graph(self):
        """
        Create a graph representation of the lattice using the HP model.
        """
        for row in range(len(self.lattice)):
            for col in range(len(self.lattice[0])):
                # Check if the position has an amino acid (non-zero)
                amino_acid_pos = self.lattice[row][col]
                if amino_acid_pos != 0:
                    # Convert the amino acid position to an index for the sequence string (1-based -> 0-based)
                    amino_acid = self.sequence[amino_acid_pos - 1]
                    # Add the node with the amino acid label ('H' or 'P')
                    self.graph.add_node((row, col), label=amino_acid)
                    
                    # Connect this node to its valid neighbors that also have amino acids
                    for neighbor in self._get_neighbors(row, col):
                        neighbor_row, neighbor_col = neighbor
                        neighbor_pos = self.lattice[neighbor_row][neighbor_col]
                        if neighbor_pos != 0:
                            neighbor_amino_acid = self.sequence[neighbor_pos - 1]
                            # Add an edge between the current node and the neighboring node
                            self.graph.add_edge((row, col), (neighbor_row, neighbor_col))

    def draw_graph(self):
        """
        Draw the graph using networkx and matplotlib with 'H' and 'P' labels.
        """
        pos = {(row, col): (col, -row) for row in range(len(self.lattice)) for col in range(len(self.lattice[0])) if self.lattice[row][col] != 0}
        
        labels = nx.get_node_attributes(self.graph, 'label')
        
        # Draw the graph
        nx.draw(self.graph, pos, with_labels=True, labels=labels, node_size=500, node_color='lightblue', font_size=12, font_weight='bold')
        plt.show()


# Example usage:
# lat object containing a 5x5 lattice with a HP sequence
lat = {
    'lattice': [
        [0, 0, 1, 0, 0],
        [0, 0, 2, 0, 0],
        [0, 3, 0, 4, 0],
        [0, 0, 5, 0, 0],
        [0, 0, 0, 0, 0]
    ],
    'sequence': 'HPHPP'
}

# Create the lattice graph
lattice_graph = LatticeHPGraph(lat)
lattice_graph.create_graph()

# Draw the graph
lattice_graph.draw_graph()
