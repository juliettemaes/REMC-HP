from typing import List, Optional, Dict
from variables import AA_DICT

class Sequence:
    """
    A class to represent a sequence of amino acids and their properties in a lattice-based protein structure.

    Attributes:
        sequence (str): The amino acid sequence.
        hp_sequence (str): The hydrophobic polar sequence.
        length (int): The length of the amino acid sequence.
        aa_coord (List[Dict[str, Optional[int]]]): A list of dictionaries containing the coordinates and properties of each amino acid.
    """

    def __init__(self, sequence: Optional[str] = None, hp_sequence: Optional[str] = None) -> None:
        """
        Initialize the Sequence object.

        Args:
            sequence (Optional[str]): The amino acid sequence.
            hp_sequence (Optional[str]): The hydrophobic polar sequence.

        Raises:
            ValueError: If an invalid amino acid is found in the sequence.
        """
        if sequence is not None and hp_sequence is None:
            for aa in sequence:
                if aa not in AA_DICT:
                    raise ValueError(f"Invalid amino acid: {aa}")
            self.sequence = sequence
            self.length = len(sequence)
            self.hp_sequence = self.HP_convert()
            self.aa_coord = self.aa_coord_generator()
        elif hp_sequence is not None and sequence is None:
            self.hp_sequence = hp_sequence
            self.sequence = hp_sequence
            self.length = len(hp_sequence)
            self.aa_coord = self.aa_coord_generator()

    def __str__(self) -> str:
        """
        Return the string representation of the sequence.

        Returns:
            str: The amino acid sequence.
        """
        return self.sequence

    def HP_convert(self) -> str:
        """
        Convert the amino acid sequence to its hydrophobic-polar model representation.

        Returns:
            str: The hydrophobic-polar sequence.
        """
        return "".join(AA_DICT[aa] for aa in self.sequence)

    def aa_coord_generator(self) -> List[Dict[str, Optional[int]]]:
        """
        Generate the initial coordinates for each amino acid in the sequence.

        Returns:
            List[Dict[str, Optional[int]]]: A list of dictionaries with initial coordinates and properties for each amino acid.
        """
        aa_coord = []
        for i, aa in enumerate(self.hp_sequence):
            dict_aa = {
                "type": aa,
                "index": i,
                "chain_index": i + 1,
                "x": None,
                "y": None
            }
            aa_coord.append(dict_aa)
        return aa_coord

    def aa_coord_update(self, index: int, x: int, y: int) -> None:
        """
        Update the coordinates of the amino acid at a specific index.

        Args:
            index (int): The index of the amino acid to update.
            x (int): The new x-coordinate.
            y (int): The new y-coordinate.
        """
        self.aa_coord[index]["x"] = x
        self.aa_coord[index]["y"] = y
