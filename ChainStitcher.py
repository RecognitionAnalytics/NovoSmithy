import py3Dmol
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.optimize import minimize
from Bio import PDB
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Atom, Residue, Chain, Model, Structure
from Bio.PDB.vectors import Vector
import sys
import os
import shutil
import subprocess
from scipy.spatial.transform import Rotation
from IPython.display import clear_output

def ClearFolder(inputFolder):
  for filename in os.listdir(inputFolder):
      file_path = os.path.join(inputFolder, filename)
      try:
          if os.path.isfile(file_path) or os.path.islink(file_path):
              os.unlink(file_path)
          elif os.path.isdir(file_path):
              shutil.rmtree(file_path)
      except Exception as e:
          print('Failed to delete %s. Reason: %s' % (file_path, e))

def load_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        while True:
            comment = file.readline().strip()
            if not comment:
                break
            sequence = file.readline().strip()
            if 'overall_confidence' in comment:
                #comment format -> ">working_04, id=1, T=0.044000000000000004, seed=37763, overall_confidence=0.2417, ligand_confidence=0.2186, seq_rec=0.2693"
                comment = 'file=' + comment[1:]
                json_data = {key: (value).strip() for key, value in [item.split('=') for item in comment.split(', ')]}
                json_data['overall_confidence']= float(json_data['overall_confidence'])
                json_data['ligand_confidence']= float(json_data['ligand_confidence'])
                json_data['seq_rec']= float(json_data['seq_rec'])
                sequences.append((json_data , sequence))
    return sequences 

class PolyGLinker:

    def __init__(self):
        # Standard amino acid geometry
        self.ca_c_length = 1.52  # Å
        self.c_n_length = 1.33  # Å
        self.n_ca_length = 1.46  # Å
        self.ca_ca_distance = 3.8  # Å
        self.ca_c_n_angle = 117.2  # degrees
        self.c_n_ca_angle = 121.7  # degrees
        self.peptide_omega = 180.0  # degrees

    def _optimize_control_points(
        self, start_pos, end_pos, start_dir, end_dir, num_residues
    ):
        """Optimize control points to achieve desired path length."""
        target_length = (num_residues - 1) * self.ca_ca_distance
        distance = np.linalg.norm(end_pos - start_pos)

        def objective(x):
            # x contains scaling factors for control points
            scale1, scale2 = x
            p0 = start_pos
            p1 = start_pos + start_dir * (distance * scale1)
            p2 = end_pos - end_dir * (distance * scale2)
            p3 = end_pos

            path_length = self._calculate_path_length((p0, p1, p2, p3))
            return (path_length - target_length) ** 2

        # Optimize scaling factors
        result = minimize(objective, [0.33, 0.33], bounds=[(0.1, 0.9), (0.1, 0.9)])
        scale1, scale2 = result.x
        scale1, scale2 = (0.33, 0.33)
        # Return optimized control points
        return (
            start_pos,
            start_pos + start_dir * (distance * scale1),
            end_pos - end_dir * (distance * scale2),
            end_pos,
        )

    def _calculate_residue_direction(self, residue):
        """Calculate the direction vector of a residue's backbone."""
        ca_pos = residue["CA"].get_coord()
        if "C" in residue:
            c_pos = residue["C"].get_coord()
            direction = c_pos - ca_pos
        else:
            n_pos = residue["N"].get_coord()
            direction = ca_pos - n_pos
        return direction / np.linalg.norm(direction)

    def calculateIdealPath(self, start_res, end_res):
        start_direction = self._calculate_residue_direction(start_res)
        end_direction = self._calculate_residue_direction(end_res)

        start_ca = start_res["CA"].get_coord()
        end_ca = end_res["CA"].get_coord()

        distance = np.linalg.norm(end_ca - start_ca)
        # Generate a test path for the Bezier curve
        test_control_points = (
            start_ca,
            start_ca + start_direction * (distance / 3),
            end_ca - end_direction * (distance / 3),
            end_ca,
        )
        test_path_length = self._calculate_path_length(test_control_points)
        num_residues = int(np.ceil(test_path_length / self.ca_ca_distance))

        # Optimize control points for target path length
        control_points = self._optimize_control_points(
            start_ca, end_ca, start_direction, end_direction, num_residues
        )

        # Generate points along optimized Bezier curve with equal arc length spacing
        # First, calculate a finely sampled path
        fine_t = np.linspace(0, 1, 1000)
        fine_path = np.array(
            [self._evaluate_bezier(control_points, ti) for ti in fine_t]
        )
        self.ca_ca_distance

        path = np.zeros((num_residues, 3))
        lastDist = 0
        currentDist = 0
        cc = 0
        for i in range(len(fine_path)):
            currentDist += np.linalg.norm(fine_path[i] - fine_path[i - 1])
            if currentDist - lastDist > self.ca_ca_distance:
                path[cc] = fine_path[i]
                cc += 1
                lastDist = currentDist

        return path, end_direction

    def generate_linker(self, segment1, segment2):
        """Generate backbone coordinates for a linker using an optimized Bezier curve."""
        start_res = segment1.c_res
        end_res = segment2.n_res

        path, end_direction = self.calculateIdealPath(start_res, end_res)

        # Create a new chain for the linker
        linker_chain = Chain.Chain("L")  # 'L' for linker
        for i in range(1, path.shape[0]):
            residue = Residue.Residue((" ", i, " "), "GLY", "")

            glyStart = path[i - 1]
            glyEnd = path[i]
            # Define atom coordinates for glycine based on position and direction
            # Calculate direction vector from previous to current position
            direction = glyEnd - glyStart
            direction = direction / np.linalg.norm(direction)

            # Calculate perpendicular vectors to create a coordinate system
            # Use cross product with arbitrary vector (0,0,1) to get perpendicular vector
            if np.allclose(direction, [0, 0, 1]) or np.allclose(direction, [0, 0, -1]):
                perp1 = np.array([1, 0, 0])
            else:
                perp1 = np.cross(direction, [0, 0, 1])
                perp1 = perp1 / np.linalg.norm(perp1)

            # Create atom positions
            ca_pos = glyStart
            c_pos = ca_pos + direction * self.ca_c_length
            n_pos = glyEnd - direction * self.n_ca_length
            o_pos = c_pos + perp1 * 1.23  # Typical C-O bond length

            # Create atoms
            ca_atom = Atom.Atom("CA", ca_pos, 0, 1.0, " ", "CA", 0, "C")
            c_atom = Atom.Atom("C", c_pos, 0, 1.0, " ", "C", 0, "C")
            n_atom = Atom.Atom("N", n_pos, 0, 1.0, " ", "N", 0, "N")
            o_atom = Atom.Atom("O", o_pos, 0, 1.0, " ", "O", 0, "O")

            # Add atoms to residue
            residue.add(n_atom)
            residue.add(ca_atom)
            residue.add(c_atom)
            residue.add(o_atom)
            # Add residue to chain
            linker_chain.add(residue)

        return linker_chain

    def _evaluate_bezier(self, control_points, t):
        """Evaluate a cubic Bezier curve at parameter t."""
        p0, p1, p2, p3 = control_points
        return (
            (1 - t) ** 3 * p0
            + 3 * (1 - t) ** 2 * t * p1
            + 3 * (1 - t) * t**2 * p2
            + t**3 * p3
        )

    def _calculate_path_length(self, control_points, num_points=100):
        """Calculate the total length of the Bezier curve."""
        t = np.linspace(0, 1, num_points)
        points = np.array([self._evaluate_bezier(control_points, ti) for ti in t])

        # Calculate total length as sum of segments
        segments = np.diff(points, axis=0)
        lengths = np.sqrt(np.sum(segments**2, axis=1))
        return np.sum(lengths)


class Segment:
    def get_named_atom(self, residue, atom_name):
        """
        Retrieves a specific atom from a residue.

        Args:
            residue (Bio.PDB.Residue.Residue): The residue to search within.
            atom_name (str): The name of the atom to find (e.g., 'N', 'CA', 'C', 'O').

        Returns:
            Bio.PDB.Atom.Atom or None: The Atom object if found, otherwise None.
        """
        if atom_name in residue:
            return residue[atom_name]
        return None

    """Represents a contiguous segment of residues. multiple segments can be made from 1 chain"""

    def __init__(self, residues, original_chain_id, segment_idx_in_chain):
        if not residues:
            raise ValueError("Segment cannot be empty")
        self.residues = list(residues)
        self.original_chain_id = original_chain_id

        first_res_id_str = (
            f"{residues[0].id[0].strip()}{residues[0].id[1]}{residues[0].id[2].strip()}"
        )
        last_res_id_str = f"{residues[-1].id[0].strip()}{residues[-1].id[1]}{residues[-1].id[2].strip()}"
        self.id = f"Seg-{original_chain_id}-{segment_idx_in_chain}({first_res_id_str}_to_{last_res_id_str})"

        self.n_res = self.residues[0]
        self.c_res = self.residues[-1]

        self.n_atom = self.get_named_atom(self.n_res, "N")
        self.c_atom = self.get_named_atom(self.c_res, "C")

        if self.n_atom is None:
            raise ValueError(
                f"Segment {self.id} missing N-terminal atom for residue {self.n_res.id}."
            )
        if self.c_atom is None:
            raise ValueError(
                f"Segment {self.id} missing C-terminal atom for residue {self.c_res.id}."
            )


class ChainStich:
    def __init__(self, model, excluded_chains=None, connection_threshold_Ang=10):
        """
        Initialize linker for connecting helix chains.

        Args:
            model: BioPython Model object containing the helices
            excluded_chains: List of chain IDs to exclude (e.g., hemes)
        """
        self.model = model
        self.excluded_chains = excluded_chains or []
        self.chain_ends = {}  # Store terminal residues for each chain
        self.connections = []  # Store pairs of chain ends to connect

        #
        # PEPTIDE_BOND_THRESHOLD: Distance (Angstroms) between C of res(i) and N of res(i+1)
        # to be considered a broken peptide bond (i.e., a new segment starts).
        self.PEPTIDE_BOND_THRESHOLD = 2.0  # Angstroms

        # CONNECTION_THRESHOLD: Maximum distance (Angstroms) between terminal atoms
        # of segments to consider them for connection.
        self.CONNECTION_THRESHOLD = connection_threshold_Ang  # Angstroms

        # Initialize structure
        self.SegmentChains()

    def displaySegments(self):
        """
        Displays each segment as a ribbon with a different color using py3Dmol.
        Creates temporary chains from segments and renders the structure.
        """
        # Create a new structure for visualization
        struct = Structure.Structure("segmented")
        model = Model.Model(0)
        struct.add(model)

        # Generate a list of distinct colors
        colors = plt.cm.rainbow(np.linspace(0, 1, len(self.segments)))
        colors = [
            f"rgb({int(r*255)},{int(g*255)},{int(b*255)})" for r, g, b, _ in colors
        ]

        # Create chains from segments
        for i, segment in enumerate(self.segments):
            chain_id = chr(65 + (i % 26))  # A-Z (looping if needed)
            chain = Chain.Chain(chain_id)

            for residue in segment.residues:
                # Create a deep copy of the residue
                new_res = Residue.Residue(residue.id, residue.resname, residue.segid)
                for atom in residue:
                    new_atom = Atom.Atom(
                        atom.name,
                        atom.coord,
                        atom.bfactor,
                        atom.occupancy,
                        atom.altloc,
                        atom.fullname,
                        atom.serial_number,
                        atom.element,
                    )
                    new_res.add(new_atom)
                chain.add(new_res)
            model.add(chain)

        # Save structure to temporary PDB file
        io = PDBIO()
        io.set_structure(struct)
        temp_pdb = "temp_segmented.pdb"
        io.save(temp_pdb)

        # Visualize with py3Dmol
        view = py3Dmol.view(width=800, height=600)
        view.addModel(open(temp_pdb, "r").read(), "pdb")

        # Color chains based on segment
        for i, segment in enumerate(self.segments):
            chain_id = chr(65 + (i % 26))
            view.setStyle({"chain": chain_id}, {"cartoon": {"color": colors[i]}})

        view.zoomTo()
        view.setBackgroundColor("white")
        return view.show()

    def SegmentChains(self):

        all_segments = []

        for chain in self.model:
            residues_in_chain = [
                res
                for res in chain.get_residues()
                if PDB.is_aa(res, standard=True) or res.id[0] != " "
            ]
            if not residues_in_chain:
                continue

            current_segment_residues = []
            segment_counter_for_chain = 0
            for i, res in enumerate(residues_in_chain):
                if not current_segment_residues:
                    current_segment_residues.append(res)
                else:
                    prev_res = current_segment_residues[-1]
                    if "C" in prev_res and "N" in res:
                        dist = np.linalg.norm(
                            prev_res["C"].get_coord() - res["N"].get_coord()
                        )
                        if dist > self.PEPTIDE_BOND_THRESHOLD:
                            all_segments.append(
                                Segment(
                                    current_segment_residues,
                                    chain.id,
                                    segment_counter_for_chain,
                                )
                            )
                            segment_counter_for_chain += 1
                            current_segment_residues = [res]
                        else:
                            current_segment_residues.append(res)
                    else:
                        all_segments.append(
                            Segment(
                                current_segment_residues,
                                chain.id,
                                segment_counter_for_chain,
                            )
                        )
                        segment_counter_for_chain += 1
                        current_segment_residues = [res]

            if current_segment_residues:
                all_segments.append(
                    Segment(
                        current_segment_residues, chain.id, segment_counter_for_chain
                    )
                )

        if not all_segments:
            print("No valid segments found in the PDB file.")
            return

        print(f"Identified {len(all_segments)} segments in total.")
        self.segments = all_segments

    def getConnections(self, segments_to_connect):
        # Find all potential connections between segments

        edges = []
        # Track which segments can form connections
        for i, seg1 in enumerate(segments_to_connect):
            bestDist = np.inf
            bestSeg2 = None
            second_bestDist = np.inf
            second_bestSeg2 = None
            for j in range(0, len(segments_to_connect)):
                # Calculate distance between C-term of seg1 and N-term of seg2
                seg2 = segments_to_connect[j]
                c_n_distance = np.linalg.norm(
                    seg1.c_atom.get_coord() - seg2.n_atom.get_coord()
                )

                # If distance is within threshold, add to potential connections
                if (
                    c_n_distance <= self.CONNECTION_THRESHOLD
                    and c_n_distance < bestDist
                ):
                    bestDist = c_n_distance
                    second_bestSeg2 = bestSeg2
                    second_bestDist = bestDist
                    bestSeg2 = j

            if bestSeg2 is not None:
                edges.append((i, bestSeg2, {"dist": bestDist}))
            if second_bestSeg2 is not None:
                edges.append((i, second_bestSeg2, {"dist": second_bestDist}))

        return edges

    def find_optimal_connections(self, edges, segments):
        """Find the permutation that maximizes the number of connected segments"""
        # Create a graph from the edges
        remaining_segments = set(range(len(segments)))
        longestPaths = []
        while len(remaining_segments) > 1:
            # Create a directed graph from the edges
            G = nx.DiGraph()
            G.add_edges_from(edges)

            # Find the largest connected component
            longestPath = nx.dag_longest_path(G, weight="dist")
            if len(longestPath) == 0:
                break
            longestPaths.append(longestPath)

            # remove the segments in longest path from remaining segments, and edges
            for i in longestPath:
                remaining_segments.remove(i)
                edges = [edge for edge in edges if edge[0] != i and edge[1] != i]

        # add the remaining segments as singletons
        for i in remaining_segments:
            longestPaths.append([i])
        return longestPaths

    def connectChains(self, chain_segments, segments_to_connect, chain_id):
        """Connect segments in a chain using PolyGLinker."""
        # Create a new chain for the connected segments
        new_chain = Chain.Chain(chain_id)
        resId = 0
        for i in range(len(chain_segments) - 1):
            residues = segments_to_connect[chain_segments[i]].residues

            for res in residues:
                newRes = res.copy()
                newRes.id = (" ", resId, " ")  # Reset residue ID
                new_chain.add(newRes)
                resId += 1

            linker = PolyGLinker()
            linkResidues = linker.generate_linker(
                segments_to_connect[chain_segments[i]],
                segments_to_connect[chain_segments[i + 1]],
            )
            for res in linkResidues:
                newRes = res.copy()
                newRes.id = (" ", resId, " ")  # Reset residue ID
                new_chain.add(newRes)
                resId += 1

        # Add the last segment to the new chain
        last_segment = segments_to_connect[chain_segments[-1]].residues
        for res in last_segment:
            newRes = res.copy()
            newRes.id = (" ", resId, " ")  # Reset residue ID
            new_chain.add(newRes)
            resId += 1
        return new_chain

    def ConnectClosest(self):
        """
        Connect segments with the closest ends, respecting protein chain directionality.
        Uses PolyGLinker to generate connecting residues between segments.



        Returns:
            A new PDB structure with connected segments
        """
        # Create a copy of the segments to avoid modifying the originals
        segments_to_connect = self.segments.copy()
        if not segments_to_connect:
            print("No segments to connect.")
            return None

        # Create a new structure for the connected protein
        connected_structure = Structure.Structure("connected_protein")
        connected_model = Model.Model(0)
        connected_structure.add(connected_model)

        potential_connections = self.getConnections(segments_to_connect)

        # Get the optimal connections
        chains = self.find_optimal_connections(
            potential_connections, segments_to_connect
        )

        print(f"Found {len(chains)} chains")

        # Create a linker for generating connecting residues

        # Organize segments by chain IDs and add them to the new structure
        chain_id_counter = 65  # Start with 'A'

        for chain_segments in chains:
            chain_id = chr(chain_id_counter)
            new_Chain = self.connectChains(
                chain_segments, segments_to_connect, chain_id
            )
            connected_model.add(new_Chain)
            chain_id_counter += 1

        return connected_model


def PDBLoader(input_pdb_path):
    """
    Load a PDB file and return the first model.

    Args:
        filename (str): Path to the PDB file.

    Returns:
        Bio.PDB.Model.Model: The loaded model.
    """
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", input_pdb_path)
    except FileNotFoundError:
        print(f"Error: Input PDB file not found at {input_pdb_path}")
    except Exception as e:
        print(f"Error parsing PDB file {input_pdb_path}: {e}")

    # Extract the first model from the structure
    return structure[0]


def PDBSaver(connected_model, output_pdb_path):
    """
    Save a PDB model to a file.

    Args:
        model (Bio.PDB.Model.Model): The model to save.
        filename (str): Path to the output PDB file.
    """
    io = PDBIO()
    io.set_structure(connected_model)
    io.save(output_pdb_path)
    print(f"Connected structure saved to {output_pdb_path}")


def PDBViewer(oldmodel):
    """
    Displays each segment as a ribbon with a different color using py3Dmol.
    Creates temporary chains from segments and renders the structure.
    """
    # Create a new structure for visualization
    struct = Structure.Structure("segmented")
    model = Model.Model(0)
    struct.add(model)

    # segments are the chains if the oldmodel
    colors = plt.cm.rainbow(np.linspace(0, 1, len(oldmodel)))
    colors = [
        f"rgb({int(r*255)},{int(g*255)},{int(b*255)})"
        for r, g, b, _ in colors
    ]
    # Create chains from segments
    for i, segment in enumerate(oldmodel):
        chain_id = chr(65 + (i % 26))  # A-Z (looping if needed)
        chain = Chain.Chain(chain_id)

        for residue in segment:
            # Create a deep copy of the residue
            new_res = Residue.Residue(residue.id, residue.resname, residue.segid)
            for atom in residue:
                new_atom = Atom.Atom(
                    atom.name,
                    atom.coord,
                    atom.bfactor,
                    atom.occupancy,
                    atom.altloc,
                    atom.fullname,
                    atom.serial_number,
                    atom.element,
                )
                new_res.add(new_atom)
            chain.add(new_res)
        model.add(chain)

    # Save structure to temporary PDB file
    io = PDBIO()
    io.set_structure(struct)
    temp_pdb = "temp_segmented.pdb"
    io.save(temp_pdb)

    # Visualize with py3Dmol
    view = py3Dmol.view(width=800, height=600)
    view.addModel(open(temp_pdb, "r").read(), "pdb")

    # Color chains based on segment
    for i, segment in enumerate(model):
        chain_id = chr(65 + (i % 26))
        view.setStyle({"chain": chain_id}, {"cartoon": {"color": colors[i]}})

    view.zoomTo()
    view.setBackgroundColor("white")
    return view.show()


def MoveChainSection(pdb_file, chain_id, residue_range, direction, distance):
    """
    Move a section of a chain in a PDB file and adjust surrounding residues to maintain CA-CA distances.

    Args:
        pdb_file (str): Path to the input PDB file.
        chain_id (str): Chain ID of the chain to modify.
        residue_range (tuple): Range of residues to move (start_residue, end_residue).
        direction (np.array): Direction vector for the movement.
        distance (float): Distance to move the section.

    Returns:
        Structure.Structure: Modified structure with the moved chain section.
    """
    # Load the PDB structure
    structure = PDBLoader(pdb_file)
    chain = structure[0][chain_id]

    # Normalize the direction vector
    direction = direction / np.linalg.norm(direction)

    # Identify residues to move
    start_res, end_res = residue_range
    residues_to_move = [
        res for res in chain.get_residues() if start_res <= res.id[1] <= end_res
    ]

    # Move the specified residues
    for residue in residues_to_move:
        for atom in residue:
            atom.coord += direction * distance

    # Adjust surrounding residues to maintain CA-CA distances
    ca_ca_distance = 3.8  # Standard CA-CA distance in Å
    all_residues = list(chain.get_residues())
    start_idx = all_residues.index(residues_to_move[0])
    end_idx = all_residues.index(residues_to_move[-1])

    # Pull residues before the range
    for i in range(start_idx - 1, -1, -1):
        current_res = all_residues[i]
        next_res = all_residues[i + 1]
        if "CA" in current_res and "CA" in next_res:
            ca_current = current_res["CA"].get_coord()
            ca_next = next_res["CA"].get_coord()
            displacement = ca_next - ca_current
            correction = (np.linalg.norm(displacement) - ca_ca_distance) * (
                displacement / np.linalg.norm(displacement)
            )
            for atom in current_res:
                atom.coord += correction

    # Pull residues after the range
    for i in range(end_idx + 1, len(all_residues)):
        current_res = all_residues[i]
        prev_res = all_residues[i - 1]
        if "CA" in current_res and "CA" in prev_res:
            ca_current = current_res["CA"].get_coord()
            ca_prev = prev_res["CA"].get_coord()
            displacement = ca_prev - ca_current
            correction = (np.linalg.norm(displacement) - ca_ca_distance) * (
                displacement / np.linalg.norm(displacement)
            )
            for atom in current_res:
                atom.coord += correction

    return structure

def CreateSequence(pdbFile, sequenceFile,addSaltBridges=False, numberAttempts=1, batch_size=1, number_batches=1):

    os.makedirs("./outputs/default/backbones", exist_ok=True)
    os.makedirs("./outputs/default/seqs", exist_ok=True)
    os.makedirs("./inputs", exist_ok=True)

    ClearFolder('/content/outputs/default/backbones')
    ClearFolder('/content/outputs/default/seqs')


    
    print(f"Loading {pdbFile}")
    cc=0  
    for seed in range(numberAttempts):
        
        attemptFile = f'./inputs/{ os.path.splitext(os.path.basename(pdbFile))[0]}_{cc}.pdb'
        shutil.copy(pdbFile, attemptFile)
        temp = np.random.randint(100)/100.0 * .1
        print(f"Running seed {seed+1} with temp {temp}")
        sseed=np.random.randint(100000)+seed

        if addSaltBridges:
            change =  np.random.randint(40)/10.0
            bias = f"A:{change/1.5:.2f},T:{change:.2f},R:{change:.2f},E:{change:.2f},Y:{change/2:.2f},G:{-1*change:.2f},S:{-1*change:.2f},K:{-1*change:.2f},D:{-1*change:.2f}"
            command = [
                "python", "./LigandMPNN/run.py",
                "--verbose", "0",
                "--model_type", "ligand_mpnn",
                "--checkpoint_ligand_mpnn", "./LigandMPNN/model_params/ligandmpnn_v_32_030_25.pt",
                "--temperature", str(temp),
                "--seed", str(sseed + 1),
                "--pdb_path", attemptFile,
                "--out_folder", "./outputs/default",
                "--batch_size", str(batch_size),
                "--number_of_batches", str(number_batches),
                "--bias_AA", bias
            ]
            subprocess.run(command, check=True)
        else:
            command = [
                "python", "./LigandMPNN/run.py",
                "--verbose", "0",
                "--model_type", "ligand_mpnn",
                "--checkpoint_ligand_mpnn", "./LigandMPNN/model_params/ligandmpnn_v_32_030_25.pt",
                "--temperature", str(temp),
                "--seed", str(sseed + 1),
                "--pdb_path", attemptFile,
                "--out_folder", "./outputs/default",
                "--batch_size", str(batch_size),
                "--number_of_batches", str(number_batches)
            ]
            subprocess.run(command, check=True)

        clear_output()
        file =(f'/content/outputs/default/seqs/{os.path.splitext(os.path.basename(attemptFile))[0]}.fa')
        generated_sequences = []
        
        sequences = load_fasta(file)
        generated_sequences.extend(sequences)
        # Sort sequences by 'overall_confidence' in json_data and select the top value
        sorted_sequences = sorted(generated_sequences, key=lambda x: x[0]['overall_confidence'], reverse=True)
        # Print out the top 15 sequences with their overall confidence
        for i, (json_data, sequence) in enumerate(sorted_sequences[:7]):
            print(f"Sequence {i+1}:")
            print(f"Overall Confidence: {json_data['overall_confidence']} {json_data['ligand_confidence']} {json_data['seq_rec']}")
            print(f"Source File: {json_data['file']}")
            parts= sequence.split(':')
            for part in parts:
                print(part) 
        print()
        print('Sorted by ligand confidence')
        print()
        sorted_sequences = sorted(generated_sequences, key=lambda x: x[0]['ligand_confidence'], reverse=True)
        for i, (json_data, sequence) in enumerate(sorted_sequences[:7]):
            print(f"Sequence {i+1}  Overall Confidence: {json_data['overall_confidence']} {json_data['ligand_confidence']} {json_data['seq_rec']}")
            parts= sequence.split(':')
            for part in parts:
                print(part)    

        print(f"outputFile: {sequenceFile}")
        with open(file, "r") as infile, open(sequenceFile, "a") as outfile:
            outfile.write(infile.read())
        print (f"Finished seed {cc} : {seed+1} with temp {temp}")
        
        cc+=1    
    return sorted_sequences



def RotateChainSection(model, chain_id, residue_range, rotation_angles, axis=None):
    """
    Rotate a section of a chain around its center of mass and adjust surrounding residues.
    
    Args:
        model: First model of the PDB file
        chain_id (str): Chain ID of the chain to modify
        residue_range (tuple): Range of residues to rotate (start_residue, end_residue)
        rotation_angles (tuple or float): Rotation angles in degrees (x, y, z) or single angle if axis is provided
        axis (np.array, optional): Custom rotation axis. If None, uses x, y, z rotation sequence.
    
    Returns:
        Structure.Structure: Modified structure with the rotated chain section
    """
     
    # Load the chain
    chain = model[chain_id]
    
    # Identify residues to rotate
    if residue_range==None:
        # If no residue range is provided, move the entire chain
        residues_to_rotate = list(chain.get_residues())
        start_res = residues_to_rotate[0].id[1]
        end_res = residues_to_rotate[-1].id[1]
    else:
        start_res, end_res = residue_range
        residues_to_rotate = [
            res for res in chain.get_residues() if start_res <= res.id[1] <= end_res
        ]
    
    # Calculate center of mass of the selected section
    com = np.zeros(3)
    atom_count = 0
    for res in residues_to_rotate:
        for atom in res:
            com += atom.coord
            atom_count += 1
    
    if atom_count > 0:
        com = com / atom_count
    
    # Create rotation matrix
    if axis is not None:
        # Single axis rotation
        axis = np.array(axis) / np.linalg.norm(axis)
        angle_rad = np.radians(rotation_angles)
        r = Rotation.from_rotvec(angle_rad * axis)
        rotation_matrix = r.as_matrix()
    else:
        # Sequential x, y, z rotations
        angles = rotation_angles if isinstance(rotation_angles, tuple) else (rotation_angles, 0, 0)
        r = Rotation.from_euler('xyz', angles, degrees=True)
        rotation_matrix = r.as_matrix()
    
    # Apply rotation to selected residues
    for res in residues_to_rotate:
        for atom in res:
            # Move to origin, rotate, then move back
            atom.coord = atom.coord - com
            atom.coord = np.dot(rotation_matrix, atom.coord)
            atom.coord = atom.coord + com
    
    # Adjust surrounding residues to maintain CA-CA distances
    ca_ca_distance = 3.8  # Standard CA-CA distance in Å
    all_residues = list(chain.get_residues())
    start_idx = all_residues.index(residues_to_rotate[0])
    end_idx = all_residues.index(residues_to_rotate[-1])
    
    # Pull residues before the range
    for i in range(start_idx - 1, -1, -1):
        current_res = all_residues[i]
        next_res = all_residues[i + 1]
        if "CA" in current_res and "CA" in next_res:
            ca_current = current_res["CA"].get_coord()
            ca_next = next_res["CA"].get_coord()
            displacement = ca_next - ca_current
            correction = (np.linalg.norm(displacement) - ca_ca_distance) * (
                displacement / np.linalg.norm(displacement)
            )
            
            # If the distance is close enough, break the loop
            if (np.linalg.norm(correction)/ca_ca_distance) < 0.01:
                break
            
            for atom in current_res:
                atom.coord += correction
    
    # Pull residues after the range
    for i in range(end_idx + 1, len(all_residues)):
        current_res = all_residues[i]
        prev_res = all_residues[i - 1]
        if "CA" in current_res and "CA" in prev_res:
            ca_current = current_res["CA"].get_coord()
            ca_prev = prev_res["CA"].get_coord()
            displacement = ca_prev - ca_current
            correction = (np.linalg.norm(displacement) - ca_ca_distance) * (
                displacement / np.linalg.norm(displacement)
            )
            
            # If the distance is close enough, break the loop
            if (np.linalg.norm(correction)/ca_ca_distance) < 0.01:
                break
                
            for atom in current_res:
                atom.coord += correction
    
    return model

def MoveChainSection(model, chain_id, residue_range, direction, distance):
    """
    Move a section of a chain in a PDB file and adjust surrounding residues to maintain CA-CA distances.

    Args:
        model : first model of the PDB file.
        chain_id (str): Chain ID of the chain to modify.
        residue_range (tuple): Range of residues to move (start_residue, end_residue).
        direction (np.array): Direction vector for the movement.
        distance (float): Distance to move the section.

    Returns:
        Structure.Structure: Modified structure with the moved chain section.
    """
    # Load the chain
    chain = model[chain_id]

    # Normalize the direction vector
    direction = direction / np.linalg.norm(direction)

    # Identify residues to move
    
    if residue_range==None:
        # If no residue range is provided, move the entire chain
        residues_to_move = list(chain.get_residues())
        start_res = residues_to_move[0].id[1]
        end_res = residues_to_move[-1].id[1]
    else:
        start_res, end_res = residue_range
        residues_to_move = [
            res for res in chain.get_residues() if start_res <= res.id[1] <= end_res
        ]

    # Move the specified residues
    for residue in residues_to_move:
        for atom in residue:
            atom.coord += direction * distance

    # Adjust surrounding residues to maintain CA-CA distances
    ca_ca_distance = 3.8  # Standard CA-CA distance in Å
    all_residues = list(chain.get_residues())
    start_idx = all_residues.index(residues_to_move[0])
    end_idx = all_residues.index(residues_to_move[-1])

    # Pull residues before the range
    for i in range(start_idx - 1, -1, -1):
        current_res = all_residues[i]
        next_res = all_residues[i + 1]
        if "CA" in current_res and "CA" in next_res:
            ca_current = current_res["CA"].get_coord()
            ca_next = next_res["CA"].get_coord()
            displacement = ca_next - ca_current
            correction = (np.linalg.norm(displacement) - ca_ca_distance) * (
                displacement / np.linalg.norm(displacement)
            )
            
            #if the distance is small enough, break the loop
            if ( np.linalg.norm( correction)/ca_ca_distance)< .01:
                break
            
            for atom in current_res:
                atom.coord += correction

    # Pull residues after the range
    for i in range(end_idx + 1, len(all_residues)):
        current_res = all_residues[i]
        prev_res = all_residues[i - 1]
        if "CA" in current_res and "CA" in prev_res:
            ca_current = current_res["CA"].get_coord()
            ca_prev = prev_res["CA"].get_coord()
            displacement = ca_prev - ca_current
            correction = (np.linalg.norm(displacement) - ca_ca_distance) * (
                displacement / np.linalg.norm(displacement)
            )
            #if the distance is small enough, break the loop
            if ( np.linalg.norm( correction)/ca_ca_distance)< .01:
                break
            for atom in current_res:
                atom.coord += correction

    return model