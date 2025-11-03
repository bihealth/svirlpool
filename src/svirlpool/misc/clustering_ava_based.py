# %%
import subprocess
import tempfile
from collections import Counter
from copy import deepcopy
from pathlib import Path
from shlex import split

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igraph import Graph
from logzero import logger as log
from stopit import threading_timeoutable as timeoutable

from ..scripts import consensus, datatypes, util


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, start_end = region.split(":")
    start, end = start_end.split("-")
    start = start.replace(",", "")
    end = end.replace(",", "")
    return (chrom, int(start), int(end))


# %%
def plot_graph(G: Graph) -> go.Figure:

    # # from the input gaph, construct a new one with filtered verteces and edges
    # new_G = Graph()
    # # add all verteces that have more than one neighbor
    # new_G.add_vertices([v for v in G.vs if len(v.neighbors()) > 1])
    # # add all edges if they connect two verteces that are in the new graph
    # new_G.add_edges([e for e in G.es if e.source in new_G.vs and e.target in new_G.vs])

    # G = new_G

    # Visualization parameters
    layout = G.layout("fr")  # Fruchterman-Reingold layout for centrality
    # edge_weights = G.es["weight"]

    # Normalize weights for coloring
    norm = mcolors.Normalize(vmin=0, vmax=1)
    cmap = plt.cm.viridis  # Choose a colormap
    # edge_colors = [mcolors.to_hex(cmap(norm(weight))) for weight in edge_weights]

    # Extract positions and prepare edge coordinates
    pos = {i: layout[i] for i in range(len(G.vs))}
    # edge_x = []
    # edge_y = []
    # weights = dict()

    fig = go.Figure()
    for edge in G.es:
        x0, y0 = pos[edge.source]
        x1, y1 = pos[edge.target]
        # edge_x.extend([x0, x1, None])
        # edge_y.extend([y0, y1, None])
        # weights[(edge.source, edge.target)] = (math.log(edge["weight"]) if edge["weight"] > 0 else 0) * 0.1
        fig.add_trace(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                line=dict(
                    width=edge["weight"], color=f'rgb({int(255.0*edge["weight"])},0,0)'
                ),
                hoverinfo="none",
                mode="lines",
            )
        )

    # Add nodes to plotly figure
    node_x = [pos[i][0] for i in range(len(G.vs))]
    node_y = [pos[i][1] for i in range(len(G.vs))]
    node_text = G.vs["name"]

    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            marker=dict(
                size=10,
                color="blue",
            ),
            text=node_text,
            textposition="top center",
            hoverinfo="text",
        )
    )
    # write the plot to a file
    return fig


def normalize_edge_weights(graph: Graph) -> Graph:
    ng = deepcopy(graph)
    max_weight = max(ng.es["weight"])
    for edge in ng.es:
        edge["weight"] = edge["weight"] / max_weight
    return ng


def plot_proto_consensuses_alignments_to_terminal(
    reconstructibles: list[datatypes.ReconstructibleSequence], terminal_width: int
) -> None:
    alignments = [rs.alignment.to_pysam() for rs in reconstructibles]
    util.display_ascii_alignments(
        alignments=alignments,
        terminal_width=terminal_width,
        min_signal_size=8,
        min_bnd_size=200,
    )


# function to make ava alignments. return true if successful
def ava_align_sequences(sequences: dict[str, SeqRecord], output: Path) -> bool:
    if len(sequences) == 0:
        raise ValueError("No sequences to align")
    # check if all sequences are of type SeqRecord
    if not all(isinstance(seq, SeqRecord) for seq in sequences.values()):
        raise ValueError("All sequences must be of type SeqRecord")
    # write SeqRecords to tmp file
    tmp_sequences = tempfile.NamedTemporaryFile(
        delete=True, suffix=".fasta", prefix="ava_align_sequences."
    )
    with open(tmp_sequences.name, "w") as f:
        for seq in sequences.values():
            SeqIO.write(seq, f, "fasta")
    util.align_reads_with_minimap(
        bamout=output,
        reads=tmp_sequences.name,
        reference=tmp_sequences.name,
        threads=1,
        tech="ava-ont",
        aln_args=" -U25,50 --sam-hit-only -D",
    )
    return True


# function to construct an overlap graph
def overlap_graph(
    alignments: list[pysam.AlignedSegment], sequences: dict[str, SeqRecord]
) -> Graph:
    dict_readname_id = {
        read_name: idx for idx, read_name in enumerate(sorted(sequences.keys()))
    }
    g = Graph()
    g.add_vertices(sorted(sequences.keys()))
    for aln in alignments:
        read_idx = dict_readname_id[aln.query_name]
        ref_idx = dict_readname_id[aln.reference_name]
        g.add_edge(
            read_idx,
            ref_idx,
            weight=abs(aln.query_alignment_end - aln.query_alignment_start),
        )
    return g


# function to construct a conflict graph
def conflict_graph(
    alignments: list[pysam.AlignedSegment],
    sequences: dict[str, SeqRecord],
    min_conflict_size: int,
) -> Graph:
    dict_readname_id = {
        read_name: idx for idx, read_name in enumerate(sorted(sequences.keys()))
    }
    g = Graph()
    g.add_vertices(sorted(sequences.keys()))
    for aln in alignments:
        read_idx = dict_readname_id[aln.query_name]
        ref_idx = dict_readname_id[aln.reference_name]
        ras: datatypes.ReadAlignmentSignals = (
            consensus.parse_ReadAlignmentSignals_from_alignment(
                alignment=aln, sampleID=0, min_signal_size=8, min_bnd_size=200
            )
        )
        edge_weight = grave_sv_signals_sum(
            read=sequences[aln.query_name],
            ref=sequences[aln.reference_name],
            ras=ras,
            max_tolerated_bnd=200,
            max_tolerated_indel=8,
        )
        if edge_weight > min_conflict_size:
            g.add_edge(read_idx, ref_idx, weight=edge_weight)
    return g


def find_non_overlapping_cliques(graph: Graph, min_clique_size: int) -> list[list[str]]:
    g = deepcopy(graph)
    cliques = []
    while len(g.vs) > 0:
        clique = g.largest_cliques()[0]
        if len(clique) < min_clique_size:
            break
        cliques.append([g.vs[c]["name"] for c in clique])
        g.delete_vertices(
            clique
        )  # changes IDs of vertices! -> need to remap the indices to the original indices
        # maybe use the name property of each vertex
    return cliques


def make_consensus_with_lamassemble(
    sequences: dict[str, SeqRecord],
    name: str,
    threads: int,
    lamassemble_mat_path: Path,
    min_coverage: int = 1,
    format="fasta",
    tmp_dir_path: Path | None = None,
) -> SeqRecord:
    if len(sequences) == 0:
        raise ValueError("No sequences to align")
    if len(sequences) == 1:
        return next(iter(sequences.values()))
    # log.debug(f"processing {name}")
    # chosen_reads is of the form [SeqRecord]
    with tempfile.TemporaryDirectory(dir=tmp_dir_path) as tdir:
        # # first write consensus sequences to tmp fasta file
        tmp_fasta = tempfile.NamedTemporaryFile(
            dir=tdir,
            suffix=".fasta",
            prefix="sequences.",
            delete=False if tmp_dir_path else True,
        )
        with open(tmp_fasta.name, "w") as f:
            for seqname, seq in sequences.items():
                SeqIO.write(seq, f, format)
        # path_reads = path_tmp_dir / f"reads.{format}"
        # SeqIO.write(chosen_reads, path_reads, format)
        cmd_consensus = split(
            f"lamassemble \
            -n {name} \
            -s {min_coverage} \
            -P {str(threads)} \
            {str(lamassemble_mat_path)} \
            {tmp_fasta.name}"
        )
        # catch output and create the seqrecord object from it
        with subprocess.Popen(cmd_consensus, stdout=subprocess.PIPE) as proc:
            output = proc.stdout.read().decode()
            # output is a fasta content, meaning the consensus sequence name starts with '>' and after the newline comes
            # the DNA sequence (that can be separated by newlines, but not necessarily)
            try:
                dna_seq = output.split("\n")[1]
                seqrec = SeqRecord(Seq(dna_seq), id=name, name=name, description="")
                # log.debug(f"lamassemble returned a consensus sequence of length {len(dna_seq)} for {name}")
            except IndexError:
                log.warning(
                    f"lamassemble did not return a consensus sequence for {name}"
                )
                seqrec = SeqRecord(Seq(""), id=name, name=name, description="")
    return seqrec


def choose_consensus_candidates(
    reads: dict[str, SeqRecord],
    consensus_candidates: dict[str, SeqRecord],
    threads: int,
    max_allowed_signal_size: int = 12,
    max_allowed_clipped_end: int = 200,
    sampleID: int | None = None,
    tmp_dir_path: Path | None = None,
    minimap_logfile: Path | None = None,
) -> dict[str, list[datatypes.ReconstructibleSequence]]:
    """From a set of consensus candidates, this function finds the consensus sequences to which the most reads align.

    Args:
        reads (dict[str,SeqRecord]): dict of reads of the form readname:SeqRecord
        consensus_candidates (dict[str,SeqRecord]): dict of consensus candidates of the form consensusname:SeqRecord
        threads (int): number of threads to use for minimap2
        sampleID (int | None, optional): is written to the ReconstructibleSequence objects if provided. If not, int(query_name.split('.')[1]) is chosen. Defaults to None.
        tmp_dir_path (Path | None, optional): Path to the temporary directory. Defaults to None.
        minimap_logfile (Path | None, optional): Path to the minimap2 logfile. Defaults to None.

    Returns:
        dict[str,list[datatypes.ReconstructibleSequence]]: _description_
    """
    tmp_reads = tempfile.NamedTemporaryFile(
        prefix="all_cutreads.",
        suffix=".fasta",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    with open(tmp_reads.name, "w") as f:
        SeqIO.write(reads.values(), f, "fasta")
    # write consensus sequences to fasta file
    tmp_consensus = tempfile.NamedTemporaryFile(
        prefix="all_consensuses.",
        suffix=".fasta",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    with open(tmp_consensus.name, "w") as f:
        SeqIO.write(consensus_candidates.values(), f, "fasta")
    # align all reads to all consensus sequences
    tmp_alignments = tempfile.NamedTemporaryFile(
        prefix="all_alignments.",
        suffix=".bam",
        dir=tmp_dir_path,
        delete=False if tmp_dir_path else True,
    )
    util.align_reads_with_minimap(
        reference=tmp_consensus.name,
        reads=tmp_reads.name,
        bamout=tmp_alignments.name,
        threads=threads,
        logfile=minimap_logfile,
        tech="map-ont",
        aln_args=" -Y --sam-hit-only -f 0.1 -U 10,25 -r200,500 -m 1000 -F 200 -H",
    )  # -H flag is use compressed homopolymers

    # iterate over all alignments and fill a dict {consensus_name:[alignments]}
    result: dict[str, list[datatypes.ReconstructibleSequence]] = {}
    for aln in pysam.AlignmentFile(tmp_alignments.name, mode="r"):
        if sampleID == None:
            sampleID = int(aln.query_name.split(".")[1])
        # get the BNDs and indels from the alignment. If there are any significant BNDs or indels, the read is not returned as reconstructible
        ras: datatypes.ReadAlignmentSignals = (
            consensus.parse_ReadAlignmentSignals_from_alignment(
                alignment=aln,
                sampleID=sampleID,
                min_signal_size=max_allowed_signal_size,
                min_bnd_size=max_allowed_clipped_end,
            )
        )
        # all sv signals in ras.SV_signals are checked for sequence complexity. If the insertion or deletion is comprised of 4 or less 2-mers,
        # then the sv signal is not counted. If the sv is not of low complexity, continue.
        sum_signal = grave_sv_signals_sum(
            ras=ras,
            read=reads[aln.query_name],
            ref=consensus_candidates[aln.reference_name],
            max_tolerated_bnd=max_allowed_clipped_end,
            max_tolerated_indel=8,
        )
        if sum_signal > max_allowed_signal_size:
            continue
        if aln.reference_name not in result:
            result[aln.reference_name] = []
        result[aln.reference_name].append(
            consensus.parse_ReconstructibleSequence(
                aln=aln,
                sampleID=sampleID,
                description=reads[aln.query_name].description,
            )
        )
    # filter result for any alignments that share the same read name
    result_filtered = {}
    for consensus_name, reconstructibles in result.items():
        result_filtered[consensus_name] = []
        counter = dict(Counter([rs.alignment.readname for rs in reconstructibles]))
        for rs in reconstructibles:
            if counter[rs.alignment.readname] == 1:
                result_filtered[consensus_name].append(rs)
    return result_filtered


def grave_sv_signals_sum(
    read: SeqRecord,
    ref: SeqRecord,
    ras: datatypes.ReadAlignmentSignals,
    max_tolerated_indel: int = 12,
    max_tolerated_bnd: int = 200,
    min_2mers: int = 4,
) -> int:
    letter_dict = util.compute_letter_dict("ACGTNXU")
    sum_signal: int = 0
    clipped_ends = [signal for signal in ras.SV_signals if signal.sv_type >= 3]
    if len(clipped_ends) > 0:
        if len(clipped_ends) > 1:
            # log.debug(f"read {read.id} has more than one clipped end")
            return len(read.seq)
        # any clipped end is only allowed within the first 100 bp of the reference or the last 100 bp
        if (
            clipped_ends[0].sv_type == 3
            and clipped_ends[0].ref_start > max_tolerated_bnd
        ):
            # log.debug(f"read {read.id} has a clipped start at {clipped_ends[0].ref_start}")
            return len(read.seq)
        if (
            clipped_ends[0].sv_type == 4
            and clipped_ends[0].ref_end < len(ref) - max_tolerated_bnd
        ):
            # log.debug(f"read {read.id} has a clipped end at {clipped_ends[0].ref_end}")
            return len(read.seq)
    for sv in ras.SV_signals:
        if sv.size > max_tolerated_indel:
            if sv.sv_type == 0:
                # insertion.
                # compute the sequence complexity of the insertion
                counter = util.kmer_counter_from_string(
                    string=read.seq[sv.read_start : sv.read_end],
                    letter_dict=letter_dict,
                    k=2,
                )
                if len(counter) >= min_2mers:
                    # log.debug(f"read {read.id} has a complex insertion: {sv.read_start}-{sv.read_end}")
                    sum_signal += sv.size
                # else:
                #     log.debug(f"read {read.id} has a simple insertion: {sv.read_start}-{sv.read_end}")
            if sv.sv_type == 1:
                # deletion
                # compute the sequence complexity of the deletion
                counter = util.kmer_counter_from_string(
                    string=ref.seq[sv.ref_start : sv.ref_end],
                    letter_dict=letter_dict,
                    k=2,
                )
                if len(counter) > 4:
                    # log.debug(f"read {read.id} has a complex deletion: {sv.ref_start}-{sv.ref_end}")
                    sum_signal += sv.size
                # else:
                # log.debug(f"read {read.id} has a simple deletion: {sv.ref_start}-{sv.ref_end}")
    return min(sum_signal, len(read.seq))


def lamassemble_proto_consensus_sequences(
    sequences: dict[str, SeqRecord],
    lamassemble_mat: Path,
    threads: int,
    min_conflict_size: int = 12,
    min_clique_size: int = 3,
    min_coverage: int = 1,
    sequences_per_assembly: int = 5,
    tmp_dir_path: Path | None = None,
) -> dict[str, SeqRecord]:

    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tdir:
        tmp_alignments = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            suffix=".bam",
            prefix="ava_alignments.",
        )
        ava_align_sequences(sequences=sequences, output=Path(tmp_alignments.name))
        ava_alignments = list(pysam.AlignmentFile(tmp_alignments.name))
        # log.debug(f"found {len(ava_alignments)} ava alignments")
        tmp_alignments.close()
        overlap_g = overlap_graph(alignments=ava_alignments, sequences=sequences)
        # plot_graph(normalize_edge_weights(overlap_g)).show()
        conflict_g = conflict_graph(
            alignments=ava_alignments,
            sequences=sequences,
            min_conflict_size=min_conflict_size,
        )
        # plot_graph(normalize_edge_weights(conflict_g)).show()
        consensus_g = overlap_g - conflict_g
        # plot_graph(normalize_edge_weights(consensus_g)).show()

        max_cliques = find_non_overlapping_cliques(
            graph=consensus_g, min_clique_size=min_clique_size
        )
        # log.debug(f"found cliques: {max_cliques}")
        dict_vertex_name_edge_sums = {
            n: v
            for n, v in zip(
                consensus_g.vs["name"],
                consensus_g.strength(weights=consensus_g.es["weight"]),
            )
        }
        sorted_clique_read_scores = [
            sorted(
                [(n, dict_vertex_name_edge_sums[n]) for n in clique],
                key=lambda x: x[1],
                reverse=True,
            )
            for clique in max_cliques
        ]

        consensus_candidates = dict()
        for i, clique in enumerate(max_cliques):
            reads = {
                name: sequences[name]
                for name, _ in sorted_clique_read_scores[i][:sequences_per_assembly]
            }
            consensus_candidates[str(i)] = make_consensus_with_lamassemble(
                sequences=reads,
                name=str(i),
                threads=threads,
                lamassemble_mat_path=lamassemble_mat,
                min_coverage=min_coverage,
                tmp_dir_path=tmp_dir_path,
            )
    return consensus_candidates


def get_unused_sequences(
    sequences: dict[str, SeqRecord],
    proto_consensuses: dict[str, list[datatypes.ReconstructibleSequence]],
) -> dict[str, SeqRecord]:
    # return all sequences that are not present in any proto consensus
    used_sequences = set()
    for reconstructible_sequences in proto_consensuses.values():
        for rs in reconstructible_sequences:
            used_sequences.add(rs.alignment.readname)
    return {name: seq for name, seq in sequences.items() if name not in used_sequences}


def polish_consensus(
    consensus_sequence: SeqRecord,
    reconstructibles: list[datatypes.ReconstructibleSequence],
    threads: int = 1,
    tmp_dir_path: Path | None = None,
) -> tuple[SeqRecord, list[datatypes.ReconstructibleSequence]] | None:
    reconstructed_reads = [
        consensus.reconstruct_ReconstructibleSequence(
            reconstructibleSequence=rc, reference=consensus_sequence.seq
        )
        for rc in reconstructibles
    ]
    read_alignments = [rc.alignment.to_pysam() for rc in reconstructibles]
    with tempfile.TemporaryDirectory(
        dir=tmp_dir_path, delete=False if tmp_dir_path else True
    ) as tdir:
        tmp_reads = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            suffix=".fasta",
            prefix="reads.",
        )
        # write reads to tmp file
        with open(tmp_reads.name, "w") as f:
            SeqIO.write(
                [
                    SeqRecord(seq=rcr.sequence, id=rcr.name, name=rcr.name)
                    for rcr in reconstructed_reads
                ],
                f,
                "fasta",
            )
        # write samf ile from alignments
        tmp_sam = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            suffix=".sam",
            prefix="alignments.",
        )
        with pysam.Samfile(
            header=read_alignments[0].header, mode="w", filename=tmp_sam.name
        ) as sf:
            for aln in read_alignments:
                sf.write(aln)
        tmp_consensus = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            suffix=".fasta",
            prefix="consensus.",
        )
        with open(tmp_consensus.name, "w") as f:
            SeqIO.write(consensus_sequence, f, "fasta")
        tmp_output = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            suffix=".fasta",
            prefix="polished.",
        )
        polished = consensus.make_consensus_with_racon(
            consensus_fasta_path=tmp_output.name,
            name=consensus_sequence.id,
            reads_fasta=tmp_reads.name,
            reference_fasta=tmp_consensus.name,
            sam_alignments=tmp_sam.name,
            threads=threads,
        )
        # debug: write all sequence names in tmp_output
        # debug_consensus = SeqIO.parse(tmp_output.name, "fasta")
        if not polished:
            return None
        # re-align the reads to the polished consensus
        # and parse the alignments to ReconstructibleSequence objects
        tmp_new_alignments = tempfile.NamedTemporaryFile(
            dir=tdir,
            delete=False if tmp_dir_path else True,
            suffix=".bam",
            prefix="new_alignments.",
        )
        descriptions_dict = {
            rc.alignment.readname: rc.description for rc in reconstructibles
        }
        # print all sequence names in tmp_reads. read the content of the file
        # debug_reads = SeqIO.parse(tmp_reads.name, "fasta")
        # for seq in debug_reads:
        #     log.debug(f"read: {seq.id}")
        util.align_reads_with_minimap(
            bamout=tmp_new_alignments.name,
            reads=tmp_reads.name,
            reference=tmp_output.name,
            aln_args=" -Y --sam-hit-only --secondary=no -f 0.1 -U 10,25 -H",
            tech="map-ont",
            threads=threads,
            logfile=None,
        )
        new_read_alignments = [
            aln for aln in pysam.AlignmentFile(tmp_new_alignments.name)
        ]
        # log.debug(f"found {len(new_read_alignments)} new read alignments")
        # for aln in new_read_alignments:
        #     log.debug(f"new alignment: {aln.query_name} reference: {aln.reference_name}, start: {aln.reference_start}, end: {aln.reference_end}, read_start: {aln.query_alignment_start}, read_end: {aln.query_alignment_end}, secondary={aln.is_secondary}, supplementary={aln.is_supplementary}")
        new_reconstructibles = [
            consensus.parse_ReconstructibleSequence(
                aln=aln, sampleID=0, description=descriptions_dict[aln.query_name]
            )
            for aln in new_read_alignments
        ]
    return polished, new_reconstructibles


def remove_redundant_clusters(
    proto_consensuses: dict[str, list[datatypes.ReconstructibleSequence]],
    chosen_candidates: dict[str, SeqRecord],
) -> list[tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]]:
    # create a set of read names of each proto consensus
    sets = {
        name: set([rs.alignment.readname for rs in proto_consensuses[name]])
        for name in proto_consensuses.keys()
    }
    # for each pair of sets, compute the overlap in percent and print it
    sets_to_delete = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            if len(sets[str(i)] - sets[str(j)]) <= 1:
                sets_to_delete.append(i)
            elif len(sets[str(j)] - sets[str(i)]) <= 1:
                sets_to_delete.append(j)
    # create tuples (name,consensus_sequence,reconstructible_reads) for each proto consensus that is not in sets_to_delete
    results: list[tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]] = []
    for name in proto_consensuses.keys():
        if name not in sets_to_delete:
            results.append((name, chosen_candidates[name], proto_consensuses[name]))
    return results


def redistribute_reads(
    non_redunant_consensuses: list[
        tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]
    ],
    proto_consensuses: dict[str, list[datatypes.ReconstructibleSequence]],
) -> list[tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]]:
    # reads need to be redistributed so that they are moved to the consensus sequences that already have most of the reads
    # if a read (reconstructible) is present in two different proto consensuses, it is removed from all proto consensuses but not from the one with most reads
    # if a read is present in only one proto consensus, it is left there
    counter_reads_per_consensus = {
        name: len(set([rc.alignment.readname for rc in reads]))
        for name, consensus, reads in non_redunant_consensuses
    }
    dict_readname_consensus = dict()
    for name, consensus, reads in non_redunant_consensuses:
        for rc in reads:
            if rc.alignment.readname not in dict_readname_consensus:
                dict_readname_consensus[rc.alignment.readname] = []
            # check if the alignment has sv signals. Only if there are no significant sv signals, the read is added to the dict_readname_consensus
            dict_readname_consensus[rc.alignment.readname].append(name)

    # iterate dict_readname_consensus. For each key,value pair, find the value with the highest count from counter_reads_per_consensus
    # remove all other values from the dict_readname_consensus
    for readname, consensuses in dict_readname_consensus.items():
        if len(consensuses) > 1:
            max_consensus = max(
                consensuses, key=lambda x: counter_reads_per_consensus[x]
            )
            dict_readname_consensus[readname] = [max_consensus]
    # for each consensusID, create a set of readnames that are in the dict_readname_consensus
    sets = {
        name: set(
            [
                readname
                for readname, consensuses in dict_readname_consensus.items()
                if consensuses[0] == name
            ]
        )
        for name in proto_consensuses.keys()
    }
    # filter the reconstructibles in the results according to sets
    for name, consensus, reads in non_redunant_consensuses:
        non_redunant_consensuses[
            non_redunant_consensuses.index((name, consensus, reads))
        ] = (
            name,
            consensus,
            [rc for rc in reads if rc.alignment.readname in sets[name]],
        )

    # filter results again. Only keep consensus sequences that have at least 2 reads
    return [
        (name, consensus, reads)
        for name, consensus, reads in non_redunant_consensuses
        if len(reads) > 1
    ]


# optionallt polish the consensus sequences with racon,
# and re-align all reconstructed reads to the polished consensus


def polish_results(
    redistributed: list[tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]],
) -> list[tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]]:
    polished_results = []
    for name, consensus, reads in redistributed:
        res = polish_consensus(
            consensus_sequence=consensus,
            reconstructibles=reads,
            threads=1,
            tmp_dir_path=None,
        )
        if res != None:
            polished_results.append((name, *res))
        else:
            log.error(f"polishing consensus {name} failed")
            raise ValueError(
                f"polishing consensus {name} failed. consensus sequence: {str(consensus.seq)} reconstructible reads: {reads}"
            )
    return polished_results


def unused_reads_after_polish(
    polished_results: list[
        tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]
    ],
    all_reads: dict[str, SeqRecord],
) -> dict[str, SeqRecord]:
    # create a set of all readnames in polished_results
    used_reads = set()
    for name, consensus, reads in polished_results:
        for rc in reads:
            used_reads.add(rc.alignment.readname)
    # return all reads that are not in used_reads
    return {name: seq for name, seq in all_reads.items() if name not in used_reads}


def consensus_object_from_polished_results(
    polished_results: list[
        tuple[str, SeqRecord, list[datatypes.ReconstructibleSequence]]
    ],
) -> dict[str, consensus_class.Consensus]:
    results = dict()
    for name, consensus, reads in polished_results:
        consensus_infos: list[datatypes.ConsensusInfo] = []
        results[name] = consensus_class.Consensus(
            ID=consensus.id,
            consensus_sequence=str(consensus.seq),
            reconstructible_reads=reads,
            crIDs=[],
            info=consensus_infos,
        )
    return results


# %%

# 1. get alignments from region
# path_alignments = Path("/data/hdd/HG002/ont-r10.32x/minimap2.ont-r10.32x.bam")
# region = parse_region("12:97,054,340-97,062,643")
# lamassemble_mat = Path("/home/vinzenz/development/LRSV-detection/tools/lamassemble/train/promethion.mat")
# sequences:dict[str, SeqRecord] = cut_reads_from_alns.cut_reads_from_alignments(alignments=path_alignments,region=region,buffer_clipped_length=200)


@timeoutable()
def consensus_while_clustering_with_lamassemble(
    sequences: dict[str, SeqRecord],
    lamassemble_mat: Path,
    threads: int = 1,
    min_conflict_size: int = 12,
    min_clique_size: int = 2,
    min_coverage: int = 1,
    sequences_per_assembly: int = 5,
    max_allowed_signal_size: int = 30,
    tmp_dir_path: Path | None = None,
) -> dict[str, consensus_class.Consensus]:
    proto_consensus_candidates = lamassemble_proto_consensus_sequences(
        sequences=sequences,
        lamassemble_mat=lamassemble_mat,
        threads=threads,
        min_conflict_size=min_conflict_size,
        min_clique_size=min_clique_size,
        min_coverage=min_coverage,
        sequences_per_assembly=sequences_per_assembly,
        tmp_dir_path=tmp_dir_path,
    )
    proto_consensuses = choose_consensus_candidates(
        reads=sequences,
        consensus_candidates=proto_consensus_candidates,
        threads=threads,
        sampleID=0,
        tmp_dir_path=tmp_dir_path,
        minimap_logfile=None,
        max_allowed_signal_size=max_allowed_signal_size,
    )
    chosen_candidates = {
        name: proto_consensus_candidates[name] for name in proto_consensuses.keys()
    }
    # unused_sequences = get_unused_sequences(sequences=sequences, proto_consensuses=proto_consensuses)
    non_redunant_consensuses = remove_redundant_clusters(
        proto_consensuses=proto_consensuses, chosen_candidates=chosen_candidates
    )
    redistributed = redistribute_reads(
        non_redunant_consensuses=non_redunant_consensuses,
        proto_consensuses=proto_consensuses,
    )
    polished_results = polish_results(redistributed)
    # unused_sequences = unused_reads_after_polish(polished_results=polished_results, all_reads=sequences)
    final_consensuses = consensus_object_from_polished_results(
        polished_results=polished_results
    )
    return final_consensuses


# %%
# align the new consensus objects to a given reference
# path_reference = Path("/home/vinzenz/development/LRSV-detection/development/reference/hs37d5/hs37d5.fa")
# tmp_consensuses = tempfile.NamedTemporaryFile(delete=True, suffix='.fasta', prefix='tmp_consensuses.')
# with open(tmp_consensuses.name, 'w') as f:
#     for name,seq,rcs in polished_results:
#         SeqIO.write(seq, f, "fasta")

# util.align_reads_with_minimap(
#     aln_args=" -U50,500",
#     bamout="/home/vinzenz/development/LRSV-detection/development/HG/test.bam",
#     reads=tmp_consensuses.name,
#     reference=path_reference,
#     tech="asm10",
#     threads=8,
# )

# %%
