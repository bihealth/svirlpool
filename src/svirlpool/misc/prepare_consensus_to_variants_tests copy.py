# %%
import subprocess
import tempfile
from pathlib import Path
from shlex import split

import pysam
from Bio import SeqIO

from . import (
    alignments_to_raf,
    consensus_lib,
    consensus_to_variants,
    crs_to_dbs,
    raf_to_signal,
    signal_to_signaldepth,
    util,
)

# %%

workdir = Path(
    "/home/vinzenz/development/LRSV-detection/tests/data/consensus_to_variants"
)

# %%
# generate data for consensus_to_variants tests
# first generate a reference, then reads with deletions, insertions

tmp_dir = tempfile.TemporaryDirectory()
tmp_dir_path = Path(tmp_dir.name)
# tmp_dir_path = workdir

reference = util.generate_sequence(size=10_000, seed=17)

readends = [10_000, 8_000, 7_000]
readstarts = [0, 1_000, 2_000]
sv_location = 4_000
sv_length = 500
margin = 1000
cut_interval = (2000, 8000)
inserted_sequence = util.generate_sequence(size=600, seed=23)

reads_deletion = [
    util.generate_read_from_rules(
        rules=[(util.delete_interval, sv_location + 10, sv_location + sv_length)],
        reference=reference,
        refStart=readstarts[0],
        refEnd=readends[0],
    ),
    util.generate_read_from_rules(
        rules=[(util.delete_interval, sv_location + 5, sv_location + sv_length)],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[2],
    ),
    util.generate_read_from_rules(
        rules=[(util.delete_interval, sv_location + 0, sv_location + sv_length)],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[1],
    ),
    util.generate_read_from_rules(
        rules=[(util.delete_interval, sv_location - 5, sv_location + sv_length)],
        reference=reference,
        refStart=readstarts[1],
        refEnd=readends[1],
    ),
    util.generate_read_from_rules(
        rules=[], reference=reference, refStart=sv_location + margin, refEnd=readends[0]
    ),
    util.generate_read_from_rules(
        rules=[], reference=reference, refStart=sv_location + margin, refEnd=readends[1]
    ),
]

reads_deletion_revcomp = [
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (util.delete_interval, sv_location + 10, sv_location + sv_length),
        ],
        reference=reference,
        refStart=readstarts[0],
        refEnd=readends[0],
    ),
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (util.delete_interval, sv_location + 5, sv_location + sv_length),
        ],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[2],
    ),
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (util.delete_interval, sv_location + 0, sv_location + sv_length),
        ],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[1],
    ),
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (util.delete_interval, sv_location - 5, sv_location + sv_length),
        ],
        reference=reference,
        refStart=readstarts[1],
        refEnd=readends[1],
    ),
    util.generate_read_from_rules(
        rules=[(util.reverse_complement,)],
        reference=reference,
        refStart=len(reference) - (sv_location - margin),
        refEnd=readends[0],
    ),
    util.generate_read_from_rules(
        rules=[(util.reverse_complement,)],
        reference=reference,
        refStart=len(reference) - (sv_location - margin),
        refEnd=readends[1],
    ),
]

reads_insertion = [
    util.generate_read_from_rules(
        rules=[
            (
                util.insertion,
                sv_location + 10,
                -1,
                inserted_sequence[20 : sv_length + 20],
            )
        ],
        reference=reference,
        refStart=readstarts[0],
        refEnd=readends[0] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[
            (
                util.insertion,
                sv_location + 5,
                -1,
                inserted_sequence[15 : sv_length + 15],
            )
        ],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[2] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[
            (util.insertion, sv_location, -1, inserted_sequence[10 : sv_length + 10])
        ],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[1] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[
            (util.insertion, sv_location - 5, -1, inserted_sequence[5 : sv_length + 25])
        ],
        reference=reference,
        refStart=readstarts[1],
        refEnd=readends[1] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[], reference=reference, refStart=sv_location + margin, refEnd=readends[0]
    ),
    util.generate_read_from_rules(
        rules=[], reference=reference, refStart=sv_location + margin, refEnd=readends[1]
    ),
]

reads_insertion_revcomp = [
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (
                util.insertion,
                sv_location + 10,
                -1,
                inserted_sequence[20 : sv_length + 20],
            ),
        ],
        reference=reference,
        refStart=readstarts[0],
        refEnd=readends[0] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (
                util.insertion,
                sv_location + 5,
                -1,
                inserted_sequence[15 : sv_length + 15],
            ),
        ],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[2] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (util.insertion, sv_location, -1, inserted_sequence[10 : sv_length + 10]),
        ],
        reference=reference,
        refStart=readstarts[2],
        refEnd=readends[1] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[
            (util.reverse_complement,),
            (
                util.insertion,
                sv_location - 5,
                -1,
                inserted_sequence[5 : sv_length + 25],
            ),
        ],
        reference=reference,
        refStart=readstarts[1],
        refEnd=readends[1] + sv_length,
    ),
    util.generate_read_from_rules(
        rules=[(util.reverse_complement,)],
        reference=reference,
        refStart=len(reference) - (sv_location - margin),
        refEnd=readends[0],
    ),
    util.generate_read_from_rules(
        rules=[(util.reverse_complement,)],
        reference=reference,
        refStart=len(reference) - (sv_location - margin),
        refEnd=readends[1],
    ),
]

# reads_paths_pairs = list(zip([reads_deletion,reads_deletion_revcomp,reads_insertion,reads_insertion_revcomp],paths_reads))

data = {
    "deletion": {
        "reads": reads_deletion,
        "path_reads_to_ref": tmp_dir_path / "deletion.reads.fasta",
    },
    "insertion": {
        "reads": reads_insertion,
        "path_reads_to_ref": tmp_dir_path / "insertion.reads.fasta",
    },
    "deletion_revcomp": {
        "reads": reads_deletion_revcomp,
        "path_reads_to_ref": tmp_dir_path / "deletion_revcomp.reads.fasta",
    },
    "insertion_revcomp": {
        "reads": reads_insertion_revcomp,
        "path_reads_to_ref": tmp_dir_path / "insertion_revcomp.reads.fasta",
    },
}

data["deletion"]["consensus"] = util.generate_read_from_rules(
    rules=[(util.delete_interval, sv_location, sv_location + sv_length)],
    reference=reference,
    refStart=readstarts[0],
    refEnd=readends[0],
)
data["deletion_revcomp"]["consensus"] = util.generate_read_from_rules(
    rules=[(util.delete_interval, sv_location, sv_location + sv_length)],
    reference=reference,
    refStart=readstarts[0],
    refEnd=readends[0],
)

data["insertion"]["consensus"] = util.generate_read_from_rules(
    rules=[(util.insertion, 4010, -1, inserted_sequence[20 : sv_length + 20])],
    reference=reference,
    refStart=readstarts[0],
    refEnd=readends[0],
)
data["insertion_revcomp"]["consensus"] = util.generate_read_from_rules(
    rules=[(util.insertion, 4010, -1, inserted_sequence[20 : sv_length + 20])],
    reference=reference,
    refStart=readstarts[0],
    refEnd=readends[0],
)


# %%
# write reads to output files
path_reference = tmp_dir_path / "reference.fasta"
util.write_sequences_to_fasta(
    seqs=[reference], path=path_reference, chrnames=True, prefix="chr"
)
util.index_reference(path=path_reference)
# load ref dict
id_to_chr_dict = util.create_ref_dict(reference=path_reference)
chr_to_id_dict = {v: k for k, v in id_to_chr_dict.items()}
# align reads to reference


for k, d in data.items():
    seqs = d["reads"]
    path = d["path_reads_to_ref"]
    d["path_aln_reads_to_ref"] = path.with_suffix(".toRef.bam")
    util.write_sequences_to_fasta(seqs=seqs, path=path, chrnames=True, prefix="read")
    util.align_reads_with_minimap(
        reference=path_reference,
        reads=path,
        bamout=d["path_aln_reads_to_ref"],
        threads=1,
    )
    # load all alignments to data
    d["alns_reads_to_ref"] = list(
        pysam.AlignmentFile(d["path_aln_reads_to_ref"]).fetch()
    )
    # create rafs
    d["rafs"] = [
        alignments_to_raf.alignment_to_alignment_segment(
            a=a, sampleID=0, min_signal_size=10, min_clipped_length=30
        )
        for a in d["alns_reads_to_ref"]
    ]
    # --- singals --- #
    # signals are also empty rows if there is no signal in a raf
    d["signals"] = [
        row
        for raf in d["rafs"]
        for row in raf_to_signal.raf_to_signals(
            dict_fai_reverse=id_to_chr_dict,
            raf=raf,
            min_signal_size=10,
            min_breakend_size=10,
            min_mapq=0,
            window_size=400,
        )
    ]
    # add depth to signals
    # write signals to tmp file
    d["path_signals"] = tmp_dir_path / f"{k}.signals"
    with open(d["path_signals"], "w") as f:
        for row in d["signals"]:
            if len(row) > 0:
                print(*row, sep="\t", file=f)
    d["path_signaldepth"] = tmp_dir_path / f"{k}.signaldepth"
    signal_to_signaldepth.add_depth_to_signal(
        bam=d["path_aln_reads_to_ref"],
        signals=d["path_signals"],
        output=d["path_signaldepth"],
        threads=1,
    )
    # load signaldepths to list of rows to d["signaldepths"]
    # d['signaldepth'] = [line.strip().split('\t') for line in open(d['path_signaldepth'],'r')]
    d["signaldepth"] = []
    with open(d["path_signaldepth"], "r") as f:
        for line in f:
            # parse line from str to mixed types
            d["signaldepth"].append(
                [
                    a(b)
                    for a, b in zip(
                        (int, int, int, str, int, int, int, str, int, int, str, int),
                        (line.strip().split("\t")),
                    )
                ]
            )
    # manually add signalstrength to signaldepths (1.5 to all entries)
    # d['path_signalstrength'] = tmp_dir_path / f"{k}.signalstrength"
    # with open(tmp_dir_path / f"{k}.signaldepth",'r') as f:
    #     with open(d['path_signalstrength'],'w') as g:
    #         # read rows in f:
    #         rows = f.readlines()
    #         for row in rows:
    #             print(row.strip(),1.5,sep='\t',file=g)
    # --- databases --- #
    # create databases
    d["path_db_reads"] = tmp_dir_path / f"{k}.reads.db"
    crs_to_dbs.create_reads_db(path_reads_db=d["path_db_reads"], timeout=300.0)
    d["path_db_intervals"] = tmp_dir_path / f"{k}.intervals.db"
    crs_to_dbs.create_intervals_db(
        path_intervals_db=d["path_db_intervals"], timeout=300.0
    )
    # add reads to reads db
    buffer_reads = [
        (a.query_name, 0, a.query_sequence, a.query_qualities)
        for a in d["alns_reads_to_ref"]
    ]
    crs_to_dbs.write_reads_to_reads_db(
        buffer_reads=buffer_reads, path_reads_database=d["path_db_reads"], timeout=300.0
    )
    buffer_intervals = [
        (
            a.query_name,
            0,
            0,
            *util.get_interval_on_read_in_region(
                a=a, start=cut_interval[0], end=cut_interval[1]
            ),
            0,
            *cut_interval,
            0,
            not a.is_reverse,
            len(a.query_sequence),
            a.query_alignment_start,
            a.query_alignment_end,
            max(a.reference_start, cut_interval[0]),
            min(a.reference_end, cut_interval[1]),
        )
        for a in d["alns_reads_to_ref"]
    ]
    crs_to_dbs.write_intervals_to_intervals_db(
        buffer_intervals=buffer_intervals,
        path_intervals_database=d["path_db_intervals"],
        timeout=300.0,
    )
    # --- cut reads --- #
    # get cut reads
    cutreads = consensus_lib.trim_reads(
        path_intervals_db=d["path_db_intervals"],
        path_reads_db=d["path_db_reads"],
        mcrID=0,
    )
    # write cutreads to fasta file
    d["path_cutreads"] = tmp_dir_path / f"{k}.cutreads.fasta"
    with open(d["path_cutreads"], "w") as f:
        for rec in cutreads[0].values():
            SeqIO.write(rec, f, "fasta")
    # --- consensus --- #
    # write fake consensus to file
    d["path_consensus"] = tmp_dir_path / f"{k}.consensus.fasta"
    util.write_sequences_to_fasta(
        seqs=[d["consensus"]],
        path=d["path_consensus"],
        chrnames=True,
        prefix="consensus.0.",
    )
    # index consensus
    util.index_reference(path=d["path_consensus"])
    # align consensus to reference
    d["consensus_to_ref"] = tmp_dir_path / f"{k}.consensus.toRef.bam"
    util.align_reads_with_minimap(
        reference=path_reference,
        reads=d["path_consensus"],
        bamout=d["consensus_to_ref"],
        threads=1,
    )
    # align cut reads to consensus
    d["cutreads_to_consensus"] = tmp_dir_path / f"{k}.cutreads.toConsensus.bam"
    util.align_reads_with_minimap(
        reference=d["path_consensus"],
        reads=d["path_cutreads"],
        bamout=d["cutreads_to_consensus"],
        threads=1,
    )
    # soft link to tmp_dir_path / f"{k}.consensus.bam"
    cmd_link = f"ln -s {d['consensus_to_ref']} {tmp_dir_path / f'{k}.consensus.bam'}"
    cmd_link_index = (
        f"ln -s {d['consensus_to_ref']}.bai {tmp_dir_path / f'{k}.consensus.bam.bai'}"
    )
    subprocess.check_call(split(cmd_link))
    subprocess.check_call(split(cmd_link_index))
    # create fake SVcandidates
    SVCs = consensus_to_variants.consensus_alignments_to_SVs(
        consensus_to_ref_alignments=d["consensus_to_ref"],
        min_clipped_length=50,
        min_signal_size=30,
    )
    assert len(SVCs) == 1
    svc = SVCs[0]
    assert svc.svSize == sv_length
    assert sv_location - 5 <= svc.readStart <= sv_location + 10
    # create signaldepth database
    d["path_db_signaldepth"] = tmp_dir_path / f"{k}.signaldepth.db"
    consensus_to_variants.create_signaldepths_db(
        path_database=d["path_db_signaldepth"], path_signaldepths=d["path_signaldepth"]
    )
    # trace back the svc
    SVCs_backtraced = consensus_to_variants.traceback_SV_candidates(
        mcrID=0,
        path_consensus_to_reference_alignments=d["consensus_to_ref"],
        path_cutreads_to_consensus_alignments=d["cutreads_to_consensus"],
        path_intervals_db=d["path_db_intervals"],
        path_signaldepths_db=d["path_db_signaldepth"],
        SVCs=SVCs,
    )
    d["svc"] = svc
    # assert the correct depth
    print(f"{k}, ", svc)
    print("=========================================================")

    # assert len(SVCs_backtraced) == 1
    # svc = SVCs_backtraced[0]
    # assert svc.attributes['depths'][0]['ref_cov_l'] == 4
    # if k.startswith('insertion'):
    #     assert svc.attributes['depths'][0]['ref_cov_r'] == 4
    # else:
    #     assert svc.attributes['depths'][0]['ref_cov_r'] == -1
    # assert svc.attributes['depths'][0]['alt_cov_l'] == 4
    # if k.startswith('insertion'):
    #     assert svc.attributes['depths'][0]['alt_cov_r'] == 4
    # else:
    #     assert svc.attributes['depths'][0]['alt_cov_r'] == -1
    # assert svc.attributes['depths'][0]['score_l'] == 1.0
    # if k.startswith('insertion'):
    #     assert svc.attributes['depths'][0]['score_r'] == 1.0
    # else:
    #     assert svc.attributes['depths'][0]['score_r'] == 0.0


# %%
# print results
for k, d in data.items():
    print(k, d["svc"])

# %%
