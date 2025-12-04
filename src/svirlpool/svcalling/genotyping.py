import attrs
import cattrs


@attrs.define
class GenotypeMeasurement:
    start_on_consensus: int
    # estimated_total_depth_start:int
    supporting_reads_start: list[str]
    end_on_consensus: int | None = None  # None for insertions and breakends
    # estimated_total_depth_end:int | None = None # None for insertions and breakends
    supporting_reads_end: list[str] | None = None  # None for insertions and breakends

    def unstructure(self):
        return cattrs.unstructure(self)
