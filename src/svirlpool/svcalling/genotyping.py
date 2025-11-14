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

    def get_depth(self, metric: str = "max") -> float:
        if metric == "mean":
            if len(self.estimated_total_depths) == 0:
                return 0.0
            return sum(self.estimated_total_depths) / len(self.estimated_total_depths)
        elif metric == "max":
            if len(self.estimated_total_depths) == 0:
                return 0.0
            return max(self.estimated_total_depths)
        elif metric == "min":
            if len(self.estimated_total_depths) == 0:
                return 0.0
            return min(self.estimated_total_depths)
        else:
            raise ValueError(f"Unknown metric: {metric}. Use 'mean', 'max' or 'min'.")

    # def add_depths(self, start_depth:int, end_depth:int | None = None) -> None:
    #     if self.estimated_total_depth_start is None:
    #         self.estimated_total_depth_start = start_depth
    #     else:
    #         self.estimated_total_depth_start += start_depth
    #     if end_depth is not None:
    #         if self.estimated_total_depth_end is None:
    #             self.estimated_total_depth_end = end_depth
    #         else:
    #             self.estimated_total_depth_end += end_depth
    #     else:
    #         self.estimated_total_depth_end = None
