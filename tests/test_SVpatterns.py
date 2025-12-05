# %%
import json
from gzip import open as gzopen
from pathlib import Path

import cattrs

from svirlpool.localassembly import SVpatterns
from svirlpool.localassembly.SVprimitives import SVprimitive

# here is how the data was saved for testing here in this script
# DEBUG START
# if SVprimitives[0].consensusID == "2.0":
#     # write function input to debugging json file with structured SVpatterns
#     debug_file_out_path = "/data/cephfs-1/work/groups/cubi/users/mayv_c/production/svirlpool/tests/data/SVpatterns/svpattern_INV_parsing.json"
#     with open(debug_file_out_path, "w") as debug_file_out:
#         json.dump(
#             {
#                 "SVprimitives": [svp.unstructure() for svp in SVprimitives]
#             },
#             debug_file_out,
#             indent=4,
#         )
# parse primitives to simple SV types (INS, DEL, BND)
# DEBUG END
# the fiel was then gzipped -> "svpattern_INV_parsing.json.gz"


def test_parse_SVprimitives_to_SVpatterns() -> None:
    # load the svPrimitives to a list from the gzipped json file
    test_data_path = (
        Path(__file__).parent / "data" / "SVpatterns" / "svpattern_INV_parsing.json.gz"
    )
    with gzopen(test_data_path, "rt") as f:
        data = json.load(f)
    SVprimitives = cattrs.structure(data["SVprimitives"], list[SVprimitive])
    # 1) test fourrelations parsing
    breakends = [svp for svp in SVprimitives if svp.sv_type > 2]
    fourrelations = SVpatterns.four_relations_of_group(group=breakends)
    # it should contain at least one SVpatterns.FOURRELATIONS.INVERSION
    assert SVpatterns.FOURRELATIONS.INVERSION in fourrelations, (
        f"Expected at least one INVERSION in fourrelations, got {fourrelations}"
    )
