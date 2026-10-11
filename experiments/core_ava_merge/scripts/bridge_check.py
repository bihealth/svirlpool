"""Run the bridging scenario of the v6 test against whatever svirlpool is on PYTHONPATH."""

import random
import sys

sys.path.insert(0, "/home/mayv_c/development/svirlpool/.claude/worktrees/agent-a7ab379d77ee1935e/tests")
import test_core_alignment_merge as t  # noqa: E402

rng = random.Random(7)
left, mid, right = t.rand_seq(1500, rng), t.rand_seq(400, rng), t.rand_seq(1500, rng)
i1, i2 = t.rand_seq(100, rng), t.rand_seq(120, rng)
child = left + i1 + mid + i2 + right
mother = left + i1 + mid[:200]
father = mid[200:] + i2 + right
group = [
    t.ins("HG002", "1.0", child, 1500, 100, 1500),
    t.ins("HG002", "1.0", child, 2000, 120, 1900),
    t.ins("HG004", "5.0", mother, 1500, 100, 1500),
    t.ins("HG003", "9.0", father, 200, 120, 1900),
]
cores = t.cores_of(
    c=("HG002", "1.0", child, False), m=("HG004", "5.0", mother, False), f=("HG003", "9.0", father, False)
)
merged, stats = t.run(group, cores, old=t.always)
print([sorted(p.samplename for p in svc.svPatterns) for svc in merged])
