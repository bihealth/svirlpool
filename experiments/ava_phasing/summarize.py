import sys, pandas as pd, numpy as np
def summ(path):
    d=pd.read_csv(path,sep="\t")
    d=d[(d.trio_n>=6)]
    # need both haplotypes among labelled reads in the container
    d=d[(d.trio_pat>=2)&(d.trio_mat>=2)]
    r=dict(n=len(d),
      acc_new=d.acc_new.mean(), acc_cur=d.acc_cur.mean(),
      perfect_new=(d.acc_new>=0.999).mean(), perfect_cur=(d.acc_cur>=0.999).mean(),
      split_new=d.split_err_new.mean(), split_cur=d.split_err_cur.mean(),
      join_new=d.join_err_new.mean(), join_cur=d.join_err_cur.mean(),
      grouped=d.grouped_frac.mean(), acc_asg=d.acc_asg.mean(), perfect_asg=(d.acc_asg>=0.999).mean(), asg_frac=d.assigned_trio_frac.mean(), groups=d.n_groups.value_counts().sort_index().to_dict())
    return r
for p in sys.argv[1:]:
    r=summ(p); print(p.split("/")[-1], " ".join(f"{k}={v:.3f}" if isinstance(v,float) else f"{k}={v}" for k,v in r.items()))
