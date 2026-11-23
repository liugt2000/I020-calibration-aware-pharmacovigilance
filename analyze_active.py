"""Active-comparator ROR + attenuation waterfall data."""
import pandas as pd, numpy as np
from pathlib import Path
HERE=Path(__file__).parent
df=pd.read_csv(HERE/"counts_active.csv",encoding="utf-8-sig")
Z=1.959963984540054
rows=[]
for _,r in df.iterrows():
    a=int(r["a"]); it=int(r["index_total"]); c=int(r["c"]); ct=int(r["comparator_total"])
    b=it-a; d=ct-c
    if min(a,b,c,d)==0:
        a,b,c,d=a+.5,b+.5,c+.5,d+.5
    ror=(a*d)/(b*c)
    se=np.sqrt(1/a+1/b+1/c+1/d)
    lo=np.exp(np.log(ror)-Z*se); hi=np.exp(np.log(ror)+Z*se)
    rows.append(dict(pair=r["pair"],index=r["index"],comparators=r["comparators"],
                     a=int(r["a"]),c=int(r["c"]),ror_active=ror,lo=lo,hi=hi,
                     lnr_active=np.log(ror)))
out=pd.DataFrame(rows)
print(out.round(3).to_string(index=False))
out.to_csv(HERE/"results_active.csv",index=False,encoding="utf-8-sig")
print("\nwrote results_active.csv")
