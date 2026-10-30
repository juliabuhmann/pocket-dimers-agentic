import pandas as pd, sys
for p in sys.argv[1:]:
    df = pd.read_parquet(p) if p.endswith('.parquet') else pd.read_csv(p, sep='\t', comment='#')
    print(p, df.shape); print(df.dtypes.to_string()); print(df.head(3).T.to_string()[:3000])
