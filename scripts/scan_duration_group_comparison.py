import pandas as pd
from scipy.stats import ttest_ind

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

df = df[df["site"].isin(["KKI", "OHSU", "SDSU"])].copy()

df["group"] = df["usable_minutes"].apply(
    lambda x: "Retained (>=5 min)" if x >= 5 else "Excluded (<5 min)"
)

variables = [
    "SRS_MOTIVATION_T",
    "FC_NAc_rACC",
    "FC_AI_a24pr",
    "DCAN_retained_mean_FD",
    "DCAN_pct_censored",
    "usable_minutes"
]

print("\n=== GROUP SIZES ===")
print(df["group"].value_counts())

print("\n=== GROUP COMPARISONS ===")

for var in variables:
    retained = df.loc[df["usable_minutes"] >= 5, var].dropna()
    excluded = df.loc[df["usable_minutes"] < 5, var].dropna()

    stat, p = ttest_ind(retained, excluded, equal_var=False)

    print(f"\n{var}")
    print(f"  Retained: {retained.mean():.4f} (SD={retained.std():.4f})")
    print(f"  Excluded: {excluded.mean():.4f} (SD={excluded.std():.4f})")
    print(f"  Welch t-test p = {p:.4f}")

print("\n=== SITE COMPOSITION ===")
print(pd.crosstab(df["site"], df["group"]))
