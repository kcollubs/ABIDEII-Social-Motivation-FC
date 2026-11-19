import pandas as pd
from statsmodels.stats.multitest import multipletests

# Exploratory reward-connectivity analyses
tests = pd.DataFrame({
    "Analysis": [
        "Excluding influential observations",
        "Three sites, >=4 minutes",
        "Three sites, >=5 minutes",
        "Lower-motion subset (FD <=0.12)",
        "Duration-group interaction",
        "Duration interaction + site-specific slopes"
    ],
    "p_uncorrected": [
        0.0831,
        0.1314,
        0.0316,
        0.5503,
        0.0413,
        0.0753
    ]
})

# Benjamini-Hochberg FDR correction
reject, p_fdr, _, _ = multipletests(
    tests["p_uncorrected"],
    alpha=0.05,
    method="fdr_bh"
)

tests["p_FDR"] = p_fdr
tests["Significant_FDR"] = reject

print("\n=== EXPLORATORY REWARD ANALYSES: BH-FDR ===\n")
print(tests.to_string(index=False, float_format="%.4f"))

print("\nNumber of tests:", len(tests))
print("FDR threshold: 0.05")
print("Significant after FDR:", sum(reject))
