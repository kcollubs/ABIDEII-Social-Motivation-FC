import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

# Keep the same three sites in both models
df = df[df["site"].isin(["KKI", "OHSU", "SDSU"])].copy()

formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC + FC_AI_a24pr "
    "+ AGE_AT_SCAN + SEX + DCAN_retained_mean_FD + C(site)"
)

baseline = (
    "SRS_MOTIVATION_T ~ AGE_AT_SCAN + SEX "
    "+ DCAN_retained_mean_FD + C(site)"
)

for threshold in [4, 5]:
    subset = df[df["usable_minutes"] >= threshold].copy()

    full = smf.ols(formula, data=subset).fit()
    base = smf.ols(baseline, data=subset).fit()

    f_stat, p_joint, df_diff = full.compare_f_test(base)

    print(f"\n=== {threshold}-MINUTE THRESHOLD ===")
    print(f"N: {int(full.nobs)}")
    print(f"R²: {full.rsquared:.4f}")
    print(f"Adjusted R²: {full.rsquared_adj:.4f}")

    print("\nReward FC (NAc–a24):")
    print(f"  b = {full.params['FC_NAc_rACC']:.4f}")
    print(f"  p = {full.pvalues['FC_NAc_rACC']:.4f}")

    print("\nSalience FC (AVI–a24pr):")
    print(f"  b = {full.params['FC_AI_a24pr']:.4f}")
    print(f"  p = {full.pvalues['FC_AI_a24pr']:.4f}")

    print("\nIncremental connectivity contribution:")
    print(f"  Change in R² = {full.rsquared - base.rsquared:.4f}")
    print(f"  Joint FC p = {p_joint:.4f}")
