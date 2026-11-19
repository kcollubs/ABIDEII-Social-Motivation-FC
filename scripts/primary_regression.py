import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC + FC_AI_a24pr "
    "+ AGE_AT_SCAN + C(SEX) + DCAN_retained_mean_FD + C(site)"
)

baseline = (
    "SRS_MOTIVATION_T ~ AGE_AT_SCAN + C(SEX) "
    "+ DCAN_retained_mean_FD + C(site)"
)

for label, subset in [
    ("Original sample", df),
    ("Mean retained FD <= 0.12 mm",
     df[df["DCAN_retained_mean_FD"] <= 0.12])
]:
    full = smf.ols(formula, data=subset).fit()
    base = smf.ols(baseline, data=subset).fit()
    f_stat, p_joint, _ = full.compare_f_test(base)

    print(f"\n=== {label} ===")
    print(f"N: {int(full.nobs)}")
    print(f"R²: {full.rsquared:.4f}")
    print(f"Adjusted R²: {full.rsquared_adj:.4f}")

    for variable in ["FC_NAc_rACC", "FC_AI_a24pr"]:
        print(
            f"{variable}: "
            f"b={full.params[variable]:.4f}, "
            f"p={full.pvalues[variable]:.4f}"
        )

    print(f"Change in R²: {full.rsquared - base.rsquared:.4f}")
    print(f"Joint FC p: {p_joint:.4f}")
