import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

baseline = smf.ols(
    "SRS_MOTIVATION_T ~ AGE_AT_SCAN + C(SEX) "
    "+ DCAN_retained_mean_FD + C(site)",
    data=df
).fit()

full = smf.ols(
    "SRS_MOTIVATION_T ~ AGE_AT_SCAN + C(SEX) "
    "+ DCAN_retained_mean_FD + C(site) "
    "+ FC_NAc_rACC + FC_AI_a24pr",
    data=df
).fit()

print("Baseline R²:", round(baseline.rsquared, 4))
print("Full R²:", round(full.rsquared, 4))
print("Change in R²:", round(full.rsquared - baseline.rsquared, 4))

print("\nNested model comparison:")
print(full.compare_f_test(baseline))
