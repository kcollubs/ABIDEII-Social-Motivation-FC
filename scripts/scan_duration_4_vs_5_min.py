import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC + FC_AI_a24pr "
    "+ AGE_AT_SCAN + SEX + DCAN_retained_mean_FD + C(site)"
)

# Original model: >= 4 minutes
original = smf.ols(formula, data=df).fit()

# Exploratory model: >= 5 minutes
df5 = df[df["usable_minutes"] >= 5].copy()
five_min = smf.ols(formula, data=df5).fit()

variables = ["FC_NAc_rACC", "FC_AI_a24pr"]

print("\n=== SAMPLE COMPARISON ===")
print(f"Original N: {int(original.nobs)}")
print(f"5-minute N: {int(five_min.nobs)}")
print(f"5-minute sites: {', '.join(sorted(df5.site.unique()))}")

print("\n=== MODEL FIT ===")
print(f"Original R²: {original.rsquared:.4f}")
print(f"5-minute R²: {five_min.rsquared:.4f}")
print(f"Original adjusted R²: {original.rsquared_adj:.4f}")
print(f"5-minute adjusted R²: {five_min.rsquared_adj:.4f}")

print("\n=== CONNECTIVITY PREDICTORS ===")

for variable in variables:
    print(f"\n{variable}")
    print(
        f"Original: b={original.params[variable]:.4f}, "
        f"p={original.pvalues[variable]:.4f}"
    )
    print(
        f"5-minute: b={five_min.params[variable]:.4f}, "
        f"p={five_min.pvalues[variable]:.4f}"
    )

# Does adding connectivity improve the 5-minute model?
baseline_formula = (
    "SRS_MOTIVATION_T ~ AGE_AT_SCAN + SEX "
    "+ DCAN_retained_mean_FD + C(site)"
)

baseline5 = smf.ols(baseline_formula, data=df5).fit()

print("\n=== INCREMENTAL CONNECTIVITY CONTRIBUTION ===")
print(f"5-minute baseline R²: {baseline5.rsquared:.4f}")
print(f"5-minute full R²: {five_min.rsquared:.4f}")
print(f"Change in R²: {five_min.rsquared - baseline5.rsquared:.4f}")

f_stat, p_value, df_diff = five_min.compare_f_test(baseline5)

print(f"Nested F-test: F={f_stat:.4f}, p={p_value:.4f}")
