import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

df = df[df["site"].isin(["KKI", "OHSU", "SDSU"])].copy()
df["five_min"] = (df["usable_minutes"] >= 5).astype(int)

# Previous interaction model
original_formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC * five_min "
    "+ FC_AI_a24pr + AGE_AT_SCAN + SEX "
    "+ DCAN_retained_mean_FD + C(site)"
)

# Add site-specific reward-connectivity slopes
adjusted_formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC * five_min "
    "+ FC_NAc_rACC:C(site) "
    "+ FC_AI_a24pr + AGE_AT_SCAN + SEX "
    "+ DCAN_retained_mean_FD + C(site)"
)

original = smf.ols(original_formula, data=df).fit()
adjusted = smf.ols(adjusted_formula, data=df).fit()

term = "FC_NAc_rACC:five_min"

print("=== DURATION INTERACTION COMPARISON ===")

for label, model in [
    ("Original interaction model", original),
    ("Site-specific slope model", adjusted)
]:
    print(f"\n{label}")
    print(f"N: {int(model.nobs)}")
    print(f"Interaction b: {model.params[term]:.4f}")
    print(f"Interaction p: {model.pvalues[term]:.4f}")
    print(
        "95% CI:",
        model.conf_int().loc[term].round(4).tolist()
    )
    print(f"R²: {model.rsquared:.4f}")
    print(f"Adjusted R²: {model.rsquared_adj:.4f}")

print("\n=== SITE-SPECIFIC REWARD SLOPES ===")

for name in adjusted.params.index:
    if "FC_NAc_rACC:C(site)" in name:
        print(
            f"{name}: b={adjusted.params[name]:.4f}, "
            f"p={adjusted.pvalues[name]:.4f}"
        )

print("\n=== MODEL COMPARISON ===")
f_stat, p_value, df_diff = adjusted.compare_f_test(original)

print(f"F = {f_stat:.4f}")
print(f"p = {p_value:.4f}")
print(f"Additional parameters = {int(df_diff)}")
