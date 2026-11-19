import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

df = df[df["site"].isin(["KKI", "OHSU", "SDSU"])].copy()
df["five_min"] = (df["usable_minutes"] >= 5).astype(int)

formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC * five_min "
    "+ FC_AI_a24pr + AGE_AT_SCAN + SEX "
    "+ DCAN_retained_mean_FD + C(site)"
)

model = smf.ols(formula, data=df).fit()

print("N:", int(model.nobs))

for term in [
    "FC_NAc_rACC",
    "five_min",
    "FC_NAc_rACC:five_min"
]:
    print(
        f"{term}: "
        f"b={model.params[term]:.4f}, "
        f"p={model.pvalues[term]:.4f}"
    )

print("\n95% CI for interaction:")
print(model.conf_int().loc["FC_NAc_rACC:five_min"].round(4))
