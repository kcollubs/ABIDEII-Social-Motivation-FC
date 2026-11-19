import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

original = pd.read_csv("results/analysis_dataset_preliminary.csv")
fc = pd.read_csv("results/fc_salience_exploratory.csv")

pathways = [
    "FC_AI_a24pr",
    "FC_AI_a32pr",
    "FC_FOP4_a24pr",
    "FC_FOP4_a32pr"
]

df = original.drop(columns=["FC_AI_a24pr"]).merge(
    fc[["participant_id", "site"] + pathways],
    on=["participant_id", "site"],
    how="left",
    validate="one_to_one"
)

df = df[df["usable_minutes"] >= 5].copy()

assert len(df) == 70, f"Unexpected N: {len(df)}"
assert df[pathways].notna().all().all(), "Missing connectivity values"

results = []

for salience in pathways:
    formula = (
        f"SRS_MOTIVATION_T ~ FC_NAc_rACC + {salience} "
        "+ AGE_AT_SCAN + SEX + DCAN_retained_mean_FD + C(site)"
    )

    model = smf.ols(formula, data=df).fit()

    results.append({
        "Pathway": salience,
        "N": int(model.nobs),
        "Salience b": model.params[salience],
        "Salience p": model.pvalues[salience],
        "Reward b": model.params["FC_NAc_rACC"],
        "Reward p": model.pvalues["FC_NAc_rACC"],
        "R2": model.rsquared
    })

out = pd.DataFrame(results)

out["Salience FDR p"] = multipletests(
    out["Salience p"], method="fdr_bh"
)[1]

print("\n=== FIVE-MINUTE SALIENCE PARCEL COMPARISON ===")
print(out.round(4).to_string(index=False))
