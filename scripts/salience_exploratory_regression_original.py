import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

original = pd.read_csv("results/analysis_dataset_preliminary.csv")
fc = pd.read_csv("results/fc_salience_exploratory.csv")

alternatives = [
    "FC_AI_a24pr",
    "FC_AI_a32pr",
    "FC_FOP4_a24pr",
    "FC_FOP4_a32pr"
]

# Merge the new connectivity measures into the original analysis sample
df = original.drop(columns=["FC_AI_a24pr"]).merge(
    fc[["participant_id", "site"] + alternatives],
    on=["participant_id", "site"],
    how="left",
    validate="one_to_one"
)

assert len(df) == 119, "Unexpected sample size"
assert df[alternatives].notna().all().all(), "Missing FC values"

results = []

for salience in alternatives:
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

# Exploratory correction across four salience tests
out["Salience FDR p"] = multipletests(
    out["Salience p"], method="fdr_bh"
)[1]

print("\n=== EXPLORATORY SALIENCE PARCEL COMPARISON ===")
print(out.round(4).to_string(index=False))
