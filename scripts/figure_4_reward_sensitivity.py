import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

# Load original analysis dataset
df = pd.read_csv("results/analysis_dataset_preliminary.csv")

os.makedirs("results/figures", exist_ok=True)

formula = (
    "SRS_MOTIVATION_T ~ FC_NAc_rACC + FC_AI_a24pr "
    "+ AGE_AT_SCAN + SEX + DCAN_retained_mean_FD + C(site)"
)

results = []

def fit_model(data, label):
    model = smf.ols(formula, data=data).fit()
    term = "FC_NAc_rACC"
    ci = model.conf_int().loc[term]

    results.append({
        "Analysis": label,
        "N": int(model.nobs),
        "b": model.params[term],
        "CI_low": ci.iloc[0],
        "CI_high": ci.iloc[1],
        "p": model.pvalues[term]
    })

# 1. Primary model
fit_model(df, "Primary analysis")

# 2. Exclude influential observations using Cook's distance
primary_model = smf.ols(formula, data=df).fit()
cooks_d = primary_model.get_influence().cooks_distance[0]
threshold = 4 / len(df)

df_no_influential = df.loc[
    cooks_d <= threshold
].copy()

fit_model(df_no_influential, "Excluding influential cases")

# 3. Three-site sample, >=4 minutes
three_sites = df[
    df["site"].isin(["KKI", "OHSU", "SDSU"])
].copy()

fit_model(three_sites, "Three sites: ≥4 min")

# 4. Three-site sample, >=5 minutes
five_min = three_sites[
    three_sites["usable_minutes"] >= 5
].copy()

fit_model(five_min, "Three sites: ≥5 min")

# 5. Lower-motion sample
low_motion = df[
    df["DCAN_retained_mean_FD"] <= 0.12
].copy()

fit_model(low_motion, "Lower motion: FD ≤0.12")

# Convert results to dataframe
res = pd.DataFrame(results)

# Forest plot
fig, ax = plt.subplots(figsize=(9, 5.2))

y = np.arange(len(res))[::-1]

for i, row in res.iterrows():
    ax.errorbar(
        row["b"],
        y[i],
        xerr=[
            [row["b"] - row["CI_low"]],
            [row["CI_high"] - row["b"]]
        ],
        fmt="o",
        color="#245A81",
        ecolor="#245A81",
        markersize=7,
        capsize=4,
        linewidth=1.8
    )

# Zero-effect reference line
ax.axvline(
    0,
    color="gray",
    linestyle="--",
    linewidth=1
)

# Labels
labels = [
    f"{row['Analysis']} (n = {row['N']})"
    for _, row in res.iterrows()
]

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=10)

ax.set_xlabel(
    "Reward-related functional connectivity coefficient (b)",
    fontsize=11
)

ax.set_title(
    "Reward Connectivity Across Sensitivity Analyses",
    fontsize=13,
    pad=15
)

ax.grid(axis="x", alpha=0.15)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

# Save publication-quality outputs
png_path = "results/figures/figure_4_reward_sensitivity_forest.png"
pdf_path = "results/figures/figure_4_reward_sensitivity_forest.pdf"

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.close()

# Print numerical results
print("\n=== FIGURE 4: REWARD FC SENSITIVITY ANALYSES ===\n")
print(
    res.to_string(
        index=False,
        float_format="%.4f"
    )
)

print("\nSaved:")
print(png_path)
print(pdf_path)
