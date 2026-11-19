import pandas as pd

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

print("SRS Social Motivation by site:")
print(
    df.groupby("site")["SRS_MOTIVATION_T"]
    .agg(["count", "mean", "std", "min", "max"])
    .round(2)
    .to_string()
)

print("\nOverall SRS Social Motivation:")
print(df["SRS_MOTIVATION_T"].describe().round(2))
