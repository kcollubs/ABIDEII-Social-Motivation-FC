import pandas as pd

df = pd.read_csv("results/analysis_dataset_preliminary.csv")

print("=== RETAINED-FRAME MEAN FD DISTRIBUTION ===")
print(df["DCAN_retained_mean_FD"].describe().round(4))

print("\n=== PARTICIPANTS RETAINED BY THRESHOLD ===")

for threshold in [0.08, 0.10, 0.12, 0.15]:
    subset = df[df["DCAN_retained_mean_FD"] <= threshold]

    print(f"\nMean FD <= {threshold:.2f} mm")
    print(f"N = {len(subset)}")
    print(subset.groupby("site").size().to_string())
