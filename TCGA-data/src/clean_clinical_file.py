import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    df = df.rename(columns={"Tobacco User": "Smoker", "Tobacco Intensity": "PackYears", "Diagnosis Age": "Age"})
    df["Smoker"] = df["Smoker"].map({"Smoker": 1, "Non-Smoker": 0})
    df.to_csv(snakemake.output[0], sep="\t")
