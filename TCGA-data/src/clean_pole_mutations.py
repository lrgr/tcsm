import pandas as pd

if __name__ == '__main__':
    pole_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    mc_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    samples = [s[0:12] for s in mc_df.index]
    # pole_exo_patients = set(pole_df.index)
    pole_exo_patients = pole_df.loc[pole_df['POLE Mutation'].isin(["p.P286R", "p.V411L"])].index
    pole_exo_df = pd.DataFrame(data=1, index=pole_exo_patients, columns=["POLEexo"])
    pole_exo_df = pole_exo_df.reindex(index=samples, fill_value=0)
    pole_exo_df.to_csv(snakemake.output[0], sep="\t")
