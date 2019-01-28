import pandas as pd

def get_hr_status(row):
    if row["BRCA1"] == 1:
        return "BRCA1"
    if row["BRCA2"] == 1:
        return "BRCA2"
    if row["RAD51C"] == 1:
        return "RAD51C"
    return "wildtype"

if __name__ == '__main__':
    bi_inactivation_df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    silencing_df = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
    samples = list(set(bi_inactivation_df.index).intersection(silencing_df.index))
    bi_inactivation_df = bi_inactivation_df.loc[samples]
    silencing_df = silencing_df.loc[samples]
    # now need to combine silencing_df with bi_inactivation_df
    genes = list(set(bi_inactivation_df.columns).intersection(silencing_df.columns))
    silencing_df = silencing_df[genes]
    bi_inactivation_df = bi_inactivation_df[genes]
    feature_df = bi_inactivation_df + silencing_df
    feature_df["HR Status"] = feature_df.apply(get_hr_status, axis=1)
    feature_df["BRCA1BRCA2"] = feature_df["BRCA1"] + feature_df["BRCA2"]
    feature_df["BRCA1BRCA2RAD51C"] = feature_df["BRCA1BRCA2"] + feature_df["RAD51C"]
    feature_df["BRCA1BRCA2RAD51C"] = feature_df["BRCA1BRCA2RAD51C"].apply(lambda x: min(x, 1))
    feature_df["BRCA1BRCA2RAD51CATM"] = feature_df["BRCA1BRCA2RAD51C"] + feature_df["ATM"]
    feature_df["BRCA1BRCA2RAD51CATM"] = feature_df["BRCA1BRCA2RAD51CATM"].apply(lambda x: min(x, 1))
    feature_df["BRCA1BRCA2RAD51CFANCM"] = feature_df["BRCA1BRCA2RAD51C"] + feature_df["FANCM"]
    feature_df["BRCA1BRCA2RAD51CFANCM"] = feature_df["BRCA1BRCA2RAD51CFANCM"].apply(lambda x: min(x, 1))
    feature_df["BRCA1BRCA2RAD51CFANCF"] = feature_df["BRCA1BRCA2RAD51C"] + feature_df["FANCF"]
    feature_df["BRCA1BRCA2RAD51CFANCF"] = feature_df["BRCA1BRCA2RAD51CFANCF"].apply(lambda x: min(x, 1))
    feature_df["HRINACTIVATIONS"] = feature_df["BRCA1"]+feature_df["BRCA2"]+ feature_df["RAD51C"] + feature_df["ATM"] + feature_df["CHEK2"] + feature_df["FANCM"] + feature_df["FANCF"]
    feature_df["HRINACTIVATIONS"] = feature_df["HRINACTIVATIONS"].apply(lambda x: min(x, 1))
    feature_df = feature_df.sort_index()
    feature_df.to_csv(snakemake.output[0], sep="\t")
