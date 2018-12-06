import pandas as pd

if __name__ == '__main__':
    output = []
    for f in snakemake.input:
        output.append(pd.read_csv(f, sep="\t", index_col=0))
    pd.concat(output).fillna(value=0).to_csv(snakemake.output[0], sep="\t")
