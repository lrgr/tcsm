import numpy as np
import pandas as pd


if __name__ == '__main__':
    lambdas = np.array([[0, -2], [0, -2], [4, -5], [0, -2]])
    df = pd.DataFrame(data=lambdas)
    df.to_csv(snakemake.output[0], sep="\t", index=False)
