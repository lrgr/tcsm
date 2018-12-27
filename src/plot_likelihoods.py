import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style('whitegrid')
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['font.family'] = 'Courier New'


if __name__ == '__main__':
    output = []
    for i in range(0, len(snakemake.input)):
        out = re.split("_|\.", snakemake.input[i])
        model = out[1]
        print(model)
        df = pd.read_csv(snakemake.input[i], sep="\t", index_col=0).transpose()
        df.index = [model]
        output.append(df)
    df = pd.concat(output)
    df.index = [m.replace("BRCA+GBM+OV+LUAD+UCEC+KIRC+HNSC+LGG+THCA+LUSC+PRAD+SKCM+COAD+STAD+LAML+ESCA+PAAD", "CancerType") for m in df.index]
    models = df.index.unique()

    plt.plot(df.loc[models[0]].transpose(), color="r", marker="o")
    plt.plot(df.loc[models[1]].transpose(), color="b", marker="o")
    plt.xlabel("Number K Signatures Used")
    plt.legend()
    plt.ylabel("Log Likelihood")
    # plt.show()
    plt.savefig(snakemake.output[0])
