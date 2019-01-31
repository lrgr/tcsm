import numpy as np
import pandas as pd

if __name__ == '__main__':
    random_state=int(snakemake.params["seed"])+int(snakemake.wildcards["iter_num"])
    np.random.seed(random_state)
    features = np.random.choice([0,1], size=int(snakemake.wildcards["total_n_samples"]), p=[.75, .25])
    df = pd.DataFrame(data={"feature1": features, "default": 1})
    df.to_csv(snakemake.output[0], sep="\t")
