import pandas as pd


# exome_signatures should be a pandas DataFrame
# exome_counts and genome_counts should both be pandas Series
def convert_exome_sigs_to_genome_sigs(exome_signatures, exome_counts, genome_counts):
    output = []
    # we want to convert the signatures, one category at a time
    for category in exome_signatures.columns:
        # first we need to find out the trinucleotide for the category
        trinucleotide = category[0] + category[2] + category[6]
        trinucleotide_exome_count = exome_counts[trinucleotide]
        trinucleotide_genome_count = genome_counts[trinucleotide]
        true_value = exome_signatures[category]/trinucleotide_exome_count
        new_sig_category = true_value*trinucleotide_genome_count
        output.append(new_sig_category)
    genome_signatures = pd.concat(output, axis=1)
    genome_signatures.index = exome_signatures.index
    return genome_signatures.div(genome_signatures.sum(axis=1), axis=0)

def test_convert_exome_sigs_to_genome_sigs():
    exome_counts = pd.Series(data={"CCC": 700, "TTT": 300})
    genome_counts = pd.Series(data={"CCC": 500, "TTT": 500})
    exome_signatures = pd.DataFrame(data={"C[C>G]C": [.875, .5], "T[T>A]T": [.125, .5]})
    genome_signatures = convert_exome_sigs_to_genome_sigs(exome_signatures, exome_counts, genome_counts)
    assert(genome_signatures.loc[0, "C[C>G]C"] == .75)
    assert(genome_signatures.loc[0, "T[T>A]T"] == .25)


if __name__ == '__main__':
    exome_counts = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)["x"]
    genome_counts = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)["x"]
    exome_signatures = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)
    genome_signatures = convert_exome_sigs_to_genome_sigs(exome_signatures, exome_counts, genome_counts)
    genome_signatures.to_csv(snakemake.output[0], sep="\t")
