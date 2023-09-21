# NOTE: 2021.02.19
#  How new sample metadata was added to sampels.tsv (issue #2, after commit 2f95b929)

import pandas

f1="config/samples.tsv"
f2="config/edited_metadata.tsv"

df1 = pandas.read_csv(f1, sep="\t", header=0, index_col=0)
df2 = pandas.read_csv(f2, sep="\t", header=0, index_col=0)

df_merged = df1.merge(
    right=df2,
    how="left",
    left_index=True,
    right_index=True
)

df_merged.to_csv("config/samples_metadata.tsv", sep="\t", na_rep='NA', header=True, index=True, index_label="Sample")
