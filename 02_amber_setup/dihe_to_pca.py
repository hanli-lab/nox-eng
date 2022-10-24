import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser()

    parser.add_argument("--mut_mna", type=str, required=True, help="path to mut mna dihedrals")
    parser.add_argument("--wt_mna", type=str, required=True, help="path to wt mna dihedrals")
    parser.add_argument("--wt_nad", type=str, required=True, help="path to wt nad dihedrals")

    return parser.parse_args()


def dihe_to_pca(
    mut_mna_df: pd.DataFrame, wt_mna_df: pd.DataFrame, wt_nad_df: pd.DataFrame
) -> np.ndarray:

    # concat all samples and pca
    df_cat = pd.concat([wt_mna_df, mut_mna_df, wt_nad_df], axis=0)
    pipe = make_pipeline(StandardScaler(), PCA(n_components=2))
    x_pca = pipe.fit_transform(df_cat)

    return x_pca


if __name__ == "__main__":

    args = parse_args()

    mut_mna_df = pd.read_csv(args.mut_mna)
    wt_mna_df = pd.read_csv(args.wt_mna)
    wt_nad_df = pd.read_csv(args.wt_nad)

    x_pca = dihe_to_pca(mut_mna_df, wt_mna_df, wt_nad_df)
    np.save('./x_pca_out.npy', x_pca)
