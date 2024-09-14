import os
from typing import List

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def _broken_stick(pca_transformer, df: pd.DataFrame):
    """
    Computes the broken stick model for a given set of eigenvalues.
    # Discussed in Diniz-Filho, J. A. F., Sant'Ana, C. E. R. d., & Bini, L. M. (1998). An Eigenvector Method for Estimating Phylogenetic Inertia. Evolution, 52(5), 1247-1262. https://doi.org/10.2307/2411294
    # Original citation Jackson, Donald A. "Stopping rules in principal components analysis: a comparison of heuristical and statistical approaches." Ecology 74.8 (1993): 2204-2214.
    Args:
        eigenvalues (numpy.ndarray): The eigenvalues obtained from PCA.
        total_variance (float): The total variance explained by all components.

    Returns:
        df with only beginning columns selected
    """

    # Allocates each eigenvalue its broken stick value, which is calculated as reversed harmonic series, and then sorted in descending order.
    eigenvalues = pca_transformer.explained_variance_
    num_components = len(eigenvalues)
    broken_stick_values = [
        1 / num_components]  # [total_variance * (1 / np.sum([1 / (i + 1) for i in range(num_components)])) * (1 / (j + 1)) for j in range(num_components)]
    for i in range(num_components)[1:]:
        broken_stick_values.append(broken_stick_values[i - 1] + (1 / (num_components - i + 1)))
    broken_stick_values = [x / num_components for x in broken_stick_values]
    broken_stick_values = sorted(broken_stick_values, reverse=True)
    total_variance = np.sum(eigenvalues)

    # the proportions of the total variance each PC is expected to account for.
    percent_expected = [x / total_variance for x in eigenvalues]

    broken_stick_mask = np.array(percent_expected) > np.array(broken_stick_values)

    num_components_to_keep = np.sum(broken_stick_mask)
    print(f"Number of components to keep: {num_components_to_keep}")
    assert num_components_to_keep == 2  # If this changes, need to update OUTPUT_PCA_COLS
    explained_variance = 0
    for x in range(num_components_to_keep):
        explained_variance += pca_transformer.explained_variance_ratio_[x]
    print(f'Explained_variance of selected components: {explained_variance}')
    X_transformed = df.iloc[:, : num_components_to_keep]

    return X_transformed


def do_PCA(data_X: pd.DataFrame, cols_to_PCA: List[str], plot=True):
    """
    Apply Principal Component Analysis (PCA) to the input data.

    :param plot: bool, whether to plot or not
    :param all_data_y:
    :param cols_to_PCA:
    :param data_X: The input data as a pandas DataFrame.
    :return: A pandas DataFrame containing the transformed data after applying PCA.
    """

    #scale data first
    standard_scaler = StandardScaler()
    standard_scaler.fit(data_X)
    X_scaled = pd.DataFrame(standard_scaler.transform(data_X),
                            index=data_X.index,
                            columns=standard_scaler.get_feature_names_out())
    # Keep same column order
    X_scaled = X_scaled[data_X.columns]


    standard_PCA = ColumnTransformer(
        transformers=[
            ("PCA_cont_vars", PCA(),
             cols_to_PCA)
        ],
        remainder='passthrough', verbose_feature_names_out=False

    )
    standard_PCA.fit(X_scaled)
    pcad_X = pd.DataFrame(standard_PCA.transform(X_scaled),
                          index=X_scaled.index,
                          columns=standard_PCA.get_feature_names_out())
    pcs = [c for c in pcad_X.columns if c not in X_scaled.columns]
    for pc in pcs:
        pcad_X = pcad_X.rename(columns={pc: pc.replace('pca', 'PC')})

    # # Rescale
    # cols_to_rescale = [c for c in df_out.columns if c.startswith('PC')]
    # df_out, s_scaler = do_scaling(pcad_X, cols_to_rescale)
    # ### Plot
    pcs = [c for c in pcad_X.columns if c not in X_scaled.columns]

    if plot:
        pcad_df_to_plot = pcad_X.copy(deep=True)[pcs]
        plot_pca(standard_PCA.named_transformers_["PCA_cont_vars"], pcad_df_to_plot)

    ## Select components with broken stick
    broken_stick_pcs = _broken_stick(standard_PCA.named_transformers_["PCA_cont_vars"], pcad_X[pcs])

    df_out = pcad_X.drop(columns=pcs)
    df_out = pd.concat([df_out, broken_stick_pcs], axis=1)
    return df_out


def plot_pca(pca_transformer, pca_df):
    pca_components_df = pd.DataFrame(pca_transformer.components_, columns=pca_transformer.feature_names_in_)
    pca_components_df.to_csv(os.path.join('outputs', 'pca_components.csv'))
    #### Bar plot of explained_variance
    plt.bar(
        [str(i) for i in range(0, len(pca_transformer.explained_variance_ratio_))],
        pca_transformer.explained_variance_ratio_
    )

    plt.plot(
        range(0, len(pca_transformer.explained_variance_ratio_)),
        np.cumsum(pca_transformer.explained_variance_ratio_),
        c='red',
        label='Cumulative Explained Variance')

    plt.legend(loc='upper left')
    plt.xlabel('Number of components')
    plt.ylabel('Explained variance')
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'pca_variance.jpg'), dpi=300)
    plt.close()

    #### Plot loadings
    import seaborn as sns

    def plot_loadings(_pca_df, var1, var2, pca_transformer):

        palette = None
        sns.scatterplot(data=_pca_df, x=var1, y=var2, palette=palette)

        useful_df = pd.DataFrame(pca_transformer.components_, columns=pca_transformer.feature_names_in_,
                                 index=_pca_df.columns)
        useful_df = useful_df.loc[[var1, var2]]

        scale_factor = 8  # Add scaling to make loadings visible
        for featr in useful_df.columns:
            x_val = useful_df[featr].loc[var1] * scale_factor
            y_val = useful_df[featr].loc[var2] * scale_factor
            plt.arrow(0, 0, x_val, y_val, color='r', alpha=0.5, shape='full', head_width=0.1, head_length=0.1)
            plt.text(x_val * 1.1, y_val * 1.1, featr, color='g', ha='center', va='center')

        plt.grid()
        plt.savefig(os.path.join('outputs', f'loadings_{var1}_{var2}.jpg'), dpi=600)
        plt.close()

    plot_loadings(pca_df, 'PC0', 'PC1', pca_transformer)
    plot_loadings(pca_df, 'PC0', 'PC2', pca_transformer)
    plot_loadings(pca_df, 'PC0', 'PC3', pca_transformer)
    plot_loadings(pca_df, 'PC1', 'PC2', pca_transformer)

    #### Matrix plot to show general groupings
    import plotly.express as px
    import kaleido # needed
    labels = {
        'PC' + str(i): f"PC {i} ({var:.1f}%)"
        for i, var in enumerate(pca_transformer.explained_variance_ratio_ * 100)
    }



    my_color_map = None
    fig = px.scatter_matrix(
        pca_df,
        labels=labels,
        dimensions=['PC0', 'PC1', 'PC2', 'PC3'],
        color_discrete_map=my_color_map
    )
    fig.update_traces(diagonal_visible=False, showupperhalf=False)
    fig.update_layout(legend=dict(title='Map', font=dict(size=20),
                                  yanchor="top",
                                  y=0.99,
                                  xanchor="left",
                                  x=0.75
                                  ))
    fig.write_image(f'outputs/pcas_matrix.jpg', scale=3)

    #### 3D plot to show general groupings
    # seaborn.set_style("whitegrid", {'axes.grid': False})
    #
    # fig = plt.figure(figsize=(6, 6))
    # ax = fig.add_subplot(111, projection='3d', )
    #
    # for s in pca_df[target_column].dropna().unique():
    #     ax.scatter(pca_df['PC0'][pca_df[target_column] == s], pca_df['PC1'][pca_df[target_column] == s], pca_df['PC2'][pca_df[target_column] == s],
    #                label=s)
    #
    # ax.legend()
    # plt.tight_layout()
    # plt.savefig(os.path.join('outputs', f'3D_pca_{target_column}.jpg'), dpi=300)
    # plt.close()
