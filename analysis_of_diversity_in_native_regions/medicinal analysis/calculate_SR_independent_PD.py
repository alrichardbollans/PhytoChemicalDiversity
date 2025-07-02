import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

from helper_functions import get_working_data_for_medicinal_species_analysis

def find_regression_model():
    y_var = 'PD'
    x_var = 'SR'
    # Step 1: Fit a linear regression model
    X = working_data[[x_var]].values  # Independent variable
    X_to_plot = working_data[x_var].values  # Independent variable (species richness)
    # scaled_data[metric] = np.log(scaled_data[metric])
    y = working_data[y_var].values  # Dependent variable (diversity)

    # fit regression model
    model = sm.OLS(y, sm.add_constant(X)).fit()

    r_squared = model.rsquared
    best_model_r_sqaured = 0
    data = [['Linear', r_squared, model.aic]]
    sns.scatterplot(x=X_to_plot, y=y, edgecolor="black", alpha=0.8)
    linear_prediction = model.predict(sm.add_constant(X))
    if r_squared > best_model_r_sqaured:
        best_model_r_sqaured = r_squared
    sns.lineplot(x=X_to_plot, y=linear_prediction, color='black', linestyle='--')
    plt.savefig(os.path.join('outputs', 'PD_SR_regression', 'linear_regression.jpg'), dpi=300)
    plt.close()

    # Fit LOESS with outlier robustness (iterations downweight outliers)
    loess_prediction = sm.nonparametric.lowess(exog=X_to_plot, endog=y, return_sorted=False)
    expected_diversity = loess_prediction
    residuals = y - expected_diversity
    # Calculate R² (coefficient of determination)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f"R² (Coefficient of Determination): {r_squared:.3f}")
    if r_squared > best_model_r_sqaured:
        best_model_r_sqaured = r_squared
    sns.scatterplot(x=X_to_plot, y=y, edgecolor="black", alpha=0.8)
    sns.lineplot(x=X_to_plot, y=expected_diversity, color='black', linestyle='--')
    plt.savefig(os.path.join('outputs', 'PD_SR_regression', 'LOESS.jpg'), dpi=300)
    plt.close()
    data.append([f'LOESS', r_squared, '?'])

    for deg in range(2, 8):
        poly = PolynomialFeatures(degree=deg)
        X_poly = poly.fit_transform(X)

        # X_poly has constant added
        model = sm.OLS(y, X_poly).fit()
        r_squared = model.rsquared
        data.append([f'Polynomial {deg}', r_squared, model.aic])
        sns.scatterplot(x=X_to_plot, y=y, edgecolor="black", alpha=0.8)
        expected_diversity = model.predict(X_poly)
        if deg == 6:
            # This has almost best R2 and very good AIC
            best_model_prediction = expected_diversity
        if r_squared > best_model_r_sqaured:
            best_model_r_sqaured = r_squared
            print(f'{deg} has best R2')
        # sns.lineplot(x=X_to_plot, y=linear_prediction, color='black', linestyle='--', label='Linear')
        sns.lineplot(x=X_to_plot, y=expected_diversity, color='black', linewidth=2.5, linestyle='--', label=f'Polynomial degree: {deg}')
        plt.xlabel('Species Richness')
        plt.ylabel('Phylogenetic Diversity')
        plt.legend()
        plt.savefig(os.path.join('outputs', 'PD_SR_regression', f'poly_{deg}_regression.jpg'), dpi=300)
        plt.close()

    df = pd.DataFrame(data, columns=['Model', 'R-squared', 'AIC'])
    df.to_csv(os.path.join('outputs', 'PD_SR_regression', 'model_comparison.csv'))
    return best_model_prediction


def main():
    best_model_prediction = find_regression_model()

    # Step 2: Predict the expected diversity based on species richness
    working_data['predicted PD'] = best_model_prediction

    # Step 3: Calculate the residuals (observed - expected)
    working_data["PD'"] = working_data['PD'] - working_data['predicted PD']
    working_data.to_csv(os.path.join('outputs', 'PD_SR_regression', 'SR-Independent PD.csv'))


if __name__ == '__main__':
    working_data = get_working_data_for_medicinal_species_analysis()
    working_data = working_data.dropna(subset=['PD', 'SR'], how='any')
    main()
