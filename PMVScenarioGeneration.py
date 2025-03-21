import math
from pythermalcomfort import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
from matplotlib.cm import get_cmap
from scipy.stats import gaussian_kde
import matplotlib.ticker as ticker


def loadAshraedbII(loc_measurements, loc_meta):
    df = pd.read_csv(loc_measurements, sep=',')
    metadata_df = pd.read_csv(loc_meta)
    return df, metadata_df


def filterDataframe(df,metadata_df):
    # Included and Excluded types
    excluded_countries = ['italy', 'portugal', 'greece']
    included_building_types = ['multifamily housing']
    selected_df = metadata_df[
        (metadata_df['region'] == 'europe') &
        (~metadata_df['country'].isin(excluded_countries))
        ]

    selected_building_ids = selected_df['building_id'].unique()
    df = df[df['building_id'].isin(selected_building_ids)]

    df_acm = df.loc[(~df['tr'].isna()) &
                    (~df['vel'].isna()) &
                    (~df['rh'].isna()) &
                    (~df['clo'].isna()) &
                    (~df['met'].isna())].copy()

    label_encoder = LabelEncoder()
    df_acm['Cluster'] = label_encoder.fit_transform(df_acm['season'])
    df_acm = df_acm[df_acm['season'] == 'winter']
    return df_acm


def generateSamples(df, n_sim, variables = ['tr', 'vel', 'rh', 'clo', 'met', 'thermal_sensation']):
    samples_list = []

    for var in variables:
        data = df[var].dropna()  # Drop NaN values for the variable

        kde = gaussian_kde(data)  # Estimate the probability density
        sampled_data = kde.resample(n_sim).flatten()

        # Apply the constraint for velocity
        if var == 'vel':
            sampled_data = sampled_data[sampled_data >= 0]
            while len(sampled_data) < n_sim:  # Re-sample until we have n_sim valid samples
                additional_samples = kde.resample(n_sim - len(sampled_data)).flatten()
                sampled_data = np.append(sampled_data, additional_samples[additional_samples >= 0])

        samples_list.append(sampled_data[:n_sim])
    return samples_list


def calculatePPD(ta_range, samples_list,n_sim, savetoCSV=False):
    tr_sample = samples_list[0]
    vel_sample = samples_list[1]
    rh_sample = samples_list[2]
    clo_sample = samples_list[3]

    met = 1.2  # met_sample[i]
    output_ta = []
    mean_PPD = []
    median_PPD = []
    for n in range(0, len(ta_range)):
        ta = ta_range[n]
        output = []
        for i in range(0, n_sim):
            clo = clo_sample[i]
            rh = rh_sample[i]
            vel = vel_sample[i]
            tr = tr_sample[i]
            output.append(pmv_ppd(ta, tr, vel, rh, met, clo)['ppd'].item())  #
        mean_PPD.append(np.nanmean(output))
        median_PPD.append(np.nanmedian(output))
        output_ta.append(output)

    ta_values = []  # To store repeated Ta values
    ppd_values = []  # To store corresponding PPD values

    for i, ta in enumerate(ta_range):
        ta_values.extend([ta] * len(output_ta[i]))  # Repeat Ta value for each PPD value
        ppd_values.extend(output_ta[i])  # Add the PPD values

    # Create a DataFrame
    df = pd.DataFrame({
        "Ta": ta_values,  # Air temperature
        "PPD": ppd_values  # Predicted percentage of dissatisfied (PPD)
    })
    median_df = pd.DataFrame({'Ta': ta_range, 'median PPD': median_PPD})
    if savetoCSV:
        median_df.to_csv('medianPPD_ta.csv')
        df.to_csv('PPD_ta.csv')
    return df, median_df,median_PPD


def plotPMVrange(df,median_PPD, ta_range,saveFigure=False):
    plt.rcParams.update({
        "font.size": 9,
        "font.family": "Times New Roman"
    })
    cmap = get_cmap('coolwarm')

    ta_plot = [item - 15 for item in ta_range]

    plt.figure(figsize=(3.5, 3.5 * 0.5))

    # Create the boxplot
    sns.boxplot(x="Ta", y="PPD", data=df, showfliers=False, color=cmap(0.1), label="Data")

    # Add the lineplot
    sns.lineplot(
        x=ta_plot,
        y=median_PPD,
        color='orange',
        marker="x",
        label="Median",
        markersize=4,
        markeredgewidth=1,
        markeredgecolor="orange"
    )

    # Adjust x-axis ticks to show every third whisker
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(3))

    # Add legend, labels, and title
    plt.legend(loc="upper right", frameon=False, borderaxespad=0.25)
    plt.xlabel("Indoor Temperature (°C)")
    plt.ylabel("PPD (%)")
    plt.ylim([0, 80])
    plt.subplots_adjust(top=0.95, bottom=0.25, left=0.2, right=0.99)
    # Show the plot
    if saveFigure:
        plt.savefig('pmv_plot.svg')

    plt.show()


def listCoefficients(median_PPD,ta_range):
    PPDcoef = []
    for n in range(0,len(median_PPD)-1):
        Dy = median_PPD[n+1] - median_PPD[n]
        Dx = 1
        m = Dy/Dx
        b = median_PPD[n] - m * ta_range[n]
        coeff = [m,b]
        PPDcoef.append(coeff)
    return PPDcoef


def ExampleSimulation(fileloc,metadataloc):
    df, metadata_df = loadAshraedbII(fileloc, metadataloc)
    df = filterDataframe(df,metadata_df)

    n_sim = 10000
    samples_list = generateSamples(df,n_sim)

    ta_range = range(15,31)
    df, median_df,median_PPD = calculatePPD(ta_range,samples_list,n_sim)
    plotPMVrange(df,median_PPD, ta_range)

if __name__ == "__main__":
    fileloc =           # path to measurements
    metadataloc =       # path to metadata
    ExampleSimulation(fileloc,metadataloc)




