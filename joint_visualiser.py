import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--counts',
                        help="file that contains transition counts.", nargs='+', type=str)
    parser.add_argument('-d', '--dates',
                        help="file that contains transition dates (same order as counts).", nargs='+', type=str)
    parser.add_argument('-f', '--families',
                        help="family names (same order as counts).", nargs='+', type=str)
    parser.add_argument('--violin',
                        help="pdf with the violin plot.", required=True, type=str)
    parser.add_argument('--scatter',
                        help="pdf with the scatter plot.", required=True, type=str)
    parser.add_argument('--heatmap',
                        help="pdf with the heatmap plot.", required=True, type=str)
    params = parser.parse_args()

    joint_ddf, joint_df = None, None

    for date_file, family in zip(params.dates, params.families):
        df = pd.read_csv(date_file, header=0, sep='\t')
        df['transition'] = df['from state'] + '>' + df['to state']
        df.sort_values(by=['transition'], inplace=True)
        df = df[df['probability'] == 1]
        df['family'] = family

        ddf = df[['transition', 'from date', 'family']].append(df[['transition', 'to date', 'family']]
                                                  .rename(columns=lambda _: _ if _ != 'to date' else 'from date'))
        ddf.columns = ['transition', 'date (Ma)', 'family']
        ddf['type'] = ['transition start'] * len(df) + ['transition end'] * len(df)

        if joint_df is None:
            joint_df = df[['transition', 'from date', 'family']]
            joint_ddf = ddf
        else:
            joint_df = joint_df.append(df[['transition', 'from date', 'family']])
            joint_ddf = joint_ddf.append(ddf)
    joint_df.columns = ['transition', 'start date (Ma)', 'family']

    rc = {'font.size': 10, 'axes.labelsize': 8, 'legend.fontsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 8,
          'ytick.labelsize': 4}
    sns.set(style="whitegrid")
    sns.set(rc=rc)
    ax = sns.scatterplot(y="transition", x="start date (Ma)", data=joint_df, hue='family', alpha=.75,
                         palette="colorblind")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(params.scatter, dpi=300)

    plt.clf()
    rc = {'font.size': 12, 'axes.labelsize': 10, 'legend.fontsize': 10, 'axes.titlesize': 10, 'xtick.labelsize': 10,
          'ytick.labelsize': 10}
    sns.set(style="whitegrid")
    sns.set(rc=rc)
    ax = sns.violinplot(y='family', x="date (Ma)", data=joint_ddf, inner='stick', hue='type', split=True, cut=0,
                        palette="colorblind")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(params.violin, dpi=300)

    plt.clf()
    rc = {'font.size': 12, 'axes.labelsize': 12, 'legend.fontsize': 12, 'axes.titlesize': 12, 'xtick.labelsize': 12,
          'ytick.labelsize': 12}
    sns.set(style="whitegrid")
    sns.set(rc=rc)
    fig, axs = plt.subplots(2, 3, figsize=(40, 20))
    row_col = ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2))
    for count_file, family, (row, col) in zip(params.counts, params.families, row_col):
        df = pd.read_csv(count_file, header=0, sep='\t')[["from state", "to state", "count"]]
        df.columns = ["from", "to", "count"]
        df = df.pivot("from", "to", "count")
        sns.heatmap(df, annot=True, cmap="YlGnBu", ax=axs[row, col])
        axs[row, col].set_title(family)
    if len(row_col) > len(axs):
        for row, col in row_col[len(params.families):]:
            axs[row, col].set_axis_off()
    plt.tight_layout()
    fig.savefig(params.heatmap, dpi=300)

