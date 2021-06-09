import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--counts',
                        help="file that contains transition counts.", required=True, type=str)
    parser.add_argument('-d', '--dates',
                        help="file that contains transition dates.", required=True, type=str)
    parser.add_argument('--fig',
                        help="pdf with visualisation.", required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.dates, header=0, sep='\t')
    df['transition'] = df['from state'] + '>' + df['to state']
    df.sort_values(by=['transition'], inplace=True)
    df = df[df['probability'] == 1]

    ddf = df[['transition', 'from date']].append(df[['transition', 'to date']]
                                              .rename(columns=lambda _: _ if _ != 'to date' else 'from date'))
    ddf.columns = ['transition', 'date (Ma)']
    ddf['type'] = ['transition start'] * len(df) + ['transition end'] * len(df)
    ddf[' '] = ''

    plt.rcParams['font.size'] = 10
    fig, axs = plt.subplots(3, 1, figsize=(15, 20), gridspec_kw={'height_ratios': [2, 1, 1]})
    axs[0].get_shared_x_axes().join(axs[0], axs[1])
    sns.set(style="whitegrid")
    sns.scatterplot(y="transition", x="date (Ma)", data=ddf, hue='type', ax=axs[0], alpha=.5, palette="colorblind")
    sns.violinplot(y=' ', x="date (Ma)", data=ddf, inner='stick', hue='type', split=True, cut=0, ax=axs[1], palette="colorblind")
    for ax in axs[:2]:
        plt.setp(ax.get_legend().get_texts(), fontsize=10)
        plt.setp(ax.get_legend().get_title(), fontsize=10)
    df = pd.read_csv(params.counts, header=0, sep='\t')[["from state", "to state", "count"]]
    df.columns = ["from", "to", "count"]
    df = df.pivot("from", "to", "count")
    sns.heatmap(df, annot=True, cmap="YlGnBu", ax=axs[2])
    plt.tight_layout()
    fig.savefig(params.fig, dpi=300)

