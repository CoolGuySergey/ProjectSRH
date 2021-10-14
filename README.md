# SRHClusterMapper: an analytics tool for the detection of assumption violations prior to tree-building

SRHClusterMapper is a web(**on-hold)/CLI tool made with the aim to clearly visualise and group sequences in an alignment into subsets where all within-subset seqeunce-pairs meet model-assumptions. The symmetry tests enable a divide-and-conquer protocol that saves resources comparing sequences that should not participate in the same phylogenetic analysis. 

![SRHClusterMapper screenshot](/images/gallery1.png)

DNA base substitution is far more complex than most popular models assume. Even the most parameter-rich of the models make the assumptions of stationary, reversible, and homogenous conditions:

I. Stationarity: The probability of obtaining a given nucleotide remains the same between different sequences

II. Reversibility: Requires stationarity. When sampling a nucleotide from the stationary distribution i.e. stationarity evolving composition, the probability of its substitution by another nucleotide is the same if sampling had happened vice versa.

III. Homogeneity: Constant rates of change assumed to cover all edges in a tree.

Or, more intuitively: when models assume SRH conditions, it is assuming that when two sequences diverge, changes accumulate symmetrically in both sequences. For a given alignment, SRHClusterMapper computes divergence matrices (shown in example below) for every sequence pair in the alignment. This is a 4x4 matrix where off-diagonal elements document base changes and forms the basis of all three symmetry tests.

![SRHClusterMapper screenshot](/images/gallery0.png)

There are two ways to use SRHClusterMapper:
1. Online at [https://srhclustermapper.com](https://srhclustermapper.com/). **on-hold

![SRHClusterMapper screenshot](/images/gallery3.png)

2. Run it from the command-line. For this only the contents of the `scripts/` directory are needed.
```
python scripts/SRHClustermapper.py --input <in_file>
```

To switch off partitioning:
```
python scripts/SRHClustermapper.py --input <in_file> --partition 0
```

To switch off correction for multiple comparisons in case of alignments with few entreis:
```
python scripts/SRHClustermapper.py --input <in_file> --alpha 0.05
```

Important:
1. The three tests and the three conditions (SRH) do not map one-to-one (see screenshot of landing page above).
2. Invalid pairs occur when a particular pair of sequences are identical or do not have overlapping non-gap sites in theM (off-diagonals sum to 0 in divergence matrix m).

REQUIREMENTS
===
* [Python 3](https://www.python.org/downloads/)
* [NumPy](http://www.numpy.org/) `pip3 install numpy`
* [SciPy](http://scipy.org/) `pip3 install scipy`
* [Seaborn](https://seaborn.pydata.org/) `pip3 install seaborn`
* [Matplotlib](https://matplotlib.org/) `pip3 install matplotlib`