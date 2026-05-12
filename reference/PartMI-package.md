# PartMI: Part Mutual Information for Direct Association Network Inference

Implements Part Mutual Information (PMI) for inferring direct
associations in biological networks. PMI extends standard Mutual
Information (MI) by conditioning on background variables, removing
indirect dependencies and reducing false-positive edges in network
reconstruction.

## Core functions

- [`mi`](mi.md) - Compute mutual information between two variables

- [`pmi`](pmi.md) - Compute part (conditional) mutual information

- [`pmi_network`](pmi_network.md) - Infer direct association network

- [`normalize_data`](normalize_data.md) - Normalize expression data

- [`discretize_data`](discretize_data.md) - Discretize continuous data

## References

Zhao et al. (2016) "Part mutual information for quantifying direct
associations in networks." PNAS 113(18): 5130-5135.

Frenzel & Pompe (2007) "Partial Mutual Information for Coupling Analysis
of Multivariate Time Series." Physical Review Letters 99(20): 204101.

Kraskov, Stögbauer & Grassberger (2004) "Estimating mutual information."
Physical Review E 69(6): 066138.

## See also

Useful links:

- <https://github.com/danymukesha/PartMI>

- Report bugs at <https://github.com/danymukesha/PartMI/issues>

## Author

**Maintainer**: Dany Mukesha <danymukesha@gmail.com>
([ORCID](https://orcid.org/0009-0001-9514-751X))

Other contributors:

- Zhao et al. (Original PMI method) \[contributor\]
