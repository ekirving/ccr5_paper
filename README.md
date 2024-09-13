# Allele frequency trajectory of the CCR5Δ32 deletion
This repository contains code for the `CLUES` selection analysis from 
[Tracing the evolutionary path of the CCR5delta32 deletion via ancient and modern genomes](https://doi.org/10.1101/2023.06.15.23290026).

![HAPI trajectories](/figure/hapi_trajectories.png?raw=true)
![Ancestral path trajectories](/figure/ancestry_trajectories.png?raw=true)

If you reuse any of this code then please cite the preprint:
> Ravn, K., Cobuccio, L., Muktupavela, R.A., Meisner, J.M., Benros, M.E., Korneliussen, T.S., Sikora, M., Willerslev, 
> E., Allentoft, M.E., Irving-Pease, E.K., Racimo, F., Rasmussen, S., 2023. Tracing the evolutionary path of the 
> CCR5delta32 deletion via ancient and modern genomes. *medRxiv* https://doi.org/10.1101/2023.06.15.23290026


## Installation
Download the code: 
```bash
git clone git@github.com:ekirving/ccr5_paper.git && cd ccr5_paper/
```

The easiest way to install all the dependencies is with the [conda package manager](https://docs.conda.io/en/latest/).

```bash
conda env create --name ccr5 --file environment.yaml
```

Then activate the environment:
```bash
conda activate ccr5
```

## Running the HAPI analyses

These commands run the `CLUES` selection tests using the deletion calls from the `HAPI` models.

```bash
# run all the CCR5 and control models
./run_hapi_analysis.sh
```

```bash
# infer the ages of the alleles
Rscript scripts/infer_allele_ages.R
```

```bash
# plot the composite figure
Rscript scripts/plot_hapi_models.R
```

## Running the ancestry path analyses

These commands run the `CLUES` selection tests using the ancestry stratified tag SNPs for the CCR5Δ32 deletion.

```bash
# run all the CCR5 and control models
./run_ancestry_analysis.sh
```

```bash
# plot the composite figure
Rscript scripts/plot_ancestry_models.R
```

## Author

Evan K. Irving-Pease, [GLOBE Institute](https://globe.ku.dk/), University of Copenhagen 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
