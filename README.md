# Allele frequency trajectory of the CCR5delta32 mutation
This repository contains code for the `CLUES` selection analysis from 
[Tracing the evolutionary path of the CCR5delta32 deletion via ancient and modern genomes](https://doi.org/10.1101/2023.06.15.23290026).

![Figure 5](./figure/ccr5_delta32_trajectory.png?raw=true)

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

## Running the code

This project contains rules for running selection tests using `CLUES` from the outputs of the different `HAPI` models.

```bash
# run all the models  
./run_clues_analysis.sh
```

```bash
# infer the ages of the mutations
Rscript scripts/infer_mutation_ages.R
```

```bash
# plot the composite figure 
Rscript scripts/ccr5_plot_models.R
```

## Author

Evan K. Irving-Pease, [GLOBE Institute](https://globe.ku.dk/), University of Copenhagen 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
