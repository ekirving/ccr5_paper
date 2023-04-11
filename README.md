# Allele frequency trajectory of the CCR5delta32 mutation
This repository contains code for the `CLUES` selection analysis from 
[The evolutionary history of the CCR5delta32 deletion as revealed by ancient and present-day genomes](https://).

![Figure 2](./figure/ccr5_delta32_trajectory.png?raw=true)

If you reuse any of this code then please cite the preprint:
> Kirstine Ravn&ast;, Leonardo Cobuccio&ast;, Rasa Audange Muktupavela&ast;, Jonas Meisner, Martin Sikora, Eske 
> Willerslev, Morten Allentoft, Evan K. Irving-Pease, Fernando Racimo, Simon Rasmussen, 2022. The evolutionary history 
> of the CCR5delta32 deletion as revealed by ancient and present-day genomes. *bioRxiv* XXX

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
./run_ccr5_analysis.sh
```

```bash
# plot the composite figure  
Rscript scripts/ccr5_plot_models.R
```

## Author

Evan K. Irving-Pease, [GLOBE Institute](https://globe.ku.dk/), University of Copenhagen 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
