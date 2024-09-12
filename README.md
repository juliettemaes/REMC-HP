# A replica exchange Monte Carlo algorithm for protein folding in the HP model

Juliette Maes
Sept. 2024

Description of this project. Put the link of the oroginal paper, from which the work is inspired. 

## Download this repository

```bash
git clone https://github.com/juliettemaes/REMC-HP.git
cd REMC-HP
```
## Install using Conda environment

Install [conda](https://docs.conda.io/en/latest/miniconda.html).

Create conda environment and install dependendies:

```bash
cd REMC-HP
conda env create -f environment.yml
```

Load conda environment:

```bash
conda activate protein-folding
```

Install dependencies:

```bash
pip install -r requirement.txt
```

## Install using venv environment

Create virtual environment:

```bash
cd REMC-HP
py -m venv protein-folding
```

Load virtual environment:

```bash
protein-folding/Scripts/activate
```

Install dependencies:

```bash
pip install -r requirement.txt
```


## Running the repo

First to run the REMC Search you have to run the main.py file found in the src folder.
That file takes two arguments, optimal_energy and hpsequence OR aasequence.

Here is an example of a command line to run the project:  

```bash
python src/run.py --hpsequence="HPHPPHHPHPPHPHHPPHPH" --optimal_energy=-9
```

The way I designed the arguments, as I said before, makes the following line not valid : 
```bash
python src/run.py --hpsequence="HPHPPHHPHPPHPHHPPHPH" --aasequence="HPHPPHHPHPPHPHHPPHPH" --optimal_energy=-9
```

