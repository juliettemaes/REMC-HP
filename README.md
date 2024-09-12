# A replica exchange Monte Carlo algorithm for protein folding in the HP model

Juliette Maes
Sept. 2024

Find the original article, from which this work is inspired, [here](https://doi.org/10.1186/1471-2105-8-342).

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

To run the REMC Search on your sequence, run the main.py file found in the src folder.
This script takes two arguments, optimal_energy and hpsequence OR aasequence.

Here is an example of a command line to run the project:  

```bash
python src/main.py --hpsequence="HPHPPHHPHPPHPHHPPHPH" --optimal_energy=-9
```

Here is an example of a INVALID command line :
```bash
python src/main.py --hpsequence="HPHPPHHPHPPHPHHPPHPH" --aasequence="AMGHICVFGEDGLKILDGEA" --optimal_energy=-9
```

