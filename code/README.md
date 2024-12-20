
#### Remarks：
* This tool passed experiment on the "Ubuntu Linux Intel i5 Xeon e7-480".
* python version : 3.6.8

## Provided Resources:
Within this open-source project, we offer essential resources for experimentation:

- **Data Folder**: Contains the necessary data files for conducting experiments.
- **categories Folder**:Contains the results of the classification of faults
- **Code Folder**: Includes the source code implementing our strategy.
- **HOMs Folder**: Contains code for generating high-order mutants.
- **baseline Folder**: Contains code to implement the baseline.
These resources are made available to facilitate reference and reproducibility of the experimental process. Developers can clone the repository to access the code and data required for utilizing our methods in fault localization experiments.

The specific execution process is as follows:
## Code Files

- **clacNew.py**: Calculate the suspicion of each entity.

- **clacRank.py**: Calculate the gradient ranking of each entity.
  
- **clacAlITypeRank.py**: Calculate the real ranking of each entity.

- **clacAlITypeRank.py**: Calculate topn.

- **clacAllTop-n.py**: Integrate all top-n.
  
- **baseline.py**: Execute the script to get the baseline results.

Execute the above files in sequence to obtain the results.

## results
All our experimental results are saved in the data folder. 
There is a README file that explains each file in detail.


