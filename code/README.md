
#### Remarks：
* This tool passed experiment on the "Ubuntu Linux 20.04.3 LTS (Focal Fossa".
* python version : 3.7.0

## Provided Resources:
Within this open-source project, we offer essential resources for experimentation:

- **Data Folder**: Contains the data files needed to perform the experiments and includes a README file that explains in detail what the individual data means.
- **categories Folder**:This folder contains the results of our classification of the faults in the dataset according to the definition of faults in the paper.
- **Code Folder**: Includes the source code implementing our strategy.
- **HOMs Folder**: Includes a script for generating high-order mutants and a high-order mutants dataset generated using this script.
- **baseline Folder**: Contains code to implement the baseline.
These resources are made available to facilitate reference and reproducibility of the experimental process. Developers can clone the repository to access the code and data required for utilizing our methods in fault localization experiments.

The specific execution process is as follows:
## Code Files
In general,after using the high-order mutant script to generate a high-order mutant, perform mutation testing, obtain mutation execution information, and execute the above script files in sequence to complete a execution process.
- **clacNew.py**: This code reads the original test, mutant execution results and some auxiliary information, calculates and sorts the suspicion of each row (or mutant) under multiple spectral methods (such as Ochiai, Dstar, Jaccard, etc.), and then outputs the results to a CSV file, thus realizing the function of calculating the suspicion of different mutants in the project.
- Command：cd to path/to/your/path,then run python clacNew.py

- **clacRank.py**: It reads the error line information from the file obtained in the previous step and converts this information into the format required for internal processing. Then, the code traverses the specified project version range, and for each version, reads the suspicion data stored in the CSV file according to various suspicion evaluation formulas (such as Dstar, Ochiai, etc.). Next, it calculates the suspicion ranking of each line of code and formats these rankings and other related information into JSON files and saves them in a specific directory structure.
- Command：cd to path/to/your/path,then run python clacRank.py
  
- **clacAlITypeRank.py**: This code implements the project suspicion assessment data (stored in JSON format) obtained from the previous step, analyzes this data, and calculates the average ranking of each error code line. It traverses different versions of multiple software projects, applies multiple suspicion assessment formulas to each version, processes and summarizes the results, and outputs these processed ranking data to a new JSON file.
- Command：cd to path/to/your/path,then run python clacAlITypeRank.py

- **clacCorrectTopn.py**: This code mainly implements reading the JSON data files generated by various suspicion formulas (such as "Dstar", "Ochiai", etc.) in a specific project (such as "Chart") from the JSON file obtained in the previous step, and calculating the distribution of fault code lines in the Top-1, Top-3, Top-5, and Top-10 suspicion rankings. Then, it saves these statistical results in JSON format to the specified directory. This process aims to evaluate the impact of different suspicion formulas on fault location efficiency and provide data support for further analysis.
- Command：cd to path/to/your/path,then run python clacCorrectTopn.py
  
- **clacAllTop-n.py**: This code mainly implements reading pre-calculated suspicion data (such as Top-1, Top-3, Top-5, Top-10 suspicion rankings) from specified directories of different projects (such as "Chart", "Time", "Lang", etc.), and evaluates and compares the effects of different formulas (such as "Dstar", "Ochiai", etc.) by merging and analyzing these data. Finally, it merges and outputs the summary results of all projects into a DataFrame for further analysis or reporting. In addition, the code exports the results to an Excel file at the final stage for more in-depth data analysis and presentation.
- Command：cd path/to/your/path,then run python clacAllTop-n.py
  
- **baseline.py**: Execute the script to get the baseline results.
- Command：cd to path/to/your/path,then run python baseline.py

## results
All our experimental results are saved in the data folder. 
There is a README file that explains each file in detail.


