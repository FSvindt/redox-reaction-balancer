# Redox Reaction Balancer

A redox reaction balancer written in Python.

The program is able to balance simple redox reactions.

To ensure correct chemical terminology, the GUI language is Danish.

## Running the Program

### Too Lazy to Install...

The app is deployed to Streamlit: https://share.streamlit.io/fsvindt/redox-reaction-balancer/redox-only/main.py

### Requirements & Installation

The program requires the following Python modules:

- Python v3.9
- Pip v21.2.4
- SymPy v1.10.1
- ChemParse v0.1.1
- Streamlit v1.8.1

It is preferred but not necessary to use a conda virtual environment.

To create a new conda environment with the necessary packages installed, run `conda env create -f environment.yml` in the project directory after cloning or downloading the project.

A generated requirements.txt file is also included in the project.

### Run using Streamlit

Make sure all the requirements have been satisfied.

Open the terminal and navigate to the location of the program.

Run the program using Streamlit by executing the command `streamlit run main.py`.

## Credits

Parts of the program are based on the following articles published by Medium in Towards Data Science and The Startup.

- How to Balance Chemical Equations in Python using Constraint Optimization (PuLP) by Rahul Banerjee: <https://towardsdatascience.com/how-to-balance-chemical-equations-in-python-using-constraint-optimization-pulp-1d7409fbe52b>

- Balancing Chemical Equations With Python by Mohammad-Ali Bandzar: <https://medium.com/swlh/balancing-chemical-equations-with-python-837518c9075b>
