# CavSim

The lid-driven cavity simulator

## Development Environment

This project is based in a `conda` environment using `conda-devenv` and `conda-lock` tools. For all basic tools instalation:
- `conda` basic package: Dor `miniconda` installation follow the instructions [here](https://conda.io/miniconda.html)
- `conda-devenv` package: Installation instructions click [here](https://conda-devenv.readthedocs.io/en/latest/installation.html) and 
- `conda-lock` package: Installation instructions click [here](https://conda.github.io/conda-lock/). 
It is important to note that `conda-devenv` and `conda-lock` MUST be installed in the `base` conda environment to allow the usage in this project.

### Setup and Environment Activation

To create and activate the conda environment called `cavsim-dev` execute the following commands into project root folder
```shell
$ conda devenv
$ conda activate cavsim-dev
```

### Updating Environment

If you add or remove dependencies, or change version requirements in the `environment.devenv.yml` file, you can update the locks with:
```shell
$ conda devenv --lock
```
This will attempt to keep the existing pins as much as possible.

To update to the latest version of one or more libraries, for all platforms:
```shell
$ conda devenv --update-locks pytest --update-locks boltons
```

If you want to update all libraries to the latest version, pass an "" (empty) value:
```shell
$ conda devenv --update-locks ""
```
All these instructions are available in the `conda-devenv` [documentation](https://conda-devenv.readthedocs.io/en/latest/usage.html#locking)
