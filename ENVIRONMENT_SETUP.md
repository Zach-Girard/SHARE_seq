## Environment setup (without Homebrew)

This document describes how to install Miniconda **without Homebrew** and create the master Conda environment for this project.

### 1. Install Miniconda

1. Open a terminal.
2. Download the latest Miniconda installer for macOS (Apple Silicon). Adjust the URL if you are on Intel.

```bash
cd ~
curl -fsSLo miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
```

3. Run the installer and follow the prompts:

```bash
bash miniconda.sh
```

Recommended choices:

- Accept the license.
- Accept the default install location.
- Allow the installer to initialize your shell (so `conda` is on your PATH).

4. When the installer finishes, **close and reopen** your terminal (or `exec zsh` / `source ~/.zshrc`) so that `conda` is available.

To verify:

```bash
conda --version
```

### 2. Create the project environment

From the project root:

```bash
cd /path/to/SHARE_seq
conda env create -f environment.yml
```

This creates an environment named `nextflow-master` (as defined in `environment.yml`).

### 3. Activate the environment

Each time you work on this project:

```bash
conda activate nextflow-master
```

You now have access to the tools defined in `environment.yml` (e.g. `nextflow`, `curl`, `pigz`, etc.).

