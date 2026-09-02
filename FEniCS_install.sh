# curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash ~/Miniconda3-latest-Linux-x86_64.sh
conda create -n fenics

conda activate fenics
conda install -c conda-forge fenics-dolfinx mpich pyvista
