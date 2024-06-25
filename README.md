# nikita_soft_spot

Required library/package
You need to first install: python, numpy, f2py, scipy ... and f90wrap

If you are using anaconda:
conda install -c conda-forge f90wrap
otherwsise: (https://github.com/jameskermode/f90wrap)
pip install f90wrap

You also need to install the softspot library for Pseudo harmonic modes
via the command: pip install softspot


go to ./soft_spot_project and do make


add the project to your python path in your bashrc:
export PYTHONPATH="...your path.../nikita_soft_spot/:$PYTHONPATH"



