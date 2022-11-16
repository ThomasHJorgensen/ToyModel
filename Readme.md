# Readme

This readme describes how to install everything required for running the `Solve Toy Model.ipynb` notebook before describing how to execute the code and how the code base is structured.

# Installation
1. Install Python: E.g. the free [Anaconda](https://www.anaconda.com/products/distribution)
2. Install packages
    1. open Anaconda Prompt (e.g. by searching)
    2. type `pip install EconModel consav`

These two packages are developed and maintained by Jeppe Druedahl for the purpose of easy implementaion of efficient (Numba and C/C++) code to solve consumption-savings models. See [consav package](https://github.com/NumEconCopenhagen/ConsumptionSaving) and [EconModel package](https://github.com/NumEconCopenhagen/EconModel) for more information.

# Cloning the code
1. open Anaconda Prompt
2. go to the location on your computer where youwant the code to be placed. This can be done by typing `cd /d PATHTOYOURFOLDER` in the prompt.
3. type `git clone https://github.com/ThomasHJorgensen/ToyModelDivorce.git`. Now the code is dowloaded into your local machine and you can do with it what you want without chaning the files in the GitHub repository.

# Compilation
The code is implemented in C++. This means that the code should be re-compiled if anything is changed. The `Solve Toy Model.ipynb` notebook does this automatically if requested (`do_compile=True`). To compile the files, Visual Studio 2022 should be installed in the path `C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Auxiliary/Build/`.

The Visual Studio compiler is free and can be downloaded for free here: [Visual Studio](https://visualstudio.microsoft.com/vs/whatsnew/).

# Execution
The code is executed from a **Jupyter Lab** notebook called `Solve Toy Model.ipynb`. This notebook executes code and illustrates results.

1. open jupyter lab by executing the `open_jupyter.bat` file. This opens up the browser. Chrome is preferred, which you can dowload [here](https://www.google.com/chrome/).
2. run the notebook. To execute the entire notebook you can press the the play button in the top. To execute individual cells (and advance to the next cell) press `shift+enter`. You can read more about the Jupyter Lab interface [here](https://jupyterlab.readthedocs.io/en/stable/user/interface.html).

# Code
The main Python code is contained in **ToyModel.py** which calls C++ code. 

**Python:** In Python the model instance (object) is initiated, memory is allocated and parameters and grids are set. Building on the `EconModel`, Python also handles all convertions to C++ that we will utilize below. This means that all parameters and grids should be set in `ToyModel.py`. The main C++ file that Python should know about is set in `cpp_filename = 'cppfuncs/solve.cpp'`.

**C++:** All the code solving and simulating the model is implemented in C++ for speed. The `.cpp` files are all stored in the `cppfuncs` folder, which the `EconModel` class assumes. All functions that should be callable from Python must have the `EXPORT` decorator in front of the function definition. All functions will then be callable from Python using `model.cpp.functionname`. See e.g. the function `solve` in `ToyModel.py`. In this example, calling `model.solve()` will automatically call the C++ implementaion. See also [EconModelNotebooks](https://github.com/NumEconCopenhagen/EconModelNotebooks) for illustrations of the link to C++.
