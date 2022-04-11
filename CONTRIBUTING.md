# Contributing to 2QAN

Contributions are welcome, and they are greatly appreciated, every little bit helps.

## Opening an issue
You can begin contributing to `2QAN` code by raising an
[issue](https://github.com/lllingoo/2QAN/issues/new), reporting a bug or
proposing a new feature request, using the labels to organize it.
Please use `2qan.about()` to document your dependencies and working environment.

## Opening a pull request
You can open a [pull request](https://github.com/lllingoo/2QAN/pulls) by pushing changes from a local branch, explaining the bug fix or new feature.

### Version control with git
git is a language that helps keeping track of the changes made. Have a look at these guidelines for getting started with [git workflow](https://www.asmeurer.com/git-workflow/).
Use short and explanatory comments to document the changes with frequent commits.

### Forking the repository
You can fork 2QAN from the github repository, so that your changes are applied with respect to the current master branch. Use the Fork button, and then use git from the command line to clone your fork of the repository locally on your machine.
```bash
(base) git clone https://github.com/your_github_username/2QAN.git
```
You can also use SSH instead of a HTTPS protocol.

### Working in a virtual environment
It is best to set up a clean environment with anaconda, to keep track of all installed applications.
```bash
(base) conda create -n myenv python=3
```
accept the configuration ([y]) and switch to the environment
```bash
(base) conda activate myenv
(myenv) conda install pip
```
Once you will finish the modifications, you can deactivate the environment with
```bash
(myenv) conda deactivate myenv
```

#### Special Note for Windows Users Using Python 3.8:
To prevent errors when running `make docs` and `make doctest`, Windows developers using Python 3.8 will also need to edit `__init__.py` in their environment's asyncio directory.
This is due to Python changing `asyncio`'s [default event loop in Windows beginning in Python 3.8](https://docs.python.org/3/library/asyncio-policy.html#asyncio.DefaultEventLoopPolicy).
The new default event loop will not support Unix-style APIs used by some dependencies.
1. Locate your environment directory (likely `C:\Users\{username}\anaconda3\envs\{your_env}`), and open `{env_dir}/Lib/asyncio/__init__.py`.
2. Add `import asyncio` to the file's import statements.
3. Find the block of code below and replace it with the provided replacement.
    * Original Code  

          if sys.platform == 'win32':  # pragma: no cover
              from .windows_events import *
              __all__ += windows_events.__all__
          else:
              from .unix_events import *  # pragma: no cover
              __all__ += unix_events.__all__
  
    * Replacement Code  

          if sys.platform == 'win32':  # pragma: no cover
              from .windows_events import *
              asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
              __all__ += windows_events.__all__
          else:
              from .unix_events import *  # pragma: no cover
              __all__ += unix_events.__all__

## Code of conduct
2QAN development abides to the [Contributors' Covenant](https://github.com/lllingoo/2QAN/blob/master/CODE_OF_CONDUCT.md).

