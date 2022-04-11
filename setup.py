# Copyright (C) 2022 Lingling Lao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from setuptools import setup, find_packages

with open("VERSION.txt", "r") as f:
    __version__ = f.read().strip()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

# save the source code in _version.py
with open("py2qan/_version.py", "r") as f:
    version_file_source = f.read()

# overwrite _version.py in the source distribution
with open("py2qan/_version.py", "w") as f:
    f.write(f"__version__ = '{__version__}'\n")

setup(
    name="py2QAN",
    version=__version__,
    install_requires=requirements,
    # packages=find_packages("py2qan"),
    # package_dir={'':'py2qan'},
    packages=find_packages(),
    include_package_data=True,
    description="2QAN is an open source compiler for two-local qubit Hamiltonian simulation algorithms.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Lingling Lao",
    author_email="laolinglingrolls@gmail.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Topic :: Scientific/Engineering",
    ],
    url = "https://github.com/lllingoo/2QAN",
    python_requires=">=3.8",
    verify=False
)
