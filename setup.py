import setuptools
from glob import glob
from os.path import basename
from os.path import splitext
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PySGM",
    version="0.1.6",
    author="Hiroyuki Goto",
    author_email="goto@catfish.dpri.kyoto-u.ac.jp",
    description="PySGM",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    # packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
