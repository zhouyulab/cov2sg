# rnatk
RNA Tool Kit

Python package rnatk>=0.4.0 for subgenome analysis


## Tested in Python3

## Installing

1. Install bx-python first, if not installed.

```
pip install bx-python>=0.8.0
pip install biopython
pip install pysam
pip install numpy cython
pip install plastid
```

2. Install from cloned repo

```
git clone https://github.com/zhouyulab/rnatk
python setup.py test # Optional
python setup.py install
```

3. Test in development
```
set PYTHONPATH=src
pytest -s tests/test_op_pemerge.py
```

4. PyPI
```
git tag v0.1.10
git push origin --tags
python setup.py bdist_wheel
```
