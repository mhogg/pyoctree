
REM Python 3.7
REM ----------

REM Activate Python 3.7
REM @call activate python37

REM Set variables for Microsoft Visual Studio
@call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64

REM Cythonize pyx file to cpp (for users without Cython installed)
del pyoctree\pyoctree.cpp
cython pyoctree\pyoctree.pyx --cplus

REM Build Python 3.7 wheel - With openmp
python setup.py sdist bdist_wheel --openmp

REM Upload entire dist folder to PyPi
REM ---------------------------------

REM Upload to PyPi
REM twine upload dist\*
