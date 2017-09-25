
REM Python 2.7
REM ----------

REM Set variables for Microsoft Visual Studio for Python 2.7 
@call "C:\Users\michaelh3\AppData\Local\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat" amd64
REM Set environ as workaround for bugs (See http://stackoverflow.com/questions/2817869/error-unable-to-find-vcvarsall-bat and http://bugs.python.org/issue23246)
SET DISTUTILS_USE_SDK=1
SET MSSdk=1

REM Cythonize pyx file to cpp (for users without Cython installed)
cython pyoctree\pyoctree.pyx --cplus

REM Build source + Python 2.7 wheel - With openmp
python setup.py sdist bdist_wheel --openmp

REM Python 3.6
REM ----------

REM Activate Python 3.6
activate python36

REM Set variables for Microsoft Visual Studio for Python 3.6 
@call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86_amd64

REM Build Python 3.6 wheel - With openmp
python setup.py bdist_wheel --openmp

REM Upload entire dist folder to PyPi
REM ---------------------------------

REM Upload to PyPi
REM twine upload dist\*
