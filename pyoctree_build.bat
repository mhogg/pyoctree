REM Set variables for Microsoft Visual Studio for Python 2.7 
REM @call "C:\Users\michaelh3\AppData\Local\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat" amd64
REM Set environ as workaround for bugs (See http://stackoverflow.com/questions/2817869/error-unable-to-find-vcvarsall-bat and http://bugs.python.org/issue23246)
REM SET DISTUTILS_USE_SDK=1
REM SET MSSdk=1

REM Set variables for Microsoft Visual Studio for Python 3.6 
@call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86_amd64

REM Test build
python setup.py sdist bdist_wheel
REM Upload to PyPi
REM python setup.py sdist bdist_wheel upload
