REM Set variables for Microsoft Visual Studio for Python 2.7 
REM "C:\Users\michaelh3\AppData\Local\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat" amd64
REM Set environ as workaround for bugs (See http://stackoverflow.com/questions/2817869/error-unable-to-find-vcvarsall-bat and http://bugs.python.org/issue23246)
SET DISTUTILS_USE_SDK=1
SET MSSdk=1
REM Test build
REM python setup.py sdist bdist_wheel 
REM Upload to PyPi
REM python setup.py sdist bdist_wheel upload
