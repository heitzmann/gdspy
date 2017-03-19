CALL "c:\Program Files (x86)\Microsoft Visual C++ Build Tools\vcbuildtools.bat" amd64
cd \gdspy
c:\WinPython-64bit-3.6.0.1Zero\python-3.6.0.amd64\python.exe setup.py test
c:\WinPython-64bit-3.6.0.1Zero\python-3.6.0.amd64\python.exe setup.py bdist_wininst
c:\WinPython-64bit-3.5.3.0Zero\python-3.5.3.amd64\python.exe setup.py test
c:\WinPython-64bit-3.5.3.0Zero\python-3.5.3.amd64\python.exe setup.py bdist_wininst

CALL "%LOCALAPPDATA%\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat" amd64
cd \gdspy
c:\WinPython-64bit-2.7.13.0Zero\python-2.7.13.amd64\python.exe setup.py test
c:\WinPython-64bit-2.7.13.0Zero\python-2.7.13.amd64\python.exe setup.py bdist_wininst
