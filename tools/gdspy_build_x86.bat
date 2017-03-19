CALL "c:\Program Files (x86)\Microsoft Visual C++ Build Tools\vcbuildtools.bat" x86

cd \gdspy

c:\WinPython-32bit-3.6.0.1Zero\python-3.6.0\python.exe setup.py test

c:\WinPython-32bit-3.6.0.1Zero\python-3.6.0\python.exe setup.py bdist_wininst

c:\WinPython-32bit-3.5.3.0Zero\python-3.5.3\python.exe setup.py test

c:\WinPython-32bit-3.5.3.0Zero\python-3.5.3\python.exe setup.py bdist_wininst



CALL "%LOCALAPPDATA%\Programs\Common\Microsoft\Visual C++ for Python\9.0\vcvarsall.bat" x86

cd \gdspy

c:\WinPython-32bit-2.7.13.0Zero\python-2.7.13\python.exe setup.py test

c:\WinPython-32bit-2.7.13.0Zero\python-2.7.13\python.exe setup.py bdist_wininst
