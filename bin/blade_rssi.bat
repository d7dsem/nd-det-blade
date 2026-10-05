@echo off
setlocal
set "ROOT=%~dp0.."
set "PY=python"
if exist "%ROOT%\.venv\Scripts\python.exe" set "PY=%ROOT%\.venv\Scripts\python.exe"
set "PYTHONPATH=%ROOT%\source;%PYTHONPATH%"
"%PY%" "%ROOT%\source\rssi_scan_blade_entry.py" %*
exit /b %ERRORLEVEL%
