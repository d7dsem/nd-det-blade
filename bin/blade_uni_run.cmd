@echo off
rem run.cmd - same as run.ps1, but works from cmd.exe and without changing ExecutionPolicy.
rem Usage: bin\run.cmd <script>[.py] [args...]
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run.ps1" %*
exit /b %ERRORLEVEL%