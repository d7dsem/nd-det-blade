# run.ps1 - wrapper: runs a Python script from ..\source with all arguments passed through.
# Usage:  .\bin\run.ps1 <script>[.py] [args...]
# Example: .\bin\run.ps1 blade_tst_entry --freq 100e6 -v "some text"
#
# No param() block on purpose: unknown "-x"/"--xyz" args land in $args untouched.

$ErrorActionPreference = 'Stop'

if ($args.Count -lt 1) {
    [Console]::Error.WriteLine('Usage: run.ps1 <script>[.py] [args...]')
    exit 2
}

$Root      = Split-Path -Parent $PSScriptRoot      # D:\Programming\blade
$SourceDir = Join-Path $Root 'source'

# Script name: ".py" extension is optional
$Name = [string]$args[0]
if (-not $Name.EndsWith('.py')) { $Name += '.py' }
$ScriptPath = Join-Path $SourceDir $Name

if (-not (Test-Path -LiteralPath $ScriptPath)) {
    [Console]::Error.WriteLine("Script not found: $ScriptPath")
    exit 2
}

# All remaining arguments go to the Python script unchanged
$Rest = @()
if ($args.Count -gt 1) { $Rest = $args[1..($args.Count - 1)] }

# Interpreter: project venv (.venv) if present, otherwise python from PATH
$VenvPython = Join-Path $Root '.venv\Scripts\python.exe'
if (Test-Path -LiteralPath $VenvPython) {
    $Python = $VenvPython
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $Python = 'python'
} else {
    [Console]::Error.WriteLine('Python not found: no .venv and no python in PATH')
    exit 2
}

# source\ on PYTHONPATH so modules there can import each other; restored afterwards
$OldPyPath = $env:PYTHONPATH
try {
    if ($OldPyPath) { $env:PYTHONPATH = "$SourceDir;$OldPyPath" }
    else            { $env:PYTHONPATH = $SourceDir }

    & $Python $ScriptPath @Rest
    $Code = $LASTEXITCODE
}
finally {
    $env:PYTHONPATH = $OldPyPath
}

exit $Code
