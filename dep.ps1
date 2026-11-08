# Installs dependencies from PowerShell by running dep.sh inside WSL (Ubuntu/Debian).
if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    Write-Error "WSL is not installed. Install it with 'wsl --install', then re-run this script."
    exit 1
}
wsl.exe --cd "$PSScriptRoot" ./dep.sh
exit $LASTEXITCODE
