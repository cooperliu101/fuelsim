param([string]$SourceDirectory = $PSScriptRoot)
$ErrorActionPreference = "Stop"
$Work = Join-Path $env:TEMP ("fuelsim_density_capacity_2025_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "density_capacity.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract.py") $Work
$SourceHash = (Get-FileHash (Join-Path $SourceDirectory "density_capacity.inp") -Algorithm SHA256).Hash
if ($SourceHash -ne (Get-FileHash (Join-Path $Work "density_capacity.inp") -Algorithm SHA256).Hash) { throw "Input changed during staging" }
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
Push-Location $Work
try {
    & "C:\SIMULIA\Commands\abq2025.bat" job=density_capacity input=density_capacity.inp cpus=1 output_precision=full interactive | Tee-Object run.log
    if ($LASTEXITCODE -ne 0 -or !(Select-String -Path density_capacity.sta -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) { throw "Native solve failed in $Work" }
    & "C:\SIMULIA\Commands\abq2025.bat" python extract.py density_capacity.odb | Tee-Object extract.log
    if ($LASTEXITCODE -ne 0 -or !(Test-Path nodes.csv) -or @(Import-Csv nodes.csv).Count -ne 4800) { throw "Extraction failed in $Work" }
    foreach ($File in @("nodes.csv", "points.csv", "density_capacity.dat", "density_capacity.msg", "density_capacity.sta", "run.log")) { Copy-Item $File $SourceDirectory }
    @("abaqus_version=2025 RELr427", "cpus=1", "input_sha256=$SourceHash", "work_directory=$Work", "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "extractor_sha256=$((Get-FileHash extract.py -Algorithm SHA256).Hash)") | Set-Content (Join-Path $SourceDirectory "provenance.txt")
    Write-Output "Completed in $Work"
} finally { Pop-Location }
