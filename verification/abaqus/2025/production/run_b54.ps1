param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b54_hex8_c3d8t_finite_thermal_load_probe"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$CapacityOutput = Join-Path $SourceDirectory "b54_hex8_c3d8t_finite_capacity.csv"
$LoadOutput = Join-Path $SourceDirectory "b54_hex8_c3d8t_finite_thermal_load.csv"
$Extractor = Join-Path $SourceDirectory "extract_b54.py"
$RunDirectory = Join-Path $env:TEMP ("fuelsim_b54_" + (Get-Date -Format "yyyyMMdd_HHmmss"))

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & "C:\SIMULIA\Commands\abq2025.bat" job=$JobName input=$InputFile interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus job failed with exit code $LASTEXITCODE" }
    & "C:\SIMULIA\Commands\abq2025.bat" python $Extractor "$JobName.odb" $CapacityOutput $LoadOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
