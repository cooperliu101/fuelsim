param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b57_hex8_c3d8t_temperature_capacity"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$Output = Join-Path $SourceDirectory "$JobName.csv"
$Extractor = Join-Path $SourceDirectory "extract_b57.py"
$RunDirectory = Join-Path $env:TEMP ("fuelsim_b57_" + (Get-Date -Format "yyyyMMdd_HHmmss"))

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & "C:\SIMULIA\Commands\abq2025.bat" job=$JobName input=$InputFile output_precision=full interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.7 probe failed with exit code $LASTEXITCODE" }
    & "C:\SIMULIA\Commands\abq2025.bat" python $Extractor "$JobName.odb" $Output
    if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.7 extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
