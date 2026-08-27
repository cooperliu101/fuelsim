param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b56_hex8_c3d8t_temperature_operator"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "b56_hex8_c3d8t_temperature_operator_nodal.csv"
$IntegrationOutput = Join-Path $SourceDirectory "b56_hex8_c3d8t_temperature_operator_integration.csv"
$Extractor = Join-Path $SourceDirectory "extract_b49.py"
$RunDirectory = Join-Path $env:TEMP ("fuelsim_b56_" + (Get-Date -Format "yyyyMMdd_HHmmss"))

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & abaqus job=$JobName input=$InputFile output_precision=full interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus job failed with exit code $LASTEXITCODE" }
    & abaqus python $Extractor "$JobName.odb" $NodalOutput $IntegrationOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
