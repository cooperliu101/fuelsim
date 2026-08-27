param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b55_hex8_c3d8t_transient_full_field"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "b55_hex8_c3d8t_transient_nodal.csv"
$IntegrationOutput = Join-Path $SourceDirectory "b55_hex8_c3d8t_transient_integration.csv"
$Extractor = Join-Path $SourceDirectory "extract_b55.py"
$RunDirectory = Join-Path $env:TEMP ("fuelsim_b55_" + (Get-Date -Format "yyyyMMdd_HHmmss"))

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & abaqus job=$JobName input=$InputFile interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus job failed with exit code $LASTEXITCODE" }
    & abaqus python $Extractor "$JobName.odb" $NodalOutput $IntegrationOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
