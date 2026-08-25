param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b52_hex8_c3d8t_thermal_contact"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "b52_hex8_c3d8t_thermal_contact_nodal.csv"
$Extractor = Join-Path $SourceDirectory "extract_b52.py"
$RunDirectory = Join-Path $env:TEMP "fuelsim_b52_abaqus"

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & abaqus job=$JobName input=$InputFile interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus job failed with exit code $LASTEXITCODE" }
    & abaqus python $Extractor "$JobName.odb" $NodalOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
