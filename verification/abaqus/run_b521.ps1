param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b521_hex8_c3d8t_thermal_contact_path"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "${JobName}_nodal.csv"
$Extractor = Join-Path $SourceDirectory "extract_b521.py"
$RunDirectory = Join-Path $env:TEMP "fuelsim_b521_abaqus"

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
