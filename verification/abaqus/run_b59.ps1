param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b59_hex8_c3d8t_multimaterial"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "${JobName}_nodal.csv"
$IntegrationOutput = Join-Path $SourceDirectory "${JobName}_integration.csv"
$Extractor = Join-Path $SourceDirectory "extract_b59.py"
$RunDirectory = Join-Path $env:TEMP ("fuelsim_b59_" + (Get-Date -Format "yyyyMMdd_HHmmss"))

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & abaqus job=$JobName input=$InputFile output_precision=full interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.9 probe failed with exit code $LASTEXITCODE" }
    & abaqus python $Extractor "$JobName.odb" $NodalOutput $IntegrationOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.9 extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
