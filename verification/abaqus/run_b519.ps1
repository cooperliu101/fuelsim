param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b519_hex8_c3d8t_finite_selective"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "${JobName}_nodal.csv"
$IntegrationOutput = Join-Path $SourceDirectory "${JobName}_integration.csv"
$Extractor = Join-Path $SourceDirectory "extract_b519.py"
$RunDirectory = Join-Path $env:TEMP ("fuelsim_b519_" + (Get-Date -Format "yyyyMMdd_HHmmss"))

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    & abaqus job=$JobName input=$InputFile output_precision=full interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.19 probe failed with exit code $LASTEXITCODE" }
    $StatusFile = Join-Path $RunDirectory "$JobName.sta"
    if (!(Test-Path $StatusFile) -or !(Select-String -Quiet -Path $StatusFile -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY")) {
        throw "Abaqus B5.19 probe did not complete successfully; inspect $RunDirectory"
    }
    & abaqus python $Extractor "$JobName.odb" $NodalOutput $IntegrationOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.19 extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
