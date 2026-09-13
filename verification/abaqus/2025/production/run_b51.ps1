param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$JobName = "b51_hex8_c3d8t_finite_heat_probe"
$InputFile = Join-Path $SourceDirectory "$JobName.inp"
$NodalOutput = Join-Path $SourceDirectory "b51_hex8_c3d8t_finite_heat_nodal.csv"
$IntegrationOutput = Join-Path $SourceDirectory "b51_hex8_c3d8t_finite_heat_integration_points.csv"
$Extractor = Join-Path $SourceDirectory "extract_b51.py"
$RunDirectory = Join-Path $env:TEMP "fuelsim_b51_abaqus"

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    Invoke-Abaqus job=$JobName input=$InputFile cpus=1 interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus job failed with exit code $LASTEXITCODE" }
    Invoke-Abaqus python $Extractor "$JobName.odb" $NodalOutput $IntegrationOutput
    if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}
