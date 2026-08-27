param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$Extractor = Join-Path $SourceDirectory "extract_b520.py"
$RunDirectory = Join-Path $env:TEMP "fuelsim_b520_abaqus"

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    foreach ($Suffix in @("clearance", "pressure")) {
        $JobName = "b520_hex8_c3d8t_gap_conductance_$Suffix"
        $InputFile = Join-Path $SourceDirectory "$JobName.inp"
        $NodalOutput = Join-Path $SourceDirectory "${JobName}_nodal.csv"
        & abaqus job=$JobName input=$InputFile interactive
        if ($LASTEXITCODE -ne 0) { throw "Abaqus job $JobName failed with exit code $LASTEXITCODE" }
        & abaqus python $Extractor "$JobName.odb" $NodalOutput
        if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction for $JobName failed with exit code $LASTEXITCODE" }
    }
}
finally {
    Pop-Location
}
