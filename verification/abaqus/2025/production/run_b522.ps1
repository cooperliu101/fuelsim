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
$Extractor = Join-Path $SourceDirectory "extract_b522.py"
$RunDirectory = Join-Path $env:TEMP "fuelsim_b522_abaqus"
$JobNames = @(
    "b522_hex8_c3d8t_faceted_thermal_contact",
    "b522_hex8_c3d8t_faceted_thermal_contact_f6",
    "b522_hex8_c3d8t_faceted_thermal_contact_f12",
    "b522_hex8_c3d8t_faceted_thermal_contact_swapped",
    "b522_hex8_c3d8t_faceted_thermal_contact_tight"
)

New-Item -ItemType Directory -Force -Path $RunDirectory | Out-Null
Push-Location $RunDirectory
try {
    foreach ($JobName in $JobNames) {
        $InputFile = Join-Path $SourceDirectory "$JobName.inp"
        $NodalOutput = Join-Path $SourceDirectory "${JobName}_nodal.csv"
        $ContactOutput = Join-Path $SourceDirectory "${JobName}_contact.csv"
        Invoke-Abaqus job=$JobName input=$InputFile cpus=1 interactive
        if ($LASTEXITCODE -ne 0) { throw "Abaqus job $JobName failed with exit code $LASTEXITCODE" }
        Invoke-Abaqus python $Extractor "$JobName.odb" $NodalOutput $ContactOutput
        if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction for $JobName failed with exit code $LASTEXITCODE" }
    }
}
finally {
    Pop-Location
}
