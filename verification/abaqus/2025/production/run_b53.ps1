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
$work = Join-Path $env:TEMP ("fuelsim_b53_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b53_hex8_c3d8t_thermal_load_probe.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b53.py") $work
Set-Location $work

Invoke-Abaqus `
    job=b53_hex8_c3d8t_thermal_load_probe `
    input=b53_hex8_c3d8t_thermal_load_probe.inp `
    output_precision=full `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.3 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b53_hex8_c3d8t_thermal_load_probe.sta") -or
    !(Select-String -Path "b53_hex8_c3d8t_thermal_load_probe.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.3 probe did not report successful completion"
}

Invoke-Abaqus `
    python extract_b53.py `
    b53_hex8_c3d8t_thermal_load_probe.odb `
    b53_hex8_c3d8t_thermal_load.csv `
    b53_hex8_c3d8t_thermal_load_integration.csv `
    b53_hex8_c3d8t_surface_geometry.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.3 extraction failed with exit code $LASTEXITCODE"
}
foreach ($artifact in @(
    "b53_hex8_c3d8t_thermal_load.csv",
    "b53_hex8_c3d8t_thermal_load_integration.csv",
    "b53_hex8_c3d8t_surface_geometry.csv"
)) {
    if (!(Test-Path $artifact)) {
        throw "Abaqus B5.3 extraction did not create $artifact"
    }
    Copy-Item $artifact $SourceDirectory
}
Write-Output "Abaqus B5.3 work directory: $work"
