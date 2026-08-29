param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$job = "b530_hex8_c3d8rt_thermal_load_probe"
$work = Join-Path $env:TEMP ("fuelsim_b530_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b530.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" job=$job input="${job}.inp" output_precision=full interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.30 thermal-load probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.30 thermal-load probe did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b530.py "${job}.odb" `
    b530_hex8_c3d8rt_thermal_load.csv `
    b530_hex8_c3d8rt_thermal_load_integration.csv `
    b530_hex8_c3d8rt_surface_geometry.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.30 thermal-load extraction failed with exit code $LASTEXITCODE"
}
foreach ($artifact in @(
    "b530_hex8_c3d8rt_thermal_load.csv",
    "b530_hex8_c3d8rt_thermal_load_integration.csv",
    "b530_hex8_c3d8rt_surface_geometry.csv"
)) {
    if (!(Test-Path $artifact)) {
        throw "Abaqus B5.30 extraction did not create $artifact"
    }
    Copy-Item $artifact $SourceDirectory
}
Write-Output "Abaqus B5.30 work directory: $work"
