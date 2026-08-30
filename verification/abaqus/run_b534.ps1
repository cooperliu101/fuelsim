param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$job = "b534_hex8_c3d8rt_finite_transient_full_field"
$work = Join-Path $env:TEMP ("fuelsim_b534_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b534.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" job=$job input="${job}.inp" output_precision=full interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.34 transient full-field probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.34 transient full-field probe did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b534.py "${job}.odb" `
    b534_hex8_c3d8rt_finite_transient_nodal.csv `
    b534_hex8_c3d8rt_finite_transient_integration.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.34 extraction failed with exit code $LASTEXITCODE"
}
foreach ($artifact in @(
    "b534_hex8_c3d8rt_finite_transient_nodal.csv",
    "b534_hex8_c3d8rt_finite_transient_integration.csv"
)) {
    if (!(Test-Path $artifact)) {
        throw "Abaqus B5.34 extraction did not create $artifact"
    }
    Copy-Item $artifact $SourceDirectory
}
Write-Output "Abaqus B5.34 work directory: $work"
