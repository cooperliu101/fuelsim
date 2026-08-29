param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory,
    [Parameter(Mandatory = $true)]
    [ValidateSet("regular", "warped", "holdout")]
    [string]$Geometry
)

$ErrorActionPreference = "Stop"
$stem = if ($Geometry -eq "regular") {
    "b528_hex8_c3d8rt_operator"
} else {
    "b528_hex8_c3d8rt_${Geometry}_operator"
}
$job = "${stem}_probe"
$work = Join-Path $env:TEMP ("fuelsim_b528_" + $Geometry + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b528.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" job=$job input="${job}.inp" output_precision=full interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.28 $Geometry probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.28 $Geometry probe did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b528.py "${job}.odb" "${stem}_nodal.csv" "${stem}_integration.csv"
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.28 $Geometry extraction failed with exit code $LASTEXITCODE"
}
foreach ($artifact in @("${stem}_nodal.csv", "${stem}_integration.csv")) {
    if (!(Test-Path $artifact)) {
        throw "Abaqus B5.28 extraction did not create $artifact"
    }
    Copy-Item $artifact $SourceDirectory
}
Write-Output "Abaqus B5.28 work directory: $work"
