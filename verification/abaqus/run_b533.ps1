param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$job = "b533_hex8_c3d8rt_finite_hourglass_probe"
$output = "b533_hex8_c3d8rt_finite_hourglass_nodal.csv"
$work = Join-Path $env:TEMP ("fuelsim_b533_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b533.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" job=$job input="${job}.inp" output_precision=full interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.33 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.33 probe did not report successful completion"
}
& "C:\SIMULIA\Commands\abaqus.bat" python extract_b533.py "${job}.odb" $output
if ($LASTEXITCODE -ne 0 -or !(Test-Path $output)) {
    throw "Abaqus B5.33 extraction failed"
}
Copy-Item $output $SourceDirectory
Write-Output "Abaqus B5.33 work directory: $work"
