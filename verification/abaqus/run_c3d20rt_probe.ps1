param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [ValidateSet("c3d20rt_thermal_probe", "c3d20rt_nonaffine_probe", "c3d20rt_small_probe")]
      [string]$JobName = "c3d20rt_thermal_probe")
$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = Join-Path $env:TEMP ("fuelsim_c3d20rt_probe_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
Write-Output "abaqus_work_directory=$Work"
Copy-Item (Join-Path $SourceDirectory "${JobName}.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_c3d20rt_probe.py") $Work
Set-Location $Work
& "C:\SIMULIA\Commands\abaqus.bat" "job=$JobName" "input=${JobName}.inp" cpus=1 output_precision=full ask_delete=OFF interactive
if ($LASTEXITCODE -ne 0) { throw "C3D20RT probe failed: $LASTEXITCODE" }
if (!(Test-Path "${JobName}.sta") -or
    !(Select-String -Path "${JobName}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "C3D20RT probe did not complete successfully"
}
Copy-Item "${JobName}.dat" (Join-Path $SourceDirectory "${JobName}_reference.dat")
Copy-Item "${JobName}.sta" (Join-Path $SourceDirectory "${JobName}_reference.sta")
& "C:\SIMULIA\Commands\abaqus.bat" python extract_c3d20rt_probe.py "${JobName}.odb" "${JobName}_nodes.csv"
if ($LASTEXITCODE -ne 0) { throw "C3D20RT reference extraction failed" }
Copy-Item "${JobName}_nodes.csv" $SourceDirectory
Write-Output "abaqus_work_directory=$Work"
