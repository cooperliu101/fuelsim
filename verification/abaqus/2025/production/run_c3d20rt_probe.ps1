param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [ValidateSet("c3d20rt_thermal_probe", "c3d20rt_nonaffine_probe", "c3d20rt_small_probe")]
      [string]$JobName = "c3d20rt_thermal_probe")
function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = Join-Path $env:TEMP ("fuelsim_c3d20rt_probe_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
Write-Output "abaqus_work_directory=$Work"
Copy-Item (Join-Path $SourceDirectory "${JobName}.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_c3d20rt_probe.py") $Work
Set-Location $Work
Invoke-Abaqus "job=$JobName" "input=${JobName}.inp" output_precision=full ask_delete=OFF cpus=1 interactive
if ($LASTEXITCODE -ne 0) { throw "C3D20RT probe failed: $LASTEXITCODE" }
if (!(Test-Path "${JobName}.sta") -or
    !(Select-String -Path "${JobName}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "C3D20RT probe did not complete successfully"
}
Copy-Item "${JobName}.dat" (Join-Path $SourceDirectory "${JobName}_reference.dat")
Copy-Item "${JobName}.sta" (Join-Path $SourceDirectory "${JobName}_reference.sta")
Invoke-Abaqus python extract_c3d20rt_probe.py "${JobName}.odb" "${JobName}_nodes.csv"
if ($LASTEXITCODE -ne 0) { throw "C3D20RT reference extraction failed" }
Copy-Item "${JobName}_nodes.csv" $SourceDirectory
Write-Output "abaqus_work_directory=$Work"
