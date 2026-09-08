param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [ValidateSet("h20_28_primary_signed_gap_probe", "h20_28_primary_signed_gap_shallow_probe", "h20_28_primary_signed_gap_32_probe")]
      [string]$Job = "h20_28_primary_signed_gap_probe")
$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = Join-Path $env:TEMP ($Job + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
Write-Output "abaqus_work_directory=$Work"
Copy-Item (Join-Path $SourceDirectory "$Job.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_h20_primary_gap_probe.py") $Work
Set-Location $Work
& "C:\SIMULIA\Commands\abaqus.bat" "job=$Job" "input=$Job.inp" cpus=1 output_precision=full ask_delete=OFF interactive
if ($LASTEXITCODE -ne 0 -or !(Select-String -Path "$Job.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) { throw "Primary gap probe failed" }
& "C:\SIMULIA\Commands\abaqus.bat" python extract_h20_primary_gap_probe.py "$Job.odb" "${Job}_gap.csv"
if ($LASTEXITCODE -ne 0 -or (Get-Content "${Job}_gap.csv").Count -ne 25) { throw "Primary gap extraction failed" }
Copy-Item "${Job}_gap.csv" $SourceDirectory
foreach ($Extension in @("sta", "msg", "dat")) { Copy-Item "$Job.$Extension" (Join-Path $SourceDirectory "${Job}_reference.$Extension") }
