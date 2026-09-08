param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [ValidateSet("h20_28_finite_disconnected_probe", "h20_28_finite_transfer_shallow_probe", "h20_28_finite_transfer_probe", "h20_28_finite_transfer_16_probe", "h20_28_finite_point_probe", "h20_28_finite_point_edge_probe", "h20_28_finite_linear_transfer_probe", "h20_28_finite_transfer_32_probe")]
      [string]$Job = "h20_28_finite_transfer_probe")
$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = Join-Path $env:TEMP ($Job + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
Write-Output "abaqus_work_directory=$Work"
Copy-Item (Join-Path $SourceDirectory ($Job + ".inp")) $Work
Copy-Item (Join-Path $SourceDirectory "extract_h20_finite_transfer.py") $Work
Set-Location $Work
& "C:\SIMULIA\Commands\abaqus.bat" "job=$Job" "input=$Job.inp" cpus=1 output_precision=full ask_delete=OFF interactive
if ($LASTEXITCODE -ne 0 -or !(Select-String -Path "$Job.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Finite transfer probe did not complete"
}
$PrimaryCount = if ($Job -eq "h20_28_finite_disconnected_probe") { 512 } elseif ($Job -eq "h20_28_finite_transfer_32_probe") { 3201 } elseif ($Job -eq "h20_28_finite_transfer_16_probe") { 833 } else { 225 }
& "C:\SIMULIA\Commands\abaqus.bat" python extract_h20_finite_transfer.py "$Job.odb" "${Job}_operator.csv" $PrimaryCount
if ($LASTEXITCODE -ne 0 -or (Get-Content "${Job}_operator.csv").Count -ne (1 + 17 * (8 + $PrimaryCount))) {
    throw "Finite transfer extraction has an unexpected row count"
}
Copy-Item "${Job}_operator.csv" $SourceDirectory
foreach ($Extension in @("sta", "msg", "dat")) {
    Copy-Item "$Job.$Extension" (Join-Path $SourceDirectory "${Job}_reference.$Extension")
}
