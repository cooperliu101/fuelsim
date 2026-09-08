param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [string]$ReuseWorkDirectory = "",
      [ValidateSet("c3d20t_finite_contact_frictionless", "c3d20rt_finite_contact_frictionless", "c3d20rt_finite_gap_matrix", "c3d20rt_gap_small_sliding_probe", "c3d20rt_gap_geometry_probe", "c3d20rt_contact_friction", "c3d20rt_finite_contact_friction", "c3d20rt_finite_matching_matrix", "c3d20rt_contact_unit_flux", "c3d20rt_contact_precision", "c3d20rt_contact_thermal_matrix", "c3d20rt_contact_matching_matrix")]
      [string]$JobName = "c3d20rt_contact_friction")
$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = $ReuseWorkDirectory
if ($Work -eq "") {
    $Work = Join-Path $env:TEMP ($JobName + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $SourceDirectory "${JobName}.inp") $Work
} elseif ((Get-FileHash (Join-Path $Work "${JobName}.inp")).Hash -ne
          (Get-FileHash (Join-Path $SourceDirectory "${JobName}.inp")).Hash) {
    throw "Existing Abaqus task uses a different input"
}
Write-Output "abaqus_work_directory=$Work"
Copy-Item (Join-Path $SourceDirectory "extract_c3d20rt_contact.py") $Work
Set-Location $Work
if ($ReuseWorkDirectory -eq "") {
    & "C:\SIMULIA\Commands\abaqus.bat" "job=$JobName" "input=${JobName}.inp" cpus=1 output_precision=full ask_delete=OFF interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus contact analysis failed" }
}
if (!(Test-Path "${JobName}.sta") -or
    !(Select-String -Path "${JobName}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus contact analysis did not complete"
}
$ContactNodes = "57,60,61,64,68,69,72,76,78,80,83,85,88"
if ($JobName -in @("c3d20rt_contact_matching_matrix", "c3d20rt_finite_matching_matrix")) { $ContactNodes = "57,60,61,64,68,69,72,76" }
& "C:\SIMULIA\Commands\abaqus.bat" python extract_c3d20rt_contact.py "${JobName}.odb" $JobName $ContactNodes
if ($LASTEXITCODE -ne 0) { throw "Abaqus contact extraction failed" }
foreach ($Suffix in @("_nodes.csv", "_contact.csv", "_points.csv")) {
    if (!(Test-Path "${JobName}${Suffix}") -or (Get-Content "${JobName}${Suffix}").Count -lt 2) {
        throw "Abaqus extraction returned an empty reference: $Suffix"
    }
    Copy-Item "${JobName}${Suffix}" $SourceDirectory
}
foreach ($Extension in @("sta", "msg", "dat")) {
    Copy-Item "${JobName}.${Extension}" (Join-Path $SourceDirectory "${JobName}_reference.${Extension}")
}
