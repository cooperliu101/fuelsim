param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [ValidateSet("small_coupled", "finite_coupled", "small_plastic", "finite_plastic", "small_creep", "finite_creep", "small_creep_hold")]
      [string]$CaseName,
      [string]$ReuseWorkDirectory = "")
function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$JobName = "c3d20rt_$CaseName"
$Work = $ReuseWorkDirectory
if ($Work -eq "") {
    $Work = Join-Path $env:TEMP ($JobName + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    Copy-Item (Join-Path $SourceDirectory "${JobName}.inp") $Work
} elseif ((Get-FileHash (Join-Path $Work "${JobName}.inp")).Hash -ne
          (Get-FileHash (Join-Path $SourceDirectory "${JobName}.inp")).Hash) {
    throw "Existing Abaqus work uses a different input file"
}
Write-Output "abaqus_work_directory=$Work"
Copy-Item (Join-Path $SourceDirectory "extract_c3d20rt_material.py") $Work
Set-Location $Work
if ($ReuseWorkDirectory -eq "") {
    Invoke-Abaqus "job=$JobName" "input=${JobName}.inp" output_precision=full ask_delete=OFF cpus=1 interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus material analysis failed" }
}
if (!(Test-Path "${JobName}.sta")) { throw "Missing Abaqus completion status" }
if (!(Select-String -Path "${JobName}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus material analysis did not complete"
}
Invoke-Abaqus python extract_c3d20rt_material.py "${JobName}.odb" $JobName
if ($LASTEXITCODE -ne 0) { throw "Abaqus material extraction failed" }
foreach ($Suffix in @("_nodes.csv", "_points.csv")) {
    if (!(Test-Path "${JobName}${Suffix}") -or (Get-Content "${JobName}${Suffix}").Count -lt 2) {
        throw "Abaqus extraction returned an empty reference: $Suffix"
    }
    Copy-Item "${JobName}${Suffix}" $SourceDirectory
}
foreach ($Extension in @("sta", "msg", "dat")) {
    Copy-Item "${JobName}.${Extension}" (Join-Path $SourceDirectory "${JobName}_reference.${Extension}")
}
