param([Parameter(Mandatory=$true)][string]$SourceDirectory, [string]$CaseName = "")
$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_thermal_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "*.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract.py") $work
Set-Location $work
$cases = @("dcax4_transient", "dcax8_transient", "dc3d8_transient", "dc3d20_transient", "dcax4_steady", "dcax8_steady", "dc3d8_steady", "dc3d20_steady", "dcax8_contact", "dcax4_capacity", "dc3d8_capacity")
if ($CaseName -ne "") { $cases = @($CaseName) }
foreach ($name in $cases) {
    & "C:\SIMULIA\Commands\abaqus.bat" job=$name input="$name.inp" cpus=1 interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus failed: $name" }
    if (!(Select-String -Path "$name.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus incomplete: $name"
    }
    & "C:\SIMULIA\Commands\abaqus.bat" python extract.py "$name.odb" "$name.nodes.csv" "$name.points.csv"
    if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $name" }
    Copy-Item "$name.nodes.csv", "$name.points.csv", "$name.sta", "$name.dat" $SourceDirectory
}
Write-Output "Abaqus thermal work directory: $work"
