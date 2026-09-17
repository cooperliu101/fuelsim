param([Parameter(Mandatory=$true)][string]$SourceDirectory,
      [string]$JobName='reference_node_probe')
$ErrorActionPreference = 'Stop'
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = Join-Path $env:TEMP ('cpeg8t_'+$JobName+'_'+(Get-Date -Format 'yyyyMMdd_HHmmss'))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory ($JobName+'.inp')) $Work
Copy-Item (Join-Path $PSScriptRoot 'extract.py') $Work
Set-Location $Work
& 'C:\SIMULIA\Commands\abaqus.bat' "job=$JobName" "input=$JobName.inp" cpus=1 output_precision=full ask_delete=OFF interactive
$SolveCode = $LASTEXITCODE
foreach ($Extension in @('sta','msg','dat','log')) {
    if (Test-Path "$JobName.$Extension") { Copy-Item "$JobName.$Extension" $SourceDirectory }
}
if ($SolveCode -ne 0) { throw "Abaqus exited with code $SolveCode" }
if (!(Test-Path "$JobName.sta")) { throw 'Abaqus did not create a status file; inspect the copied dat file' }
if (!(Select-String -Path "$JobName.sta" -Pattern 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' -Quiet)) {
    throw 'Abaqus did not complete the prescribed step'
}
& 'C:\SIMULIA\Commands\abaqus.bat' python extract.py "$JobName.odb"
if ($LASTEXITCODE -ne 0) { throw 'Native field extraction failed' }
Copy-Item ($JobName+'_fields.json') $SourceDirectory
Write-Output "Abaqus artifacts: $Work"
