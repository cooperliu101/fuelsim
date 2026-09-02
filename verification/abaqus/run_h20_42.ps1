param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory,
    [Parameter(Mandatory = $true)]
    [string]$ProbeDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "h20_42_hex20_pressure_recovery_probe"
$Manifest = "h20_42_hex20_pressure_recovery_manifest.csv"
$Result = "h20_42_hex20_pressure_recovery.csv"
$Work = Join-Path $env:TEMP ("fuelsim_h20_42_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $ProbeDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $ProbeDirectory $Manifest) $Work
Copy-Item (Join-Path $SourceDirectory "extract_h20_42.py") $Work
Set-Location $Work

& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$JobName" `
    "input=$JobName.inp" `
    cpus=1 `
    output_precision=full `
    ask_delete=OFF `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.42 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.42 did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_h20_42.py `
    "$JobName.odb" `
    $Manifest `
    $Result
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.42 extraction failed with exit code $LASTEXITCODE"
}
Copy-Item $Result $ProbeDirectory
Write-Output "Abaqus H20.42 work directory: $Work"
