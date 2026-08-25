param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b50_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b50_hex8_c3d8t_capacity_probe.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b50.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=b50_hex8_c3d8t_capacity_probe `
    input=b50_hex8_c3d8t_capacity_probe.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.0 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b50_hex8_c3d8t_capacity_probe.sta") -or
    !(Select-String -Path "b50_hex8_c3d8t_capacity_probe.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.0 probe did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_b50.py `
    b50_hex8_c3d8t_capacity_probe.odb `
    b50_hex8_c3d8t_capacity.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.0 extraction failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b50_hex8_c3d8t_capacity.csv")) {
    throw "Abaqus B5.0 extraction did not create the reference CSV file"
}
Copy-Item b50_hex8_c3d8t_capacity.csv $SourceDirectory
Write-Output "Abaqus B5.0 work directory: $work"
