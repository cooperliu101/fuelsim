param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b49_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b49_hex8_c3d8t_operator_probe.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b49.py") $work
Set-Location $work

Invoke-Abaqus `
    job=b49_hex8_c3d8t_operator_probe `
    input=b49_hex8_c3d8t_operator_probe.inp `
    output_precision=full `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.9 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b49_hex8_c3d8t_operator_probe.sta") -or
    !(Select-String -Path "b49_hex8_c3d8t_operator_probe.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.9 probe did not report successful completion"
}

Invoke-Abaqus `
    python extract_b49.py `
    b49_hex8_c3d8t_operator_probe.odb `
    b49_hex8_c3d8t_operator_nodal.csv `
    b49_hex8_c3d8t_operator_integration_points.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.9 extraction failed with exit code $LASTEXITCODE"
}
foreach ($artifact in @(
    "b49_hex8_c3d8t_operator_nodal.csv",
    "b49_hex8_c3d8t_operator_integration_points.csv"
)) {
    if (!(Test-Path $artifact)) {
        throw "Abaqus B4.9 extraction did not create $artifact"
    }
    Copy-Item $artifact $SourceDirectory
}
Write-Output "Abaqus B4.9 work directory: $work"
