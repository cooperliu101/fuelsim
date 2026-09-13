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
$job = "b531_hex8_c3d8rt_transient_full_field"
$work = Join-Path $env:TEMP ("fuelsim_b531_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b531.py") $work
Set-Location $work

Invoke-Abaqus job=$job input="${job}.inp" output_precision=full cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.31 transient full-field probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.31 transient full-field probe did not report successful completion"
}

Invoke-Abaqus python extract_b531.py "${job}.odb" `
    b531_hex8_c3d8rt_transient_nodal.csv `
    b531_hex8_c3d8rt_transient_integration.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.31 transient full-field extraction failed with exit code $LASTEXITCODE"
}
foreach ($artifact in @(
    "b531_hex8_c3d8rt_transient_nodal.csv",
    "b531_hex8_c3d8rt_transient_integration.csv"
)) {
    if (!(Test-Path $artifact)) {
        throw "Abaqus B5.31 extraction did not create $artifact"
    }
    Copy-Item $artifact $SourceDirectory
}
Write-Output "Abaqus B5.31 work directory: $work"
