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
$work = Join-Path $env:TEMP ("fuelsim_b38_probe_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b38_hex8_sts_operator_probe.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b38.py") $work
Set-Location $work

Invoke-Abaqus `
    job=b38_hex8_sts_operator_probe `
    input=b38_hex8_sts_operator_probe.inp `
    output_precision=full `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B3.8 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b38_hex8_sts_operator_probe.sta") -or
    !(Select-String -Path "b38_hex8_sts_operator_probe.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B3.8 probe did not report successful completion"
}

Invoke-Abaqus `
    python extract_b38.py `
    b38_hex8_sts_operator_probe.odb `
    b38_hex8_sts_operator.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B3.8 probe extraction failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b38_hex8_sts_operator.csv")) {
    throw "Abaqus B3.8 probe extraction did not create the reference CSV file"
}

$summary = Get-Content b38_hex8_sts_operator_probe.dat |
    Where-Object {
        $_ -match "PENALTY CONSTRAINT ENFORCEMENT" -or
        $_ -match "SURFACE TO SURFACE WITH THICKNESS" -or
        $_ -match "CONSTRAINT POSITION IS AT NODE" -or
        $_ -match "SUPPLEMENTARY CONSTRAINTS" -or
        $_ -match "NUMBER OF INTERNAL ELEMENTS GENERATED FOR CONTACT"
    } |
    ForEach-Object { $_.Trim() } |
    Select-Object -Unique
[System.IO.File]::WriteAllText(
    (Join-Path $work "b38_hex8_sts_contact_summary.txt"),
    (($summary -join "`n") + "`n"),
    (New-Object System.Text.UTF8Encoding($false)))

Copy-Item b38_hex8_sts_operator.csv $SourceDirectory
Copy-Item b38_hex8_sts_contact_summary.txt $SourceDirectory
Write-Output "Abaqus B3.8 probe work directory: $work"
