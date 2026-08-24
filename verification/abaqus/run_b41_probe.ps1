param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b41_probe_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b41_hex8_sts_finite_sliding_operator_probe.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b38.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=b41_hex8_sts_finite_sliding_operator_probe `
    input=b41_hex8_sts_finite_sliding_operator_probe.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.1 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b41_hex8_sts_finite_sliding_operator_probe.sta") -or
    !(Select-String -Path "b41_hex8_sts_finite_sliding_operator_probe.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.1 probe did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_b38.py `
    b41_hex8_sts_finite_sliding_operator_probe.odb `
    b41_hex8_sts_finite_sliding_operator.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.1 probe extraction failed with exit code $LASTEXITCODE"
}

$summary = Get-Content b41_hex8_sts_finite_sliding_operator_probe.dat |
    Where-Object {
        $_ -match "PENALTY CONSTRAINT ENFORCEMENT" -or
        $_ -match "SURFACE TO SURFACE WITH THICKNESS" -or
        $_ -match "CONSTRAINT POSITION IS AT" -or
        $_ -match "SUPPLEMENTARY CONSTRAINTS" -or
        $_ -match "NUMBER OF INTERNAL ELEMENTS GENERATED FOR CONTACT"
    } |
    ForEach-Object { $_.Trim() } |
    Select-Object -Unique
[System.IO.File]::WriteAllText(
    (Join-Path $work "b41_hex8_sts_finite_sliding_contact_summary.txt"),
    (($summary -join "`n") + "`n"),
    (New-Object System.Text.UTF8Encoding($false)))

Copy-Item b41_hex8_sts_finite_sliding_operator.csv $SourceDirectory
Copy-Item b41_hex8_sts_finite_sliding_contact_summary.txt $SourceDirectory
Write-Output "Abaqus B4.1 probe work directory: $work"
