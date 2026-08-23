param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_26_probe_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_26_hex20_sts_operator_probe.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_26.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_26_hex20_sts_operator_probe `
    input=h20_26_hex20_sts_operator_probe.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.26 probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_26_hex20_sts_operator_probe.sta") -or
    !(Select-String -Path "h20_26_hex20_sts_operator_probe.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.26 probe did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_26.py `
    h20_26_hex20_sts_operator_probe.odb `
    h20_26_hex20_sts_operator.csv `
    h20_26_hex20_sts_history.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.26 probe extraction failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_26_hex20_sts_operator.csv") -or
    !(Test-Path "h20_26_hex20_sts_history.csv")) {
    throw "Abaqus H20.26 probe extraction did not create both reference CSV files"
}
$summary = Get-Content h20_26_hex20_sts_operator_probe.dat |
    Where-Object {
        $_ -match "PENALTY CONSTRAINT ENFORCEMENT" -or
        $_ -match "SURFACE TO SURFACE WITH THICKNESS" -or
        $_ -match "CONSTRAINT POSITION IS AT NODE" -or
        $_ -match "SUPPLEMENTARY CONSTRAINTS" -or
        $_ -match "NUMBER OF INTERNAL ELEMENTS GENERATED FOR CONTACT"
    } |
    ForEach-Object { $_.Trim() } |
    Select-Object -Unique
if ($summary.Count -ne 5) {
    throw "Abaqus H20.26 contact summary did not contain the five expected identification lines"
}
[System.IO.File]::WriteAllText(
    (Join-Path $work "h20_26_hex20_sts_contact_summary.txt"),
    (($summary -join "`n") + "`n"),
    (New-Object System.Text.UTF8Encoding($false)))

Copy-Item h20_26_hex20_sts_operator.csv $SourceDirectory
Copy-Item h20_26_hex20_sts_history.csv $SourceDirectory
Copy-Item h20_26_hex20_sts_contact_summary.txt $SourceDirectory
Write-Output "Abaqus H20.26 probe work directory: $work"
