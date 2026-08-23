param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_27_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_24_hex20_nonmatching_surface_contact.inp") $work
Copy-Item (Join-Path $SourceDirectory "generate_h20_27.py") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_27.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    python generate_h20_27.py `
    h20_24_hex20_nonmatching_surface_contact.inp `
    h20_27_hex20_nonmatching_sts_operator.inp
if (!(Test-Path "h20_27_hex20_nonmatching_sts_operator.inp")) {
    throw "Abaqus H20.27 input generation failed"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_27_hex20_nonmatching_sts_operator `
    input=h20_27_hex20_nonmatching_sts_operator.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.27 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_27_hex20_nonmatching_sts_operator.sta") -or
    !(Select-String -Path "h20_27_hex20_nonmatching_sts_operator.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.27 solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_27.py `
    h20_27_hex20_nonmatching_sts_operator.odb `
    h20_27_hex20_nonmatching_sts_operator.csv `
    h20_27_hex20_nonmatching_sts_history.csv
if (!(Test-Path "h20_27_hex20_nonmatching_sts_operator.csv") -or
    !(Test-Path "h20_27_hex20_nonmatching_sts_history.csv")) {
    throw "Abaqus H20.27 extraction did not create both CSV files"
}

$summary = Get-Content h20_27_hex20_nonmatching_sts_operator.dat |
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
    throw "Abaqus H20.27 contact summary did not contain the five expected identification lines"
}
[System.IO.File]::WriteAllText(
    (Join-Path $work "h20_27_hex20_nonmatching_sts_contact_summary.txt"),
    (($summary -join "`n") + "`n"),
    (New-Object System.Text.UTF8Encoding($false)))

Copy-Item h20_27_hex20_nonmatching_sts_operator.inp $SourceDirectory
Copy-Item h20_27_hex20_nonmatching_sts_operator.csv $SourceDirectory
Copy-Item h20_27_hex20_nonmatching_sts_history.csv $SourceDirectory
Copy-Item h20_27_hex20_nonmatching_sts_contact_summary.txt $SourceDirectory
Write-Output "Abaqus H20.27 work directory: $work"
