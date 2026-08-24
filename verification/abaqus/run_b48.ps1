param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b48_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b48_hex8_multi_contact.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b48.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=b48_hex8_multi_contact `
    input=b48_hex8_multi_contact.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.8 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b48_hex8_multi_contact.sta") -or
    !(Select-String -Path "b48_hex8_multi_contact.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.8 did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_b48.py `
    b48_hex8_multi_contact.odb `
    b48_hex8_multi_contact_pair_a.csv `
    b48_hex8_multi_contact_pair_b.csv `
    b48_hex8_multi_contact_reaction.csv `
    b48_hex8_multi_contact_energy.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.8 extraction failed with exit code $LASTEXITCODE"
}
foreach ($result in @(
        "b48_hex8_multi_contact_pair_a.csv",
        "b48_hex8_multi_contact_pair_b.csv",
        "b48_hex8_multi_contact_reaction.csv",
        "b48_hex8_multi_contact_energy.csv")) {
    if (!(Test-Path $result)) {
        throw "Abaqus B4.8 extraction did not create $result"
    }
    Copy-Item $result $SourceDirectory
}
Write-Output "Abaqus B4.8 work directory: $work"
