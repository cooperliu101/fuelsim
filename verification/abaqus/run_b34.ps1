param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b34_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b34_hex8_sliding_contact.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b34.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=b34_hex8_sliding_contact `
    input=b34_hex8_sliding_contact.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B3.4 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b34_hex8_sliding_contact.sta") -or
    !(Select-String -Path "b34_hex8_sliding_contact.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B3.4 did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_b34.py `
    b34_hex8_sliding_contact.odb `
    b34_hex8_sliding_contact_nodes.csv `
    b34_hex8_sliding_contact_contact.csv `
    b34_hex8_sliding_contact_reaction.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B3.4 extraction failed with exit code $LASTEXITCODE"
}
foreach ($result in @(
    "b34_hex8_sliding_contact_nodes.csv",
    "b34_hex8_sliding_contact_contact.csv",
    "b34_hex8_sliding_contact_reaction.csv")) {
    if (!(Test-Path $result)) {
        throw "Abaqus B3.4 extraction did not create $result"
    }
    Copy-Item $result $SourceDirectory
}
Write-Output "Abaqus B3.4 work directory: $work"
