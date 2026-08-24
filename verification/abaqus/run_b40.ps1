param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b40_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b40_hex8_sts_friction.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b40.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=b40_hex8_sts_friction `
    input=b40_hex8_sts_friction.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.0 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "b40_hex8_sts_friction.sta") -or
    !(Select-String -Path "b40_hex8_sts_friction.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.0 did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_b40.py `
    b40_hex8_sts_friction.odb `
    b40_hex8_sts_friction_contact.csv `
    b40_hex8_sts_friction_reaction.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.0 extraction failed with exit code $LASTEXITCODE"
}
foreach ($result in @("b40_hex8_sts_friction_contact.csv", "b40_hex8_sts_friction_reaction.csv")) {
    if (!(Test-Path $result)) {
        throw "Abaqus B4.0 extraction did not create $result"
    }
    Copy-Item $result $SourceDirectory
}
Write-Output "Abaqus B4.0 work directory: $work"
